#!/usr/bin/env python
__author__ = 'jiri1kolar'
#Generate seeds, call Neper and extract foam structure

import logging
import vtk
from vtk.util.colors import *
import math
import argparse
import common
import numpy as np
import foam_geoextractor

#load loggingger
class foam_model:
    def __init__(self):
        self.th=1e-6
        #3D geometry
        #primitives
        self.primVertices=[]
        self.primEdges=[]
        self.primObj=[]
        #status
        self.prim_status=False
        self.primSource_status = False
        self.sample_status=False
        self.render_status=False
        self.voxsample_status=False
        # outline visualisation
        outline = vtk.vtkOutlineFilter()
        points = vtk.vtkPoints()
        points.InsertNextPoint([.0, .0, .0])
        points.InsertNextPoint([1., 1., 1.])
        box = vtk.vtkPolyData()
        box.SetPoints(points)
        outline.SetInputData(box)
        outlineMapper = vtk.vtkPolyDataMapper()
        outlineMapper.SetInputConnection(outline.GetOutputPort())
        outlineActor = vtk.vtkActor()
        outlineActor.SetMapper(outlineMapper)
        outlineActor.GetProperty().SetColor([0, 0, 0])
        self.outlineActor = outlineActor

    def initFromGeoFile(self,geofile):
        v, e, s, vol = foam_geoextractor.loadGeoFoamFile(geofile)
        self.vertexList = v
        self.edgeList = e
        self.surfaceList = s
        self.volumeList = vol

    def initFromList(self,v,e,s,vol):
        self.vertexList = v
        self.edgeList = e
        self.surfaceList = s
        self.volumeList = vol

    def saveAsGeo(self,filename):
        foam_geoextractor.writeAsGeo(self.vertexList, self.edgeList, self.surfaceList, self.volumeList, filename)

    def createPrimitivesSource(self,imax):
        logging.info("Creating primitive objects...")
        self.primSource_status = True
        rVertex = 0.015
        rEdge = 0.01
        dSurface=0.001
        i = 0
        self.resSource=10
        sceneSource=vtk.vtkAppendPolyData()
        for vert in self.vertexList:
            source = vtk.vtkSphereSource()
            source.SetCenter(vert)
            source.SetRadius(rVertex)
            sceneSource.AddInputConnection(source.GetOutputPort())
            i += 1
            if i > imax:
                break
        i=0
        if imax==-1:
            imax=1e9

        for edge in self.edgeList:
            v0 = self.vertexList[edge[0] - 1]
            v1 = self.vertexList[edge[1] - 1]
            source = self.makeConeEdgeSource(v0, v1, rEdge)
            sceneSource.AddInputConnection(source.GetOutputPort())
            [v0c, v1c, dv] = self.clipEdge(edge)
            if (np.linalg.norm(dv) > self.th):
                source = self.makeConeEdgeSource(v0c, v1c, rEdge)
                sceneSource.AddInputConnection(source.GetOutputPort())
            i += 1
            if i > imax:
                break
        i=0
        for surface in self.surfaceList:
            dv=np.array([.0,.0,.0])
            source = self.makeFaceSource(surface,dSurface,dv)
            sceneSource.AddInputConnection(source.GetOutputPort())
            dv=self.getSurfaceClip(surface)
            if (np.linalg.norm(dv)>self.th):
                source = self.makeFaceSource(surface, dSurface, dv)
                sceneSource.AddInputConnection(source.GetOutputPort())
            i += 1
            if i > imax:
                break

        sceneSource.Update()
        self.sceneSource=sceneSource



    def createPrimitivesImplicit(self,imax):
        #create primitives at each vertex, edge
        logging.info("Creating primitive objects from implicit functions...")
        self.prim_status=True
        rVertex=0.025
        rEdge=0.024
        dSurface=0.015
        i=0
        if imax==-1:
            imax=1e9

        for vert in self.vertexList:
            obj=vtk.vtkSphere();
            obj.SetCenter(vert)
            obj.SetRadius(rVertex)
            self.primVertices.append(obj)
            i+=1
            if i>imax:
                break
        i = 0
        for edge in self.edgeList:
            v0 = self.vertexList[edge[0] - 1]
            v1 = self.vertexList[edge[1] - 1]
            #obj=self.makeEdgeImplicit(v0,v1,rEdge)
            #self.primEdges.append(obj)
            obj=self.makeConeEdgeImplicit(v0,v1,rEdge)
            self.primEdges.append(obj)
            [v0c,v1c,dv]=self.clipEdge(edge)
            if (np.linalg.norm(dv)>self.th):
                #obj = self.makeEdgeImplicit(v0c, v1c, rEdge)
                obj = self.makeConeEdgeImplicit(v0, v1, rEdge)
                #self.primEdges.append(obj)
            i += 1
            if i > imax:
                break
        i=0
        for surface in self.surfaceList:
            dv = np.array([.0, .0, .0])
            obj = self.makeSurfaceImplicit(surface, dSurface, dv)
            self.primObj.append(obj)
            dv = self.getSurfaceClip(surface)
            if (np.linalg.norm(dv) > self.th):
                obj = self.makeSurfaceImplicit(surface, dSurface, dv)
                self.primObj.append(obj)
            i+=1
            if i > imax:
                break

        scene=vtk.vtkImplicitBoolean()
        scene.SetOperationTypeToUnion()
        for obj in self.primVertices:
            scene.AddFunction(obj)
        for obj in self.primEdges:
            scene.AddFunction(obj)
        for obj in self.primObj:
            scene.AddFunction(obj)
        self.scene=scene

    def sampleModelSource(self):
        if not self.primSource_status:
            logging.error("There are no primitive structures to sample.")
            return
        self.sample_status=True

        logging.info("Sampling model of implicit functions...")

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputConnection(self.sceneSource.GetOutputPort())
        for i in range(3):
            for j in [-1,1]:
                n=np.zeros(3)
                n[i]=j
                c=(np.ones(3)-j*np.ones(3))/2
                cp=vtk.vtkPlane()
                cp.SetOrigin(c)
                cp.SetNormal(n)
                mapper.AddClippingPlane(cp)
        mapper.ScalarVisibilityOff()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(slate_grey_dark)
        #triangulation for output
        triangleFilter=vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(self.sceneSource.GetOutputPort())
        triangleFilter.Update()
        # save
        self.outputData=triangleFilter
        self.mapper=mapper
        self.actor=actor


    def sampleModelImplicit(self):
        if not self.prim_status:
            logging.error("There are no primitive structures to sample.")
            return
        self.sample_status = True

        logging.info("Sampling model of implicit functions...")

        res = 200  # resolution
        surf = vtk.vtkContourFilter()
        sample = vtk.vtkSampleFunction()
        sample.SetImplicitFunction(self.scene)
        #sample.SetModelBounds(-1, 2, -1, 2, -1, 2)
        sample.SetModelBounds(0, 1, 0, 1, 0, 1)
        sample.SetSampleDimensions(res, res, res)
        sample.ComputeNormalsOff()
        surf.SetInputConnection(sample.GetOutputPort())
        surf.SetValue(0, 0.0)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(surf.GetOutputPort())
        mapper.ScalarVisibilityOff()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        # actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetColor(alizarin_crimson)
        # save
        self.mapper = mapper
        self.surf = surf
        self.outputData=surf
        self.actor = actor



    def sampleVoxels(self,res):
        if not self.prim_status:
            logging.error("There are no primitive structures to sample.")
            return
        self.voxsample_status=True
        logging.info("Sampling voxel model...")

        sample = vtk.vtkSampleFunction()
        sample.SetImplicitFunction(self.scene)
        sample.SetModelBounds(0, 1, 0, 1, 0, 1)
        sample.SetSampleDimensions(res, res, res)
        sample.ComputeNormalsOff()
        extract=vtk.vtkExtractGeometry()
        extract.SetInputConnection(sample.GetOutputPort())
        extract.SetImplicitFunction(self.scene)

        vox = vtk.vtkVoxelModeller()
        vox.SetSampleDimensions(res, res, res)
        vox.SetModelBounds(0, 1, 0, 1, 0, 1)
        vox.SetScalarTypeToFloat()
        vox.SetMaximumDistance(.1)
        vox.SetInputConnection(extract.GetOutputPort())
        #self.vox=extract
        dataMapper = vtk.vtkDataSetMapper()
        dataMapper.SetInputConnection(extract.GetOutputPort())
        dataActor = vtk.vtkActor()
        dataActor.SetMapper(dataMapper)
        #save
        self.vox = vox
        self.voxActor = dataActor

    def sampleVoxelsSource(self,res):
        if not self.primSource_status:
            logging.error("There are no primitive structures to sample.")
            return
        self.voxsample_status=True
        logging.info("Sampling model to voxels...")

        vox = vtk.vtkVoxelModeller()
        vox.SetSampleDimensions(res, res, res)
        vox.SetModelBounds(0, 1, 0, 1, 0, 1)
        vox.SetScalarTypeToFloat()
        vox.SetMaximumDistance(.1)
        vox.SetInputConnection(self.sceneSource.GetOutputPort())
        #save voxels
        self.vox = vox



    def renderModel(self):
        # Create the usual rendering stuff
        if not self.sample_status:
            logging.error("There are no sampled data to render.")
            return
        logging.info("Rendering model...")
        self.render_status=True
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        # Add the actors to the renderer, set the background and size

        ren.AddActor(self.actor)
        ren.AddActor(self.outlineActor)
        ren.SetBackground(1, 1, 1)
        renWin.SetSize(1500, 1400)
        ren.ResetCamera()
        ren.GetActiveCamera().Roll(90)
        ren.GetActiveCamera().Dolly(1.5)
        ren.ResetCameraClippingRange()

        iren.Initialize()
        renWin.Render()
        iren.Start()

    def renderVox(self):
        # Create the usual rendering stuff
        if not self.voxsample_status:
            logging.error("There are no voxel data to render.")
            return
        logging.info("Rendering voxel model...")
        self.render_status = True
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        # Add the actors to the renderer, set the background and size

        ren.AddActor(self.voxActor)
        ren.AddActor(self.outlineActor)
        ren.SetBackground(1, 1, 1)
        renWin.SetSize(1500, 1400)
        ren.ResetCamera()
        ren.GetActiveCamera().Roll(90)
        ren.GetActiveCamera().Dolly(1.5)
        ren.ResetCameraClippingRange()

        iren.Initialize()
        renWin.Render()
        iren.Start()

    '''
        Building functions
    '''
    def makeFaceSource(self,edgeList,d,dv):
        # compute geometric center of surface
        c = np.array([.0, .0, .0])
        points = vtk.vtkPoints()
        normals = vtk.vtkPoints()
        for edgeId in edgeList:
            edge = self.edgeList[abs(edgeId) - 1]
            v0 = self.vertexList[edge[0] - 1]
            if (edgeId < 0):
                v0 = self.vertexList[edge[1] - 1]
            c += v0
        # comupte normal to the facet
        p0 = self.vertexList[edge[0] - 1]
        p1 = self.vertexList[edge[1] - 1]
        c = c / len(edgeList)
        v0 = p0 - c
        v1 = p1 - c
        n = np.cross(v0, v1)
        n = n / np.linalg.norm(n)
        points=vtk.vtkPoints()
        poly=vtk.vtkCellArray()
        poly.InsertNextCell(len(edgeList))
        i=0
        for edgeId in edgeList:
            edge=self.edgeList[abs(edgeId)-1]
            v0 = self.vertexList[edge[0] - 1]
            if (edgeId<0):
                v0 = self.vertexList[edge[1] - 1]
            v1=v0+dv
            points.InsertNextPoint(v1-0.5*d*n)
            poly.InsertCellPoint(i)
            i+=1
        profile=vtk.vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(poly)
        extrude=vtk.vtkLinearExtrusionFilter()
        extrude.SetInputData(profile)
        #extrude.SetExtrusionTypeToVectorExtrusion()
        extrude.SetExtrusionTypeToNormalExtrusion()
        extrude.SetVector(n)
        extrude.SetScaleFactor(d)
        extrude.Update()
        normals=vtk.vtkPolyDataNormals()
        normals.SetInputConnection(extrude.GetOutputPort())
        normals.SetFeatureAngle(60)
        return normals

    def makeEdgeSource(self,v0,v1,rEdge):
        source = vtk.vtkLineSource()
        c = (v0 + v1) / 2.0
        n = (v1 - v0)
        n = n / np.linalg.norm(n)
        source.SetPoint1(v0)
        source.SetPoint2(v1)
        filter=vtk.vtkTubeFilter()
        filter.SetRadius(rEdge)
        filter.AddInputConnection(source.GetOutputPort())
        return filter

    def makeConeEdgeSource(self,v0,v1,rEdge):
        c = (v0 + v1) / 2.0
        n = (v1 - v0)
        n = n / np.linalg.norm(n)
        ang=100
        cone0 = vtk.vtkConeSource()
        cone0.SetCenter(c)
        cone0.SetDirection(n)
        cone0.SetHeight(np.linalg.norm(v1-v0))
        cone0.SetRadius(rEdge)
        cone0.SetResolution(self.resSource)
        cone1 = vtk.vtkConeSource()
        cone1.SetCenter(c)
        cone1.SetDirection(-n)
        cone1.SetRadius(rEdge)
        cone1.SetHeight(np.linalg.norm(v1-v0))
        cone1.SetResolution(self.resSource)
        filter=vtk.vtkAppendPolyData()
        filter.AddInputConnection(cone0.GetOutputPort())
        filter.AddInputConnection(cone1.GetOutputPort())
        return filter


    def makeSurfaceImplicit(self,edgeList,d,dv):
        #compute geometric center of surface
        radMax=1.0
        c=np.array([.0,.0,.0])
        i=0
        for edgeId in edgeList:
            edge = self.edgeList[abs(edgeId) - 1]
            v0 = self.vertexList[edge[0] - 1]
            if (edgeId < 0):
                v0 = self.vertexList[edge[1] - 1]
            c+=v0
            i+=1
        # comupte normal to the facet
        p0 = self.vertexList[edge[0] - 1]
        p1 = self.vertexList[edge[1] - 1]
        c=c/i
        v0=p0-c
        v1=p1-c
        n=np.cross(v0,v1)
        n=n/np.linalg.norm(n)
        c+=dv
        #add planes
        pl0=vtk.vtkPlane()
        pl0.SetOrigin(c-d/2*n)
        pl0.SetNormal(n)
        pl1 = vtk.vtkPlane()
        pl1.SetOrigin(c + d / 2 * n)
        pl1.SetNormal(-n)
        sp=vtk.vtkSphere()
        sp.SetCenter(c)
        sp.SetRadius(radMax)
        #add edges as surfaces
        obj = vtk.vtkImplicitBoolean()
        obj.SetOperationTypeToDifference()
        obj.AddFunction(sp)
        obj.AddFunction(pl0)
        obj.AddFunction(pl1)
        for edgeId in edgeList:
            edge=self.edgeList[abs(edgeId)-1]
            tv0 = self.vertexList[edge[0] - 1]
            tv1 = self.vertexList[edge[1] - 1]
            if (edgeId < 0):
                tv0 = self.vertexList[edge[1] - 1]
                tv1=self.vertexList[edge[0]-1]
            v0=tv0+dv
            v1=tv1+dv
            cEdge=(v0+v1)/2
            ve0=v1-v0
            ne=np.cross(ve0,n)
            ve1=cEdge-c
            ve1=ve1/np.linalg.norm(ve1)
            ne=ne/np.linalg.norm(ne)
            ve2 = ve1 + ne
            ple = vtk.vtkPlane()
            if (np.linalg.norm(ve2)>1):
                ple.SetNormal(-ne)
            else:
                ple.SetNormal(ne)
            ple.SetOrigin(cEdge)
            obj.AddFunction(ple)
        return obj

    def makeEdgeImplicit(self,v0,v1,rEdge):
        cylinder = vtk.vtkCylinder()
        c = (v0 + v1) / 2.0
        n = (v1 - v0)
        n = n / np.linalg.norm(n)
        cylinder.SetCenter(c)
        cylinder.SetRadius(rEdge)
        cylinder.SetAxis(n)
        p0 = vtk.vtkPlane()
        p0.SetOrigin(v0)
        p0.SetNormal(-n)
        p1 = vtk.vtkPlane()
        p1.SetOrigin(v1)
        p1.SetNormal(n)
        obj = vtk.vtkImplicitBoolean()
        obj.SetOperationTypeToIntersection()
        obj.AddFunction(cylinder)
        obj.AddFunction(p0)
        obj.AddFunction(p1)
        return obj

    def makeConeEdgeImplicit(self,v0,v1,rEdge):
        c = (v0 + v1) / 2.0
        n = (v1 - v0)
        n = n / np.linalg.norm(n)
        dir=np.array([1.,0.,0.])
        axis=np.cross(dir,-n)
        angle=-math.acos(np.dot(dir,-n))/math.pi*180.0
        ang = math.atan(rEdge/np.linalg.norm(v1-v0))/np.pi*180
        #set cones
        #cone0
        cone0 = vtk.vtkCone()
        cone0.SetAngle(ang)
        t0 = vtk.vtkTransform()
        t0.RotateWXYZ(angle,axis)
        t0.Translate(-v0)
        cone0.SetTransform(t0)
        #cone1
        cone1 = vtk.vtkCone()
        cone1.SetAngle(ang)
        t1=vtk.vtkTransform()
        t1.RotateWXYZ(angle,axis)
        t1.Translate(-v1)
        cone1.SetTransform(t1)
        #combine
        cone=vtk.vtk.vtkImplicitBoolean()
        cone.SetOperationTypeToUnion()
        cone.AddFunction(cone0)
        cone.AddFunction(cone1)
        p0 = vtk.vtkPlane()
        p0.SetOrigin(v0)
        p0.SetNormal(-n)
        p1 = vtk.vtkPlane()
        p1.SetOrigin(v1)
        p1.SetNormal(n)
        obj = vtk.vtkImplicitBoolean()
        obj.SetOperationTypeToIntersection()
        obj.AddFunction(cone)
        obj.AddFunction(p0)
        obj.AddFunction(p1)
        return obj

    def getRotAngle(self,dir0,dir1,i):
        dir=dir1+np.zeros(3)
        dir[i]=0
        cosa=np.dot(dir0,dir1)/np.linalg.norm(dir0)/np.linalg.norm(dir)
        return math.atan(dir0[0]/dir1[1])/math.pi*180

    def getSurfaceClip(self,edgeList):
        dx=.0
        dy=.0
        dz=.0
        for edgeId in edgeList:
            edge = self.edgeList[abs(edgeId) - 1]
            v0 = self.vertexList[edge[0] - 1]
            if (edgeId < 0):
                v0 = self.vertexList[edge[1] - 1]
            if (abs(dx)<self.th):
                dx+=self.getVertexClip(v0[0])
            if (abs(dy) < self.th):
                dy+=self.getVertexClip(v0[1])
            if (abs(dz) < self.th):
                dz+=self.getVertexClip(v0[2])
        return np.array([dx,dy,dz])

    def clipEdge(self,edge):
        v0 = self.vertexList[edge[0] - 1]
        v1 = self.vertexList[edge[1] - 1]
        dx = self.getEdgeClip(v0[0],v1[0])
        dy = self.getEdgeClip(v0[1], v1[1])
        dz = self.getEdgeClip(v0[2], v1[2])
        dv=np.array([dx,dy,dz])
        return [v0+dv,v1+dv,dv]


    def getEdgeClip(self,x0,x1):
        dx0=self.getVertexClip(x0)
        dx1=self.getVertexClip(x1)
        if (abs(dx0)>self.th):
            return dx0
        return dx1

    def getVertexClip(self,x0):
        if x0<0:
            return 1.0
        elif x0>1:
            return -1.0
        return 0

    '''
        Writers to files
    '''
    def saveAsSTL(self,fileName):
        if not self.sample_status:
            logging.error("There are no sampled data to write to file.")
            return
        logging.info("Saving model to %s in stl format...",fileName)
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName(fileName)
        stlWriter.SetInputConnection(self.outputData.GetOutputPort())
        stlWriter.Write()

    def saveAsPLY(self,fileName):
        if not self.sample_status:
            logging.error("There are no sampled data to write to file.")
            return
        logging.info("Saving model to %s in ply format...", fileName)
        plyWriter = vtk.vtkPLYWriter()
        plyWriter.SetFileName(fileName)
        plyWriter.SetInputConnection(self.outputData.GetOutputPort())
        plyWriter.Write()

    def saveAsVox(self,fileName):
        if not self.voxsample_status:
            logging.error("There are no sampled voxel data to write to file.")
            return
        logging.info("Saving model as voxel 3D image to file %s in vtk format...",fileName)
        voxWriter=vtk.vtkDataSetWriter()
        voxWriter.SetFileName(fileName)
        voxWriter.SetInputConnection(self.vox.GetOutputPort())
        voxWriter.Write()

def main():
    common.init_logging()
    if args.input_file is None or args.input_file == "":
        logging.error("Specify valid geo input-file.")
        return 0
    model = foam_model()
    model.initFromGeoFile(args.input_file)
    # Modeling
    if args.object_method == 'implicit':
        model.createPrimitivesImplicit(args.limit_num_obj)
        model.sampleModelImplicit()
        if (args.save_vox is not None):
            model.sampleVoxels(args.voxel_resolution)
    elif args.object_method == 'source':
        model.createPrimitivesSource(args.limit_num_obj)
        model.sampleModelSource()
        if (args.save_vox is not None):
            model.sampleVoxelsSource(args.voxel_resolution)
    # Output
    if args.save_vox is not None:
        model.saveAsVox(args.save_vox)
    if args.save_stl is not None:
        model.saveAsSTL(args.save_stl)
    if args.display:
        model.renderModel()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        help="Input geo file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-d", "--dimension",
                        help="Dimension of structure",
                        type=int,
                        default=3)
    parser.add_argument("--object-method",
                        help="Object method for structure modelling.",
                        choices=['implicit', 'source'],
                        default="implicit")
    parser.add_argument("-r", "--voxel-resolution",
                        help="Resolution for voxel output",
                        type=int,
                        default=100)
    parser.add_argument("--limit-num-obj",
                        help="Limit maximal number of rendered objects of each type: vertex,edge,surface.",
                        type=int,
                        default=-1)
    parser.add_argument("-disp", "--display",
                        help="Display results",
                        action='store_true')
    parser.add_argument("--save-vox",
                        help="Save voxel vtk format",
                        metavar='FILE',
                        type=str)
    parser.add_argument("--save-stl",
                        help="Save in stl format",
                        metavar='FILE',
                        type=str)
    args = parser.parse_args()
    main()