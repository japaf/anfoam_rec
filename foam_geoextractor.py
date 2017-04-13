#!/usr/bin/env python
__author__ = 'jiri'

import os
import os.path
import logging
import numpy as np
import re
import time
import random
import math
import sys
import json
import packing_alg2d
import packing_alg3d
import common

#load logger
log=logging.getLogger('anfoamRec.GeoExtractor')


def loadGeoFoam2DAndTransform(geofile,nE,nS):
    """
    :param geofile: input geo file
    :param nE: number of macro cells
    :param nS: number of spheres per cell
    :return:
    """
    vertices, edges, surfaces, volumes = loadGeoFoamFile(geofile+".geo")
    edgesToDelete=[]
    cells2d=[]
    #identify and merge subcells that form cell
    k=0
    maxId=0
    maxvId=0
    for i in range(0,nE):
        surf=[]
        for j in range(0,(2*nS+1)):
            for edge in surfaces[k]:
                #for each edge go through surf
                removeAtIndx=-1
                for l in range(0,len(surf)):
                    if (surf[l]==edge or surf[l]==-edge):
                        removeAtIndx=l
                        break;
                if (removeAtIndx>-1):
                    edgesToDelete.append(abs(surf[removeAtIndx]))
                    surf.pop(removeAtIndx)
                else:
                    surf.append(edge)
                    if(abs(edge)>maxId):
                        maxId=abs(edge)
            k+=1
        #arrange surf
        surfAranged=[]
        surfAranged.append(surf[0])
        for j in range(1,len(surf)):
            eId=surfAranged[j-1]
            if (eId > 0):
                vId = edges[eId - 1][1]
            else:  # eId<0
                vId = edges[-eId - 1][0]
            if (edges[abs(eId) - 1][0]>maxvId):
                maxvId=edges[abs(eId) - 1][0]
            if (edges[abs(eId) - 1][1]>maxvId):
                maxvId=edges[abs(eId) - 1][1]
            for l in range(1,len(surf)):
                if (surf[l]==eId):
                    continue
                eeId=surf[l]
                if (edges[abs(eeId) - 1][0]==vId):
                    surfAranged.append(abs(eeId))
                    break;
                if (edges[abs(eeId) - 1][1]==vId):
                    surfAranged.append(-abs(eeId))
                    break;
        cells2d.append(surfAranged)
    #correct indices of edges
    edgeMapping=list(range(maxId))
    shift=0
    j=0
    edgesToDelete.sort()
    for i in range(maxId):
       si=i+shift
       if (j>len(edgesToDelete)-1 or i+1<edgesToDelete[j]):
            edgeMapping[i]=si+1
       else:
            j+=1
            shift-=1
            edgeMapping[i]=-1
    finvertices = list(range(maxvId))
    for i in range(maxvId):
        finvertices[i]=vertices[i]
    finedges = []
    for i in range(maxId):
        if (edgeMapping[i]==-1):
            continue
        finedges.append(edges[i])
    for i in range(len(cells2d)):
        for j in range(len(cells2d[i])):
            edge2d=abs(cells2d[i][j])
            cells2d[i][j]=edgeMapping[abs(edge2d)-1]
    writeAsGeo(finvertices, finedges, cells2d,[], geofile+"_mod.geo")

def loadGeoFoamFile(geofile):
    vertices = []
    edges = []
    surfaces = []
    volumes = []
    ####   Reading the GEO file
    #geofile_path = os.path.join(mypath, geofile)
    geofile_stream = open(geofile, "r")
    lines = geofile_stream.readlines()
    for i in range(0, len(lines)):
        currentline = lines[i].split("{")
        idString = "None"
        values = []
        idNum = 0
        if len(currentline) >= 1:
            idLine = currentline[0].split("(")
            if len(idLine) > 1:
                idStringL = idLine[0].split(' ')
                idStringL = list(filter(None,idStringL))
                idString = idStringL[0]
                if len(idStringL) > 1:
                    for j in range(1, len(idStringL)):
                        idString = idString + ' ' + idStringL[j]
                idNumS = idLine[1].split(")")[0]
                idNum = np.int32(idNumS)
        if len(currentline)>=2:
            values = currentline[1].split("}")[0].split(",")
            valuesArr = np.array(values)
        if (idString == "Point"):
            valuesArr=valuesArr[:3]
            vertices.append(valuesArr.astype(np.float))
        elif (idString == "Line"):
            edges.append(valuesArr.astype(np.int32))
        elif (idString == "Line Loop"):
            surfaces.append(valuesArr.astype(np.int32))
        elif (idString == "Surface Loop"):
            if len(valuesArr)>1:
                volumes.append(valuesArr.astype(np.int32))
            else:
                volumes.append([])
                log.error("Neper failed to generate all cells")
        else:
            continue
    geofile_stream.close()
    return [vertices,edges,surfaces,volumes]

def loadGeoFoamFileAsDict(geofile):
    vertices = {}
    edges = {}
    surfaces = {}
    volumes = {}
    ####   Reading the GEO file
    #geofile_path = os.path.join(mypath, geofile)
    geofile_stream = open(geofile, "r")
    lines = geofile_stream.readlines()
    for i in range(0, len(lines)):
        currentline = lines[i].split("{")
        idString="None"
        values=[]
        idNum=0
        if len(currentline)>=1:
            idLine=currentline[0].split("(")
            if len(idLine)>1:
                idStringL = idLine[0].split(' ')
                idStringL = list(filter(None, idStringL))
                idString=idStringL[0]
                if len(idStringL)>1:
                    for j in range(1,len(idStringL)):
                        idString=idString+' '+idStringL[j]
                idNumS = idLine[1].split(")")[0]
                idNum=np.int32(idNumS)
        if len(currentline)>=2:
            values = currentline[1].split("}")[0].split(",")
            valuesArr = np.array(values)
        if (idString == "Point"):
            valuesArr=valuesArr[:3]
            vertices[idNum]=[valuesArr.astype(np.float),-1]
        elif (idString == "Line"):
            edges[idNum]=[valuesArr.astype(np.int32),-1]
        elif (idString == "Line Loop"):
            surfaces[idNum]=[valuesArr.astype(np.int32),-1]
        elif (idString == "Surface Loop"):
            if len(valuesArr)>1:
                volumes[idNum]=[valuesArr.astype(np.int32),0]
            else:
                volumes[idNum]=[]
                print("Neper failed to generate all cells")
        else:
            continue
    geofile_stream.close()
    return [vertices,edges,surfaces,volumes]

def processNeperGeoFoam3D(geofile,nE,nS,unify,clip):
    """
       :param geofile: input geo file
       :param nE: number of macro cells
       :param nS: number of spheres per cell
       :return:
       """
    # load geo file
    vertices,edges,surfaces,volumes=loadGeoFoamFile(geofile+".geo")
    if clip:
        centers=[]
        centers_file = os.path.join(mypath, 'centers.txt')
        fc_stream = open(centers_file, "r")
        lines = fc_stream.readlines()
        for i in range(len(lines)):
            values = lines[i].split()
            valuesArr = np.array(values)
            centers.append(valuesArr.astype(np.float))
        fc_stream.close()
        clippedGeoFile,clipId=clipCells(centers,vertices,edges,surfaces,volumes,nE,nS)
        dvertices, dedges, dsurfaces, dvolumes = loadGeoFoamFileAsDict(clippedGeoFile)
        volKeys=list(dvolumes.keys())
        j=0
        for i in range(len(dvolumes)-len(clipId),len(dvolumes)):
            #co kdyz bude jinak usporadane?
            oldKey=volKeys[i]
            dvolumes[clipId[j]]=dvolumes[oldKey]
            del(dvolumes[oldKey])
            j+=1

        vertices, edges, surfaces, volumes = orderStructure(dvertices, dedges, dsurfaces, dvolumes)
        volumes = addOrientation(vertices, edges, surfaces, volumes)

    if unify:
        surfaces, volumes=unifyVolumes(surfaces, volumes,nE,nS)
        #edges, surfaces, volumes = unifySurfaces(edges, surfaces, volumes)
        #surfaces=usurfaces
        #volumes=uvolumes
    writeAsGeo(vertices, edges, surfaces, volumes, geofile+"_mod.geo")

def orderStructure(dvertices, dedges, dsurfaces, dvolumes):

    vertices = []
    edges = []
    surfaces = []
    volumes = []

    #extract only used entities
    markEntity(dvolumes,dsurfaces)
    markEntity(dsurfaces,dedges)
    markEntity(dedges,dvertices)
    # create mapping
    assignNumbers(dvolumes)
    assignNumbers(dsurfaces)
    assignNumbers(dedges)
    assignNumbers(dvertices)

    # read mapping and copy to list
    volumes=copyToListViaMapping(dvolumes,dsurfaces)
    surfaces=copyToListViaMapping(dsurfaces,dedges)
    edges=copyToListViaMapping(dedges,dvertices)
    vertices=copyItems(dvertices)

    return [vertices,edges,surfaces,volumes]

def markEntity(dict0,dict1):
    #mark all entities from dict0 in dict1
    for k in sorted(dict0.keys()):
        v=dict0[k]
        if v[1]==0:
            list0=v[0]
            for item0 in list0:
                itemKey=abs(item0)
                dict1[itemKey][1]=0


def assignNumbers(dict0):
    i = 1
    for k in sorted(dict0.keys()):
        if dict0[k][1] == 0:
            dict0[k][1] = i
            i += 1

def copyToListViaMapping(dict0,dict1):
    result=[]
    for k in sorted(dict0.keys()):
        v=dict0[k]
        list0=v[0]
        if (v[1]>=0):
            list1 = []
            for item0 in list0:
                itemKey=abs(item0)
                mappedKey=dict1[itemKey][1]*np.sign(item0)
                list1.append(mappedKey)
            result.append(list1)
    return result

def copyItems(dict0):
    result=[]
    for k in sorted(dict0.keys()):
        v=dict0[k]
        if(v[1]>=0):
            result.append(v[0])
    return result

def addOrientation(vertices,edges,surfaces,volumes):
    for i in range(len(volumes)):
        surfL=volumes[i]
        cp = getCellCenter(vertices, edges, surfaces, surfL)
        for j in range(len(surfL)):
            surfId=abs(surfL[j])
            edgeL=surfaces[surfId-1]
            v0=getVector(edgeL[0],vertices,edges)
            v1=getVector(edgeL[1],vertices,edges)
            v2=getVector(edgeL[2], vertices, edges)
            n01=np.cross(v0,v1)
            n02=np.cross(v1,v2)
            if (np.linalg.norm(n01)>np.linalg.norm(n02)):
                n=n01
            else:
                n=n02
            n=n/np.linalg.norm(n)
            cf=[0.0,0.0,0.0]
            for k in range(len(edgeL)):
                cf+=getVertex(edgeL[k],vertices,edges)
            cf=cf/len(edgeL)
            cv=cp-cf
            cvn=np.dot(cv/np.linalg.norm(cv),n)
            if (cvn>0):
                surfL[j]=-surfId
            else:
                surfL[j] = surfId
        volumes[i]=surfL
    return volumes

def getCellCenter(vertices,edges,surfaces,surfL):
    cp = [0.0, 0.0, 0.0]
    verts=[]
    vertsId=[]
    for i in range(len(surfL)):
        surfId = abs(surfL[i])
        edgeL = surfaces[surfId - 1]
        for j in range(len(edgeL)):
            edge=edges[abs(edgeL[j])-1]
            v0Id=edge[0]
            v1Id=edge[1]
            v0 = vertices[edge[0] - 1]
            v1=vertices[edge[1]-1]
            if v0Id not in vertsId:
                vertsId.append(v0Id)
                verts.append(v0)
            if v1Id not in vertsId:
                vertsId.append(v1Id)
                verts.append(v1)
    for v in verts:
        cp+=v
    return cp/len(verts)

def getVector(edgeId,vertices,edges):
    vertL=edges[abs(edgeId)-1]
    v0=vertices[vertL[0]-1]
    v1 = vertices[vertL[1] - 1]
    if (edgeId>0):
        return v1-v0;
    return v0-v1;

def getVertex(edgeId,vertices,edges):
    startPoint=0
    if (edgeId<0):
        startPoint=1
    vertL=edges[abs(edgeId)-1]
    return vertices[vertL[startPoint]-1]

def unifyVolumes(surfaces,volumes,nE,nS):
    log.info("Unifying volumes..")
    usurfaces=[]
    surfacesToDel=[]
    surfMapping=[]
    uvolumes=[]
    k=0
    for i in range(nE):
        vol=[]
        for j in range((2*nS+1)):
            for l in range(len(volumes[k])):
                addSurf = True
                surf0=volumes[k][l]
                for m in range(len(vol)):
                    if abs(vol[m])==abs(surf0):
                        addSurf=False
                        surfacesToDel.append(abs(surf0))
                        vol.pop(m)
                        break
                if addSurf:
                    vol.append(surf0)
            k+=1
        uvolumes.append(vol)
    k=0
    for i in range(len(surfaces)):
        delSurf=False
        for surf in surfacesToDel:
            if (surf==i+1):
                delSurf=True
                break
        if not delSurf:
            k+=1
            surfMapping.append(k)
            usurfaces.append(surfaces[i])
        else:
            surfMapping.append(-1)
    for i in range(len(uvolumes)):
        for j in range(len(uvolumes[i])):
            surf=uvolumes[i][j]
            uvolumes[i][j]=surfMapping[abs(surf)-1]*np.sign(surf)
    return usurfaces,uvolumes

def unifyEntities(idU,listL,listU):
    #upper dim
    #lower dim
    edgesToDel=[]
    surfL=listU[idU[0] - 1]
    for i in range(len(idU)):
        uL0=listU[idU[i]-1]
        for j in range(i+1,len(idU)):
            uL1=listU[idU[j]-1]
            cL=getCommonItems(uL0,uL1)
            if (len(cL)>0):
                surfL=unifyCycles(cL[0],surfL,uL1)
                edgesToDel.append(cL[0])
            if (len(cL)>1):
                log.error('More than 1 common subEntity!')
    return [surfL,edgesToDel]

def unifyCycles(edgeId,cycle0,cycle1):
    cycle=[]
    for e0 in cycle0:
        if (abs(e0)==edgeId):
            #add second cycle
            #find index of id
            i1=getIndexOfItem(edgeId,cycle1)
            if (i1>-1):
                dir=1
                if (e0*cycle1[i1]>0):
                    dir=-1
                for k in range(len(cycle1)-1):
                    i1 = (i1 + dir) % len(cycle1)
                    cycle.append(dir*cycle1[i1])
        else:
            cycle.append(e0)
    return cycle

def getIndexOfItem(item,list):
    for i in range(len(list)):
        if (abs(list[i])==item):
            return i
    return -1

def unifySurfaces(edges,surfaces,volumes):
    log.info("Unifying surfaces...")

    for i in range(len(volumes)):
        surfL0 = volumes[i]
        for j in range(i+1,len(volumes)):
            surfL1=volumes[j]
            cSurf=getCommonItems(surfL0,surfL1)
            if (len(cSurf)>1):
                newSurf,edgesToDel0=unifyEntities(cSurf,edges,surfaces)
                minSurfId=findMinItem(cSurf)[0]
                cSurf.remove(minSurfId)
                surfaces[minSurfId-1]=newSurf
                for surfId in cSurf:
                    surfaces[surfId-1]=[]
                for edgeId in edgesToDel0:
                    edges[edgeId-1]=[]
                volumes[i]= updateEntity([],cSurf,volumes[i])
                volumes[j] = updateEntity([], cSurf, volumes[j])
    #renumber edges,surfaces
    surfMapping=[0]*len(surfaces)
    edgeMapping=[0]*len(edges)
    is0=1
    ie0=1
    for i in range(len(surfaces)):
        if(len(surfaces[i])>0):
            surfMapping[i]=is0
            is0+=1
    for i in range(len(edges)):
        if(len(edges[i])>0):
            edgeMapping[i]=ie0
            ie0+=1
    #correct
    edges=copyNonEmptyItems(edges)
    surfaces = copyNonEmptyItems(surfaces)
    reMap(surfaces,edgeMapping)
    reMap(volumes,surfMapping)
    return edges, surfaces, volumes

def copyNonEmptyItems(list0):
    list1=[]
    for item in list0:
        if (len(item)>0):
            list1.append(item)
    return list1

def reMap(entities,mapping):
    for i in range(len(entities)):
        list0=reMapList(entities[i],mapping)
        entities[i]=list0
    return entities
def reMapList(list0,mapping):
    for i in range(len(list0)):
        id=list0[i]
        mapId=mapping[abs(id)-1]*np.sign(id)
        list0[i]=mapId
    return list0

def updateEntity(addList,delList,entity):
    for item in delList:
        for i in range(len(entity)):
            if (abs(entity[i])==item):
                entity.pop(i)
                break
    for item in addList:
        entity.append(item)
    return entity


def findMinItem(list0):
    min=sys.maxsize
    imin=0
    for i in range(len(list0)):
        if abs(list0[i])<min:
            min=list0[i]
            imin=i
    return [min,i]

def getCommonItems(list0,list1):
    common=[]
    for item0 in list0:
        for item1 in list1:
            if (abs(item0)==abs(item1)):
                common.append(abs(item0))
    return common

def clipCells(centers,vertices,edges,surfaces,volumes,nE,nS):
    clipId=[]
    clipDir=[]
    for i in range(nE):
        ccId=i*(2*nS+1)
        for j in range(2*nS+1):
            volId = ccId + j
            pos0=centers[volId]
            volL=volumes[volId]
            clip=False
            for fId in range(len(volL)):
                if clip:
                    break;
                ffId=abs(volL[fId])-1
                surfL=surfaces[ffId]
                for sId in range(len(surfL)):
                    if clip:
                        break;
                    eeId=abs(surfL[sId])-1
                    edgeL=edges[eeId]
                    for vId in range(2):
                        vvId=abs(edgeL[vId])-1
                        vert=vertices[vvId]
                        pos1=np.array([vert[0],vert[1],vert[2]])
                        dpos=np.abs(pos0-pos1)
                        if (np.linalg.norm(pos0-pos1)>0.5):
                            clip=True
                            clipId.append(int(volId+1))
                            clipDir.append(np.round(pos0-pos1))
                            #print('Clip, vol:', volId+1,'dir: ',np.round(dpos),'dist: ',np.linalg.norm(pos0-pos1))
                            break;
    #generate geo file
    fg = open('clipFoam.geo', 'w')
    fg.write('Include \"foam.geo\";\n')
    if (len(clipId) >= 1):
        for i in range(len(clipId)):
            dr=clipDir[i]
            fg.write('Translate {{ {0:f}, {1:f}, {2:f} }} {{ Duplicata {{ Volume {{ {3:d} }} ; }} }}\n'.format(dr[0], dr[1], dr[2],clipId[i]))
        fg.write('Delete { Volume {')
        for i in range(len(clipId)-1):
            fg.write('{0:d}, '.format(clipId[i]))
        fg.write('{0:d} }};}}\n'.format(clipId[len(clipId)-1]))
    fg.write('Delete Physicals;\nCoherence;\n')
    fg.close()
    log.info('Running gmsh...')
    cmd_gmsh = "gmsh clipFoam.geo -tol 1e-9 -0"
    os.system(cmd_gmsh)
    return "clipFoam.geo_unrolled",clipId




def writeAsGeo(vertices,edges,surfaces,volumes,filename):
    out = open(filename, 'w')
    for i in range(len(vertices)):
        out.write("Point ({0:d}) = {{{1:2.10f},{2:2.10f},{3:2.10f}}};\n".format(i+1,vertices[i][0],vertices[i][1],vertices[i][2]))
    for i in range(len(edges)):
        out.write("Line ({0:d}) = {{{1:d},{2:d}}};\n".format(i+1,edges[i][0],edges[i][1]))
    for i in range(len(surfaces)):
        out.write("Line Loop ({0:d}) = {{".format(i+1))
        for j in range(len(surfaces[i])-1):
            out.write("{0:d},".format(surfaces[i][j]))
        out.write("{0:d}".format(surfaces[i][len(surfaces[i]) - 1]))
        out.write("};\n")
        out.write("Plane Surface ({0:d}) = {{{0:d}}}; Physical Surface ({0:d}) = {{{0:d}}};\n".format(i + 1))
    k=0
    for i in range(len(volumes)):
        out.write("Surface Loop ({0:d}) = {{".format(k+1))
        if len(volumes[i])<=1:
            continue
        for j in range(len(volumes[i])-1):
            out.write("{0:d},".format(volumes[i][j]))
        out.write("{0:d}".format(volumes[i][len(volumes[i]) - 1]))
        out.write("};\n")
        out.write("Volume ({0:d}) = {{{0:d}}};\n".format(k + 1))
        k+=1
    out.close()