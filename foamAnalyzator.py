__author__ = 'jiri'

import math

import numpy as np
from scipy.spatial import ConvexHull

import geoExtractor

def computeCellSize(filename):
    #evaluate max size of the cell in directions of axis x,y,z
    vertices, edges, surfaces, volumes= geoExtractor.loadGeoFoamFile(filename)
    cellSize=[]
    for vol in volumes:
        min = np.array([1e9, 1e9, 1e9])
        max = -1*min
        for surf in vol:
            surfId=abs(surf)-1
            surfL=surfaces[surfId]
            for edge in surfL:
                edgeId=abs(edge)-1
                edgeL=edges[edgeId]
                vId=edgeL[0]-1
                v=vertices[vId]
                for i in range(3):
                    if v[i]>max[i]:
                        max[i]=v[i]
                for i in range(3):
                    if v[i]<min[i]:
                        min[i]=v[i]
        cellSize.append(max-min)
    return cellSize

def analyzeAnizotropyXZXY(cellSizes):
    #evaluate cell anizotropy: everage x/z and x/y size
    Rxy=0
    Rxz=0
    dRxy=0
    dRxz=0
    Raveg=0
    if len(cellSizes)>0:
        for size in cellSizes:
            Rxy+=size[0]/size[1]
            Rxz+=size[0]/size[2]
            dRxy+=(size[0]/size[1])**2
            dRxz+=(size[0]/size[2])**2
        Rxy=Rxy/len(cellSizes)
        Rxz=Rxz/len(cellSizes)
        dRxy=math.sqrt(dRxy/len(cellSizes)-Rxy*Rxy)
        dRxz=math.sqrt(dRxz/len(cellSizes)-Rxz*Rxz)
        Raveg=(Rxy+Rxz)/2
        dRaveg=(dRxy+dRxz)/2
    return [Raveg,dRaveg,Rxy,dRxy,Rxz,dRxz]

def directionBasedSize(filename,nEl,nS):
    #evaluate direction based cell anizotropy
    #uses vector between side subcells of each ellipsoidical cell
    #needs original geo file with ununified cells!
    #nS: total number of spheres per ellipsoid:
    vertices, edges, surfaces, volumes = geoExtractor.loadGeoFoamFile(filename)
    dirCell=[]
    vertByEls=[]
    for i in range(nEl):
        iL=i*nS+nS-1
        iR=i*nS+nS-2
        volL=volumes[iL]
        volR=volumes[iR]
        vertList=[]
        rightAdded=False
        #computeAllvectors from left cell to right cell
        dirVect=np.array([0.0,0.0,0.0])
        for surfL in volL:
            surfIdL=abs(surfL)-1
            edgeListL=surfaces[surfIdL]
            for edge in edgeListL:
                edgeIdL=abs(edge)-1
                vertListL=edges[edgeIdL]
                vIdL=vertListL[0]-1
                vL=vertices[vIdL]
                vertList.append(vL)
                #test with all vertices from the other
                for surfR in volR:
                    surfIdR = abs(surfR) - 1
                    edgeListR = surfaces[surfIdR]
                    for edge in edgeListR:
                        edgeIdR = abs(edge) - 1
                        vertListR = edges[edgeIdR]
                        vIdR = vertListR[0] - 1
                        vR = vertices[vIdR]
                        if not rightAdded:
                            vertList.append(vR)
                        dirVect+=vR-vL
                rightAdded=True
        vertByEls.append(vertList)
        #average direction for cell
        dirCell.append(dirVect/np.linalg.norm(dirVect))
    k=0
    l=0
    cellSize=[]
    for ii in range(nEl): #cell
        min = np.array([1e9, 1e9, 1e9])
        max = -1 * min
        for jj in range(nS):   #subcell
            vol=volumes[k]
            for surf in vol:
                surfId=abs(surf)-1
                surfL=surfaces[surfId]
                for edge in surfL:
                    edgeId=abs(edge)-1
                    edgeL=edges[edgeId]
                    vId=edgeL[0]-1
                    v=vertices[vId]
                    xbase=dirCell[l]
                    [ybase,zbase]=getOrtoBase(xbase)
                    vImg=np.array([np.dot(v,xbase),np.dot(v,ybase),np.dot(v,zbase)])
                    for i in range(3):
                        if vImg[i] > max[i]:
                            max[i] = vImg[i]
                    for i in range(3):
                        if vImg[i] < min[i]:
                            min[i] = vImg[i]
            k+=1
        cellSize.append(max - min)
        l+=1
    return cellSize

def getVolumes(filename,nEl,nS):
    # compute volume for each cell from its subcells
    # needs original geo file with ununified cells!
    # nS: total number of spheres per ellipsoid:
    vertices, edges, surfaces, volumes = geoExtractor.loadGeoFoamFile(filename)

    cellVol = []
    k=0
    for ii in range(nEl):  # cell
        scellVol=0
        for jj in range(nS):  # subcell
            vol = volumes[k]
            vertIdList=[]
            vertList=[]
            for surf in vol:
                surfId = abs(surf) - 1
                surfL = surfaces[surfId]
                for edge in surfL:
                    edgeId = abs(edge) - 1
                    edgeL = edges[edgeId]
                    vId = edgeL[0] - 1
                    addVertex=True
                    for vId0 in vertIdList:
                        if vId0==vId:
                            addVertex=False
                    if addVertex:
                        vertIdList.append(vId)
            for vId0 in vertIdList:
                vertList.append(vertices[vId0])
            scellVol+=ConvexHull(vertList).volume
            k += 1
        cellVol.append(scellVol)
    return cellVol

def getHistogram(nbins,valList):
    min=1e9
    max=-min
    for val in valList:
        if val>max:
            max=val
        if val<min:
            min=val
    h=(max-min)/nbins
    hist=np.zeros(nbins)
    for val in valList:
        i=int(math.floor((val-min)/h))
        if (i>nbins-1):
            i=nbins-1
        hist[i]+=1
    return [hist,min,max,h]

def getOrtoBase(dir):
    #return 2 ortogonal vector to dir and to each other
    #suppose that dir point mainly to axis x
    a=dir[0]
    b=dir[1]
    c=dir[2]
    y1=1.0
    z1=0.1
    z2=1.0
    x1=(-b*y1-c*z1)/a
    x2=(b*z1-c)/(a-b*x1)
    y2=-z1-x1*x2
    ybase=np.array([x1,y1,z1])
    ybase=ybase/np.linalg.norm(ybase)
    zbase=np.array([x2,y2,z2])
    zbase=zbase/np.linalg.norm(zbase)
    return [ybase,zbase]