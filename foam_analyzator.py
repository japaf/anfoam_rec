#!/usr/bin/env python

import logging
import argparse
import subprocess
import json
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

import foam_geoextractor
import common

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"



'''
    #analyze resulting foam
    cellSizes = fa.computeCellSize(filename+"_mod")
    R=fa.analyzeAnizotropyXZXY(cellSizes)
    print("Average R: {0:f} +- {1:f}".format(R[0],R[1]))
    cellSizesDir=fa.directionBasedSize(filename+"_ununified_mod",nEl,2*nS+1)
    Rdir=fa.analyzeAnizotropyXZXY(cellSizesDir)
    print("Average direction based R: {0:f} +- {1:f}".format(Rdir[0],Rdir[1]))
    volumes=fa.getVolumes(filename+"_ununified_mod",nEl,2*nS+1)
    sumVol=0

    eqDiameter=[]
    x=[]
    eqDmax=0
    eqDmin=1
    for vol in volumes:
        sumVol+=vol
    for vol in volumes:
        d=vol*3/4/np.pi/sumVol
        d=math.pow(d,1.0/3)*2
        eqDiameter.append(d)
    # make histogram
    hist,min,max,h=fa.getHistogram(10,eqDiameter)
    volDist_file = os.path.join(mypath, 'eqDiameter_hist.txt')
    volDistFrac_file = os.path.join(mypath, 'eqDiameter_hist_frac.txt')
    volDist_out = open(volDist_file, 'w')
    volDistFrac_out = open(volDistFrac_file, 'w')
    i=0
    volDist_out.write("{0:f}\t{1:f}\t{2:f}\t{3:d}\n".format(min, min, min, 0))
    volDistFrac_out.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\n".format(min, min, min, 0))
    for val in hist:
        l=i*h+min
        i+=1
        r=i*h+min
        m=(r+l)/2
        x.append(m)
        volDist_out.write("{0:f}\t{1:f}\t{2:f}\t{3:d}\n".format(m,l,r,int(val)))
        volDistFrac_out.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\n".format(m, l, r, float(val)/nEl))
    volDist_out.write("{0:f}\t{1:f}\t{2:f}\t{3:d}\n".format(max, max, max, 0))
    volDistFrac_out.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\n".format(max, max, max, 0))
    hist=hist/len(eqDiameter)
    plt.plot(x,hist)
    plt.show()
'''




def computeCellSize(filename):
    #evaluate max size of the cell in directions of axis x,y,z
    vertices, edges, surfaces, volumes= foam_geoextractor.loadGeoFoamFile(filename)
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
    vertices, edges, surfaces, volumes = foam_geoextractor.loadGeoFoamFile(filename)
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
    vertices, edges, surfaces, volumes = foam_geoextractor.loadGeoFoamFile(filename)

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

def execute_se_api_analyzator(file_name,nvolpercell):
    output_file=file_name[:-4]+'.json'
    logging.info('Executing se_api analyzator ...')
    thread = subprocess.Popen(
        ['se_api', '-i',file_name,
        '--all-union', nvolpercell.__str__(),
         '-a',output_file])
    thread.wait()
    logging.info('Executing se_api analyzator ... done')
    return output_file

def analyze_cell_types(cell_types):
    _known_types=[[0,12,0],[2,8,2],[1,10,2],
                  [2,8,3],[3,6,4],[0,12,2],
                  [1,10,3],[2,8,4],[3,6,5],
                  [4,4,6],[0,12,3],[1,10,4],
                  [2,8,5],[0,12,4]]
    results={}
    for ct0 in cell_types:
        for ct1 in _known_types:
            if abs(ct0[0]-ct1[0])<2 and abs(ct0[1]-ct1[1])<2 and abs(ct0[2]-ct1[2])<2:
                key=ct1.__str__()
                if (key not in results):
                    results[key]=[]
                results[key].append(ct0)

    for key,value in results.items():
        print(key,len(value)/len(cell_types)*100)

def analyze_edge_length(length_data,bins):
    min=np.nanmin(length_data)
    max=np.nanmax(length_data)
    print(min,max)
    h=(max-min)/bins
    xv=np.zeros(bins)
    for i in range(bins):
        xv[i]=(2*i+1)*h/2
    hist=np.zeros(bins)
    for length in length_data:
        i=(int)(np.floor((length-min)/h))
        if i>=bins: i=bins-1
        hist[i]+=1
    hist=hist/len(length_data)
    plt.xlabel('length')
    plt.ylabel('partial fraction')
    plt.title('Edge length distribution')
    plt.grid(True)
    plt.plot(xv,hist)
    plt.show()

def analyze_angle(angle_data):
    for vertex in angle_data:
        avg=(np.average(vertex))
        vrc=(np.std(vertex))
        print(avg,vrc)




def main():
    common.init_logging()
    json_foam=execute_se_api_analyzator(args.input_file,args.nvolpercell)
    print(json_foam)
    data=json.load(open(json_foam,'r'))
    if args.all:
        analyze_cell_types(data['cell-types'])
        print('avg facet per cell: ',data['avg-facet-per-cell'])
        print('avg edges per facet: ', data['avg-edges-per-facet'])
        print('vertex-angle: ',data['vertex-angle'])
        #analyze_angle(data['angles-by-vertex'])
        analyze_edge_length(data['length-by-edges'], 10)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        required=True,
                        metavar='FILE',
                        type=str,
                        help="Input geofile for analysis")
    parser.add_argument("-a","--all",
                        action='store_true',
                        help="Analyze all foam features.")
    parser.add_argument("-c","--cell-types",
                        action='store_true',
                        help="Analyze only cell types")
    parser.add_argument("-e", "--edge-length",
                        action='store_true',
                        help="Analyze only edge length distribution")
    parser.add_argument("-n","--nvolpercell",
                        default="1",
                        type=int,
                        help="Number of volumes per one cell, they will be merged.")


    args = parser.parse_args()
    main()