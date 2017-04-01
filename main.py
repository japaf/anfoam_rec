#!/usr/bin/env python
__author__ = 'jiri'
import json
import logging
import os
import numpy as np
import NeperAnFoam
import foam3Dmodel

import const
#logger setting
#create file handler
log = logging.getLogger('anfoamRec')
log.setLevel(logging.INFO)
# create file handler which logs even debug messages
fh=logging.FileHandler('anfoamRec.log')
fh.setLevel(logging.INFO)
## create console handler
ch=logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
log.addHandler(fh)
log.addHandler(ch)



data={}
mypath=const.work_path
with open('input.json') as data_file:
    data = json.load(data_file)
locals().update(data) # Creates variables from dictionary
#change dir to output
os.chdir(mypath)


if tesselation:

    #MUA = 0  # angle
    SIGMAA = SIGMAA*np.pi / 180
    nS = 4
    #ps=ellipsoidPacking.randSpherePack(0.15,0.01,1,1,1,30,True)
    #sps1=ellipsoidPacking.sphereGrow(sps,1.05)
    #els = ellipsoidPacking.randEllipsoidPack(MUX, MUY, MUZ, MUA, SIGMAX, SIGMAY, SIGMAZ, SIGMAA, 1, 1, 1, nS, 10, True)
    #els2= ellipsoidPacking.ellipsoidGrow(els,1.05)
    model=foam3Dmodel.foamModel()
    #model.initFromGeoFile('foam_dmpRemoved.geo')
    #model.initFromGeoFile('foam_modSimple.geo')# foam_modSimple.geo
    model.initFromGeoFile('foamBig2_final.geo')# foamBig2_final.geo
    model.createPrimitivesSource(10000)
    model.sampleModelSource()
    #model.sampleVoxelsSource(100)

    #Implicit modelling
    #model.createPrimitivesImplicit(10000)
    #model.sampleModelImplicit()
    #model.sampleVoxels(200)
    #Output
    #model.saveAsVox('foam_voxelBig.vtk')
    #model.saveAsSTL('foam_renderedBig.stl')
    model.renderModel()

    #nEl = NeperAnFoam.laguerre3dperiodic(MUX, MUY, MUZ, MUA, SIGMAX, SIGMAY, SIGMAZ, SIGMAA, nS, 0, True, filename)


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
