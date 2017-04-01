#!/usr/bin/env python
__author__ = 'jiri'
#Generate seeds, call Neper and extract foam structure

import logging
import math
import os
import os.path

import numpy as np
import ellipsePacking
import ellipsoidPacking
import geoExtractor
import const
mypath=const.work_path

#load logger
log=logging.getLogger('anfoamRec.Neper')

def laguerre2d(MUX,MUY,SIGMAX,SIGMAY,NumOfCells,filename):
    a=0 #angle
    SIGMAa=np.pi/6
    nCirc=5
    els= ellipsePacking.approxEllipsePack(MUX, MUY, a, SIGMAX, SIGMAY, SIGMAa, 1, 1, nCirc)
    #copy box to make consistent periodic surrounding and prevetn cells from random clipping
    centers_file = os.path.join(mypath, 'centers.txt')
    rads_file = os.path.join(mypath, 'rads.txt')
    log_file = os.path.join(mypath, 'log.txt')
    fc = open(centers_file, 'w')
    fr = open(rads_file, 'w')
    flog = open(log_file, 'w')
    hx = 1
    hy=1
    lcy = 1 * 5
    lcx = 1 * 5
    seeds=0
    for ix in [0,-1,1]:
        for iy in [0,-1,1]:
            for el in els:
                for cir in el.circles:
                    pr=cir.r
                    px=2*hx+el.pos[0]+cir.relpos[0]+ix*hx
                    py=2*hy+el.pos[1]+cir.relpos[1]+iy*hy
                    fc.write('{0:f}\t{1:f}\t{2:f}\n'.format(px, py, 0.0))
                    fr.write('{0:f}\n'.format(pr))
                    seeds+=1
                    if ix == 0 and iy == 0:
                        flog.write('{0:f}\t{1:f}\t{2:f}\n'.format(px - 2*hx, py - 2*hy, 0.0))
                if ix == 0 and iy == 0:
                    drawEllipse(flog,el)
    fc.close()
    fr.close()
    flog.close()
    commandTessellation = "neper -T -n {0:d} -dim 2 -domain 'square({1:f},{2:f})' -periodicity 0 -morphooptiini coo:file\(centers.txt\),weight:file\(rads.txt\) -o {3:s} -format geo,tess,vtk".format(seeds, lcx, lcy,filename)
    # #commandTessellation = "neper -T -n {0:d} -dim 2 -domain 'square({1:f},{2:f})' -morphooptiini coo:file\(centers.txt\),weight:file\(rads.txt\) -o {3:s} -format geo,tess,vtk".format(NumOfCells, lcx, lcy,filename)
    #
    nEl=len(els)
    os.system(commandTessellation)
    geoExtractor.loadGeoFoam2DAndTransform(filename, nEl, nCirc)

def drawEllipse(flog,el):
    px = el.pos[0] + el.circles[0].relpos[0]
    py= el.pos[1] + el.circles[0].relpos[1]
    flog.write('{0:f}\t{1:f}\t{2:f}\n'.format(px, py, 0.0))
    for j in range(100):
        tx = 0.01 * j * np.float(el.ax)
        ty = math.sqrt(-el.by ** 2 / el.ax ** 2 * tx ** 2 + el.by ** 2)
        for ix in [-1,1]:
            for iy in [-1,1]:
                xc=ix*tx
                yc=iy*ty
                x=xc*math.cos(el.ang)-yc*math.sin(el.ang)
                y=xc*math.sin(el.ang)+yc*math.cos(el.ang)
                flog.write('{0:f}\t{1:f}\t{2:f}\n'.format(x+px, y+py, 0.0))

def laguerre3dperiodic(MUX,MUY,MUZ,MUA,SIGMAX,SIGMAY,SIGMAZ,SIGMAA,nS,nEls,sat,filename):
    #els = ellipsoidPacking.approxEllipsoidPack(MUX, MUY, MUZ, MUA, SIGMAX, SIGMAY, SIGMAZ, SIGMAA, 1, 1, 1, nS)

    els = ellipsoidPacking.randEllipsoidPack(MUX, MUY, MUZ, MUA, SIGMAX, SIGMAY, SIGMAZ, SIGMAA, 1, 1, 1, nS, nEls, sat)
    els = ellipsoidPacking.ellipsoidGrow(els, 1.05)
    centers_file = os.path.join(mypath, 'centers.txt')
    rads_file = os.path.join(mypath, 'rads.txt')
    log_file = os.path.join(mypath, 'log3D.txt')
    fc = open(centers_file, 'w')
    fr = open(rads_file, 'w')
    flog = open(log_file, 'w')
    hx = 1
    hy = 1
    hz = 1
    lcy = 1
    lcx = 1
    lcz = 1
    seeds = 0
    origseeds = 0

    for el in els:
        for sp in el.spheres:
            pr = sp.r
            px = el.pos[0] + sp.relpos[0]
            py = el.pos[1] + sp.relpos[1]
            pz = el.pos[2] + sp.relpos[2]
            fc.write('{0:f}\t{1:f}\t{2:f}\n'.format(px, py, pz))
            fr.write('{0:f}\n'.format(pr))
            seeds += 1
    origseeds=seeds
    fc.close()
    fr.close()
    flog.close()
    log.info('Running neper...')
    commandTessellation = "neper -T -n {0:d} -domain 'cube({1:f},{2:f},{3:f})' -periodicity 1 -morphooptiini coo:file\(centers.txt\),weight:file\(rads.txt\) -o {4:s} -format geo,tess,vtk".format(
        seeds, lcx, lcy, lcz, filename)
    commandNeperVisualisation = "neper --rcfile none -V foam.tess -cameraangle 40 -imagesize 800:800 -datacellcol id -datacelltrs 0.5 \\\n -showcell \"id<={0:d}\" -showedge \"cell_shown||(domtype==1)\" \\\n -print selectedfoam".format(
        origseeds)
    nEl = len(els)
    os.system(commandTessellation)

    #os.system(commandNeperVisualisation)
    #geoExtractor.processNeperGeoFoam3D(filename, nEl, nS, False,False)
