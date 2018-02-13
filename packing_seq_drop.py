#!/usr/bin/env python
__author__ = 'jiri1kolar'


import logging
import os.path
import numpy as np
import random
import common
import time
import packing_alg3d as pa
"""
Ellipsoid packing algorithm - Random sequential droping alg
"""
def denseEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,sp_per_el,num_cell):
    logging.info('Generating ellipsoid packing (hexab)...')
    ncell=int(1.0/MUY/2)
    ncellx=int(1.0/MUX/2)
    space = 1.0 / (ncell)
    spacex = 1.0 / (ncellx)
    els = []
    ge = pa.GeneratorEllipsoid(MUX, MUY, MUZ, MUA, SX, SY, SZ, SA, sp_per_el)
    centerx=[]
    centerz=[]
    centery=[]
    els_vol=0
    for i in range(ncellx):
        for j in range(ncell):
            for k in range(ncell):
                if (k % 2 == 0):
                    if (j % 2 == 0):
                        centerx.append (spacex / 4 + spacex * i)
                    else:
                        centerx.append (spacex * 3 / 4 + spacex * i)
                    centery.append (space / 4 + space * j)
                else:
                    if (j % 2 == 0):
                        centerx.append (spacex * 3 / 4 + spacex * i)
                    else:
                        centerx.append (spacex / 4 + spacex * i)
                    centery.append (space * 3 / 4 + space * j)
                centerz.append (space / 4 + space * k);
    for x,y,z in zip(centerx,centery,centerz):
        e0=ge.getEl()
        e0.pos=np.array([x,y,z])
        els_vol += 4 / 3 * np.pi * e0.ax * e0.by * e0.cz
        els.append(e0)
    logging.info('Volume fraction: %.3f Number of ellipsoids: %d', els_vol, len(els))
    return els


def dropEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,sp_per_el,num_cell):
    xmax = ymax = zmax = 1.0
    if num_cell == -1:  # generate all
        num_cell = 1e9
    logging.info('Generating ellipsoid packing (dropping)...')
    els = []
    random.seed()
    ge = pa.GeneratorEllipsoid(MUX, MUY, MUZ, MUA, SX, SY, SZ, SA, sp_per_el)
    j = 0
    els_vol = 0
    ntrials=10
    dx = 0.1
    maxx=-1
    jprev=-1
    timeout = time.time() + 6000
    while j < num_cell:
        if jprev==j:
            break
        jprev=j
        if time.time() > timeout:
            break
        e0 = ge.getEl()
        min_pos=np.array([2,0,0])
        add_new_el=False
        #drop
        for i in range(ntrials):
            e0.pos = np.array([0, ymax * random.random(), zmax * random.random()])
            #drop ellipsoid
            pos=rise_ellipsoid(e0,els,dx)
            if pos[0]<min_pos[0]:
                min_pos=pos.copy()
        if not intersect(e0,els):
            e0.pos=min_pos.copy()
            if e0.pos[0]>maxx:
                maxx=e0.pos[0]
            els.append(e0)
            els_vol += 4 / 3 * np.pi * e0.ax * e0.by * e0.cz
            j += 1
        print(j, maxx)
    logging.info('Volume fraction: %.3f Number of ellipsoids: %d', els_vol, len(els))

    return els

def rise_ellipsoid(e0,els,dx):
    for exp in range(0,5):
        cdx=dx/2**exp
        while (e0.pos[0]<1):
            if intersect(e0,els):
                e0.pos[0] += cdx
            else:
                break
        if e0.pos[0]-cdx>0:
            e0.pos[0] -= cdx
    e0.pos[0] += cdx
    return e0.pos

def intersect(e0,els):
    for el in els:
        if (e0.testFastInter(el)):
            if (e0.testSphereInter(el)):
                return True
    return False