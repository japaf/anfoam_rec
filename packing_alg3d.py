#!/usr/bin/env python
__author__ = 'jiri'

import os
import logging
import os.path
import numpy as np
import re
import time
import random
import math
import sys
import json
import common


#load logger
log=logging.getLogger('anfoamRec.ellipsoidPacking')

class Sphere:
    def __init__(self,vr,pos):
        self.r=vr
        self.relpos=pos

class Ellipsoid:
    def __init__(self,vax,vby,vcz,vdir,a1,a2,n):
        self.pos=np.array([0.0,0.0,0.0])
        self.ax=vax
        self.by=vby #by=cz
        self.cz=vcz
        self.dir=vdir   #direction vector
        self.ang1=a1
        self.ang2=a2
        self.spheres=[]
        if n%2==1:
            self.spheres.append(Sphere(self.by,self.pos))
        self.aproxBySpheres((int)(n//2))
    def aproxBySpheres(self,n):
        #aprox the ellipsoid by 2n+1 spheres
        if n>0:
            xi=0
            rn = self.by ** 2 / self.ax
            xn = self.ax - rn
            h = xn / n
            for j in range(n):
                xi = xi + h
                if j==n-1:
                    ri=rn
                    xi=xn
                else:
                    ri = math.sqrt(self.by ** 2 * (1 - xi ** 2 / (self.ax ** 2 - self.by ** 2)))
                #rotated coor
                dirx=np.array([1.0,0.0,0.0])
                diry=np.array([0.0,1.0,0.0])
                dirz=np.array([0.0,0.0,1.0])
                rxi=xi*np.dot(self.dir,dirx)
                ryi=xi*np.dot(self.dir,diry)
                rzi=xi*np.dot(self.dir,dirz)
                s1=Sphere(ri, np.array([rxi,ryi,rzi]))
                s2=Sphere(ri, np.array([-rxi, -ryi,-rzi]))
                self.spheres.append(s1)
                self.spheres.append(s2)
    def testFastInter(self,el1):
        if (getDist(self.pos,el1.pos)<=(el1.ax+self.ax)):
            return True
        return False
        # cosa1=np.dot(v,self.rot)/np.linalg.norm(v)/np.linalg.norm(self.rot)
        # cosa2 = np.dot(v, el1.rot) / np.linalg.norm(v) / np.linalg.norm(el1.rot)
        # d1=self.getRectDist(cosa1)
        # d2=el1.getRectDist(cosa2)
        # if (maxd<d1+d2):
        #     return True
        # else:
        #     return False

    def testSphereInter(self,el1):
        for s1 in self.spheres:
            p1 = self.pos + s1.relpos
            for s2 in el1.spheres:
                p2=el1.pos+s2.relpos
                if (getDist(p1,p2)<=(s1.r+s2.r)):
                    return True
        return False

#assuming clipping along faces of the cube
def getDist(vX,vY):
    dv=np.abs(vX-vY)
    dist=np.linalg.norm(dv)
    distRev=np.linalg.norm(dv-np.round(dv))
    return min(dist,distRev)

class GeneratorEllipsoid():
    def __init__(self,vMUX,vMUY,vMUZ,vMUA,vSX,vSY,vSZ,vSA,vn):
        self.MUX=vMUX
        self.MUY = vMUY
        self.MUZ = vMUZ
        self.MUA = vMUA
        self.SX=vSX
        self.SY=vSY
        self.SZ=vSZ
        self.SA=vSA
        self.n=vn
    def getEl(self):
        ax = abs(np.random.normal(self.MUX, self.SX))
        by = abs(np.random.normal(self.MUY, self.SY))
        cz = abs(np.random.normal(self.MUZ, self.SZ))
        alfa = abs(np.random.normal(self.MUA, self.SA))
        beta = abs(np.random.normal(self.MUA, self.SA))
        cosa = math.cos(alfa)
        sina = math.sin(alfa)
        cosb = math.cos(beta)
        sinb = math.sin(beta)
        rot= np.array([cosb*cosa,cosb*sina,sinb])
        return Ellipsoid(ax,by,cz,rot,alfa,beta,self.n)
    def getEl_proporcional(self):
        ax = abs(np.random.normal(self.MUX, self.SX))
        axmuxproportion=ax/self.MUX
        by = self.MUY*axmuxproportion
        cz = by
        alfa = abs(np.random.normal(self.MUA, self.SA))
        beta = abs(np.random.normal(self.MUA, self.SA))
        cosa = math.cos(alfa)
        sina = math.sin(alfa)
        cosb = math.cos(beta)
        sinb = math.sin(beta)
        rot = np.array([cosb * cosa, cosb * sina, sinb])
        return Ellipsoid(ax, by, cz, rot, alfa, beta, self.n)

"""
Ellipsoid packing algorithm - Random sequential adsorption (RSA)
"""
def _randEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,sp_per_el,num_cell):
    '''
    :param MUX:
    :param MUY:
    :param MUA:
    :param SX:
    :param SY:
    :param SA:
    :param xmax:
    :param ymax:
    :return:
    '''
    xmax= ymax= zmax =1.0
    sf=1.0
    if num_cell==-1: #generate non scaled cells
        num_cell=1e9
    else:
        #estimate scale parameter to make sure that specified number of cells can fit into the box
        el_vol=4 / 3 * np.pi * MUX * MUY * MUZ
        max_vol=0.35
        sf=np.power(max_vol/num_cell/el_vol,1.0/3)
    els=[]
    random.seed()
    ge=GeneratorEllipsoid(MUX*sf,MUY*sf,MUZ*sf,MUA*sf,SX*sf,SY*sf,SZ*sf,SA*sf,sp_per_el)
    j = 0
    els_vol = 0
    timeout = time.time() + 4.0+num_cell/50
    while j<num_cell:
        if time.time()>timeout:
            break
        e0 = ge.getEl_proporcional() #consistent also for one sphere per ellipsoids= sphere packing
        for i in range(50):
            e0.pos = np.array([xmax * random.random(), ymax * random.random(), zmax * random.random()])
            inter=False
            for el in els:
                if (e0.testFastInter(el)):
                    if (e0.testSphereInter(el)):
                        inter = True
                        break
            if not inter:
                els.append(e0)
                els_vol += 4 / 3 * np.pi * e0.ax * e0.by * e0.cz
                j+=1
                break
    #mayavi_visu.show_els(els)
    return els,els_vol

def randEllipsoidPack(MUX,MUYZ,MUA,SX,SA,sp_per_el,num_cell):
    els,els_vol=_randEllipsoidPack(MUX, MUYZ, MUYZ, MUA, SX, SX, SX, SA, sp_per_el, num_cell)
    log.info('Volume fraction: %.3f Number of ellipsoids: %d', els_vol, len(els))
    return els

def randEllipsoidPack(MUX,MUYZ,MUA,SX,SYZ,SA,sp_per_el,num_cell):
    els,els_vol=_randEllipsoidPack(MUX, MUYZ, MUYZ, MUA, SX, SYZ, SYZ, SA, sp_per_el, num_cell)
    log.info('Volume fraction: %.3f Number of ellipsoids: %d', els_vol, len(els))
    return els

"""
Sphere packing algorithm - Random sequential adsorption (RSA)
"""
def randSpherePack(MUR,SR,num_cells):
    '''

    :param MUR: average radius
    :param SR: deviation of radius
    :param xmax: x size of the box
    :param ymax: y size of the box
    :param zmax: z size of the box
    :param nSpheres: number of spheres
    :return:
    '''
    els, els_vol = _randEllipsoidPack(MUR, MUR, MUR, 0, SR, SR, SR, 0, 1, num_cells)
    log.info('Volume fraction: %.3f Number of spheres: %d', els_vol, len(els))
    return els

def sphereVol(sps):
    vol=0
    for sp in sps:
        vol += 4 / 3 * np.pi * sp.r ** 3
    return vol

def ellipsoidVol(els):
    vol=0
    for el in els:
        vol+=4 / 3 * np.pi * el.ax * el.by * el.cz
    return vol

def _ellipsoidGrow(els,gf):
    timeout = time.time()  + 4.0+len(els)/50
    while True:
        if time.time()>timeout:
            break
        for i in range(len(els)):
            ax=els[i].ax*gf
            by=els[i].by*gf
            cz=els[i].cz*gf
            e0=Ellipsoid(ax,by,cz,els[i].dir,els[i].ang1,els[i].ang2,len(els[i].spheres))
            e0.pos=els[i].pos
            inter=False
            for j in range(len(els)):
                if ( i != j ):
                    if (e0.testFastInter(els[j])):
                        if (e0.testSphereInter(els[j])):
                            inter = True
                            break
            if not inter:
                els[i]=e0
    logging.info('Inceased volume fraction: %.3f',ellipsoidVol(els))
    #mayavi_visu.show_els(els)
    return els

def sphereGrow(els,gf):
    log.info('Increasing radius of spheres...')
    els=_ellipsoidGrow(els, gf)
    return els

def ellipsoidGrow(els,gf):
    log.info('Increasing radius of ellipsoids...')
    els = _ellipsoidGrow(els, gf)
    return els

'''
Not implemetned
'''

def approxEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,xmax,ymax,zmax,nCirc):
    '''

    :param MUX:
    :param MUY:
    :param MUA:
    :param SX:
    :param SY:
    :param SA:
    :param xmax:
    :param ymax:
    :return:
    '''
    els=[]
    random.seed()
    ge=GeneratorEllipsoid(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,nCirc)
    minFound = False
    z=0
    while z<zmax-MUZ*2.1:
        y = 0
        while y<ymax-MUY*2.1:
            x = 0
            while x<xmax-MUX*2.1:
                distMin=1e9
                e0 = ge.getEl()
                xs = x
                ys = y
                zs = z
                for i in range(50):
                    inter = False
                    e0.pos=np.array([xs+e0.ax*random.random(),ys+e0.by*random.random(),zs+e0.cz*random.random()])
                    dist=0
                    weight=0
                    for el in els:
                        currDist=np.linalg.norm(el.pos-e0.pos)
                        if currDist<e0.ax*4:
                            dist+=currDist
                            weight+=1
                        if (e0.testFastInter(el)):
                            if (e0.testSphereInter(el)):
                                inter=True
                    if not inter and len(els)>0:
                        dist=dist/weight
                        if (dist<distMin):
                            minFound=True
                            distMin=dist
                            minpos=e0.pos
                            eMin=e0
                if minFound:
                    minFound=False
                    eMin.pos=minpos
                    els.append(eMin)
                if len(els)==0:
                    els.append(e0)
                x+=e0.ax*1.5
            y+=MUY*1.5
            print("Ellipsoids placed: ",len(els))
        z+=MUZ*1.5
    #mayavi_visu.show_els(els)
    return els