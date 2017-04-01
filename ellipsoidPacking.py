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
import mayavi_visu
import const
mypath=const.work_path

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
        self.spheres.append(Sphere(self.by,self.pos))
        self.aproxBySpheres(n)
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
        #cz = abs(np.random.normal(self.MUZ, self.SZ))
        cz = by #####
        alfa = abs(np.random.normal(self.MUA, self.SA))
        beta = abs(np.random.normal(self.MUA, self.SA))
        cosa = math.cos(alfa)
        sina = math.sin(alfa)
        cosb = math.cos(beta)
        sinb = math.sin(beta)
        rot= np.array([cosb*cosa,cosb*sina,sinb])
        return Ellipsoid(ax,by,cz,rot,alfa,beta,self.n)

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


"""
Ellipsoid packing algorithm - Random sequential adsorption (RSA)
"""
def randEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,xmax,ymax,zmax,nCirc,nEls,sat):
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
    log.info('Generating ellipsoid packing...')
    els=[]
    random.seed()
    ge=GeneratorEllipsoid(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,nCirc)
    z=0
    j = 0
    els_vol = 0
    timeout = time.time() + 2
    while j<=nEls or sat:
        if time.time()>timeout:
            log.warning('TIMEOUT: Can not place all ellipsoids. Placed: %d',j)
            break
        e0 = ge.getEl()
        inter = False
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
    print('Volume fraction: ', els_vol, '\tNumber of ellipsoids: ', len(els))
    #mayavi_visu.show_els(els)
    return els


"""
Sphere packing algorithm - Random sequential adsorption (RSA)
"""
def randSpherePack(MUR,SR,xmax,ymax,zmax,nSpheres,sat):
    '''

    :param MUR: average radius
    :param SR: deviation of radius
    :param xmax: x size of the box
    :param ymax: y size of the box
    :param zmax: z size of the box
    :param nSpheres: number of spheres
    :return:
    '''
    log.info('Generating sphere packing...')
    sps = [] #spheres
    random.seed()
    rad=abs((np.random.normal(MUR,SR,nSpheres)))
    j=0
    sphere_vol = 0
    timeout = time.time() + 5
    while j<=nSpheres or sat:
        if time.time()>timeout:
            log.warning('TIMEOUT: Can not place all spheres. Placed: %d',j)
            break
        pos0=np.array([xmax*random.random(),ymax*random.random(),zmax*random.random()])
        inter=False
        for sp in sps:
            if (getDist(pos0,sp.relpos)<(rad[j]+sp.r)):
                inter=True
                break
        if not inter:
            sps.append(Sphere(rad[j],pos0))
            sphere_vol += 4 / 3 * np.pi * rad[j] ** 3
            j+=1
    #EdgeCubeSize=[math.ceil(MAXcenters[0]-Mincenters[0]),math.ceil(MAXcenters[1]-Mincenters[1]),math.ceil(MAXcenters[2]-Mincenters[2])]
    #EdgeRVESize=int(3.0*max(EdgeCubeSize)) #For NEPER: Size of edge of RVE
    print('Volume fraction: ',sphere_vol,'\tNumber of spheres: ',len(sps))
    #mayavi_visu.show_sps(sps)
    return sps

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

def sphereGrow(sps,gf):
    log.info('Increasing radius of spheres...')
    timeout = time.time() + 5
    while True:
        if time.time()>timeout:
            break
        for i in range(len(sps)):
            r0=sps[i].r*gf
            inter=False
            for j in range(len(sps)):
                if ( i != j ) and (getDist(sps[i].relpos, sps[j].relpos) < (r0 + sps[j].r)):
                    inter = True
                    break
            if not inter:
                sps[i].r=r0
    print(sphereVol(sps))
    #mayavi_visu.show_sps(sps)
    return sps

def ellipsoidGrow(els,gf):
    log.info('Increasing radius of ellipsoids...')
    timeout = time.time() + 5
    while True:
        if time.time()>timeout:
            break
        for i in range(len(els)):
            ax=els[i].ax*gf
            by=els[i].by*gf
            cz=els[i].cz*gf
            e0=Ellipsoid(ax,by,cz,els[i].dir,els[i].ang1,els[i].ang2,(len(els[i].spheres)-1)//2)
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
    print(ellipsoidVol(els))
    #mayavi_visu.show_els(els)
    return els
