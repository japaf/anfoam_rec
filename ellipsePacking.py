__author__ = 'jiri'

import os
import os.path
import numpy as np
import re
import time
import random
import math
import sys
import json
import const
mypath=const.work_path

class Circ:
    def __init__(self,vr,pos):
        self.r=vr
        self.relpos=pos

class Ellipse:
    def __init__(self,vax,vby,vang,n):
        self.pos=np.array([0.0,0.0])
        self.ax=vax
        self.by=vby
        self.ang=vang
        self.rot=np.array([math.cos(self.ang),math.sin(self.ang)])
        c=Circ(vby,np.array([0.0,0.0]))
        self.circles=[]
        self.circles.append(c)
        self.leftCircle = c
        self.bottomCircle= c
        self.topCicle=c
        self.rightCircle=c
        self.aproxByCircles(n)
    def aproxByCircles(self,n):
        #aprox the ellipse by 2n+1 circles

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
            ryi=xi*math.sin(self.ang)
            rxi=xi*math.cos(self.ang)
            maxbot=self.bottomCircle.relpos[1]-self.bottomCircle.r
            maxtop=self.topCicle.relpos[1]+self.topCicle.r
            c1=Circ(ri, np.array([rxi,ryi]))
            c2=Circ(ri, np.array([-rxi, -ryi]))
            if ryi-ri<maxbot:
                self.bottomCircle=c1
            elif -ryi-ri<maxbot:
                self.bottomCircle=c2
            if ryi+ri>maxtop:
                self.topCicle=c1
            elif -ryi+ri>maxtop:
                self.topCicle=c2
            self.circles.append(c1)
            self.circles.append(c2)
            self.rightCircle=c1
            self.leftCircle=c2
    def testRectInter(self,el1):
        pos1=self.pos
        pos2=el1.pos
        v=pos2-pos1
        maxd = np.linalg.norm(v)
        if (maxd>el1.ax+self.ax):
            return False
        return True
        # cosa1=np.dot(v,self.rot)/np.linalg.norm(v)/np.linalg.norm(self.rot)
        # cosa2 = np.dot(v, el1.rot) / np.linalg.norm(v) / np.linalg.norm(el1.rot)
        # d1=self.getRectDist(cosa1)
        # d2=el1.getRectDist(cosa2)
        # if (maxd<d1+d2):
        #     return True
        # else:
        #     return False

    def testCircInter(self,el1):
        for c1 in self.circles:
            p1 = self.pos + c1.relpos
            for c2 in el1.circles:
                p2=el1.pos+c2.relpos
                v=p2-p1
                maxd=np.linalg.norm(v)
                if (maxd<=c1.r+c2.r):
                    return True
        return False

    def getRectDist(self,cosa):
        if (abs(cosa)>=math.sqrt(2)/2):
            return self.ax/abs(cosa)
        else:
            return self.by/math.sqrt(1-cosa*cosa)


class GeneratorEllipse():
    MUX=1
    MUY=1
    MUA=0
    SX=0
    SY=0
    SA=0
    n=3
    def __init__(self,vMUX,vMUY,vMUA,vSX,vSY,vSA,vn):
        self.MUX=vMUX
        self.MUY = vMUY
        self.MUA = vMUA
        self.SX=vSX
        self.SY=vSY
        self.SA=vSA
        self.n=vn
    def getEl(self):
        ax = abs(np.random.normal(self.MUX, self.SX))
        by = abs(np.random.normal(self.MUY, self.SY))
        a = abs(np.random.normal(self.MUA, self.SA))
        return Ellipse(ax,by,a,self.n)
def randomEllipsePack(MUX,MUY,MUA,SX,SY,SA,xmax,ymax,nCirc):
    els=[]
    random.seed()
    ge = GeneratorEllipse(MUX, MUY, MUA, SX, SY, SA, nCirc)
    minFound = False
    for i in range(50):
        distMin=1e9
        for j in range(50):
            inter = False
            e0 = ge.getEl()
            e0.pos = np.array([(xmax-e0.ax)* random.random(), (ymax-e0.by) * random.random()])
            dist = 0
            for el in els:
                dist += np.linalg.norm(el.pos - e0.pos)
                if (e0.testCircInter(el)):
                    inter = True
            if not inter and len(els) > 0:
                dist = dist / len(els)
                if (dist < distMin):
                    minFound = True
                    distMin = dist
                    eMin = e0
        if minFound:
            minFound = False
            els.append(eMin)
        if len(els) == 0:
            els.append(e0)
    print("Ellipses placed: ", len(els))
    return els

def approxEllipsePack(MUX,MUY,MUA,SX,SY,SA,xmax,ymax,nCirc):
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
    ge=GeneratorEllipse(MUX,MUY,MUA,SX,SY,SA,nCirc)
    #placeBottom(els,ge,xmax)
    #placeLeft(els,ge,ymax)
    minFound=False
    y=0
    while y<ymax-MUY*2:
        x = 0
        while x<xmax-MUX*2:
            distMin=1e9
            for i in range(10):
                e0 = ge.getEl()
                xc = x + e0.ax / 4
                yc = y + e0.by / 4
                inter = False
                e0.pos=np.array([xc+e0.ax*random.random(),yc+e0.by*random.random()])
                dist=0
                for el in els:
                    dist+=np.linalg.norm(el.pos-e0.pos)
                    if (e0.testCircInter(el)):
                        inter=True
                if not inter and len(els)>0:
                    dist=dist/len(els)
                    if (dist<distMin):
                        minFound=True
                        distMin=dist
                        eMin=e0
            if minFound:
                minFound=False
                els.append(eMin)
            if len(els)==0:
                els.append(e0)
            x+=math.cos(e0.ang)*e0.ax*3/4
        y+=MUY*1.5
        print("Ellipses placed: ",len(els))
    return els

def placeBottom(els,ge,xmax):
    #place bottom layer
    n=int(math.ceil(1/ge.MUX/2))
    x=0
    ef=ge.getEl()
    ef.pos=np.array([ef.ax/2,ef.by/2])
    els.append(ef)
    for i in range(n):
        e0=els[len(els)-1]
        e1=ge.getEl()
        rc=e0.rightCircle
        lc=e1.leftCircle



def placeLeft(els,ge,ymax):
    #place left column
    n=1