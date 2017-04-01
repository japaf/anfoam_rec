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
import ellipsoidPacking
from numpy import pi, sin, cos, mgrid
#from mayavi import mlab

def show_els(els):
    dphi, dtheta = pi / 50.0, pi / 50.0
    [phi, theta] = mgrid[0:pi + dphi * 1.5:dphi, 0:2 * pi + dtheta * 1.5:dtheta]
    msh=[]
    for el in els:
        pos=el.pos
        for s in el.spheres:
            r=s.r
            relpos=s.relpos
            x = r * sin(phi) * cos(theta)
            y = r * cos(phi)
            z = r * sin(phi) * sin(theta)
            x,y,z=x+pos[0]+relpos[0],y+pos[1]+relpos[1],z+pos[2]+relpos[2]
            #msh.append(mlab.mesh(x, y, z))

    m=1
    ax=[0,0,0,0,m,m,m,m]
    ay=[0,m,0,m,0,m,0,m]
    az=[0,0,m,m,0,0,m,m]
    size=[0.1]*8
    #cube=mlab.points3d(ax,ay,az,size)
    #mlab.outline()
    #mlab.show()

def show_sps(sps):
    dphi, dtheta = pi / 50.0, pi / 50.0
    [phi, theta] = mgrid[0:pi + dphi * 1.5:dphi, 0:2 * pi + dtheta * 1.5:dtheta]
    msh = []
    for sp in sps:
        r = sp.r
        pos = sp.relpos
        x = r * sin(phi) * cos(theta)
        y = r * cos(phi)
        z = r * sin(phi) * sin(theta)
        x, y, z = x + pos[0], y + pos[1], z + pos[2]
        #msh.append(mlab.mesh(x, y, z))
    m = 1
    ax = [0, 0, 0, 0, m, m, m, m]
    ay = [0, m, 0, m, 0, m, 0, m]
    az = [0, 0, m, m, 0, 0, m, m]
    size = [0.1] * 8
    #cube = mlab.points3d(ax, ay, az, size)
    #mlab.outline()
    #mlab.show()


def rotAroundAxis(pos,alfa,id):
    cosa=math.cos(alfa)
    sina=math.sin(alfa)
    id1=(id+1)%3
    id2=(id+2)%3
    pos[id1]=pos[id1]*cosa-pos[id2]*sina
    pos[id2]=pos[id1]*sina+pos[id2]*cosa
    return pos