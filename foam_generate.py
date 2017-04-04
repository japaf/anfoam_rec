#!/usr/bin/env python

import logging
import packing_alg3d
import subprocess
import os

__name__ = 'anfoam_rec.foam_generate'
__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger
log=logging.getLogger(__name__)

def generate_3d(output_file,gen_options):
    MUX=gen_options['MUX']
    MUY=gen_options['MUY']
    MUZ=gen_options['MUZ']
    MUA=gen_options['MUA']
    SX=gen_options['SIGMAX']
    SY=gen_options['SIGMAY']
    SZ=gen_options['SIGMAZ']
    SA=gen_options['SIGMAA']
    sp_per_el=gen_options['spheres-per-ellipsoid']
    num_cell=gen_options['num-cell']
    ellipsoids=[]
    if gen_options['packing-alg']=='ellipsoid-random':
        ellipsoids = packing_alg3d.randEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,sp_per_el,num_cell)
        if gen_options['grow-alg']:
            ellipsoids = packing_alg3d.ellipsoidGrow(ellipsoids, 1.05)

    elif gen_options['packing-alg']=='sphere-random':
        MUR=gen_options['MUR']
        SR=gen_options['SIGMAR']
        ellipsoids = packing_alg3d.randEllipsoidPack(MUR,MUR,MUR,0,SR,SR,SR,0,1,num_cell)
        if gen_options['grow-alg']:
            ellipsoids = packing_alg3d.ellipsoidGrow(ellipsoids, 1.05)
    center_file='centers.txt'
    rads_file='rads.txt'
    fc = open(center_file, 'w')
    fr = open(rads_file, 'w')
    seeds = 0
    for el in ellipsoids:
        for sp in el.spheres:
            pr = sp.r
            px = el.pos[0] + sp.relpos[0]
            py = el.pos[1] + sp.relpos[1]
            pz = el.pos[2] + sp.relpos[2]
            fc.write('{0:f}\t{1:f}\t{2:f}\n'.format(px, py, pz))
            fr.write('{0:f}\n'.format(pr))
            seeds += 1
    fc.close()
    fr.close()
    #os.remove(center_file)
    #os.remove(rads_file)
    log.info('Executing neper ...')
    thread = subprocess.Popen(
        ['neper', '-T', '-n', seeds.__str__(),
         '-domain', 'cube(1.0,1.0,1.0)',
         '-periodicity', 'all',
         '-morphooptiini', 'coo:file({0:s}),weight:file({1:s})'.format(center_file, rads_file),
         '-o', output_file,
         '-format', 'geo,tess,vtk'])
    thread.wait()
    log.info('Executing neper ... done')