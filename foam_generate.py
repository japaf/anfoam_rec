#!/usr/bin/env python

import logging
import packing_alg3d
import packing_seq_drop
import subprocess
import os
import argparse
import json
import common
import numpy as np
#local
import lib.periodicBox

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger

def generate_3d_neper(output_file,opt):
    MUX=opt['MUX']
    MUY=opt['MUY']
    MUZ=opt['MUZ']
    MUA=opt['MUA']
    SX=opt['SIGMAX']
    SY=opt['SIGMAY']
    SZ=opt['SIGMAZ']
    SA=opt['SIGMAA']*np.pi/180 #radians
    sp_per_el=opt['spheres-per-ellipsoid']
    num_cells=opt['num-cells']
    ellipsoids=[]

    if opt['packing-alg']=='elrand':
        logging.info('Generating ellipsoid packing (adsorption)...')
        ellipsoids = packing_alg3d.randEllipsoidPack(MUX,MUY,MUZ,MUA,SX,SY,SZ,SA,sp_per_el,num_cells)
        if opt['grow-alg']:
            ellipsoids = packing_alg3d.ellipsoidGrow(ellipsoids, 1.05)

    elif opt['packing-alg']=='sprand':
        if 'MUR' not in opt:
            MUR=opt['MUX']
        else:
            MUR=opt['MUR']
        if 'SIGMAR' not in opt:
            SR=opt['SIGMAX']
        else:
            SR = opt['SIGMAR']
        logging.info('Generating sphere packing (adsorption)...')
        ellipsoids = packing_alg3d.randEllipsoidPack(MUR,MUR,MUR,0,SR,SR,SR,0,1,num_cells)
        if opt['grow-alg']:
            ellipsoids = packing_alg3d.ellipsoidGrow(ellipsoids, 1.05)
    elif opt['packing-alg'] == 'elseqdrop':
        ellipsoids=packing_seq_drop.denseEllipsoidPack(MUX, MUY, MUZ, MUA, SX, SY, SZ, SA, sp_per_el, num_cells)

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
    exec_neper(output_file,seeds,center_file,rads_file)
    os.remove(center_file)
    os.remove(rads_file)

def exec_neper(output_file,seeds,center_file,rads_file):
    logging.info('Executing neper ...')
    thread = subprocess.Popen(
        ['neper', '-T', '-n', seeds.__str__(),
         '-domain', 'cube(1.0,1.0,1.0)',
         '-periodicity', 'all',
         '-morphooptiini', 'coo:file({0:s}),weight:file({1:s})'.format(center_file, rads_file),
         '-o', output_file,
         '-format', 'geo'])
    logging.info('Executing neper .. DONE')
    thread.wait()

    logging.info('Executing neper ... done')

def main():
    common.init_logging()
    configuration = set()
    if args.config_file != "" and args.config_file is not None:
        logging.info('Loading config file: %s', args.config_file)
        with open(args.config_file) as data_file:
            configuration = json.load(data_file)
    output_file = args.output_file
    if 'packing-alg' not in configuration:
        configuration['packing-alg']=args.packing_alg
    if 'grow-alg' not in configuration:
        configuration['grow-alg']=args.grow_alg

    if args.grow_alg:
        configuration['grow-alg']=args.grow_alg
    generate_3d_neper(output_file, configuration)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config-file",
                        help="Json configuration file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-i", "--input-file",
                        help="Input geo file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-o", "--output-file",
                        help="Output geo file from Neper",
                        metavar='FILE',
                        type=str,
                        default="foam.geo")
    parser.add_argument("-p", "--packing-alg",
                        help="Define packing algorithm",
                        choices=['sprand', 'elrand','elseqdrop'],
                        default='elrand')
    parser.add_argument("-g", "--grow-alg",
                        help="Turn on grow algorihm",
                        action='store_true')
    args = parser.parse_args()
    main()