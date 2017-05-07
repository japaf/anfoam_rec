#!/usr/bin/env python

import logging
import packing_alg3d
import packing_seq_drop
import packing_extalg3d
import subprocess
import os
import argparse
import json
import common
import numpy as np
#local
import lib.periodicBox
import foam_convert

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger

def generate_3d_neper(output_file,opt):
    num_cells=opt['num-cells']
    ellipsoids=[]
    if 'MUR' not in opt:
        MUR = opt['MUX']
    else:
        MUR = opt['MUR']
    if 'SIGMAR' not in opt:
        SR = opt['SIGMAX']
    else:
        SR = opt['SIGMAR']

    if opt['packing-alg']=='elrand':
        MUX = opt['MUX']
        MUYZ = opt['MUYZ']
        MUA = opt['MUA']
        SX = opt['SIGMAX']
        SA = opt['SIGMAA'] * np.pi / 180  # radians
        sp_per_el = opt['spheres-per-ellipsoid']
        logging.info('Generating ellipsoid packing (RSP)...')
        ellipsoids = packing_alg3d.randEllipsoidPack(MUX,MUYZ,MUA,SX,SA,sp_per_el,num_cells)
        if opt['grow-alg']:
            ellipsoids = packing_alg3d.ellipsoidGrow(ellipsoids, 1.05)

    elif opt['packing-alg']=='sprand':
        logging.info('Generating sphere packing (RSP)...')
        ellipsoids = packing_alg3d.randSpherePack(MUR,SR,num_cells)
        if opt['grow-alg']:
            ellipsoids = packing_alg3d.sphereGrow(ellipsoids, 1.05)
    #elif opt['packing-alg'] == 'elseqdrop':
    #    ellipsoids=packing_seq_drop.denseEllipsoidPack(MUX, MUY, MUZ, MUA, SX, SY, SZ, SA, sp_per_el, num_cells)
    elif opt['packing-alg'] == 'densesp':
        logging.info('Generating close sphere packing (RCP)...')
        ellipsoids=packing_extalg3d.dense_sphere_pack(MUR,SR)

    center_file=output_file+'_centers.txt'
    weights_file=output_file+'_weights.txt'
    fc = open(center_file, 'w')
    fw = open(weights_file, 'w')
    seeds = 0
    for el in ellipsoids:
        for sp in el.spheres:
            pr = sp.r#**2 #weight in Laguerre tesselation r=sqrt(w) but in neper NOT?
            px = el.pos[0] + sp.relpos[0]
            py = el.pos[1] + sp.relpos[1]
            pz = el.pos[2] + sp.relpos[2]
            fc.write('{0:f}\t{1:f}\t{2:f}\n'.format(px, py, pz))
            fw.write('{0:f}\n'.format(pr))
            seeds += 1
    fc.close()
    fw.close()
    exec_neper(output_file,seeds,center_file,weights_file)
    #os.remove(center_file)
    #os.remove(rads_file)

def generate_3d_neper_and_fe(output_file,opt):
    for i in range(10):
        generate_3d_neper(output_file, opt)
        if os.path.isfile(output_file + '.fe'):
            os.remove(output_file + '.fe')
        successful_conversion=foam_convert.geo2fe(output_file+'.geo',output_file,opt['spheres-per-ellipsoid'])
        #test evolver file
        successful_generation=False
        if successful_conversion:
            successful_generation=foam_convert.test_evolver(output_file + '.fe')
        if not successful_generation or not successful_conversion:
            logging.info("An error occured during foam generation or conversion to .fe file. Repeating generation...")
        else:
            logging.info("Foam structure was succesfully generated.")
            return True
    logging.info("ERROR: Generation of structure failed!")
    return False

def exec_neper(output_file,seeds,center_file,rads_file):
    logging.info('Executing neper ...')
    thread = subprocess.Popen(
        ['neper', '-T', '-n', seeds.__str__(),
         '-domain', 'cube(1.0,1.0,1.0)',
         '-periodicity', 'all',
         '-morphooptiini', 'coo:file({0:s}),weight:file({1:s})'.format(center_file, rads_file),
         '-o', output_file,
         '-format', 'geo'])
    thread.wait()
    logging.info('Executing neper ... done')

def main():
    common.init_logging()
    configuration = {}
    if args.config_file != "" and args.config_file is not None:
        logging.info('Loading config file: %s', args.config_file)
        with open(args.config_file) as data_file:
            configuration_all = json.load(data_file)
            configuration=configuration_all['generate']
    if args.output_file=="":
        output_file=configuration_all['output-file']
    else:
        output_file = args.output_file
    if 'packing-alg' not in configuration:
        configuration['packing-alg']=args.packing_alg
    if 'grow-alg' not in configuration:
        configuration['grow-alg']=args.grow_alg

    if args.grow_alg:
        configuration['grow-alg']=args.grow_alg
    common.create_parent_directory(output_file)
    generate_3d_neper(output_file, configuration)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config-file",
                        help="Json configuration file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-o", "--output-file",
                        help="Output geo file from Neper",
                        metavar='FILE',
                        type=str,
                        default="")
    parser.add_argument("-p", "--packing-alg",
                        help="Define packing algorithm",
                        choices=['sprand', 'elrand','elseqdrop','densesp'],
                        default='densesp')
    parser.add_argument("-g", "--grow-alg",
                        help="Turn on grow algorihm",
                        action='store_true')
    args = parser.parse_args()
    main()