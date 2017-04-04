#!/usr/bin/env python

import json
import logging
import os
import numpy as np
import common
import foam_generate
import foam_model
import argparse
import distutils.spawn as ds

__name__ = 'anfoam_rec'
__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

#logger setting
#set logger
common.init_logging()
log = logging.getLogger(__name__)

def check_exec(exec):
    return ds.find_executable(exec)

def main():
    exec={'neper':'neper',
          'evolver': 'evolver',
          'se_api': 'se_api'}
    neper_path=check_exec(exec['neper'])
    se_api_path=check_exec(exec['evolver'])
    evolver_path=check_exec(exec['se_api'])
    if neper_path is None:
        log.error("neper not installed")
        return 0
    if se_api_path is None:
        log.error("se_api not installed")
        return 0
    if evolver_path is None:
        log.error("evolver not installed")
        return 0
    configuration=set()
    if args.config_file != "" and args.config_file is not None:
        log.info('Loading config file: %s', args.config_file)
        with open(args.config_file) as data_file:
            configuration = json.load(data_file)
        #if ('work-dir' in configuration):
        #    log.info("Setting working dir to %s",configuration['work-dir'])
        #    os.chdir(configuration['work-dir'])
    if args.module=='generate':
        gen_options={'packing-alg':args.packing_alg,
                     'grow-alg':args.grow_alg,
                     'MUX': 0.08,
                     'MUY': 0.06,
                     'MUZ': 0.06,
                     'MUA': 0.0,
                     'SIGMAX': 0.003,
                     'SIGMAY': 0.002,
                     'SIGMAZ': 0.002,
                     'SIGMAA': 10,
                     'num-cell':27,
                     'spheres-per-ellipsoid': 4,
                     'MUR': 0.08,
                     'SIGMAR' : 0.003
                     }
        if 'MUX' in configuration:
            gen_options['MUX']=configuration['MUX']
            gen_options['MUR'] = configuration['MUX']
        if 'MUY' in configuration:
            gen_options['MUY']=configuration['MUY']
        if 'MUZ' in configuration:
            gen_options['MUZ']=configuration['MUZ']
        if 'MUR' in configuration:
            gen_options['MUR']=configuration['MUR']
        if 'SIGMAX' in configuration:
            gen_options['SIGMAX']=configuration['SIGMAX']
            gen_options['SIGMAR'] = configuration['SIGMAX']
        if 'SIGMAY' in configuration:
            gen_options['SIGMAY']=configuration['SIGMAY']
        if 'SIGMAZ' in configuration:
            gen_options['SIGMAZ']=configuration['SIGMAZ']
        if 'SIGMAA' in configuration:
            gen_options['SIGMAA']=configuration['SIGMAA']
        if 'SIGMAR' in configuration:
            gen_options['SIGMAR']=configuration['SIGMAR']
        foam_generate.generate_3d(args.output_file,gen_options)

        # ps=ellipsoidPacking.randSpherePack(0.15,0.01,1,1,1,30,True)
        # sps1=ellipsoidPacking.sphereGrow(sps,1.05)
        # els = ellipsoidPacking.randEllipsoidPack(MUX, MUY, MUZ, MUA, SIGMAX, SIGMAY, SIGMAZ, SIGMAA, 1, 1, 1, nS, 10, True)
        # els2= ellipsoidPacking.ellipsoidGrow(els,1.05)
        # nEl = NeperAnFoam.laguerre3dperiodic(MUX, MUY, MUZ, MUA, SIGMAX, SIGMAY, SIGMAZ, SIGMAA, nS, 0, True, filename)


    #locals().update(data)  # Creates variables from dictionary
    # change dir to output
    #os.chdir(mypath)

    #memory omezeni?
    #
    if args.module=='model':
        if args.input_file is None or args.input_file=="":
            log.error("Specify valid geo input-file.")
            return 0
        model = foam_model.foam_model()
        model.initFromGeoFile(args.input_file)
        #Modeling
        if args.object_method=='implicit':
            model.createPrimitivesImplicit(args.limit_num_obj)
            model.sampleModelImplicit()
            if (args.save_vox is not None):
                model.sampleVoxels(args.voxel_resolution)
        elif args.object_method=='source':
            model.createPrimitivesSource(args.limit_num_obj)
            model.sampleModelSource()
            if (args.save_vox is not None):
                model.sampleVoxelsSource(args.voxel_resolution)

        #Output
        if args.save_vox is not None:
                model.saveAsVox(args.save_vox)
        if args.save_stl is not None:
                model.saveAsSTL(args.save_stl)
        if args.display:
            model.renderModel()

if __name__ == 'anfoam_rec':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--module",
                        required=True,
                        help="Module:\n\tgenerate - generate foam with neper and packing alg.\n\trelax - relax foam in evolver.\n\tmodel - generate 3D foam model",
                        choices=['generate', 'relax', 'model'],
                        default="generate")
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
                        choices=['sphere-random', 'ellipsoid-random'],
                        default='sphere-random')
    parser.add_argument("-g", "--grow-alg",
                        help="Turn on grow algorihm",
                        action='store_true')
    parser.add_argument("-d", "--dimension",
                        help="Dimension of structure",
                        type=int,
                        default=3)
    parser.add_argument("--object-method",
                        help="Object method for structure modelling.",
                        choices=['implicit', 'source'],
                        default="implicit")
    parser.add_argument("-r", "--voxel-resolution",
                        help="Resolution for voxel output",
                        type=int,
                        default=100)
    parser.add_argument("--limit-num-obj",
                        help="Limit maximal number of rendered objects of each type: vertex,edge,surface.",
                        type=int,
                        default=-1)
    parser.add_argument("-disp", "--display",
                        help="Display results",
                        action='store_true')
    parser.add_argument("--save-vox",
                        help="Save voxel vtk format",
                        metavar='FILE',
                        type=str)
    parser.add_argument("--save-stl",
                        help="Save in stl format",
                        metavar='FILE',
                        type=str)
    args = parser.parse_args()
    main()