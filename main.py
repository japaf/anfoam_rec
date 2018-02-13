#!/usr/bin/env python

import json
import logging
import os
import numpy as np
import common
import foam_generate
import foam_relax
import foam_analyze
import foam_convert

import argparse


__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

'''
Main file
can handle
    foam generation
    foam relaxation
    foam analysis

Input config file should include all informations for other modules see examples/config.json
'''
def main():
    common.init_logging()
    with open(args.config_file,'r') as stream:
        config=json.load(stream)
    output_file = config['output-file']
    file_name=output_file.split('/')[-1]
    parent_directory_path=common.getparrentdir(output_file)
    common.create_parent_directory(parent_directory_path)



    if args.generate:
        generate_config = config['generate']
        successful_generation=foam_generate.generate_3d_neper_and_fe(output_file, generate_config)
        if not successful_generation:
            return False

    if args.relax_dry or args.relax_strut or args.relax_por or args.relax_all:
        if args.anisotropic_relax is not None:
            output_file_anis=output_file+'_anis'
            ax,by,cz=args.anisotropic_relax
            foam_relax.prepare_anisotropic_relax_fe(ax,by,cz,output_file+'.fe',output_file_anis+'.fe')
            output_file=output_file_anis
        relax_config = config['relax']
        relax_files = foam_relax.init_file_option(output_file+'_relaxed')
        relax_config.update(relax_files)
        if args.diffusion:
            relax_config['convert-vtk-to-txt'] = True
        if args.relax_all:
            foam_relax.relax_all(output_file+'.fe',relax_config)
        else:
            if args.relax_dry:
                foam_relax.relax_dry_foam(output_file + '.fe', relax_config)
            if args.relax_strut:
                foam_relax.optimize_strut_content(relax_config['relaxed-dry'], relax_config)
            if args.relax_por:
                foam_relax.relax_porosity(relax_config)


    if args.analyze:
        analyze_config = config['analyze']
        if args.recursive:
            json_files=common.find_files_recursively(parent_directory_path,file_name+'*.json')
        else:
            json_files=common.find_files(parent_directory_path,output_file+'*.json')
        for json_file in json_files:
            logging.info('Analyzing file:')
            print(json_file)
            foam_analyze.analyze_and_plot(json_file,'auto',analyze_config)
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config-file",
                        help="Json configuration file",
                        metavar='FILE',
                        type=str,
                        required=True)
    parser.add_argument("-g", "--generate",
                        action='store_true',
                        help="Generate foam according to settings in config file [output: geo,fe]")
    parser.add_argument("-d", "--relax-dry",
                        action='store_true',
                        help="Relax dry-foam, according to settings in config file [output: fe]")
    parser.add_argument("-s", "--relax-strut",
                        action='store_true',
                        help="Create wet-foam structure, optimize strut content, according to settings in config file [output: fe,stl]")
    parser.add_argument("-p", "--relax-por",
                        action='store_true',
                        help="Optimize foam porosity, according to settings in config file [output: vtk(vox)]")
    parser.add_argument("-r", "--relax-all",
                        action='store_true',
                        help="Relax dry foam, prepare wet foam, optimize strut content and porosity according to settings in config file [output fe,stl,vtk(vox)]")
    parser.add_argument("--diffusion",
                        action='store_true',
                        help="Create suitable file for diffusion simulation after porosity optimization[output txt]")
    parser.add_argument("-a", "--analyze",
                        action='store_true',
                        help="Analyze foam according to settings in config file")
    parser.add_argument("--recursive",
                        action='store_true',
                        help="Look for json files recursive in all subfolders.")
    parser.add_argument("--anisotropic-relax",
                        nargs='*',
                        help="Relax foam with set anisotropy ratio, supply three, space separated values: 1.5 1 1")
    args = parser.parse_args()
    main()