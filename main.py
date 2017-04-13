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

Input config file should include all informations for other modules se examples for example config.json
'''
def main():
    common.init_logging()
    with open(args.config_file,'r') as stream:
        config=json.load(stream)
    output_file = config['output-file']
    file_name=output_file.split('/')[-1]
    parent_directory_path=output_file[:-len(file_name)]
    common.create_parent_directory(parent_directory_path)

    if args.generate:
        generate_config = config['generate']
        foam_generate.generate_3d_neper(output_file, generate_config)
        foam_convert.geo2fe(output_file+'.geo',output_file,generate_config['spheres-per-ellipsoid'])
    if args.relax:
        relax_config = config['relax']
        relax_files = foam_relax.init_file_option(output_file+'_relaxed')
        relax_config.update(relax_files)
        foam_relax.relax(output_file+'.fe',relax_config)
    if args.analyze:
        analyze_config = config['analyze']
        if args.recursive:
            json_files=common.find_files_recursively(parent_directory_path,file_name+'*.json')
        else:
            json_files=common.find_files(parent_directory_path,output_file+'*.json')
        for json_file in json_files:
            foam_analyze.analyze(json_file,'auto',analyze_config)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--generate",
                        action='store_true',
                        help="Generate foam according to settings in config file")
    parser.add_argument("-r", "--relax",
                        action='store_true',
                        help="Relax foam according to settings in config file")
    parser.add_argument("-a", "--analyze",
                        action='store_true',
                        help="Analyze foam according to settings in config file")
    parser.add_argument("--recursive",
                        action='store_true',
                        help="Look for json files recursive in all subfolders.")
    parser.add_argument("-c", "--config-file",
                        help="Json configuration file",
                        metavar='FILE',
                        type=str,
                        required=True)
    args = parser.parse_args()
    main()