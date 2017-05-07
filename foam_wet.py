#!/usr/bin/env python

import logging
import packing_alg3d
import subprocess
import os
import argparse
import json
import common
#local

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger

def main():
    common.init_logging()
    configuration = set()
    if args.config_file != "" and args.config_file is not None:
        logging.info('Loading config file: %s', args.config_file)
        with open(args.config_file) as data_file:
            configuration = json.load(data_file)
    output_file = args.output_file
    if ('work-dir' in configuration):
        output_file = configuration['work-dir'] + '/' + output_file
    if 'packing-alg' not in configuration:
        configuration['packing-alg']=args.packing_alg
    if 'grow-alg' not in configuration:
        configuration['grow-alg']=args.grow_alg

    if args.grow_alg:
        configuration['grow-alg']=args.grow_alg
    #generate_3d_neper(output_file, configuration)

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
                        choices=['sphere-random', 'ellipsoid-random'],
                        default='ellipsoid-random')
    parser.add_argument("-g", "--grow-alg",
                        help="Turn on grow algorihm",
                        action='store_true')
    args = parser.parse_args()
    main()