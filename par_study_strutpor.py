#!/usr/bin/env python

'''
This script runs several random generation and annealing procedure, process and report results
'''

import json
import logging
import os
import argparse
# local
import common
import foam_relax
import foam_analyze


__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

por_range = [0.92,0.93,0.94,0.95,0.96,0.97]
strut_range = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
por_range=[0.92]
strut_range=[0.7]

def main():
    common.init_logging()
    logging.info("Running parametric study for various porosity and strut content:")
    print('Porosity range:', por_range)
    print('Strut range:', strut_range)
    with open(args.output_file, 'a') as outstream:
        outstream.write('#Porosity strut-content resolution\n')
        for porosity in por_range:
            for strut in strut_range:
                relax_config = {"porosity": porosity, "strut-content": strut}
                relax_config['relax-cmd'] = args.relax_cmd_file
                relax_files = foam_relax.init_file_option(args.input_file + '_')
                relax_config.update(relax_files)
                foam_relax.optimize_strut_content(relax_config['relaxed-dry'], relax_config)
                parent_directory_path=common.getparrentdir(args.input_file)
                input_file=common.find_files_recursively(parent_directory_path, '*.wet.json')
                print(input_file)
                total_surface=foam_analyze.evolver_json_wet(input_file[0], None)['total-surface'][0]
                [final_porosity, final_resolution] = foam_relax.relax_porosity(relax_config)
                #[final_porosity,final_resolution]= [porosity,0] #compute only surfaces
                #if os.path.isfile(relax_config['vox-out']):
                #    os.remove(relax_config['vox-out'].__str__())
                outstream.write('{0:f}, {1:f}, {2:d}, {3:f}\n'.format(final_porosity, strut, final_resolution,total_surface))
                outstream.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        required=True,
                        help="Input relaxed .fe file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-o", "--output-file",
                        required=True,
                        help="Output results of analysis",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-f", "--relax-cmd-file",
                        required=True,
                        help="Relax cmd file for Surface Evolver with relax_dry, relax_wet procedure",
                        metavar='FILE',
                        type=str)
    args = parser.parse_args()
    main()
