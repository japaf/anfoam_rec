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

"""
Run parametric study, compare different methods for creation of anisotropic foam structure 
"""

def main():
    common.init_logging()
    with open(args.config_file, 'r') as stream:
        config = json.load(stream)
    output_file = config['output-file']
    file_name = output_file.split('/')[-1]
    parent_directory_path = common.getparrentdir(output_file)
    common.create_parent_directory(parent_directory_path)
    zeros = "000"
    logging.info("Testing: %s", parent_directory_path)
    if args.generate or args.relax_strut or args.relax_porosity or args.relax_dry or args.diffusion:
        for i in range(args.nrepetition):
            current_dir = zeros[len(str(i + 1)):] + str(i + 1)
            logging.info("Run %s", current_dir)
            current_output_file = os.path.join(parent_directory_path, current_dir, file_name)
            common.create_parent_directory(current_output_file)
            if args.generate:
                generate_config = config['generate']
                foam_generate.generate_3d_neper_and_fe(current_output_file, generate_config)
            if args.relax_strut or args.relax_porosity or args.relax_dry or args.diffusion:
                if args.anisotropic_relax is not None:
                    output_file_anis = current_output_file + '_anis'
                    ax, by, cz = args.anisotropic_relax
                    foam_relax.prepare_anisotropic_relax_fe(ax, by, cz, current_output_file + '.fe', output_file_anis + '.fe')
                    current_output_file = output_file_anis
                relax_config = config['relax']
                if args.relax_cmd_file!="":
                    relax_config['relax-cmd'] = args.relax_cmd_file
                if args.scale_vector_before is not None:
                    relax_config['scale-vector-before']=args.scale_vector_before
                if args.scale_vector_after is not None:
                    relax_config['scale-vector-after']=args.scale_vector_after
                relax_files = foam_relax.init_file_option(current_output_file + '_relaxed')
                relax_config.update(relax_files)
                if args.relax_dry:
                    foam_relax.relax_dry_foam(current_output_file + '.fe', relax_config)
                if args.relax_strut:
                    if not os.path.isfile(relax_config['relaxed-dry']):
                        foam_relax.relax_dry_foam(current_output_file + '.fe', relax_config)
                    foam_relax.optimize_strut_content(relax_config['relaxed-dry'], relax_config)
                if args.diffusion:
                    relax_config['convert-vtk-to-txt']=True
                if args.const_res:
                    relax_config['const-res'] = True
                    relax_config['resolution'] = 362
                if args.relax_porosity:
                    foam_relax.relax_porosity(relax_config)


    if args.analyze:
        # load results
        analyze_config = config['analyze']

        files = {}
        relax_cmd_file_name = args.relax_cmd_file.split('/')[-1][:-4]
        files['before-dry'] = common.find_files_recursively(parent_directory_path, '*_before.dry.json')
        files['anis-before-dry'] = common.find_files_recursively(parent_directory_path, '*_anis_relaxed_an_before.dry.json')
        files['after-dry'] = common.find_files_recursively(parent_directory_path, '*_after.dry.json')
        files['anis-after-dry'] = common.find_files_recursively(parent_directory_path, '*_anis_relaxed_an_after.dry.json')
        #filter _anis_relaxed* files from  after-dry and before-dry
        files['before-dry']=[file for file in files['before-dry'] if file not in files['anis-before-dry']]
        files['after-dry']=[file for file in files['after-dry'] if file not in files['anis-after-dry']]
        files['after-wet'] = common.find_files_recursively(parent_directory_path, '*.wet.json')
        for key in ['before-dry','anis-before-dry', 'after-dry','anis-after-dry', 'after-wet']:  #
            if len(files[key]) > 0:
                analyze_config['n-repetition'] = len(files[key])
                logging.info('Analyzis of %s json files...', key)
                plot_results(files[key], analyze_config, parent_directory_path + relax_cmd_file_name + key)
            else:
                logging.info('No %s files for analysis in %s... skipped', key, parent_directory_path)


def plot_results(json_files, analyze_config, output_file_name_generic):
    collected_resutls = collect_results(json_files, analyze_config)
    resultsToPlot = foam_analyze.analyze_from_data(collected_resutls, analyze_config)
    analyze_config['parent-directory'] = output_file_name_generic
    foam_analyze.plot_results(resultsToPlot, analyze_config)
    save_data_for_gnuplot(resultsToPlot,output_file_name_generic)

def save_data_for_gnuplot(results,output_file_name):
    keys=['cell-size','edgesperfacet','facetspercell','angles','edge-length']
    for key in keys:
        if key in results:
            with open(output_file_name+"_"+key+".dat",'w') as stream:
                [hist, min_val, max_val, domain_values, range_values] = results[key]['hist']
                for i in range(len(domain_values)):
                    stream.write("{0:f}, {1:f}\n".format(domain_values[i],range_values[i]))


def collect_results(json_files, analyze_config):
    cres = {}
    for file in json_files:
        logging.info("Collecting file %s", file)
        results = foam_analyze.analyze(file, 'auto', analyze_config)
        for key in ['cell-types', 'data-angles', 'data-cells', 'length-by-edges','total-surface']:
            add_arr_todict(key, results, cres)
    return cres


def add_arr_todict(key, dict0, dict1):
    if key in dict0:
        if key not in dict1:
            dict1[key] = []
        dict1[key] = dict1[key] + dict0[key]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config-file",
                        required=True,
                        help="Json configuration file for parametric study",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-f", "--relax-cmd-file",
                        help="Relax cmd file for Surface Evolver with relax_dry, relax_wet procedure",
                        metavar='FILE',
                        default="",
                        type=str)
    parser.add_argument("-g", "--generate",
                        action='store_true',
                        help="Generate foam according to settings in config file")
    parser.add_argument("-d", "--relax-dry",
                        action='store_true',
                        help="Relax only dry foam structure according to settings in config file")
    parser.add_argument("-p", "--relax-porosity",
                        action='store_true',
                        help="Optimize porosity according to settings in config file [output stl,vox]")
    parser.add_argument("--const-res",
                        action='store_true',
                        help="With constant resolution [output stl,vox]")
    parser.add_argument("--diffusion",
                        action='store_true',
                        help="Create suitable txt file for diffusion simulation")
    parser.add_argument("-s", "--relax-strut",
                        action='store_true',
                        help="Relax foam and optimize strut concent, according to settings in config file [output: stl]")
    parser.add_argument("-a", "--analyze",
                        action='store_true',
                        help="Analyze foam according to settings in config file")
    parser.add_argument("-n", "--nrepetition",
                        help="Number of repetition",
                        type=int,
                        default=10)
    parser.add_argument("--anisotropic-relax",
                        nargs='*',
                        help="Relax foam with set anisotropy ratio, supply three, space separated values: 1.5 1 1")
    parser.add_argument("--scale-vector-before",
                        nargs='*',
                        help="Linear scale of axis before relaxation, supply three, space separated values: 1.5 1 1")
    parser.add_argument("--scale-vector-after",
                        nargs='*',
                        help="Linear scale of axis after relaxation, supply three, space separated values: 1.5 1 1")
    args = parser.parse_args()
    main()
