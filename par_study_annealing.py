#!/usr/bin/env python

'''
This script runs several random generation and annealing procedure, process and report results
'''

import json
import logging
import os
import argparse
#local
import common
import foam_generate
import foam_relax
import foam_analyze

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

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
    if args.generate or args.relax_strut or args.relax_porosity or args.relax:
        for i in range(args.nrepetition):
            current_dir = zeros[len(str(i + 1)):] + str(i + 1)
            logging.info("Run %s", current_dir)
            current_output_file = os.path.join(parent_directory_path, current_dir, file_name)
            common.create_parent_directory(current_output_file)
            if args.generate:
                generate_config = config['generate']
                foam_generate.generate_3d_neper_and_fe(current_output_file, generate_config)
            if args.relax_strut or args.relax_porosity or args.relax:
                relax_config = config['relax']
                relax_config['relax-cmd'] = args.relax_cmd_file
                relax_files = foam_relax.init_file_option(current_output_file + '_')
                relax_config.update(relax_files)
                if args.relax:
                    foam_relax.relax_dry_foam(current_output_file + '.fe', relax_config)
                if args.relax_strut:
                    if not os.path.isfile(relax_config['relaxed-dry']):
                        foam_relax.relax_dry_foam(current_output_file + '.fe', relax_config)
                    foam_relax.optimize_strut_content(relax_config['relaxed-dry'], relax_config)
                if args.relax_porosity:
                    foam_relax.relax_porosity(relax_config)

    if args.analyze:
        # load results
        analyze_config = config['analyze']
        files = {}
        relax_cmd_file_name = args.relax_cmd_file.split('/')[-1][:-4]
        files['before-dry'] = common.find_files_recursively(parent_directory_path, '*before.dry.json')
        files['after-dry'] = common.find_files_recursively(parent_directory_path, '*after.dry.json')
        files['after-wet'] = common.find_files_recursively(parent_directory_path, '*.wet.json')
        for key in ['before-dry', 'after-dry', 'after-wet']:  #
            if len(files[key]) > 0:
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

def save_data_for_gnuplot(results,outpu_file_name):
    keys=['cell-size','edgesperfacet','facetspercell','angles','edge-length']
    for key in keys:
        if key in results:
            with open(outpu_file_name+"_"+key+".dat",'w') as stream:
                [hist, min_val, max_val, domain_values, range_values] = results[key]['hist']
                for i in range(len(domain_values)):
                    stream.write("{0:f}, {1:f}\n".format(domain_values[i],range_values[i]))


def collect_results(json_files, analyze_config):
    cres = {}
    for file in json_files:
        logging.info("Collecting file %s", file)
        results = foam_analyze.analyze(file, 'auto', analyze_config)
        for key in ['cell-types', 'data-angles', 'data-cells', 'length-by-edges']:
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
                        help="Json configuration file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-f", "--relax-cmd-file",
                        required=True,
                        help="Relax cmd file for Surface Evolver with relax_dry, relax_wet procedure",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-g", "--generate",
                        action='store_true',
                        help="Generate foam according to settings in config file")
    parser.add_argument("-r", "--relax",
                        action='store_true',
                        help="Relax only dry foam structure according to settings in config file")
    parser.add_argument("-p", "--relax-porosity",
                        action='store_true',
                        help="Optimize porosity according to settings in config file [output stl,vox]")
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
    args = parser.parse_args()
    main()
