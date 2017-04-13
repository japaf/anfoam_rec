#!/usr/bin/env python

import logging
import packing_alg3d
import subprocess
import os
import argparse
import json
import common
from scipy.optimize import minimize_scalar as minimize_scalar
#local
import lib.periodicBox

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

def stlbox2ply(stl_file):
    #return output file
    logging.info('Executing meshconv ...')
    thread = subprocess.Popen(
        ['meshconv',stl_file,
         '-c','ply']) #,stdout=subprocess.DEVNULL)
    logging.info('Executing meshconv .. done')
    thread.wait()
    return stl_file[:-4]+'.ply'

def ply2vox(ply_file,resolution):
    logging.info('Executing binvox ...')
    execute_binvox(ply_file, resolution)
    logging.info('Executing binvox .. done')
    return ply_file[:-4] + '.vtk'

def stl2stlbox(stl_file,stl_box):
    logging.info("Moving to periodic box...")
    lib.periodicBox(stl_file, stl_box, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0)
    return stl_box

def stl2vox(stl_file,resolution):
    stl_box=stl_file[:-4] + 'Box.stl'
    stl_box=stl2stlbox(stl_file,stl_box)
    logging.info("Converting .stl to .ply ...")
    print(stl_box)
    ply_file=stlbox2ply(stl_box)
    logging.info("Generating voxel file ...")
    vox_file=ply2vox(ply_file,resolution)
    #logging.info("Generating ascii voxel file")
    #asciivox_file = vox_file[:4] + 'ascii.vtk'
    #lib.vtkbin2vtkascii(vox_file,asciivox_file, [0,0,0], [1.0/resolution,1.0/resolution,1.0/resolution])
    return vox_file

def geo2fe(file_name, output_file, nvolpercell):
    logging.info('Executing se_api ...')
    thread = subprocess.Popen(
        ['se_api',
         '-i',file_name,
         '-o',output_file,
        '--all-union', nvolpercell.__str__()])
    thread.wait()
    logging.info('Executing se_api ... done')
    return output_file

def dry2wet(fe_dry,fe_wet,geo_out,spread):
    tmpcmd='tmp.cmd'
    with open(tmpcmd,'w') as cmd_stream:
        cmd_stream.write('read "se_cmd/foamface.cmd"\n')
        cmd_stream.write('read "se_cmd/wetfoam2.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2geo.cmd"\n')
        cmd_stream.write('read "se_cmd/an_foam.cmd"\n')
        cmd_stream.write('read "se_cmd/relax.cmd"\n')
        cmd_stream.write('connected\n')
        if geo_out is not None:
            cmd_stream.write(' an_dry_json >>> "{0:s}"\n'.format(geo_out[:-4] + 'before.dry.json'))
        cmd_stream.write('relax\n')
        cmd_stream.write('spread:= {0:f}\n'.format(spread))
        cmd_stream.write('wetfoam >>> "{0:s}"\n'.format(fe_wet))
        if geo_out is not None:
            cmd_stream.write(' an_dry_json >>> "{0:s}"\n'.format(geo_out[:-4] + 'after.dry.json'))
        cmd_stream.write('detorus\n')
        if geo_out is not None:
            cmd_stream.write(' geo >>> "{0:s}"\n'.format(geo_out))
        cmd_stream.write('q\n')
        cmd_stream.write('q\n')
    execute_evolver(fe_dry,tmpcmd)
    #os.remove(tmpcmd)

def wet2stlgeo(fe_wet,geo_out,stl_out,iter):
    tmp_cmd = 'tmp.cmd'
    foam_stat='foam_stat.json'
    with open(tmp_cmd, 'w') as cmd_stream:
        cmd_stream.write('read "se_cmd/foamface.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2stl.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2geo.cmd"\n')
        cmd_stream.write('read "se_cmd/an_foam.cmd"\n')
        cmd_stream.write('read "se_cmd/strut_content.cmd"\n')
        cmd_stream.write('connected\n')
        cmd_stream.write('g 10;u;r;{{u;g 10}}{0:d}\n'.format(iter))
        cmd_stream.write('print_stat >>> "{0:s}"\n'.format(foam_stat))
        if geo_out is not None:
            cmd_stream.write(' an_wet_json >>> "{0:s}"\n'.format(geo_out[:-4] + '.wet.json'))
        cmd_stream.write('detorus\n')
        if geo_out is not None:
            cmd_stream.write(' geo >>> "{0:s}"\n'.format(geo_out))
        if stl_out is not None:
            cmd_stream.write(' stl >>> "{0:s}"\n'.format(stl_out))
        cmd_stream.write('q\n')
        cmd_stream.write('q\n')
    execute_evolver(fe_wet, tmp_cmd)
    strut_volume=0.0
    with open(foam_stat,'r') as stream:
        stat=json.load(stream)
        strut_volume=stat['volume']
    os.remove(tmp_cmd)
    os.remove(foam_stat)
    return strut_volume

'''
    Execution program
'''


def execute_evolver(fe_file,cmd_file):
    thread = subprocess.Popen(
        ['evolver',
         '-f', cmd_file,
         fe_file],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    thread.wait()
    warning_count=0
    error_count=0
    for line in thread.stderr:
        if "WARNING" in line.__str__():
            warning_count+=1
        if "ERROR" in line.__str__():
            error_count+=1
    if warning_count>0 or error_count>0:
        logging.error("Surface evolver reported %d warnings and %d errors.",warning_count,error_count)

def execute_binvox(ply_file,resolution):
    thread = subprocess.Popen(
        ['binvox',
         '-e', '-d', resolution.__str__(),
         '-t', 'vtk', ply_file], stdout=subprocess.PIPE)
    thread.wait()
    vtk_file=""
    for line in thread.stdout:
        if "counted" in line.__str__():
            numbers = [int(s) for s in line.__str__()[:-3].split() if s.isdigit()]
        if "VoxelFile::write_file" in line.__str__():
            vtk_file = line.__str__()[:-4].split('(')[-1]
    [foam_porosity,vtk_out]=execute_vox_fill(vtk_file,vtk_file)
    return [foam_porosity,vtk_out]

def execute_vox_fill(vtk_in,vtk_out):
    # fill holes in structure
    thread = subprocess.Popen(
        ['vox_fill',
         '-i', vtk_in,
         '-o', vtk_out,
         '--radius', '5',
         '--step','5'], stdout=subprocess.PIPE)
    thread.wait()
    foam_porosity=0.0
    for line in thread.stdout:
        if "Foam porosity" in line.__str__():
            foam_porosity = float(line.__str__()[:-4].split()[-1])
    return [foam_porosity,vtk_out]

def main():
    common.init_logging()
    if args.convert=='geo2fe':
        geo2fe(args.input_file, args.output_file, args.vol_per_cell)
    elif args.convert=='stl2vox':
        vox_file=stl2vox(args.input_file,args.resolution)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        required=True,
                        help="Input file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-o", "--output-file",
                        help="Output file",
                        metavar='FILE',
                        type=str,
                        default="foam_mod")
    parser.add_argument("-n", "--vol-per-cell",
                        help="Number of volume per cell",
                        type=int,
                        default=1)
    parser.add_argument("-r", "--resolution",
                        help="Resolution for voxel output",
                        type=int,
                        default=200)
    parser.add_argument("-c", "--convert",
                        help="Define output format",
                        choices=['geo2fe', 'stl2vox'],
                        required=True)
    args = parser.parse_args()
    main()
