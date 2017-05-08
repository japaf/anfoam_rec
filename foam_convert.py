#!/usr/bin/env python

import logging
import packing_alg3d
import subprocess
import os
import argparse
import json
import common
import time
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
        '--all-union', nvolpercell.__str__()],stderr=subprocess.DEVNULL,stdout=subprocess.DEVNULL)
    thread.wait()
    if not os.path.isfile(output_file+'.fe'):
        return False
    return True
'''
    Execution program
'''
def test_evolver(fe_file):
    tmp_cmd = 'tmp_' + time.time().__str__() + '.cmd'
    with open(tmp_cmd, 'w') as cmd_stream:
        cmd_stream.write('q\n')
        cmd_stream.write('q\n')
    report=execute_evolver_test(fe_file,tmp_cmd)
    os.remove(tmp_cmd)
    success=False
    if report['errors']==0 and report ['warnings']==0:
        success=True
    return success

def execute_evolver(fe_file,cmd_file):
    thread = subprocess.Popen(
        ['evolver',
         '-x', # exit if an error occur
         '-f', cmd_file,
         fe_file],stderr=subprocess.PIPE,stdout=subprocess.DEVNULL) # stdout=subprocess.DEVNULL
    thread.wait()
    warning_count=0
    error_count=0
    for line in thread.stderr:
        print(line)
        if "WARNING" in line.__str__() or "warning" in line.__str__():
            warning_count+=1
        elif "ERROR" in line.__str__() or "error" in line.__str__():
            print(line)
            error_count+=1
    if warning_count>0 or error_count>0:
        logging.error("Surface evolver reported %d warnings and %d errors.",warning_count,error_count)
    if error_count>0:
        return False
    return True

def execute_evolver_test(fe_file,cmd_file):
    thread = subprocess.Popen(
        ['evolver',
         '-x', # exit if an error occur
         '-f', cmd_file,
         fe_file],stderr=subprocess.PIPE,stdout=subprocess.DEVNULL) # stdout=subprocess.DEVNULL
    thread.wait()
    warning_count=0
    error_count=0
    for line in thread.stderr:
        print(line)
        if "WARNING" in line.__str__() or "warning" in line.__str__():
            warning_count+=1
        elif "ERROR" in line.__str__() or "error" in line.__str__():
            print(line)
            error_count+=1
    if warning_count>0 or error_count>0:
        logging.error("Surface evolver reported %d warnings and %d errors.",warning_count,error_count)
    return {'errors':error_count,'warnings':warning_count}

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
    radius=int(3*resolution//100)
    [foam_porosity,vtk_out]=execute_vox_fill(vtk_file,vtk_file,radius,radius//4)
    return [foam_porosity,vtk_out]

def execute_vox_fill(vtk_in,vtk_out,radius,step):
    # fill holes in structure
    if step<1:
        step=1
    if radius<1:
        radius=1
    thread = subprocess.Popen(
        ['vox_fill',
         '-i', vtk_in,
         '-o', vtk_out,
         '--radius', str(radius),
         '--step',str(step)], stdout=subprocess.PIPE)
    thread.wait()
    foam_porosity=0.0
    found_cells=0
    for line in thread.stdout:
        if "Found" in line.__str__():
            found_cells=int(line.__str__().split()[-2])
            print("Found cells: ",found_cells)
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
