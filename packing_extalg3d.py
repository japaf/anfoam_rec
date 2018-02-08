#!/usr/bin/env python
__author__ = 'jiri1kolar'
'''
Use external program for packing generation,
Lubachevsky-Stillinger algorithm and The Force-Biased Algorithm
'''

# !/bin/env python
import struct
import numpy as np
import subprocess
import logging
import os
import argparse
import common
import packing_alg3d

def execute_externalPackingGenerator(conf_file,opt):
    write_config_file(conf_file, opt)
    logging.info("Executing PackingGeneration.exe ...")
    successful_generation = _execute_externalPackingGenerator('-fba')
    if successful_generation:
        opt['gen-start'] = 0
        write_config_file(conf_file, opt)
        successful_generation = _execute_externalPackingGenerator('-ls')
    return successful_generation

def _execute_externalPackingGenerator(mode):
    warning_count=0
    thread = subprocess.Popen(
        ['PackingGeneration.exe',
         mode,  # exit if an error occur
         ],stderr=subprocess.PIPE,stdout=subprocess.PIPE)  # stderr=subprocess.PIPE stdout=subprocess.DEVNULL
    for line in thread.stdout:
        if "WARNING" in line.__str__() or "warning" in line.__str__():
            #print(line)
            warning_count+=1
        if warning_count>5:
            logging.error('Unsucessful run of external PackingGenerator, algorithm is not converging.')
            thread.kill()
            return False
    thread.wait()
    return True


def write_config_file(file_name, opt):
    with open(file_name, 'w') as os:
        os.write('Particles count: {0:d}\n'.format(opt['particles-count']))
        ps = opt['packing-size']
        os.write('Packing size: {0:d} {1:d} {2:d}\n'.format(ps[0], ps[1], ps[2]))
        os.write('Generation start: {0:d}\n'.format(opt['gen-start']))
        os.write('Seed: {0:d}\n'.format(opt['seed']))
        os.write('Steps to write: {0:d}\n'.format(opt['steps-to-write']))
        os.write('Boundaries mode: {0:d}\n'.format(opt['boundaries-mode']))
        os.write('Contraction rate: {0:f}\n'.format(opt['contraction-rate']))
        os.write('Generation mode: {0:d}\n'.format(opt['gen-mode']))


def load_scaling_factor(file_name):
    teor_por = 0
    final_por = 0
    with open(file_name, 'r') as ins:
        for line in ins:
            arr = line.strip().split()
            if arr[0] == 'Theoretical':
                teor_por = float(arr[2])
            elif arr[0] == 'Final':
                final_por = float(arr[2])
    scaling_factor = ((1 - final_por) / (1 - teor_por)) ** (1.0 / 3)
    return scaling_factor




def dense_sphere_pack(MUR,SR,n_particles):

    diameters_file='diameters.txt'
    conf_file='generation.conf'
    info_file='packing.nfo'
    output_data='packing.xyzd'
    #clean files
    for filef in [diameters_file,conf_file,info_file,output_data]:
        if os.path.isfile(filef):
            os.remove(filef)
    if n_particles==-1:
        n_particles=27
    elif n_particles>74 and n_particles<400:
        logging.warning('Number of cells may cause problems during generation, 3rd party software have some bugs in range of 75-400.')
    volume = 0
    radius_arr=[]
    for i in range(n_particles):
        radius=abs(np.random.normal(MUR, SR))
        radius_arr.append(radius)
        volume += 4.0 / 3 * (radius) ** 3 * np.pi
    scaling=0.5/volume #make sure that radius wont be too big

    with open(diameters_file,'w') as stream:
        for i in range(n_particles):
            stream.write("{0:f}\n".format(2*radius_arr[i]*np.power(scaling,1.0/3)))

    logging.info ('Number of spheres: %d',n_particles)

    opt = {'particles-count': n_particles,
           'packing-size': [1, 1, 1],
           'gen-start': 1,
           'seed': 341,
           'steps-to-write': 1000,
           'boundaries-mode': 1,
           'contraction-rate': 1e-3,
           'gen-mode': 1}
    successful_generation=execute_externalPackingGenerator(conf_file,opt)
    els = None
    if successful_generation:
        scaling = load_scaling_factor(info_file)
        els_vol = 0
        els = []
        with open(output_data, "rb") as f:
            byte = f.read(8)
            while byte != b"":
                pos = np.array([0.0, 0.0, 0.0, 0.0])
                for i in range(4):
                    num = struct.unpack('<d', byte)
                    pos[i] = num[0]
                    byte = f.read(8)
                pos[3] = scaling * pos[3]
                els_vol += 4.0 / 3 * (pos[3] / 2) ** 3 * np.pi
                radius=pos[3] / 2
                el=packing_alg3d.Ellipsoid(radius,radius,radius,0,0,0,1)
                el.pos=pos[:-1]
                els.append(el)
            logging.info('Final volume fraction: %.3f', els_vol)
    return els