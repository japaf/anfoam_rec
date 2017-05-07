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
import packing_alg3d
import logging
import os
import argparse


def execute_externalPackingGenerator(mode):
    thread = subprocess.Popen(
        ['PackingGeneration.exe',
         mode,  # exit if an error occur
         ])  # stderr=subprocess.PIPE stdout=subprocess.DEVNULL
    thread.wait()


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




def dense_sphere_pack(MUR,SR):
    volume=0
    diameters_file='diameters.txt'
    conf_file='generation.conf'
    info_file='packing.nfo'
    output_data='packing.xyzd'
    n_particles=0
    with open(diameters_file,'w') as stream:
        for i in range(1000):
            radius=abs(np.random.normal(MUR, SR))
            volume+=4.0 / 3 * (radius) ** 3 * np.pi
            stream.write("{0:f}".format(2*radius))
            n_particles+=1
            if volume>0.65:
                break
    print (volume,n_particles)

    opt = {'particles-count': n_particles,
           'packing-size': [1, 1, 1],
           'gen-start': 1,
           'seed': 341,
           'steps-to-write': 1000,
           'boundaries-mode': 1,
           'contraction-rate': 1e-3,
           'gen-mode': 1}
    write_config_file(conf_file, opt)
    execute_externalPackingGenerator('-fba')
    opt['gen-start'] = 0
    write_config_file(conf_file, opt)
    execute_externalPackingGenerator('-ls')
    scaling = load_scaling_factor(info_file)
    print('Scaling::',scaling)
    els=[]
    els_vol = 0
    with open(output_data, "rb") as f:
        particles = []
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
        logging.info('Volume fraction: %.3f Number of spheres: %d', els_vol, len(els))
    return els