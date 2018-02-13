#!/usr/bin/env python
__author__ = 'jiri1kolar'
'''
Use external program for packing generation,
Lubachevsky-Stillinger algorithm and The Force-Biased Algorithm
'''

# !/bin/env python
import struct
import json
import numpy as np
import subprocess
import logging
import os
import argparse
import common
import packing_alg3d

def execute_externalSpherePackingGenerator(conf_file, opt):
    write_config_file(conf_file, opt)
    logging.info("Executing PackingGeneration.exe ...")
    successful_generation = _execute_externalSpherePackingGenerator('-fba')
    if successful_generation:
        opt['gen-start'] = 0
        write_config_file(conf_file, opt)
        successful_generation = _execute_externalSpherePackingGenerator('-ls')
    return successful_generation

def _execute_externalSpherePackingGenerator(mode):
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
    successful_generation=execute_externalSpherePackingGenerator(conf_file, opt)
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

def execute_externalEllipsoidPackingGenerator(config_file):
    warning_count=0
    thread = subprocess.Popen(
        ['ellipsoids-openmp',
         config_file, 'dense_els.log'  # exit if an error occur
         ])#,stderr=subprocess.PIPE,stdout=subprocess.PIPE)  # stderr=subprocess.PIPE stdout=subprocess.DEVNULL
    # for line in thread.stdout:
    #     if "WARNING" in line.__str__() or "warning" in line.__str__():
    #         #print(line)
    #         warning_count+=1
    #     if warning_count>5:
    #         logging.error('Unsucessful run of external PackingGenerator, algorithm is not converging.')
    #         thread.kill()
    #         return False
    thread.wait()
    return True

def dense_ellipsoid_pack(MUX,MUYZ,sp_per_el,num_cell):
    def rot_mat(a,axis):
        s=np.sin(a)
        c=np.cos(a)
        mat=np.eye(3)
        i1=(axis+1)%3
        i2=(axis+2)%3
        mat[i1,i1]=mat[i2,i2]=c
        mat[i1,i2]=-s
        mat[i2,i1]=s
        return mat
    els=[]
    #Prepare config_file for external app ellipsoids-openmp
    config_file='eldense_config.dat'

    config_string="""	1.1  BASIC CONTROL DATA
=true if coordinates are to be saved    (         Lreslt_x)            0
	2.1  PRINCIPAL SYSTEM PARAMETERS
Number of particles                     (         No_parts)          {0}
Number of species                       (       No_species)            1
Force scaling factor                    (          Epsilon)          0.1
Rotation scaling factor                 (          Eps_rot)         0.01
Diameter increasing factor              (        Diam_incr)         0.01
Number of cells in x direction          (       No_cells_x)            1
Number of cells in y direction          (       No_cells_y)            1
Number of cells in z direction          (       No_cells_z)            1
Contraction rate of outer diameter      (             Ntau)       102400
Initial packing density                 (            Pnom0)          0.8
Maximal number of steps                 (        Max_steps)     10000000
 = true if constant volume run          (          Leq_vol)            1
Number of parts (species 0)             (           Number)          {0}
Diameter of parts (species 0)           (         Diameter)          {1}
Diameter of parts (species 0)           (         Diameter)          {2}
Diameter of parts (species 0)           (         Diameter)          {2}
	4.1  PRINTOUT CONTROL
Number of lines per page                (        Npage_len)           56
number of steps between printouts       (      Nprint_step)          100
number of steps between coord. storage  (       Nrslt_step)      1000000
number of steps between rotations       (        Nrot_step)            1""".format(num_cell,10*MUX,10*MUYZ)

    with open(config_file,'w') as stream:
        stream.write(config_string)
    execute_externalEllipsoidPackingGenerator(config_file)
    with open('ellipsoid.json') as data_file:
        els_data = json.load(data_file)
        els_vol = 0
        for el in els_data[:-1]:
            a, b, c = np.array(el['axis'])
            els_vol += 4 / 3 * np.pi * a*b*c
            rotx,roty,rotz=np.array(el['rotation'])/180.0*np.pi
            dir=rot_mat(rotz,2).dot(rot_mat(roty,1).dot(rot_mat(rotx,0).dot(np.array([1.0,0,0]))))
            dir=dir/np.linalg.norm(dir)
            el0 = packing_alg3d.Ellipsoid(a, b, c, dir, 0, 0, sp_per_el)
            el0.pos = np.array(el['position'])
            els.append(el0)
        logging.info('Final volume fraction: %.3f', els_vol)
    return els