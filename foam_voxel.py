#!/usr/bin/env python
__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

import logging
import subprocess
import os
import foam_geoextractor
import vtk
import numpy as np
import common
import argparse
import queue

def load_voxel_vtk(file_name):
    r=vtk.vtkStructuredPointsReader()
    r.SetFileName(file_name)
    r.Update()
    data=vtk.vtkStructuredPoints()
    data.SetDimensions(r.GetOutput().GetDimensions())
    data.SetExtent(r.GetOutput().GetExtent())
    data.SetScalarType(r.GetOutput().GetScalarType(), r.GetOutputInformation(0))
    data.SetNumberOfScalarComponents(r.GetOutput().GetNumberOfScalarComponents(),r.GetOutputInformation(0))
    data.GetPointData().SetScalars(r.GetOutput().GetPointData().GetScalars())
    dim=data.GetDimensions()
    array=np.zeros(dim,dtype=int)
    logging.info("Loading file ...")
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                if data.GetScalarComponentAsDouble(i,j,k,0)>0:
                    array[i, j, k] = 1
                else:
                    array[i, j, k] = 0
    return array,data

def fill_holes_vtk(array):
    dim=len(array)
    logging.info('Filling holes')
    farr = np.ones((dim, dim, dim),dtype=int)

    step=10
    nfill=0
    for i in range(0,dim,step):
        for j in range(0,dim,step):
            for k in range(0,dim,step):
                if array[i, j, k]==0:
                    isfoamcell=is_foam_cell(array,[i,j,k],step//2)
                    if isfoamcell:
                        nfill+=1
                        fill_area(farr,array,[i,j,k])
                else:
                    continue
    logging.info('Filling report: found %d cells',nfill)
    return farr

def fill_area(arr_fill0,arr_fill1,c):
    dim=len(arr_fill0)
    if (dim!=len(arr_fill1)):
        logging.error("Arrays have different dimension!")
        return False
    arr_fill1[c[0], c[1], c[2]] = 1
    arr_fill0[c[0], c[1], c[2]] = 0
    q = queue.Queue()
    q.put(c)
    while (q.qsize() > 0):
        p = q.get()
        for j in [0, 1, 2]:
            for dx in [-1, 1]:
                p1 = [p[0], p[1], p[2]]
                p1[j] = (p1[j] + dx) % dim
                if (arr_fill1[p1[0], p1[1], p1[2]] == 0 and arr_fill0[p1[0], p1[1], p1[2]] == 1):
                    arr_fill0[p1[0], p1[1], p1[2]] = 0
                    arr_fill1[p1[0], p1[1], p1[2]] = 1
                    q.put(p1)
    return True

def is_foam_cell(array,center,radius):
    dim=len(array)
    i0=center[0]-radius
    j0 = center[1] - radius
    k0 = center[2] - radius
    for i in range(-radius, 1, radius):
        for j in range(-radius, 1, radius):
            for k in range(-radius,1,radius):
                i1=(i0+i)%dim
                j1=(j0+j)%dim
                k1=(k0+k)%dim
                if array[i1,j1,k1]==1:
                    return False
    return True

def write(array,file):
    logging.info("Writing to file")
    dim=len(array)
    with open(file,'w') as stream:
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    stream.write(array[i, j, k].__str__())
                stream.write('\n')
            stream.write('\n')
            stream.write('\n')
            stream.write('\n')

def write_to_vtk(array,data,file):
    logging.info('Writing vtk output to file %s ...',file)
    dim=data.GetDimensions()
    setvox=0
    totvox=0
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                if array[i,j,k]>0:
                    data.SetScalarComponentFromDouble(i,j,k,0,128)
                    setvox += 1
                else:
                    data.SetScalarComponentFromDouble(i, j, k, 0, 0)
                totvox+=1
    print(setvox,totvox)
    logging.info("Foam porosity: %.3f",1-setvox/totvox)
    w = vtk.vtkStructuredPointsWriter()
    w.SetInputData(data)
    w.SetFileName(file)
    w.Write()

def main():
    common.init_logging()
    arr,data=load_voxel_vtk(args.input_file)
    #write(arr, data, 'foam_orig.vtk')
    farr=fill_holes_vtk(arr)
    write_to_vtk(farr,data,args.output_file)
    #    vox_file=stl2vox(args.input_file,args.resolution)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        required=True,
                        help="Input file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-o", "--output-file",
                        help="Input file",
                        metavar='FILE',
                        type=str,
                        default="foam_tmp.vtk")
    args = parser.parse_args()
    main()
