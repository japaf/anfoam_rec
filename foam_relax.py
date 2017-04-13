#!/usr/bin/env python

import logging
import subprocess
import os
import common
import argparse

#local
import foam_convert

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger

def relax(fe_file,opt):
    [strut_content, opt_spread]=optimize_strut_content(fe_file,opt)
    stl_box=foam_convert.stl2stlbox(opt['stl-out'],opt['stl-box-out'])
    ply_file=foam_convert.stlbox2ply(stl_box)
    [porosity, resolution]=optimize_porosity(ply_file,opt)
    strut_volume=strut_content2volume(strut_content,opt)
    opt['porosity']=porosity
    # correction if porosity optimization failed than it shows real strut content
    strut_content=strut_volume2content(strut_volume,opt)
    logging.info("Relaxing completed:")
    logging.info("\tPorosity: %.3f",porosity)
    logging.info("\tStrut content: %.3f",strut_content)
    logging.info("\tResolution: %d voxels",resolution)

def optimize_strut_content(fe_file,opt):
    logging.info('Optimizing strut content: Target %.3f ...',opt['strut-content'])
    opt_spread=newton_method_strut_cont(fe_file,opt)
    foam_convert.dry2wet(fe_file, opt['wet-out'], opt['geo-out'], opt_spread)
    strut_volume=foam_convert.wet2stlgeo(opt['wet-out'], opt['geo-out'], opt['stl-out'], 2)
    strut_content=strut_volume2content(strut_volume,opt)
    logging.info('Resulting foam has strut content: %.3f with spread: %.3f',strut_content,opt_spread)
    return [strut_content,opt_spread]

def optimize_porosity(ply_file, opt):
    logging.info('Optimizing porosity: Target %.3f ...',opt['porosity'])
    resolution = newton_method_por(ply_file, opt['porosity'])
    porosity,opt['vtk-out'] = foam_convert.execute_binvox(ply_file, resolution)
    logging.info('Resulting foam has voxel porosity: %.3f with %d resolution.', porosity,resolution)
    return [porosity,resolution]

def newton_method_strut_cont(fe_dry,opt):
    target_strut_volume=strut_content2volume(opt['strut-content'],opt)
    tmp_wet='wettmp.fe'
    x=0.2
    minx=0.02
    maxx=0.6
    dx=0.01
    th=1e-5
    print("Strut content\tError\tSpread parameter")
    for i in range(10):
        res1 = obj_fn_strut(fe_dry, tmp_wet, None, None, x+dx,2)-target_strut_volume
        strut_volume = obj_fn_strut(fe_dry, tmp_wet, None, None, x,2)
        res0 = strut_volume - target_strut_volume
        df = (res1 - res0) / dx
        newx = x - res0 / df
        if abs(x-newx)<th:
            os.remove(tmp_wet)
            return newx
        print(strut_volume2content(strut_volume,opt), strut_volume2content(res0,opt),x)
        if newx>maxx: newx=maxx
        if newx < minx: newx = minx
        x=newx
    logging.error("Optimum was not found!")
    return x

def strut_volume2content(strut_volume,opt):
    return strut_volume/(1-opt['porosity'])

def strut_content2volume(strut_content,opt):
    return strut_content*(1-opt['porosity'])

def obj_fn_strut(fe_dry,fe_wet,geo,stl,strut_volume,iter):
    foam_convert.dry2wet(fe_dry, fe_wet,geo, strut_volume)
    return foam_convert.wet2stlgeo(fe_wet, geo, stl,iter)

def newton_method_por(ply_file,target_por):
    x=150
    minx=50
    maxx=500
    dx=50
    print("Porosity\tError\tResolution")
    for i in range(10):
        res1 = obj_fn_por(ply_file,x+dx)-target_por
        por=obj_fn_por(ply_file, x)
        res0 = por - target_por
        df=(res1-res0)/dx
        newx=int(x-res0/df)
        if abs(x-newx)<2: return x
        print(por, res0, x)
        if newx>maxx: newx=maxx
        if newx < minx: newx = minx
        x=newx
    logging.error("Optimum was not found!")
    return x

def obj_fn_por(ply_file,resolution):
    porosity,vtk_file= foam_convert.execute_binvox(ply_file, resolution)
    os.remove(vtk_file)
    return porosity

def init_file_option(output_file):
    opt = {'stl-out': output_file + '.stl',
           'stl-box-out': output_file + 'Box.stl',
           'ply-out': output_file + '.ply',
           'geo-out': output_file + '.geo',
           'wet-out': output_file + '.fe',
           'vox-out': output_file + '.vtk'}
    return opt

def main():
    common.init_logging()
    opt=init_file_option(args.output_file)
    opt['porosity']=args.porosity
    opt['strut-content']=args.strut_content
    relax(args.input_file,opt)
    #elif args.convert=='evolver':
    #    vox_file=stl2vox(args.input_file,args.resolution)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        required=True,
                        help="Input Surface evolver *.fe file",
                        metavar='FILE',
                        type=str)
    parser.add_argument("-o", "--output-file",
                        help="Output file - without extension",
                        metavar='FILE',
                        type=str,
                        default='results/res_foam')
    parser.add_argument("-p", "--porosity",
                        help="Foam porosity",
                        type=float,
                        default=0.95)
    parser.add_argument("-s", "--strut-content",
                        help="Volume of struts in foam",
                        type=float,
                        default=0.6)
    args = parser.parse_args()
    main()
