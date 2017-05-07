#!/usr/bin/env python

import logging
import os
import common
import argparse
import json
#local
import foam_convert
import time

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger

def relax(fe_file,opt):
    relax_dry_foam(fe_file, opt)
    optimize_strut_content(opt['relaxed-dry'], opt)
    relax_porosity(opt)

def relax_dry_foam(fe_file,opt):
    relaxdryfoam_savefe(fe_file, opt['relaxed-dry'], opt['analyze-out'], opt['relax-cmd'])

def relax_porosity(opt):
    ply_file = foam_convert.stlbox2ply(opt['stl-box-out'])
    [porosity, resolution] = optimize_porosity(ply_file, opt)
    return [porosity,resolution]

def optimize_strut_content(fe_file,opt):
    logging.info('Optimizing strut content: Target %.3f ...',opt['strut-content'])
    opt_spread=newton_method_strut_cont(fe_file,opt)
    savewet(fe_file, opt['wet-out'], opt_spread)
    strut_volume=relaxwetfoam_savestlandgeo(opt['wet-out'], opt['analyze-out'], opt['stl-out'], opt['relax-cmd'])
    strut_content=strut_volume2content(strut_volume,opt)
    stl_box = foam_convert.stl2stlbox(opt['stl-out'], opt['stl-box-out'])
    logging.info('Resulting foam has strut content: %.3f with spread: %.3f', strut_content, opt_spread)
    logging.info("\tSaved to: %s", stl_box)
    return [strut_content, opt_spread]

def optimize_porosity(ply_file, opt):
    logging.info('Optimizing porosity: Target %.3f ...',opt['porosity'])
    resolution = newton_method_por(ply_file, opt['porosity'])
    porosity,opt['vtk-out'] = foam_convert.execute_binvox(ply_file, resolution)
    strut_volume = strut_content2volume(opt['strut-content'], opt)
    # correction if porosity optimization failed than it shows real strut content
    opt['porosity'] = porosity
    strut_content = strut_volume2content(strut_volume, opt)
    logging.info('Resulting foam has voxel porosity: %.3f with %d resolution.', porosity, resolution)
    logging.info("\tReal strut content: %.3f", strut_content)
    logging.info("\tResolution: %d voxels", resolution)

    return [porosity,resolution]

def relaxdryfoam_savefe(fe_dry,fe_dry_relaxed,analyze_file,relax_cmd_path):
    logging.info('Relaxing dry foam: %s',fe_dry)
    logging.info('\t used relax-cmd-file: %s', relax_cmd_path)
    tmpcmd = 'tmp_' + time.time().__str__() + '.cmd'
    with open(tmpcmd, 'w') as cmd_stream:
        cmd_stream.write('read "se_cmd/foamface.cmd"\n')
        cmd_stream.write('read "se_cmd/distribution.cmd"\n')
        cmd_stream.write('read "se_cmd/wetfoam2.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2geo.cmd"\n')
        cmd_stream.write('read "se_cmd/an_foam.cmd"\n')
        cmd_stream.write('read "' + relax_cmd_path + '"\n')
        cmd_stream.write('connected\n')
        if analyze_file is not None:
            cmd_stream.write(' an_dry_json >>> "{0:s}"\n'.format(analyze_file + '_before.dry.json'))
        cmd_stream.write('relax_dry\n')
        cmd_stream.write('dump "{0:s}"\n'.format(fe_dry_relaxed))
        if analyze_file is not None:
            cmd_stream.write(' an_dry_json >>> "{0:s}"\n'.format(analyze_file + '_after.dry.json'))
        cmd_stream.write('q\n')
        cmd_stream.write('q\n')
    success = foam_convert.execute_evolver(fe_dry, tmpcmd)
    os.remove(tmpcmd)
    if not success:
        logging.info("Fatal error: Relax procedure for dryfoam failed!")
        exit()

def savewet(fe_dry,fe_wet,spread):
    tmpcmd='tmp_'+time.time().__str__()+'.cmd'
    with open(tmpcmd,'w') as cmd_stream:
        cmd_stream.write('read "se_cmd/foamface.cmd"\n')
        cmd_stream.write('read "se_cmd/distribution.cmd"\n')
        cmd_stream.write('read "se_cmd/wetfoam2.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2geo.cmd"\n')
        cmd_stream.write('read "se_cmd/an_foam.cmd"\n')
        cmd_stream.write('connected\n')
        cmd_stream.write('spread:= {0:f}\n'.format(spread))
        cmd_stream.write('wetfoam >>> "{0:s}"\n'.format(fe_wet))
        cmd_stream.write('q\n')
        cmd_stream.write('q\n')
    success=foam_convert.execute_evolver(fe_dry,tmpcmd)
    os.remove(tmpcmd)
    if not success:
        logging.info("Fatal error: Relax procedure for dryfoam failed!")
        exit()

def relaxwetfoam_savestlandgeo(fe_wet,analyze_file,stl_out,relax_cmd_path):
    tmp_cmd = 'tmp_'+time.time().__str__()+'.cmd'
    foam_stat=fe_wet+'foam_stat.json'
    with open(tmp_cmd, 'w') as cmd_stream:
        cmd_stream.write('read "se_cmd/foamface.cmd"\n')
        cmd_stream.write('read "se_cmd/distribution.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2stl.cmd"\n')
        cmd_stream.write('read "se_cmd/fe2geo.cmd"\n')
        cmd_stream.write('read "se_cmd/an_foam.cmd"\n')
        cmd_stream.write('read "se_cmd/strut_content.cmd"\n')
        cmd_stream.write('read "'+relax_cmd_path+'"\n')
        cmd_stream.write('connected\n')
        cmd_stream.write('relax_wet\n')
        cmd_stream.write('print_stat >>> "{0:s}"\n'.format(foam_stat))
        if analyze_file is not None:
            cmd_stream.write(' an_wet_json >>> "{0:s}"\n'.format(analyze_file + '.wet.json'))
        cmd_stream.write('detorus\n')
        if stl_out is not None:
            cmd_stream.write(' stl >>> "{0:s}"\n'.format(stl_out))
        cmd_stream.write('q\n')
        cmd_stream.write('q\n')
    success=foam_convert.execute_evolver(fe_wet, tmp_cmd)
    os.remove(tmp_cmd)
    if not success:
        logging.info("Fatal error: Relax procedure for wetfoam failed!")
        exit()
    strut_volume=0.0
    with open(foam_stat,'r') as stream:
        stat=json.load(stream)
        strut_volume=stat['volume']
    os.remove(foam_stat)
    return strut_volume


def newton_method_strut_cont(fe_dry,opt):
    target_strut_volume=strut_content2volume(opt['strut-content'],opt)
    tmp_wet='wettmp_'+time.time().__str__()+'.fe'
    x=0.2
    minx=0.02
    maxx=0.6
    dx=0.01
    th=1e-5
    res0 = obj_fn_strut(fe_dry, tmp_wet, None, None, x - dx, opt) - target_strut_volume
    print("Strut content\tError\tSpread parameter")
    for i in range(10):
        strut_volume = obj_fn_strut(fe_dry, tmp_wet, None, None, x,opt)
        res1 = strut_volume - target_strut_volume
        df = (res1 - res0) / dx
        res0=res1
        newx = x - res0 / df
        dx=newx-x
        print(strut_volume2content(strut_volume, opt), strut_volume2content(res0, opt), x)
        if abs(x-newx)<th:
            os.remove(tmp_wet)
            return newx
        if newx>maxx: newx=maxx
        if newx < minx: newx = minx
        x=newx
    logging.error("Optimum was not found!")
    return x

def strut_volume2content(strut_volume,opt):
    return strut_volume/(1-opt['porosity'])

def strut_content2volume(strut_content,opt):
    return strut_content*(1-opt['porosity'])

def obj_fn_strut(fe_dry,fe_wet,analysis_file,stl,spread,opt):
    savewet(fe_dry, fe_wet, spread)
    return relaxwetfoam_savestlandgeo(fe_wet, analysis_file, stl,opt['relax-cmd'])

def newton_method_por(ply_file,target_por):
    x=250
    minx=50
    maxx=600
    dx=50
    res0 = obj_fn_por(ply_file, x - dx) - target_por
    print("Porosity\tError\tResolution")
    for i in range(10):
        por=obj_fn_por(ply_file, x)
        res1 = por - target_por
        df=(res1-res0)/dx
        res0=res1
        if abs(df) < 1e-5: return x
        newx=int(x-res0/df)
        dx=newx-x
        print(por, res0, x)
        if abs(x-newx)<2: return x
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
            'relaxed-dry': output_file+'dryrelaxed.fe',
            'analyze-out': output_file+'_an',
            'wet-out': output_file + '.fe',
            'vox-out': output_file + '.vtk'}
    return opt

def main():
    common.init_logging()
    opt=init_file_option(args.output_file)
    opt['porosity']=args.porosity
    opt['strut-content']=args.strut_content
    opt['relax-cmd']=args.relax_cmd_file
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
    parser.add_argument("-f", "--relax-cmd-file",
                        required=True,
                        help="Relax cmd file for Surface Evolver with relax_dry, relax_wet procedure",
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
