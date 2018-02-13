#!/usr/bin/env python
"""
Analyze foam structure, can process json output file from evolver defining structure.
"""

import logging
import os
import argparse
import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.spatial import ConvexHull

# import foam_geoextractor
import common

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"


#Matzke experiment
#n4 n5 n6 n7
_known_types_full = [
[3,6,2,0],
[0,12,0,0],
[2,8,2,0],
[3,6,3,0],
[4,4,4,0],
[1,10,2,0],
[2,8,3,0],
[3,6,4,0],
[4,4,5,0],
[5,2,6,0],
[3,7,2,1],
[0,12,2,0],
[1,10,3,0],
[2,8,4,0],
[3,6,5,0],
[4,4,6,0],
[2,9,2,1],
[3,7,3,1],
[4,5,4,1],
[0,12,3,0],
[1,10,4,0],
[2,8,5,0],
[3,6,6,0],
[4,4,7,0],
[1,11,2,1],
[2,9,3,1],
[3,7,4,1],
[4,5,5,1],
[5,3,6,1],
[0,12,4,0],
[1,10,5,0],
[2,8,6,0],
[3,6,7,0],
[2,9,4,1],
[3,7,5,1],
[2,8,7,0]]
_known_types_values_full = [2,50,15,7,1,118,19,36,4,1,1,39,73,64,17,
                            12,4,8,1,21,35,24,7,2,2,6,5,2,2,10,4,1,
                            1,3,1,2]


_known_types = [
    [0, 12, 0], [2, 8, 2], [1, 10, 2],
    [2, 8, 3], [3, 6, 4], [0, 12, 2],
    [1, 10, 3], [2, 8, 4], [3, 6, 5],
    [4, 4, 6], [0, 12, 3], [1, 10, 4],
                    [2, 8, 5], [0, 12, 4]]


def analyze_cell_anisotropy(cells):
    average_ratio_xy = 0.0
    average_ratio_xz = 0.0
    deviation_ratio_xy = 0.0
    deviation_ratio_xz = 0.0
    average_spatial_vector=np.array([0.0,0.0,0.0])
    for cell in cells:
        numpy_points = np.array(cell['vertices'])
        unwrapped_points = unwrap_points(numpy_points)
        [ratio_xy, ratio_xz] = distance_ratio_xyz(unwrapped_points)
        spatial_vector=spatial_orientation(unwrapped_points)
        average_spatial_vector+=spatial_vector
        average_ratio_xy += ratio_xy
        deviation_ratio_xy += ratio_xy ** 2
        average_ratio_xz += ratio_xz
        deviation_ratio_xz += ratio_xz ** 2
    average_ratio_xy /= len(cells)
    deviation_ratio_xy = np.sqrt(deviation_ratio_xy / len(cells) - average_ratio_xy ** 2)
    average_ratio_xz /= len(cells)
    deviation_ratio_xz = np.sqrt(deviation_ratio_xz / len(cells) - average_ratio_xz ** 2)
    average_spatial_vector=average_spatial_vector/len(cells)
    average_spatial_vector=average_spatial_vector/np.linalg.norm(average_spatial_vector)
    avg_anis=np.array([1.0, 1.0/average_ratio_xy, 1.0/average_ratio_xz])
    avg_anis=avg_anis/np.min(avg_anis)
    anisotropy_results={
        'avg-xy':average_ratio_xy,
        'dev-xy':deviation_ratio_xy,
        'avg-xz':average_ratio_xz,
        'dev-xz':deviation_ratio_xz,
        'spatial-vector':average_spatial_vector,
        'avg-anis': avg_anis
    }
    return anisotropy_results

def plot_cell_anisotropy(anisotropy_results,opt):
    average_ratio_xy=anisotropy_results['avg-xy']
    deviation_ratio_xy=anisotropy_results['dev-xy']
    average_ratio_xz=anisotropy_results['avg-xz']
    deviation_ratio_xz=anisotropy_results['dev-xz']
    average_spatial_vector=anisotropy_results['spatial-vector']
    with open(os.path.join(opt['parent-directory'] + '_anisotropy.txt'),'w') as stream:
        stream.write("XY, {0:f}, {1:f}\n".format(average_ratio_xy, deviation_ratio_xy))
        stream.write("XZ, {0:f}, {1:f}\n".format(average_ratio_xz, deviation_ratio_xz))
        stream.write("dir, {0:f}, {1:f}, {2:f}\n".format(average_spatial_vector[0],average_spatial_vector[1],average_spatial_vector[2]))
    print("XY: ",average_ratio_xy, deviation_ratio_xy)
    print("XZ: ",average_ratio_xz, deviation_ratio_xz)
    print("dir: ",average_spatial_vector)
    print("axis: ",anisotropy_results['avg-anis'])

def unwrap_points(points):
    '''
    unwrap points from unit cube
    suppose:
        all points are inside unit cube
        maximum cell diameter is less than 0.5
    :param points: array of points [[x,y,z],...]
    :return:
    '''
    mincoor = np.array([2.0, 2.0, 2.0])
    maxcoor = np.array([-1.0, -1.0, -1.0])
    for point in points:
        for k in range(3):
            if maxcoor[k] < point[k]: maxcoor[k] = point[k]
            if mincoor[k] > point[k]: mincoor[k] = point[k]
    maxdist = np.abs(maxcoor - mincoor)
    shift = np.array([0.0, 0.0, 0.0])
    clip_threshold = 0.8
    for k in range(3):
        if maxdist[k] > clip_threshold: shift[k] = 1.0
    unwrapped_points = []
    for point in points:
        unwrapped_point = point.copy()
        for k in range(3):
            if unwrapped_point[k] < 0.5: unwrapped_point[k] += shift[k]
        unwrapped_points.append(unwrapped_point)
    return unwrapped_points

def spatial_orientation(cell_points):
    zero=1e-6
    max_distance=0.0
    max_dist_vector=np.array([0.0,0.0,0.0])
    dist_matrix=distance_matrix(cell_points,cell_points)
    for i in range(len(cell_points)):
        for j in range(i+1,len(cell_points)):
            if max_distance<dist_matrix[i,j]:
                max_distance=dist_matrix[i,j]
                max_dist_vector=cell_points[j]-cell_points[i]
    max_dist_vector=max_dist_vector/np.linalg.norm(max_dist_vector)
    abs_max_val=np.nanmax(abs(max_dist_vector))
    max_val=np.nanmax(max_dist_vector)
    if abs(abs_max_val-max_val)>zero:
        max_dist_vector=-max_dist_vector
    return max_dist_vector

def distance_ratio_xyz(cell_points):
    # evaluate ration of maximum x,y,z distance: x/z x/y size
    mincoor = np.array([1e9, 1e9, 1e9])
    maxcoor = np.array([-1e9, -1e9, -1e9])
    for point in cell_points:
        for k in range(3):
            if maxcoor[k] < point[k]: maxcoor[k] = point[k]
            if mincoor[k] > point[k]: mincoor[k] = point[k]
    distance = np.abs(maxcoor - mincoor)
    ratio_xy = distance[0] / distance[1]
    ratio_xz = distance[0] / distance[2]
    return [ratio_xy, ratio_xz]


def analyze_equivalent_diameter(cells, bins):
    eq_diameters = []
    for cell in cells:
        eq_diameters.append(get_equivalent_diameter(cell['volume']))
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(eq_diameters, bins,0,np.max(eq_diameters))
    avg=np.average(eq_diameters)
    dev=0
    for val in eq_diameters:
        dev+=val**2
    dev=np.sqrt(dev/len(eq_diameters)-avg**2)
    cell_size_results={
        'hist': [hist, min_val, max_val, domain_values, range_values],
        'avg': avg,
        'dev': dev
    }
    return cell_size_results

def plot_equivalent_diameter(cell_size_results,opt):
    [hist, min_val, max_val, domain_values, range_values]=cell_size_results['hist']
    avg=cell_size_results['avg']
    print ("d_eq = ",avg," +- ",cell_size_results['dev'])
    plt.clf()
    plt.xlabel('equivalent diameter [-]')
    plt.ylabel('partial fraction')
    plt.title('Cell size distribution')
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.plot([avg, avg], [0, np.max(range_values)], 'r')
    plt.annotate(str(avg), xy=(avg + avg / 10, np.max(range_values) / 2))
    save_figure(os.path.join(opt['parent-directory']+'_cellsize_dist.png'))

def get_equivalent_diameter(cell_volume):
    # compute equivalent diameter of a sphere with same volume
    return (cell_volume * 6 / np.pi) ** (1 / 3)


def get_histogram(values, bins,min_val,max_val):
    h = (max_val - min_val) / bins
    domain_values = np.zeros(bins + 2)
    domain_values[0] = min_val
    domain_values[bins + 1] = max_val
    for i in range(bins):
        domain_values[i + 1] = min_val + (2 * i + 1) * h / 2
    hist = np.zeros(bins)
    for val in values:
        i = (int)(np.floor((val - min_val) / h))
        if i < 0:
            i = 0
        if i >= bins:
            i = bins - 1
        hist[i] += 1
    hist = hist / len(values)
    range_values = np.zeros(bins + 2)
    range_values[1:-1] = hist
    return [hist, min_val, max_val, domain_values, range_values]


def getOrtoBase(dir):
    # return 2 ortogonal vector to dir and to each other
    # suppose that dir point mainly to axis x
    a = dir[0]
    b = dir[1]
    c = dir[2]
    y1 = 1.0
    z1 = 0.1
    z2 = 1.0
    x1 = (-b * y1 - c * z1) / a
    x2 = (b * z1 - c) / (a - b * x1)
    y2 = -z1 - x1 * x2
    ybase = np.array([x1, y1, z1])
    ybase = ybase / np.linalg.norm(ybase)
    zbase = np.array([x2, y2, z2])
    zbase = zbase / np.linalg.norm(zbase)
    return [ybase, zbase]


def execute_se_api_analyzator(file_name, nvolpercell):
    output_file = file_name[:-4] + '.json'
    logging.info('Executing se_api analyzator ...')
    thread = subprocess.Popen(
        ['se_api', '-i', file_name,
         '--all-union', nvolpercell.__str__(),
         '-a', output_file])
    thread.wait()
    logging.info('Executing se_api analyzator ... done')
    return output_file

def analyze_angles(angles_data):
    average_angle=0
    variance_angle=0
    angles_list=[]
    nangles=0
    for vertex_angles in angles_data:
        for angle in vertex_angles:
            average_angle+=angle
            variance_angle+=angle**2
            angles_list.append(angle)
            nangles+=1
    average_angle=average_angle/nangles
    variance_angle=variance_angle/nangles-average_angle**2
    mina = 0
    maxa = 180
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(angles_list, (maxa - mina) // 2, mina, maxa)
    angle_results = {
        'hist': [hist, min_val, max_val, domain_values, range_values],
        'avg': average_angle,
        'dev': np.sqrt(variance_angle)
    }
    return angle_results

def plot_angles(angle_results,opt):
    mina = 0
    maxa = 180
    [hist, min_val, max_val, domain_values, range_values]=angle_results['hist']
    average_angle=angle_results['avg']
    deviation_angle=angle_results['dev']
    plt.clf()
    plt.xlabel('angle')
    plt.ylabel('partial fraction')
    plt.title('Distribution of angles between edges')
    plt.xticks(np.arange(mina, maxa + 1, 30.0))
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.plot([average_angle, average_angle], [0, np.max(range_values)], 'r')
    plt.annotate(str(average_angle) + ' +- ' + str(deviation_angle),
                 xy=(average_angle + 5, np.max(range_values) / 2))
    save_figure(os.path.join(opt['parent-directory'] + '_angles_dist.png'))
    print('a_avg = ', average_angle, " +- ", deviation_angle)

def analyze_cell_types(cell_types,opt):
    #n3 n4 n5
    results = {}
    for ct0 in cell_types:
        for ct1 in _known_types_full:
            match=True
            for i in range(4):
                if abs(ct0[i+1] - ct1[i]) != 0:
                    match=False
                    break
            if match and np.sum(ct1)==np.sum(ct0): #cell match and not other faces
                key = ct1.__str__()
                if (key not in results):
                    results[key] = []
                results[key].append(ct0)

    bar_labels=[]
    range_values=[]
    total_percent=0
    for ct1 in _known_types_full:
        key = ct1.__str__()
        if key in results:
            value=results[key]
            percent=len(value) / len(cell_types) * 100
            total_percent+=percent
            range_values.append(percent)
        else:
            range_values.append(0)
        bar_labels.append(key.__str__())
    bar_labels.append('different')
    range_values.append(100-total_percent)
    N = len(range_values)

    ind = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars
    range_values_matzke=[]
    for val in _known_types_values_full:
        range_values_matzke.append(val/6.) #percent out of total
    range_values_matzke.append(0)
    fig, ax = plt.subplots(figsize=(20, 10))
    rects1 = ax.bar(ind, range_values, width, color='g')
    rect_matzke = ax.bar(ind+width,range_values_matzke , width, color='r')
    ax.legend((rects1[0], rect_matzke[0]), ('Model', 'Matzke'))
    # add some text for labels, title and axes ticks
    ax.set_ylabel('partial fraction [%]')
    ax.set_title('Distribution of cell types')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(bar_labels)
    for rect in rects1:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                '%.1f' % int(height),
                ha='center', va='bottom')
    for rect in rect_matzke:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                '%.1f' % int(height),
                ha='center', va='bottom')
    save_figure(os.path.join(opt['parent-directory'] + '_celltypes_dist.png'))


def facetspercell_distribution(cell_types):
    nfacets=0
    facetspercell=[]
    minf=4.5
    maxf=22.5
    for cell in cell_types:
        nfacets+=np.sum(cell)
        facetspercell.append(np.sum(cell))
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(facetspercell, int(maxf-minf),minf,maxf)
    avg=np.average(facetspercell)
    facetspercell_results={
        'hist': [hist, min_val, max_val, domain_values, range_values],
        'avg': avg,
        'minf': minf,
        'maxf': maxf
    }
    return facetspercell_results

def plot_facetspercell(facetspercell_results,opt):
    [hist, min_val, max_val, domain_values, range_values]=facetspercell_results['hist']
    avg=facetspercell_results['avg']
    minf=facetspercell_results['minf']
    maxf=facetspercell_results['maxf']
    plt.clf()
    plt.xlabel('f')
    plt.ylabel('partial fraction')
    plt.title('Distribution of cells with f faces')
    plt.xticks(np.arange(minf, maxf + 1, 1.0))
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.plot([avg, avg], [0, np.max(range_values)], 'r')
    plt.annotate(str(avg), xy=(avg + avg / 20, np.max(range_values) / 2))
    save_figure(os.path.join(opt['parent-directory'] + '_f_dist.png'))
    print('f_avg = ', avg)

def edgesperfacet_distribution(cell_types):
    nedges=0
    nfacets=0
    minn=2.5
    maxn=10.5
    edgesperfacet=[]
    for cell in cell_types:
        for i in range(len(cell)):
            nedges+=(i+3)*cell[i]
            nfacets+=cell[i]
            for j in range(cell[i]):
                edgesperfacet.append((i+3))
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(edgesperfacet, int(maxn-minn), minn,maxn)
    avg=np.average(edgesperfacet)
    edgesperfacet_results = {
        'hist': [hist, min_val, max_val, domain_values, range_values],
        'avg': avg,
        'minn': minn,
        'maxn': maxn
    }
    return edgesperfacet_results

def plot_edgesperfacet(edgesperfacet_results,opt):
    [hist, min_val, max_val, domain_values, range_values] = edgesperfacet_results['hist']
    avg = edgesperfacet_results['avg']
    minn = edgesperfacet_results['minn']
    maxn = edgesperfacet_results['maxn']
    plt.clf()
    plt.xlabel('f')
    plt.ylabel('partial fraction')
    plt.title('Distribution of faces with n edges')
    plt.xticks(np.arange(minn, maxn + 1, 1.0))
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.plot([avg, avg], [0, np.max(range_values)], 'r')
    plt.annotate(str(avg), xy=(avg + avg / 20, np.max(range_values) / 2))
    # plt.show()
    save_figure(os.path.join(opt['parent-directory'] + '_n_dist.png'))
    print ('n_avg = ', avg)

def analyze_edge_length(length_data, bins):
    mine=np.min(length_data)
    maxe=np.max(length_data)
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(length_data, bins,mine,maxe)
    avg=np.average(length_data)

    edge_results={
        'hist':[hist, min_val, max_val, domain_values, range_values],
        'avg': avg
    }
    return edge_results

def plot_edge_length(edge_results,opt):
    [hist, min_val, max_val, domain_values, range_values]=edge_results['hist']
    avg=edge_results['avg']
    print("l_avg = ", avg)
    plt.clf()
    plt.xlabel('length')
    plt.ylabel('partial fraction')
    plt.title('Edge length distribution')
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.plot([avg, avg], [0, np.max(range_values)], 'r')
    plt.annotate(str(avg), xy=(avg + avg / 20, np.max(range_values) / 2))
    save_figure(os.path.join(opt['parent-directory'] + '_edgelength_dist.png'))

def analyze_angle(angle_data):
    for vertex in angle_data:
        avg = (np.average(vertex))
        vrc = (np.std(vertex))
        print(avg, vrc)

def analyze_surface_area(facet_area_list,border_surfaces):
    total_surface=0
    for i in range(len(facet_area_list)-border_surfaces):
            total_surface+=facet_area_list[i]
    return [total_surface]

def plot_total_surface(surface_area_results,opt):
    outfile=opt['parent-directory'] + '_total_surface.dat'
    total_surface=np.average(surface_area_results)
    with open(outfile,'w') as stream:
        stream.write(total_surface.__str__())

def se_api_json_an(json_file, opt):
    data = json.load(open(json_file, 'r'))
    print('avg facet per cell: ', data['avg-facet-per-cell'])
    print('avg edges per facet: ', data['avg-edges-per-facet'])
    print('vertex-angle: ', data['vertex-angle'])
    # analyze_angle(data['angles-by-vertex'])
    results={}
    results['cell-types']=data['cell-types']
    results['edge-length']=analyze_edge_length(data['length-by-edges'], opt['bins-edge'])
    return results

def evolver_json_dry(json_file, opt):
    data = json.load(open(json_file, 'r'))
    results={}
    # each value is counted twice, because we have list of edges per each volume unit
    cell_types = []
    for item in data['cell-types']:
        # crop the array to hold only 3,4,5... sided faces
        cell_types.append(item[2:])
    length_by_edges = extract_data(data['cells'], 'edge-length')
    results['cell-types']=cell_types
    results['data-angles']=data['angles']
    results['data-cells']=data['cells']
    results['length-by-edges']=length_by_edges
    results['facetspercell'] =facetspercell_distribution(cell_types)
    results['edgesperfacet'] =edgesperfacet_distribution(cell_types)
    results['edge-length']=analyze_edge_length(length_by_edges, opt['bins-edge'])
    results['cell-anisotropy']=analyze_cell_anisotropy(data['cells'])
    results['cell-size'] =analyze_equivalent_diameter(data['cells'], opt['bins-cell'])
    results['angles'] =analyze_angles(data['angles'])
    facet_area_list=[]
    if 'facet-area' not in data:
        for cell in data['cells']:
            facet_area_list=facet_area_list+cell['facet-area']
        for i in range(len(facet_area_list)):
            facet_area_list[i]=facet_area_list[i]/2.0
        results['total-surface'] = analyze_surface_area(facet_area_list, 0)
    else:
        results['total-surface'] = analyze_surface_area(data['facet-area'],0)
    return results

def analyze_from_data(result_data,opt):
    results={}
    if 'cell-types' in result_data:
        results['cell-types']=result_data['cell-types']
        results['facetspercell'] = facetspercell_distribution(result_data['cell-types'])
        results['edgesperfacet'] = edgesperfacet_distribution(result_data['cell-types'])
    if 'length-by-edges' in result_data:
        results['edge-length'] = analyze_edge_length(result_data['length-by-edges'], opt['bins-edge'])
    if 'data-cells' in result_data:
        results['cell-anisotropy'] = analyze_cell_anisotropy(result_data['data-cells'])
        results['cell-size'] = analyze_equivalent_diameter(result_data['data-cells'], opt['bins-cell'])
    if 'data-angles' in result_data:
        results['angles'] = analyze_angles(result_data['data-angles'])
    if 'total-surface' in result_data:
        results['total-surface']=result_data['total-surface']
    return results

def evolver_json_wet(json_file, opt):
    data = json.load(open(json_file, 'r'))
    border_body = data['border-body']
    cells = data['cells'][:-1]# + data['cells'][border_body:]
    results = {}
    if opt is not None:
        results['data-cells'] = cells
        results['cell-anisotropy'] = analyze_cell_anisotropy(cells)  # skip border body
        results['cell-size'] = analyze_equivalent_diameter(cells,opt['bins-cell'])
        for i in range(len(data['facet-area'])):
            data['facet-area'][i]=data['facet-area'][i]/2.0
    results['total-surface']= analyze_surface_area(data['facet-area'],len(cells))
    return results

def plot_results(results,opt):
    plt.subplots(figsize=(20, 10))
    if 'cell-size' in results:
        plot_equivalent_diameter(results['cell-size'],opt)
    if 'edgesperfacet' in results:
        plot_edgesperfacet(results['edgesperfacet'],opt)
    if 'facetspercell' in results:
        plot_facetspercell(results['facetspercell'],opt)
    if 'angles' in results:
        plot_angles(results['angles'],opt)
    if 'edge-length' in results:
        plot_edge_length(results['edge-length'],opt)
    if 'cell-anisotropy' in results:
        plot_cell_anisotropy(results['cell-anisotropy'],opt)
    if 'cell-types' in results:
        analyze_cell_types(results['cell-types'], opt)
    if 'total-surface' in results:
        plot_total_surface(results['total-surface'],opt)

def extract_data(list, item_name):
    # exract all data from dictionary list {{'item_name':value},{..}} => [value,...]
    res = []
    for item in list:
        res += item[item_name]
    return res

def analyze(input_file,format,opt):
    if format == 'auto':
        if input_file.endswith('.geo'):
            format = 'geo'
        else:
            if input_file.endswith('.wet.json'):
                format='wetjson'
            else:
                format = 'dryjson'
    #parent_directory=os.path.join(common.getparrentdir(input_file),)
    results=None
    opt['parent-directory'] = input_file
    if format == 'geo':
        json_foam = execute_se_api_analyzator(input_file, opt['spheres-per-ellipsoid'])
        results =se_api_json_an(json_foam, opt)
    elif format == 'dryjson':
        results=evolver_json_dry(input_file, opt)
    elif format == 'wetjson':
        results=evolver_json_wet(input_file, opt)
    return results

def analyze_and_plot(input_file,format,opt):
    results=analyze(input_file, format, opt)
    if results is not None:
        opt['parent-directory'] = input_file
        plot_results(results, opt)

def save_figure(figure_name):
    plt.savefig(figure_name, bbox_inches='tight')

def main():
    common.init_logging()
    opt = {'cell-types': args.cell_types,
           'edge-length': args.edge_length,
           'all': args.all,
           'anisotropy': args.anisotropy,
           'bins-cell': args.binscell,
           'bins-edge': args.binsedge,
           'spheres-per-ellipsoid': args.nvolpercell}
    results=analyze(args.input_file,args.format,opt)
    if opt['anisotropy']:
        print('Foam anisotropy analysis:',results['cell-anisotropy']['avg-anis'])
    print('Total surface area:',results['total-surface'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",
                        required=True,
                        metavar='FILE',
                        type=str,
                        help="Input file for analysis")
    parser.add_argument("-f", "--format",
                        help="Define input format, by default file format is determined from the extension",
                        choices=['geo', 'wetjson', 'dryjson', 'jsonseapi','auto'],
                        default='auto')
    parser.add_argument("--all",
                        action='store_true',
                        help="Analyze all foam features.")
    parser.add_argument("-a", "--anisotropy",
                        action='store_true',
                        help="Analyze anisotropy of foam.")
    parser.add_argument("-c", "--cell-types",
                        action='store_true',
                        help="Analyze only cell types")
    parser.add_argument("-e", "--edge-length",
                        action='store_true',
                        help="Analyze only edge length distribution")
    parser.add_argument("--binsedge",
                        default="50",
                        type=int,
                        help="Number of bins for edge dist histogram.")
    parser.add_argument("--binscell",
                        default="12",
                        type=int,
                        help="Number of bins for cell-size dist histogram.")
    parser.add_argument("-n", "--nvolpercell",
                        default="1",
                        type=int,
                        help="Number of volumes per one cell, they will be merged.")

    args = parser.parse_args()
    main()

    '''
        #analyze resulting foam
        cellSizes = fa.computeCellSize(filename+"_mod")
        R=fa.analyzeAnizotropyXZXY(cellSizes)
        print("Average R: {0:f} +- {1:f}".format(R[0],R[1]))
        cellSizesDir=fa.directionBasedSize(filename+"_ununified_mod",nEl,2*nS+1)
        Rdir=fa.analyzeAnizotropyXZXY(cellSizesDir)
        print("Average direction based R: {0:f} +- {1:f}".format(Rdir[0],Rdir[1]))
        volumes=fa.getVolumes(filename+"_ununified_mod",nEl,2*nS+1)
        sumVol=0

        eqDiameter=[]
        x=[]
        eqDmax=0
        eqDmin=1
        for vol in volumes:
            sumVol+=vol
        for vol in volumes:
            d=vol*3/4/np.pi/sumVol
            d=math.pow(d,1.0/3)*2
            eqDiameter.append(d)
        # make histogram
        hist,min,max,h=fa.getHistogram(10,eqDiameter)
        volDist_file = os.path.join(mypath, 'eqDiameter_hist.txt')
        volDistFrac_file = os.path.join(mypath, 'eqDiameter_hist_frac.txt')
        volDist_out = open(volDist_file, 'w')
        volDistFrac_out = open(volDistFrac_file, 'w')
        i=0
        volDist_out.write("{0:f}\t{1:f}\t{2:f}\t{3:d}\n".format(min, min, min, 0))
        volDistFrac_out.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\n".format(min, min, min, 0))
        for val in hist:
            l=i*h+min
            i+=1
            r=i*h+min
            m=(r+l)/2
            x.append(m)
            volDist_out.write("{0:f}\t{1:f}\t{2:f}\t{3:d}\n".format(m,l,r,int(val)))
            volDistFrac_out.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\n".format(m, l, r, float(val)/nEl))
        volDist_out.write("{0:f}\t{1:f}\t{2:f}\t{3:d}\n".format(max, max, max, 0))
        volDistFrac_out.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\n".format(max, max, max, 0))
        hist=hist/len(eqDiameter)
        plt.plot(x,hist)
        plt.show()

    def computeCellSize(filename):
        #evaluate max size of the cell in directions of axis x,y,z
        vertices, edges, surfaces, volumes= foam_geoextractor.loadGeoFoamFile(filename)
        cellSize=[]
        for vol in volumes:
            min = np.array([1e9, 1e9, 1e9])
            max = -1*min
            for surf in vol:
                surfId=abs(surf)-1
                surfL=surfaces[surfId]
                for edge in surfL:
                    edgeId=abs(edge)-1
                    edgeL=edges[edgeId]
                    vId=edgeL[0]-1
                    v=vertices[vId]
                    for i in range(3):
                        if v[i]>max[i]:
                            max[i]=v[i]
                    for i in range(3):
                        if v[i]<min[i]:
                            min[i]=v[i]
            cellSize.append(max-min)
        return cellSize

    def analyzeAnizotropyXZXY(cellSizes):
        #evaluate cell anizotropy: everage x/z and x/y size
        Rxy=0
        Rxz=0
        dRxy=0
        dRxz=0
        Raveg=0
        if len(cellSizes)>0:
            for size in cellSizes:
                Rxy+=size[0]/size[1]
                Rxz+=size[0]/size[2]
                dRxy+=(size[0]/size[1])**2
                dRxz+=(size[0]/size[2])**2
            Rxy=Rxy/len(cellSizes)
            Rxz=Rxz/len(cellSizes)
            dRxy=math.sqrt(dRxy/len(cellSizes)-Rxy*Rxy)
            dRxz=math.sqrt(dRxz/len(cellSizes)-Rxz*Rxz)
            Raveg=(Rxy+Rxz)/2
            dRaveg=(dRxy+dRxz)/2
        return [Raveg,dRaveg,Rxy,dRxy,Rxz,dRxz]

    def directionBasedSize(filename,nEl,nS):
        #evaluate direction based cell anizotropy
        #uses vector between side subcells of each ellipsoidical cell
        #needs original geo file with ununified cells!
        #nS: total number of spheres per ellipsoid:
        vertices, edges, surfaces, volumes = foam_geoextractor.loadGeoFoamFile(filename)
        dirCell=[]
        vertByEls=[]
        for i in range(nEl):
            iL=i*nS+nS-1
            iR=i*nS+nS-2
            volL=volumes[iL]
            volR=volumes[iR]
            vertList=[]
            rightAdded=False
            #computeAllvectors from left cell to right cell
            dirVect=np.array([0.0,0.0,0.0])
            for surfL in volL:
                surfIdL=abs(surfL)-1
                edgeListL=surfaces[surfIdL]
                for edge in edgeListL:
                    edgeIdL=abs(edge)-1
                    vertListL=edges[edgeIdL]
                    vIdL=vertListL[0]-1
                    vL=vertices[vIdL]
                    vertList.append(vL)
                    #test with all vertices from the other
                    for surfR in volR:
                        surfIdR = abs(surfR) - 1
                        edgeListR = surfaces[surfIdR]
                        for edge in edgeListR:
                            edgeIdR = abs(edge) - 1
                            vertListR = edges[edgeIdR]
                            vIdR = vertListR[0] - 1
                            vR = vertices[vIdR]
                            if not rightAdded:
                                vertList.append(vR)
                            dirVect+=vR-vL
                    rightAdded=True
            vertByEls.append(vertList)
            #average direction for cell
            dirCell.append(dirVect/np.linalg.norm(dirVect))
        k=0
        l=0
        cellSize=[]
        for ii in range(nEl): #cell
            min = np.array([1e9, 1e9, 1e9])
            max = -1 * min
            for jj in range(nS):   #subcell
                vol=volumes[k]
                for surf in vol:
                    surfId=abs(surf)-1
                    surfL=surfaces[surfId]
                    for edge in surfL:
                        edgeId=abs(edge)-1
                        edgeL=edges[edgeId]
                        vId=edgeL[0]-1
                        v=vertices[vId]
                        xbase=dirCell[l]
                        [ybase,zbase]=getOrtoBase(xbase)
                        vImg=np.array([np.dot(v,xbase),np.dot(v,ybase),np.dot(v,zbase)])
                        for i in range(3):
                            if vImg[i] > max[i]:
                                max[i] = vImg[i]
                        for i in range(3):
                            if vImg[i] < min[i]:
                                min[i] = vImg[i]
                k+=1
            cellSize.append(max - min)
            l+=1
        return cellSize

    def getVolumes(filename,nEl,nS):
        # compute volume for each cell from its subcells
        # needs original geo file with ununified cells!
        # nS: total number of spheres per ellipsoid:
        vertices, edges, surfaces, volumes = foam_geoextractor.loadGeoFoamFile(filename)

        cellVol = []
        k=0
        for ii in range(nEl):  # cell
            scellVol=0
            for jj in range(nS):  # subcell
                vol = volumes[k]
                vertIdList=[]
                vertList=[]
                for surf in vol:
                    surfId = abs(surf) - 1
                    surfL = surfaces[surfId]
                    for edge in surfL:
                        edgeId = abs(edge) - 1
                        edgeL = edges[edgeId]
                        vId = edgeL[0] - 1
                        addVertex=True
                        for vId0 in vertIdList:
                            if vId0==vId:
                                addVertex=False
                        if addVertex:
                            vertIdList.append(vId)
                for vId0 in vertIdList:
                    vertList.append(vertices[vId0])
                scellVol+=ConvexHull(vertList).volume
                k += 1
            cellVol.append(scellVol)
        return cellVol
    '''
