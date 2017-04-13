#!/usr/bin/env python

import logging
import argparse
import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
# import pylab
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.spatial import ConvexHull

# import foam_geoextractor
import common

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"


def analyze_cell_anisotropy(cells):
    average_ratio_xy = 0.0
    average_ratio_xz = 0.0
    variance_ratio_xy = 0.0
    variance_ratio_xz = 0.0
    average_spatial_vector=np.array([0.0,0.0,0.0])
    for cell in cells:
        numpy_points = np.array(cell['vertices'])
        unwrapped_points = unwrap_points(numpy_points)
        [ratio_xy, ratio_xz] = distance_ratio_xyz(unwrapped_points)
        spatial_vector=spatial_orientation(unwrapped_points)
        average_spatial_vector+=spatial_vector
        average_ratio_xy += ratio_xy
        variance_ratio_xy += ratio_xy ** 2
        average_ratio_xz += ratio_xz
        variance_ratio_xz += ratio_xz ** 2
    average_ratio_xy /= len(cells)
    variance_ratio_xy = np.sqrt(variance_ratio_xy / len(cells) - average_ratio_xy ** 2)
    average_ratio_xz /= len(cells)
    variance_ratio_xz = np.sqrt(variance_ratio_xz / len(cells) - average_ratio_xz ** 2)
    average_spatial_vector=average_spatial_vector/len(cells)
    average_spatial_vector=average_spatial_vector/np.linalg.norm(average_spatial_vector)
    print(average_ratio_xy, variance_ratio_xy)
    print(average_ratio_xz, variance_ratio_xz)
    print(average_spatial_vector)

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
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(eq_diameters, bins)
    print(eq_diameters)
    plt.xlabel('equivalent diameter [-]')
    plt.ylabel('partial fraction')
    plt.title('Cell size distribution')
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.show()


def get_equivalent_diameter(cell_volume):
    # compute equivalent diameter of a sphere with same volume
    return (cell_volume * 6 / np.pi) ** (1 / 3)


def get_histogram(values, bins):
    min_val = 0.0
    max_val = np.nanmax(values)
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


def analyze_cell_types(cell_types):
    logging.info("Cell type analysis: ")
    _known_types = [[0, 12, 0], [2, 8, 2], [1, 10, 2],
                    [2, 8, 3], [3, 6, 4], [0, 12, 2],
                    [1, 10, 3], [2, 8, 4], [3, 6, 5],
                    [4, 4, 6], [0, 12, 3], [1, 10, 4],
                    [2, 8, 5], [0, 12, 4]]
    results = {}
    for ct0 in cell_types:
        for ct1 in _known_types:
            if abs(ct0[0] - ct1[0]) < 2 and abs(ct0[1] - ct1[1]) < 2 and abs(ct0[2] - ct1[2]) < 2:
                key = ct1.__str__()
                if (key not in results):
                    results[key] = []
                results[key].append(ct0)

    for key, value in results.items():
        print(key, len(value) / len(cell_types) * 100)


def analyze_edge_length(length_data, bins):
    logging.info("Edge length analysis: ")
    [hist, min_val, max_val, domain_values, range_values] = get_histogram(length_data, bins)
    plt.xlabel('length')
    plt.ylabel('partial fraction')
    plt.title('Edge length distribution')
    plt.grid(True)
    plt.plot(domain_values, range_values)
    plt.show()


def analyze_angle(angle_data):
    for vertex in angle_data:
        avg = (np.average(vertex))
        vrc = (np.std(vertex))
        print(avg, vrc)


def se_api_json_an(json_file, opt):
    data = json.load(open(json_file, 'r'))
    print('avg facet per cell: ', data['avg-facet-per-cell'])
    print('avg edges per facet: ', data['avg-edges-per-facet'])
    print('vertex-angle: ', data['vertex-angle'])
    # analyze_angle(data['angles-by-vertex'])
    if opt['all'] or opt['cell-types']:
        analyze_cell_types(data['cell-types'])
    if opt['all'] or opt['edge-length']:
        analyze_edge_length(data['length-by-edges'], opt['bins'])
    return data


def evolver_json_dry(json_file, opt):
    data = json.load(open(json_file, 'r'))
    # each value is counted twice, because we have list of edges per each volume unit
    if opt['all'] or opt['cell-types']:
        cell_types = []
        for item in data['cell-types']:
            # crop the array to hold only 3,4,5... sided faces
            cell_types.append(item[2:])
        analyze_cell_types(cell_types)
    if opt['all'] or opt['edge-length']:
        length_by_edges = extract_data(data['cells'], 'edge-length')
        analyze_edge_length(length_by_edges, opt['bins'])
    if opt['all'] or opt['anisotropy']:
        analyze_cell_anisotropy(data['cells'])
    analyze_equivalent_diameter(data['cells'], opt['bins'])

def evolver_json_wet(json_file, opt):
    data = json.load(open(json_file, 'r'))
    border_body = data['border-body']
    cells = data['cells'][:border_body - 1] + data['cells'][border_body:]
    if opt['all'] or opt['anisotropy']:
        analyze_cell_anisotropy(cells)  # skip border body
    analyze_equivalent_diameter(cells,opt['bins'])


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
    if format == 'geo':
        json_foam = execute_se_api_analyzator(input_file, nvolpercell)
        se_api_json_an(json_foam, opt)
    elif format == 'dryjson':
        evolver_json_dry(input_file, opt)
    elif format == 'wetjson':
        evolver_json_wet(input_file, opt)

def main():
    common.init_logging()
    opt = {'cell-types': args.cell_types,
           'edge-length': args.edge_length,
           'all': args.all,
           'anisotropy': args.anisotropy,
           'bins': args.bins}
    analyze(args.input_file,args.format,opt)

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
    parser.add_argument("-b", "--bins",
                        default="11",
                        type=int,
                        help="Number of bins for histogram.")
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
