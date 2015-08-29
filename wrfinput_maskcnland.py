#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# mask region out of china as sea
# Author: Hui Zheng

import os.path
import numpy as np
import matplotlib.path as mpath
import netCDF4 as nc

def inchina(lat, lon):
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    POLYGONDIR_CANDIDATE = [os.path.expandvars('$HOME/soft/cnpolygon'),
                            '/Volumes/Works/soft/cnpolygon']
    for POLYGONDIR in POLYGONDIR_CANDIDATE:
        if os.path.isdir(POLYGONDIR):
            break
    path_ml = mpath.Path(np.loadtxt(os.path.join(POLYGONDIR, 'china_mainland.txt')))
    path_hn = mpath.Path(np.loadtxt(os.path.join(POLYGONDIR, 'china_hainan.txt')))
    path_tw = mpath.Path(np.loadtxt(os.path.join(POLYGONDIR, 'china_taiwan.txt')))
    locs = np.column_stack((lon.ravel(), lat.ravel()))
    in_china = np.logical_or(np.logical_or(path_ml.contains_points(locs),
                                           path_hn.contains_points(locs)),
                             path_tw.contains_points(locs)).reshape(lat.shape)
    return in_china

def main(wrfinput):
    with nc.Dataset(wrfinput, 'r+') as ncf:
        latxy = ncf.variables['XLAT'][0,:,:]
        lonxy = ncf.variables['XLONG'][0,:,:]
        inchn = inchina(latxy, lonxy)
        land = ncf.variables['LANDMASK'][0,:,:].astype('b')
        landinchn = np.logical_and(land, inchn)
        landinchn = landinchn.astype('f')
        ncf.variables['LANDMASK'][0,:,:] = landinchn[:]
    return

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mask region out of China as sea.')
    parser.add_argument('wrfinput', type=str,
                        help='wrfinput file')
    args = parser.parse_args()
    main(args.wrfinput)
