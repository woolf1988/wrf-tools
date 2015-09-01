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
                            os.path.expandvars('$WORK/soft/cnpolygon')]
    for POLYGONDIR in POLYGONDIR_CANDIDATE:
        if os.path.isdir(POLYGONDIR):
            break
    path_ml = mpath.Path(np.loadtxt(os.path.join(POLYGONDIR, 'china_mainland_hi.txt')))
    path_hn = mpath.Path(np.loadtxt(os.path.join(POLYGONDIR, 'china_hainan_hi.txt')))
    path_tw = mpath.Path(np.loadtxt(os.path.join(POLYGONDIR, 'china_taiwan_hi.txt')))
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
        iswater = getattr(ncf, 'ISWATER')
        isoilwater = getattr(ncf, 'ISOILWATER')

        land = ncf.variables['XLAND'][0,:,:].astype('i')
        landinchn = np.logical_and(land==1, inchn)
        nlandinchn = np.logical_not(landinchn)

        land[nlandinchn] = 2.0
        ncf.variables['XLAND'][0,:,:] = land[:]

        landmask = ncf.variables['LANDMASK'][0,:,:]
        landmask[nlandinchn] = 0.0
        ncf.variables['LANDMASK'][0,:,:] = landmask[:]

        lu_index = ncf.variables['LU_INDEX'][0,:,:]
        lu_index[nlandinchn] = iswater
        ncf.variables['LU_INDEX'][0,:,:] = lu_index

        isltyp = ncf.variables['ISLTYP'][0,:,:]
        isltyp[nlandinchn] = isoilwater
        ncf.variables['ISLTYP'][0,:,:] = isltyp[:]

        ivgtyp = ncf.variables['IVGTYP'][0,:,:]
        ivgtyp[nlandinchn] = iswater
        ncf.variables['IVGTYP'][0,:,:] = ivgtyp[:]

        hgt = ncf.variables['HGT'][0,:,:]
        hgt[nlandinchn] = 0.0
        ncf.variables['HGT'][0,:,:] = hgt

        tmn = ncf.variables['TMN'][0,:,:]
        sst = ncf.variables['SST'][0,:,:]
        tmn[nlandinchn] = sst[nlandinchn]
    return

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mask region out of China as sea.')
    parser.add_argument('wrfinput', type=str,
                        help='wrfinput file')
    args = parser.parse_args()
    main(args.wrfinput)
