
import os, sys
import copy
import numpy as np
import healpy as hp
import itertools

from scipy.stats import norm

import ephem

import astropy.coordinates
from astropy.time import Time, TimeDelta
import astropy.units as u

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.path

def read_skymap(params,is3D=False):

    filename = params["skymap"]
    map_struct = {}

    if is3D:
        healpix_data = hp.read_map(filename, field=(0,1,2,3), verbose=False)

        distmu_data = healpix_data[1]
        distsigma_data = healpix_data[2]
        prob_data = healpix_data[0]
        norm_data = healpix_data[3]

        map_struct["distmu"] = distmu_data / params["DScale"]
        map_struct["distsigma"] = distsigma_data / params["DScale"]
        map_struct["prob"] = prob_data
        map_struct["distnorm"] = norm_data

    else:
        prob_data = hp.read_map(filename, field=0, verbose=False)
        prob_data = prob_data / np.sum(prob_data)

        map_struct["prob"] = prob_data

    nside = hp.pixelfunc.get_nside(prob_data)
    nside = params["nside"]
    map_struct["prob"] = hp.ud_grade(map_struct["prob"],nside,power=-2)
    if is3D:
        map_struct["distmu"] = hp.ud_grade(map_struct["distmu"],nside,power=-2) 
        map_struct["distsigma"] = hp.ud_grade(map_struct["distsigma"],nside,power=-2) 
        map_struct["distnorm"] = hp.ud_grade(map_struct["distnorm"],nside,power=-2) 

    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5*np.pi - theta)

    map_struct["ra"] = ra
    map_struct["dec"] = dec

    sort_idx = np.argsort(prob_data)[::-1]
    csm = np.empty(len(prob_data))
    csm[sort_idx] = np.cumsum(prob_data[sort_idx])

    map_struct["cumprob"] = csm

    pixarea = hp.nside2pixarea(nside)
    pixarea_deg2 = hp.nside2pixarea(nside, degrees=True)

    map_struct["nside"] = nside
    map_struct["npix"] = npix
    map_struct["pixarea"] = pixarea
    map_struct["pixarea_deg2"] = pixarea_deg2

    return map_struct   

def read_transients(params, map_struct):

    nside = params["nside"]
    prob_data_sorted = np.sort(map_struct["prob"])[::-1]
    prob_data_indexes = np.argsort(map_struct["prob"])[::-1]
    prob_data_cumsum = np.cumsum(prob_data_sorted)

    lines = [line.rstrip('\n') for line in open(params["transientsFile"])]
    lines = lines[1:]
    lines = filter(None,lines)

    transients_struct = {}
    transients_struct["data"] = np.empty((0,8))
    transients_struct["name"] = []
    transients_struct["classification"] = []

    for line in lines:
        line = line.replace('"','')
        lineSplit = line.split(",")
        if lineSplit[0] == "NULL":
            continue

        name = lineSplit[0]
        ra = float(lineSplit[2])
        dec = float(lineSplit[3])
        classification = lineSplit[6]

        if not lineSplit[11] == "":
            mag = float(lineSplit[11])
        else:
            mag = -1

        if not lineSplit[13] == "":
            z = float(lineSplit[13])
        else:
            z = -1
        if not lineSplit[14] == "":
            mpc = float(lineSplit[14])
        else:
            mpc = -1        

        if not lineSplit[15] == "":
            discovery_mjd = Time(lineSplit[15], format='isot').mjd
        else:
            discovery_mjd = -1

        if not lineSplit[16] == "":
            nondetection_mjd = Time(lineSplit[16], format='isot').mjd
        else:
            #nondetection_mjd = -1
            nondetection_mjd = discovery_mjd

        event_mjd = Time(params["gpstime"], format='gps', scale='utc').mjd
        if (event_mjd < nondetection_mjd) or (event_mjd > discovery_mjd + params["dt"]):
            continue

        ipix = hp.ang2pix(nside, theta=ra, phi=dec, lonlat=True)
        cumprob = prob_data_cumsum[ipix]

        transients_struct["data"] = np.append(transients_struct["data"],np.array([[ra,dec,mag,discovery_mjd,nondetection_mjd,z,mpc,cumprob]]),axis=0)
        transients_struct["name"].append(name)
        transients_struct["classification"].append(classification)

    return transients_struct

