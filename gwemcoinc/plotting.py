
import os, sys, copy
import numpy as np
import healpy as hp

from scipy.stats import norm

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import matplotlib.pyplot as plt

def add_edges():

    hp.graticule(verbose=False)
    plt.grid(True)
    lons = np.arange(-150.0,180,30.0)
    lats = np.zeros(lons.shape)
    for lon, lat in zip(lons,lats):
        hp.projtext(lon,lat,"%.0f"%lon,lonlat=True)
    lats = np.arange(-60.0,90,30.0)
    lons = np.zeros(lons.shape)
    for lon, lat in zip(lons,lats):
        hp.projtext(lon,lat,"%.0f"%lat,lonlat=True)

def skymap(params,map_struct):

    unit='Gravitational-wave probability'
    cbar=False

    lons = np.arange(-150.0,180,30.0)
    lats = np.zeros(lons.shape)

    plotName = os.path.join(params["outputDir"],'prob.pdf')
    hp.mollview(map_struct["prob"],title='',unit=unit,cbar=cbar)
    add_edges()
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

    if "distmu" in map_struct:
        plotName = os.path.join(params["outputDir"],'dist.pdf')
        hp.mollview(map_struct["distmu"],unit='Distance [Mpc]',min=0.0,max=100.0)
        add_edges()
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

def transients(params, map_struct, transients_struct):

    unit='Gravitational-wave probability'
    cbar=False

    ra = transients_struct["data"][:,0]
    dec = transients_struct["data"][:,1]

    plotName = os.path.join(params["outputDir"],'transients.pdf')
    plt.figure()
    hp.mollview(map_struct["prob"],unit=unit,cbar=cbar)
    hp.projplot(ra, dec, 'wx', lonlat=True, coord='G')
    add_edges()
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

