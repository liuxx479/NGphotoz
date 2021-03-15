import numpy as np
from matplotlib import pyplot as plt
from pylab import *
from scipy.stats import norm
from scipy.stats import uniform
from astropy.io import fits
from astropy import units as u
import os
from lenstools import ConvergenceMap

root = '/global/u1/j/jialiu/NGphotoz/'
dir_storage = root+'NGphotoz_scratch/'
dir_cosmos = dir_storage+'Cosmo_maps/'

## constants
map_side_deg = 10*u.degree
map_pix = 7745
pixel_angular_side = map_side_deg / map_pix

### parameters
sigma_e=0.26
Nbin=20
theta_g_arr = [1,5,10]
l_edges = logspace(log10(40),log10(4000),Nbin+1)
sigma_bin_edges = linspace(-5,5,Nbin+1)

### derived parameters
N_theta_g = len(theta_g_arr)
ngal_tomo = loadtxt(root+'ngal_tomo.txt')
sigma_kappa_arr = loadtxt(root+'sigma_kappa.txt')
kappa_bin_edges = kron(sigma_kappa_arr,sigma_bin_edges).reshape(3,5,-1)
sigma_pix_arr =[ (sigma_e / (pixel_angular_side * sqrt(ingal / u.arcmin**2))).decompose().value
                for ingal in ngal_tomo]


### tomo runs from 1-5, cone run from 1-5, cosmo run from 0-24, 25 is the fiducial model
cosmos = [ '%02d_%s'%(i, j) for i in range(25) for j in ['a','f']]
cosmos += ['fid_a', 'fid_f']
cosmo_dir = '/global/cscratch1/sd/jialiu/desc-sprint-raytracing/Cosmo_maps/'
cosmo_fn_gen = lambda cosmo, tomo, cone: cosmo_dir+cosmo+'/kappa_LSST-SRD_tomo%i_cone%i.fits'%(tomo, cone)

def map_stats (cosmo_tomo_cone):
    '''for fits file fn, generate ps, peaks, minima, pdf, MFs
    fn: input file name, including full path
    tomo=1, 2,..5: int, for tomographic bins
    cone=1, 2,..5: int, for light cones'''
    cosmo, tomo, cone = cosmo_tomo_cone
    fn = cosmo_fn_gen(cosmo, tomo, cone)
    print (fn)
    imap = fits.open(fn)[0].data ## open the file
    ### add noise
    ### generate random see, such that it is the same for all cosmology
    ### but different for tomo and cone
    if cosmo[-1]=='a':
        iseed=int(cone*100+tomo)
    else: ##'f' starts with a different seed from the a cosmology
        iseed=int(1000+cone*100+tomo)
    seed(iseed)
    noise_map = np.random.normal(loc=0.0, scale=sigma_pix_arr[tomo-1], size=(map_pix, map_pix))
    kappa_map = ConvergenceMap(data=imap+noise_map, angle=map_side_deg)
        
    ### compute stats
    ## 3 smoothing
    ## 9 cols: ell, ps, kappa, peak, minima, pdf, v0, v1, v2
    ps_noiseless=ConvergenceMap(data=imap, angle=map_side_deg).powerSpectrum(l_edges)
    ps_unsmoothed=kappa_map.powerSpectrum(l_edges) ## power spectrum should be computed on unsmoothed maps

    s=0
    for theta_g in theta_g_arr:        
        out_fn = dir_cosmos+cosmo+'_tomo%i_cone%i_s%i.npy'%(tomo, cone, theta_g)
        if os.path.isfile(out_fn): ## skip the calculatoin if the file is already there
            print (out_fn,'exist; skip computation.\n')
            continue
        print (out_fn,'does NOT exist; compute NG stats..\n')
        imap = kappa_map.smooth(theta_g*u.arcmin)
        out=zeros(shape=(11, Nbin)) 
        kappa_bins = kappa_bin_edges[s][tomo-1] 
        ps=imap.powerSpectrum(l_edges)
        peak=imap.peakCount(kappa_bins)
        minima = ConvergenceMap(data=-imap.data, angle=map_side_deg).peakCount(kappa_bins)
        pdf=imap.pdf(kappa_bins)
        mfs=imap.minkowskiFunctionals(kappa_bins)
        out[0] = ps[0]
        out[1] = ps_noiseless[1]
        out[2] = ps_unsmoothed[1]
        out[3] = ps[1]
        out[4] = peak[0]
        out[5] = peak[1]
        out[6] = minima[1][::-1]
        out[7] = pdf[1]
        out[8] = mfs[1]
        out[9] = mfs[2]
        out[10] = mfs[3]        
        save(out_fn, out)   
        s+=1
    #return out

cosmo_tomo_cone_arr = [[cosmo, tomo, cone] 
                       for cosmo in cosmos 
                       for tomo in range(1,6)
                       for cone in range(1,6)]

## test on one map
map_stats (cosmo_tomo_cone_arr[1])