import numpy as np
from matplotlib import pyplot as plt
from pylab import *
from scipy.stats import norm
from scipy.stats import uniform
from astropy.io import fits
from astropy import units as u
import os
from lenstools import ConvergenceMap
from emcee.utils import MPIPool 

root = '/global/u1/j/jialiu/NGphotoz/'
dir_storage = root+'NGphotoz_scratch/'
dir_cosmos = dir_storage+'Cosmo_maps/'
dir_cov = dir_storage+'Cov_maps/'

## constants
map_side_deg = 10*u.degree
map_pix = 7745
pixel_angular_side = map_side_deg / map_pix

### parameters
sigma_e = 0.26
sigma_kappa = sqrt(2)*sigma_e
Nbin=10 ## using 10 bins just to be safe
theta_g_arr = [1,5,10]
l_edges = logspace(log10(40),log10(4000),Nbin+1)
sigma_bin_edges = linspace(-5,5,Nbin+1)
neff = 37.0 ## chang+2013 
ngal = neff 
## chang+2013, n=46, raw number for fiducial case

### derived parameters
N_theta_g = len(theta_g_arr)
ngal_tomo = loadtxt(root+'ngal_tomo.txt') * ngal ## use neff instead of ngal
sigma_kappa_arr = loadtxt(root+'sigma_kappa.txt')
kappa_bin_edges = kron(sigma_kappa_arr,sigma_bin_edges).reshape(3,5,-1)
sigma_pix_arr =[ (sigma_kappa / (pixel_angular_side * sqrt(ingal / u.arcmin**2))).decompose().value
                for ingal in ngal_tomo]


### tomo runs from 1-5, cone run from 1-5, cosmo run from 0-24, 25 is the fiducial model
cosmos = [ '%02d_%s'%(i, j) for i in range(25) for j in ['a','f']]
cosmos += ['fid_a', 'fid_f']
cosmo_dir = '/global/cscratch1/sd/jialiu/desc-sprint-raytracing/Cosmo_maps/'
cosmo_fn_gen = lambda cosmo, tomo, cone: cosmo_dir+cosmo+'/kappa_LSST-SRD_tomo%i_cone%i.fits'%(tomo, cone)

### cov runs from 74 to 199
cov_fn_gen = lambda tomo, cone: '/global/cscratch1/sd/jialiu/desc-sprint-raytracing/Cov_maps/kappa_LSST-SRD_tomo%i_LOS%i.fits'%(tomo, cone)
#/global/cscratch1/sd/jialiu/desc-sprint-raytracing/Cov_maps/kappa_LSST-SRD_tomo4_LOS74.fits

def map_stats (cosmo_tomo_cone):
    '''for fits file fn, generate ps, peaks, minima, pdf, MFs
    fn: input file name, including full path
    tomo=1, 2,..5: int, for tomographic bins
    cone=1, 2,..5: int, for light cones'''
    cosmo, tomo, cone = cosmo_tomo_cone
    
    ### generate random see, such that it is the same for all cosmology
    ### but different for tomo and cone
    if cosmo=='cov':
        iseed=int(10000+cone*10+tomo)
        out_dir = dir_cov
        fn = cov_fn_gen(tomo, cone)
        
#     elif cosmo=='biased':
#         iseed=20000        
        #/global/cscratch1/sd/jialiu/desc-sprint-raytracing/biased_maps/kappa_LSST-SRD_tomo5_pz_zbias0.0015_simgaz0.04_outlier0.15.txt_LOS_cone1.fits
        
    else: ## all cosmologies
        if cosmo[-1]=='a':
            iseed=int(cone*10+tomo) ## cone goes from 1 to 25, so 10 to 250
        else cosmo[-1]=='f': ##'f' starts with a different seed from the a cosmology
            iseed=int(1000+cone*10+tomo)
        out_dir = dir_cosmos
        fn = cosmo_fn_gen(cosmo, tomo, cone)
        print (fn)
    
    out_fn_arr = [out_dir+cosmo+'_tomo%i_cone%i_s%i.npy'%(tomo, cone, theta_g) 
                  for theta_g in theta_g_arr]
    if np.prod(array([os.path.isfile(out_fn) for out_fn in out_fn_arr])):
        print (out_fn,'exist; skip computation.\n')
        return 0 ### all files already exist, no need to process
    
    ########## map operations
    imap = fits.open(fn)[0].data ## open the file
    
    ### add noise
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
        out_fn = out_dir+cosmo+'_tomo%i_cone%i_s%i.npy'%(tomo, cone, theta_g)
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

cov_tomo_con_arr = [['cov', tomo, cone] 
                       for tomo in range(1,6)
                       for cone in range(74,200)]
## test on one map
# map_stats (cosmo_tomo_cone_arr[1])

pool=MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

out=pool.map(map_stats, cosmo_tomo_cone_arr+cov_tomo_con_arr)
pool.close()
print ('DONE-DONE-DONE')