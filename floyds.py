######################################
# Script to simulate LCOGT FLOYDS results
# Katie O'Neill, Liam Kirkpatrick
#
#KO 3/29/16: Divided into two preliminary functions, read_spectra and bin_spectra
#KO/LK 3/29/16: Met study hall to discuss progress/next steps 
#KO/LK 4/4/16: Met study hall to regrid_spectra
#KO/LK 4/5/16: Met study hall to revise regrid_spectra
#KO 4/7/16: Tried to fix regrid_spectra, failed. Started get_logg
#KO/LK 4/12/16: added poisson noise to bin_spectra
#KO 4/13/16: started work on get_values, combined get_values and get_logg
#KO 4/19/16: started get_snr
#KO 4/20/16: added comments, cleaned up
######################################


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import resample
from scipy import interpolate
import pdb
from scipy.constants import constants as c
import matplotlib.patches as mpatches
import glob as glob
import math

# read spectrum and plot (no frills version)
def read_spectrum(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
                 wavefile='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits', plot=True):
     spec,spech = fits.getdata(specfile,header=True)
     wave,waveh = fits.getdata(wavefile,header=True)
     if plot:
         plt.ion()
         plt.figure(1)
         plt.clf()
         plt.plot(wave,spec)
         plt.xlabel('Wavelength')
         plt.ylabel('Flux')
                        
     return spec, wave
     
# read spectrum and plot, limited to FLOYDS wavelength values, with poisson noise added
def bin_spectrum(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
                 wavefile='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits', R=550):
    spec = fits.getdata(specfile)
    wave = fits.getdata(wavefile)
    inds, = np.where((wave >= 5400) & (wave <= 10000))
    plt.figure(3)
    plt.clf()
    plt.plot(wave[inds],spec[inds])
    dl = np.median(wave[inds])/R
    num = np.int(np.round((np.max(wave[inds])-np.min(wave[inds]))/dl))
    wave_resamp = []
    spec_resamp = []
    for i in range(num):
        try:
            bin, = np.where( (wave >= np.min(wave[inds])+ dl*i) &
                      (wave < np.min(wave[inds]) + dl*(i+1)) )
            wave_resamp = np.append(wave_resamp,np.mean(wave[bin]))
            spec_resamp = np.append(spec_resamp,np.mean(spec[bin]))
        except:
            print 'Skipping iteration '+str(i)
    s = spec_resamp
    scaled_s = (1000*s)/(np.median(s))
    noisy_spec = np.random.poisson(scaled_s)
    plt.clf()
    plt.figure(2)
    plt.plot(wave_resamp,noisy_spec,'r-')
    plt.plot(wave_resamp,scaled_s,'g-')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    red_patch = mpatches.Patch(color='R', label = 'Noise')
    green_patch = mpatches.Patch(color='G', label = 'No noise')
    plt.legend(handles=[red_patch,green_patch], loc=4)
    plt.show()
    plt.xlim(np.min(wave[inds]),np.max(wave[inds])) 
    
    return wave_resamp, noisy_spec

# read spectrum and plot, resampled into logspace
def regrid_spectrum(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
                 wavefile='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'):
    spec = fits.getdata(specfile)
    wave = fits.getdata(wavefile)
    inds, = np.where((wave >= 5400) & (wave <= 10000))
    startwave = 5400
    stopwave = 10000   
    lnwave = np.linspace(np.log(startwave),np.log(stopwave),len(wave))
    wave_logspace = np.exp(lnwave)
    wave_interpolate = interpolate.interp1d(wave, spec)
    wave_final = wave_interpolate(wave_logspace)
    plt.clf()    
    plt.ion()
    plt.figure(6)
    plt.plot(wave_final,spec,'r-')
    plt.xlim(np.min(wave[inds]),np.max(wave[inds])) 
    plt.show()
    
    return wave_logspace, spec

'''
# will eventually return the value of logg based on mass and radius of star
def get_logg(mass,radius):
    G = c.G
    #M = 
    # need to figure out how to input M
    R = float(specfile[9:-44])
    logg = np.log((G*M)/(R^2)
'''
    
# returns values of Teff, Radius, and Metallicity from PHOENIX file
    # should eventually read the header of the file
def get_values(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'):
    temp = float(specfile[3:-48])
    radius = float(specfile[9:-44])
    metal = float(specfile[14:-39])
    
    return temp, radius, metal
    
# calculates SNR based on given Teff and magnitude
def get_snr(teff=700,mag=15):
    B = 1.89 * (10 ** -6)
    # this value is incorrect but placehold for now
    T = math.sqrt(teff)
    M = 10 ** (-0.4 * mag)
    SNR = B * T * M
    # this returns an obviously false result
    
    return SNR