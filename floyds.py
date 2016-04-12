import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import resample
from scipy import interpolate
import pdb
from scipy.constants import constants as c
import matplotlib.patches as mpatches
import glob as glob


#KO 3/29: Divided into two preliminary functions, read_spectra and bin_spectra
#KO/LK 3/29: Met study hall to discuss progress/next steps 
#KO/LK 4/4: Met study hall to regrid_spectra
#KO/LK 4/5: Met study hall to revise regrid_spectra
#KO 4/7: Tried to fix regrid_spectra, failed. Started get_logg
#KO/LK 4/12: added poisson noise to bin_spectra
#KO 4/13: started work on get_values

def read_spectra(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
                 wavefile='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits', plot=True):
     spec,spech = fits.getdata(specfile,header=True)
     wave,waveh = fits.getdata(wavefile,header=True)
     if plot:
         plt.ion()
         plt.figure(1)
         plt.clf()
         plt.plot(wave,spec)
                        
     return spec, wave
     
def get_values(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
                 wavefile='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'):
    temp = float(specfile[3:-48])
    radius = float(specfile[9:-44])
    
    return temp, radius
     
def bin_spectra(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
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

def regrid_spectra(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
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
def get_logg(mass,radius):
    G = c.G
    #M = 
    #R = 
    logg = np.log((G*M)/(R^2)

# Load synthetic spectrum
#specfile = 'lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
#spec,spech = fits.getdata(specfile,header=True)

# Load corresponding wavelength file
#wavefile = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
#wave,waveh = fits.getdata(wavefile,header=True)

# Plot the entire spectrum
#plt.ion()
#plt.figure(1)
#plt.clf()
#plt.plot(wave,spec)

# Spectral grasp of FLOYDS
#inds, = np.where((wave >= 5400) & (wave <= 10000))
#plt.figure(2)
#plt.clf()
#plt.plot(wave[inds],spec[inds])

# Floyds specifications
#R = 550 # Resolution = lambda over delta lambda
#dl = np.median(wave[inds])/R
# Number of resolution elements across spectrum
#num = np.int(np.round((np.max(wave[inds])-np.min(wave[inds]))/dl))

# Resample spectrum to the resolution of Floyds
#wave_resamp = []
#spec_resamp = []
#for i in range(num):
    #try:
     #   bin, = np.where( (wave >= np.min(wave[inds])+ dl*i) &
      #                (wave < np.min(wave[inds]) + dl*(i+1)) )
       # wave_resamp = np.append(wave_resamp,np.mean(wave[bin]))
        #spec_resamp = np.append(spec_resamp,np.mean(spec[bin]))
    #except:
     #   print 'Skipping iteration '+str(i)

# Overplot the resampled spectrum
#plt.figure(2)
#plt.plot(wave_resamp,spec_resamp,'r-')
#plt.xlim(np.min(wave[inds]),np.max(wave[inds]))
'''
