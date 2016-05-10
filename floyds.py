####################################################################################
# Script to simulate LCOGT FLOYDS results
# Katie O'Neill, Liam Kirkpatrick (with some help from jswift)
#
# KO 3/29/16: Divided into two preliminary functions, read_spectra and bin_spectra
# KO/LK 3/29/16: Met study hall to discuss progress/next steps 
# KO/LK 4/4/16: Met study hall to regrid_spectra
# KO/LK 4/5/16: Met study hall to revise regrid_spectra
# KO 4/7/16: Tried to fix regrid_spectra, failed. Started get_logg
# KO/LK 4/12/16: added poisson noise to bin_spectra
# KO 4/13/16: started work on get_values, combined get_values and get_logg
# KO 4/19/16: started get_snr
# KO 4/20/16: added comments, cleaned up
# jswift 4/24/16: general clean up, separated add_noise function from rebin_spectrum,
#                 import homegrown constants (in utils) instead of scipy constants,
#                 smooth_spectrum, and flatten_spectrum added as key elements in
#                 preparation for autocorrelation
# KO/LK 5/5/16: Started cross_correlate
# jswift 5/6/16: Bug fixes. There was a bug involving numerical error in regrid_spec
# KO 5/8/16: Started work on velocity x-axis conversion in cross_correlate
####################################################################################


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import resample
from scipy import interpolate
import pdb
import constants as c
import matplotlib.patches as mpatches
import glob as glob
import math
from scipy.signal import savgol_filter
from scipy.signal import correlate

def read_spectrum(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
                  wavefile='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits', plot=False):
     '''
     Read spectrum and plot (no frills version)
     '''
     spec,spech = fits.getdata(specfile,header=True)
     wave,waveh = fits.getdata(wavefile,header=True)
     if plot:
         plt.ion()
         plt.figure(1)
         plt.clf()
         plt.plot(wave,spec)
         plt.xlabel('Wavelength')
         plt.ylabel('Flux')
                        
     return wave, spec

######################################################################
def bin_spectrum(wave,spec, R=550,plot=False):
     '''
     Read spectrum and plot, limited to FLOYDS wavelength values
     (js: you will want noise addition to be a separate procedure)

     '''

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
 
     if plot:
          plt.ion()
          plt.figure(2)
          plt.clf()
          plt.plot(wave,spec,'k-')
          plt.plot(wave_resamp,spec_resamp,'r-')
          plt.xlabel('Wavelength')
          plt.ylabel('Flux')
          black_patch = mpatches.Patch(color='K', label = 'Original')
          blue_patch = mpatches.Patch(color='R', label = 'Re-binned')
          plt.legend(handles=[black_patch,blue_patch], loc=4)
          plt.xlim(np.min(wave[inds]),np.max(wave[inds])) 
          plt.show()
     
     return wave_resamp, spec_resamp



######################################################################
def add_noise(wave,spec,SNR=50.0,plot=False):
     '''
     Add poisson noise to a binned spectrum
     '''

     scale = SNR**2
     scaled_s = (scale*spec)/(np.median(spec))
     noisy_spec = np.random.poisson(scaled_s)

     if plot:
          plt.ion()
          plt.figure(3)
          plt.clf()
          plt.plot(wave,scaled_s,'k-')
          plt.plot(wave,noisy_spec,'r-')
          plt.xlabel('Wavelength')
          plt.ylabel('Flux')
          black_patch = mpatches.Patch(color='K', label = 'Original')
          red_patch = mpatches.Patch(color='R', label = 'Noisy')
          plt.legend(handles=[black_patch,red_patch], loc=4)
          plt.xlim(np.min(wave),np.max(wave)) 
          plt.show()

     return wave, noisy_spec
          


######################################################################
def regrid_spectrum(wave,spec,plot=False):
     """
     Resample spectrum into logspace
     """

     startwave = np.min(wave)
     stopwave = np.max(wave)
     lnwave = np.linspace(np.log(startwave),np.log(stopwave),len(wave))
     wave_logspace = np.exp(lnwave)

     # Fix for numerical errors in going to log and normal space
     if np.min(wave_logspace) < startwave:
          wave_logspace[0] = startwave
     if np.max(wave_logspace) > stopwave:
          wave_logspace[-1] = stopwave

     wave_interpolate = interpolate.interp1d(wave, spec)
     
     wave_final = wave_interpolate(wave_logspace)
     
     if plot:
          plt.ion()
          plt.clf()    
          plt.figure(4)
          plt.plot(wave_final,spec,'r-')
          plt.xlim(np.min(wave[inds]),np.max(wave[inds])) 
          
     return wave_logspace, wave_final


######################################################################
def smooth_spectrum(spec,window=45,polyorder=3):
     '''
     Smooth spectrum so that it may be "flattened" in preparation for
     cross correlation
     '''
     if window % 2 == 0:
          print('Filter window must be odd!')
          return None

     smooth_spec = savgol_filter(spec,window,polyorder)

     return smooth_spec

######################################################################
def flatten_spec(spec,window=45,polyorder=3,plot=False):
     '''
     Flatten given spectrum using a Savitzky-Golay filter
     '''

     smooth_spec = smooth_spectrum(spec,window=window,polyorder=polyorder)

     flat_spec = spec/smooth_spec - 1

     if plot:
          plt.ion()
          plt.figure(5)
          plt.clf()
          plt.plot(spec/np.median(spec)-1,'k-')
          plt.plot(flat_spec,'r-')
          plt.axhline(y=0,color='g',linestyle='--')
          plt.xlabel('Pixel')
          plt.ylabel('Flux')
          black_patch = mpatches.Patch(color='K', label = 'Original')
          red_patch = mpatches.Patch(color='R', label = 'Flattened')
          plt.legend(handles=[black_patch,red_patch], loc=4)
          plt.show()
     
     return flat_spec



######################################################################
def get_logg(mass,radius):
     '''
     Return the value of logg based on mass and radius of star (in solar units)
     '''
     G = c.G
     M = c.Msun * mass 
     R = c.Rsun * radius

     return np.log10((G*M)/(R**2))

                  

######################################################################
def get_values(specfile='lte03800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'):                   
     """
     Returns values of temp, logg, and metallicity from PHOENIX file
     should eventually read the header of the file

     """

     temp = float(specfile[3:-48])
     logg = float(specfile[9:-44])
     metal = float(specfile[14:-39])
    
     return temp, logg, metal


######################################################################
def get_snr(time=700,mag=15):
     """ 
     Calculates expected FLOYDS SNR based on given time and magnitude
     """
     B = 1889.9
     T = np.sqrt(time)
     M = 10**(-0.2 * mag)
     SNR = B * T * M

     return SNR

######################################################################
def get_inttime(mag=15,SNR=20):
     """
     Calculates integration time needed to acheive a specified SNR for a 
     star of a given a magnitude.
     """
     return 10**(0.4*mag)*(SNR/1889.9)**2
     
######################################################################
def cross_correlate(SNR=10.0):
     wave,spec = read_spectrum()
     wave_resamp, spec_resamp = bin_spectrum(wave,spec)

     #With noise
     wave_resamp, noisy_spec = add_noise(wave_resamp, spec_resamp, SNR=SNR)    
     flat_spec_noise = flatten_spec(noisy_spec)
     wave_logspace_noise, flat_spec_log_noise = regrid_spectrum(wave_resamp, flat_spec_noise)

     #Without noise
     flat_spec_no_noise = flatten_spec(spec_resamp)
     wave_logspace_no_noise, flat_spec_log = regrid_spectrum(wave_resamp,flat_spec_no_noise)
    
     #Plot with noise and without noise, both flattened and regridded into logspace
     plt.ion()
     plt.figure(1)
     plt.clf()
     plt.plot(wave_logspace_noise,flat_spec_log_noise)
     plt.plot(wave_logspace_no_noise,flat_spec_log)
     
     #Cross correlate
     #Need to check if this is the right order for np.correlate
     cor = correlate(flat_spec_log_noise,flat_spec_log, mode='full' )
     plt.ion()
     plt.figure(2)
     plt.clf()
     plt.plot(cor)
     
     #Cross correlate with velocity on x-axis
     delta_lnwave = np.log(wave_logspace_noise)
     #take median of difference and convert to meters from angstroms (check unit!)
     diff = (np.diff(delta_lnwave))/(10e10)
     # set value speed light in m/s
     c = 2.99792458e8
     # create vector for n from -329 to 329
     len = len(cor)
     r = np.array(range(len))
     n_init = ((len-1)/2.0)
     n = r - n_init
     
'''  
     n =    
     v = c * n * diff
     plt.ion()
     plt.figure(3)
     plt.clf()
     plt.plot(v, cor)
     plt.xlabel('Velocity m/s')
''' 
     