import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate

lam = np.linspace(100,200,1000)
refspec1 = 10 - np.exp(-(lam-150.0)**2/(2*2**2))
refspec2 = 10 - np.exp(-(lam-120.0)**2/(2*2**2))

spec = refspec1 + np.random.normal(0,0.5,len(refspec1))
    
plt.ion()
plt.figure(1)
plt.clf()
plt.plot(lam,spec)
plt.plot(lam,refspec2,'r-')
plt.ylim(0,12)


xcor = correlate(spec-10,refspec2-10,mode='full')

# Check lag
lag = np.linspace(-999,999,1999)
plt.figure(2)
plt.clf()
plt.plot(lag,xcor)

# Peak lag is at ~ 100 PIXELS
# Wavelength scale is ~0.1 per pixel.
# Wavelength shift of ~ 10
