from PyAstronomy import pyasl
import datetime as dt
from astroquery.simbad import Simbad
from astropysics.coords import AngularCoordinate as angcor
import constants as const
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil
import glob
import fnmatch
import os
import pdb
import math
#
# Formatting conversions:


def dec2sex(deci):
	""" 
	Converts a Decimal number (in hours or degrees) to Sexagesimal.
	
	Parameters
	----------
	deci : float
		A decimal number to be converted to Sexagismal.
	
	Returns
	-------
	hd : int
		hours or degrees
	m : int
		minutes or arcminutes
	s : float
		seconds or arcseconds
	
	"""
	(hfrac, hd) = math.modf(deci)
	(min_frac, m) = math.modf(hfrac * 60)
	s = min_frac * 60.
	return (int(hd), int(m), s)


def sex2dec(hd, min, sec):
	""" 
	Converts a Sexagesimal number to a Decimal number.
	
	Parameters
	----------
	hd : int
		hours or degrees
	m : int
		minutes or arcminutes
	s : float
		seconds or arcseconds
	
	Returns
	-------
	hd : float
		A decimal number
	
	"""
	return float(hd) + min/60.0 + sec/3600.0

def recursive_glob(rootdir='.', pattern='*'):
    return [os.path.join(rootdir, filename)
            for rootdir, dirnames, filenames in os.walk(rootdir)
            for filename in filenames
            if fnmatch.fnmatch(filename, pattern)]


def get_kic_file(kic):

    files = get_files()
    for f in files:
        kicmatch,path = get_kic(f)
        if kicmatch == str(kic):
            return f
        else:
            pass

def get_files():
    
    files = recursive_glob('../outdata', '*[0-9].out')

    fout = []
    for f in files:
        pdir = f.split('/')[-2]
        if pdir == 'Refine':
            fout.append(f)
        
    return fout


def get_kic(file):
    fname = file.split('/')[-1]
    kic = fname.split('.')[0]

    path = '/Users/jonswift/Astronomy/EBs/outdata/'+kic+'/Refine/'

    return kic,path



def get_info(file,secondary=False):

    fname = file.split('/')[-1]
    kic = fname.split('.')[0]

    path = '/Users/jonswift/Astronomy/EBs/outdata/'+kic+'/Refine/'
    
    if kic == '10268903':
        rast  = '19 18 51.13'
        decst = '47 20 36.60'
    elif kic == '12106934':
        rast  = '19 16 41.439'
        decst = '50 37 50.01'
    elif kic == '1575690':
        rast  = '19 27 31.658'
        decst = '37 11 20.80'
    elif kic == '2437060':
        rast  = '19 20 47.882'
	decst = '+37 46 37.16'
    elif kic == '4454219':
        rast  = '19 14 15.153'	
	decst = '39 32 29.93'
    elif kic == '5769943':
        rast  = '18 55 00.093'	
	decst = '41 05 09.64'
    elif kic == '7287391':
        rast  = '19 33 51.254'
	decst = '42 53 02.61'
    elif kic == '7938883':
        rast  = '18 48 46.342'
	decst = '43 43 59.59'
    elif kic == '8264097':
        rast  = '20 03 29.352'
	decst = '44 08 39.54'
    elif kic == '8937019':
        rast  = '18 59 09.713'	
	decst = '45 13 48.38'
    elif kic == '8949316':
        rast  = '19 22 04.551'	
	decst = '45 14 07.33'
    else:
        eb = Simbad.query_object('KIC '+kic)
        rast = eb["RA"][0]
        decst = eb["DEC"][0]

    RAdeg = angcor(rast,sghms=True).d
    DECdeg = angcor(decst).d

    M1 = 0.5 * const.Msun
    M2 = 0.4 * const.Msun

    data = np.loadtxt(file,delimiter=',')
    period = data[1]
    if secondary != False:
        t0 = data[9]
        tdur = data[13]/24.0
        tag = '_sec'
    else:
        t0 = data[3]
        tdur = data[7]/24.0
        tag = '_prim'

    sma = ((period*24.*3600.)**2 * const.G * (M1+M2)/
           (4.0*np.pi**2))**(1.0/3.0) / const.AU



    # Object data must contain
    # ra	Right ascension of object [deg]
    # dec	Declination of object [deg]
    # T0	Time reference point (HJD)
    # orbPer	Orbital period [d]
    # orbInc	Orbital inclination [deg]
    # SMA	Semi-major axis [AU]
    # RpJ	Planetary radius [Jovian radii]
    # RsSun	Stellar Radius [solar]
    # Tdur	OPTIONAL, Transit duration [d]
    
    object = {
        "plName": kic+tag,
        "ra": RAdeg,
        "dec": DECdeg,
        "T0": t0,
        "orbPer": period,
        "orbInc": 90.0,
        "SMA": sma,
        "RpJ": 5.0,
        "RsSun": 0.5,
        "Tdur": tdur  }
 
    return object



def etime(kicin,date=[2014,9,10,0],ndays=60,secondary=False,observatory='Caltech',minalt=30.0,
          twilight='nautical',buffer=1.0,write=True):

    if observatory != None:
        otag = '_'+observatory
    else:
        latstr = '+34 08 10.0'
        lonstr = '-118 07 34.5'
        alt = 100.0
        otag = ''

    if observatory == 'Lowell':
        latstr = '+35 12 10.0'
        lonstr = '-111 39 52.0'
        alt = 2210.0
    if observatory == 'Hopkins':
        latstr = '+31 40 49.35'
        lonstr = '-110 52 44.60'
        alt = 2210.0
    if observatory == 'Caltech':
        latstr = '+34 08 10.0'
        lonstr = '-118 07 34.5'
        alt = 100.0
    if observatory == 'Keck':
        latstr = '19 49 34.9'
        lonstr =  '-155 28 30.04'
        alt = 4145.0

    lat = angcor(latstr).d
    lon = angcor(lonstr).d + 360

    d = dt.datetime(date[0],date[1],date[2],date[3])
    jd = pyasl.jdcnv(d)
    
    file = get_kic_file(kicin)

    object = get_info(file,secondary=secondary)

    kic,path = get_kic(file)

    if secondary != False:
        tag = '_sec'
    else:
        tag = '_prim'

    fn = path+kic+tag+otag+'_obs.dat'
    rmfile = glob.glob(fn)

    buffer /= 24.0

    if rmfile:
        os.remove(rmfile[0])
    rmfile = []
    if write == True:
        dat = pyasl.transitTimes(jd, jd+ndays, object, lon=lon, lat=lat, alt=alt, \
                                 obsOffset=buffer,minAltitude=minalt, \
                                 showTwilight=twilight,fileOutput=fn)
    else: 
        dat = pyasl.transitTimes(jd, jd+ndays, object, lon=lon, lat=lat, alt=alt, \
                                 obsOffset=buffer,minAltitude=minalt, \
                                 showTwilight=twilight)
    if dat:
        if write == True:
            pyasl.transitVisibilityPlot(dat, markTransit=True,print2file=True)
            if observatory != None:
                fold = 'transVis-'+object["plName"]+'.png'
                rmfile = glob.glob(path+fold)
                if rmfile:
                        os.remove(rmfile[0])
                fn = path+object["plName"]+otag+'.png'
                os.rename(fold,fn)
            else:
                fn = object["plName"]+'.png'
                rmfile = glob.glob(path+fn)
                if rmfile:
                    os.remove(rmfile[0])
                shutil.move(fn,path)
        else:
            pyasl.transitVisibilityPlot(dat, markTransit=True,print2file=False)
            
    else:
        print "No eclipses found!"


def do_etimes(start=None,stop=None,date=[2015,9,19,0],ndays=1,
              observatory='Hopkins',minalt=30.0,
              twilight='nautical',buffer=1.0,write=True):

    files = get_files()

    n = len(files)
    if stop == None:
        stop = n
    if start == None:
        start = 0

    for i in np.arange(start,stop):
        f = files[i]
        kic,path = get_kic(f)
        print "Starting "+f
        print "Doing primary..."
        etime(kic,date=date,ndays=ndays,observatory=observatory,
              minalt=minalt,twilight=twilight,
              buffer=buffer,write=write)
        print "...now doing secondary"
        etime(kic,secondary=True,date=date,ndays=ndays,observatory=observatory,
              minalt=minalt,twilight=twilight,buffer=buffer,write=write)

    return

    

# This procedure is acting weird!! (works sometimes, but not other times).
def obsblock(kic,date=[2014,6,18,0]):
    print " "
    str1 = "%2d/%2d/%4d" % (date[1],date[2],date[0])
    print 'KIC '+str(kic)+'; UT '+str1

    file = get_kic_file(kic)
    object1 = get_info(file,secondary=False)
    object2 = get_info(file,secondary=True)
    
    kic,path = get_kic(file)

    d = dt.datetime(date[0],date[1],date[2],date[3])
    jd = pyasl.jdcnv(d)


### Primary ###
    # time after mid-eclipse at JD
    t1 = np.mod(jd - object1["T0"],object1["orbPer"])

    # phase of eclipse at JD
    p1 = np.mod(jd - object1["T0"],object1["orbPer"])/object1["orbPer"]

    # Time until ingress
    t1ing  = (object1["orbPer"]-object1["Tdur"]/2.0)-t1
            
    # Time until egress
    t1egr  = (object1["orbPer"]+object1["Tdur"]/2.0)-t1

    if (t1ing >= 0 and t1ing <= 0.75) or (t1egr >= 0 and t1egr <= 0.75):
            print "Potential conflict with primary eclipse:"

   # Does an ingress happen in the "night" window?
    if t1ing >= 0 and t1ing <= 0.75:

        jd1ing = t1ing + jd

        y,m,d,h = pyasl.daycnv(jd1ing)
        
        str1 = "%2d/%2d/%4d" % (m,d,y)
        hour = dec2sex(h) 
        stri = "%2d:%2d:%.2f" % hour

        # Observe before the start of ingress
        print '... observe before UT '+stri

    else: pass
        #print 'No primary eclipse conflict'

    if t1egr >= 0 and t1egr <= 0.75:
        jd1egr = t1egr + jd
        
        y,m,d,h = pyasl.daycnv(jd1egr)
        
        str1 = "%2d/%2d/%4d" % (m,d,y)
        hour = dec2sex(h) 
        stre = "%2d:%2d:%.2f" % hour
        print '... observe after  UT '+stre
    else: pass

### Secondary ###
    # time after mid-eclipse at JD
    t2 = np.mod(jd - object2["T0"],object2["orbPer"])

    # phase of eclipse at JD
    p2 = np.mod(jd - object2["T0"],object2["orbPer"])/object2["orbPer"]

    # Time until ingress
    t2ing  = (object2["orbPer"]-object2["Tdur"]/2.0)-t2
    # Time until egress
    t2egr  = (object2["orbPer"]+object2["Tdur"]/2.0)-t2
    
    if (t2ing >= 0 and t2ing <= 0.75) or (t2egr >= 0 and t2egr <= 0.75):
            print "Potential conflict with seconary eclipse:"

    # Does an ingress happen in the "night" window?
    if t2ing >= 0 and t2ing <= 0.75:
        jd2ing = t2ing + jd

        y,m,d,h = pyasl.daycnv(jd2ing)
        
        str2 = "%2d/%2d/%4d" % (m,d,y)
        hour = dec2sex(h) 
        stri = "%2d:%2d:%.2f" % hour

        # Observe before the start of ingress
        print '... observe before UT '+stri

    else: pass
        #print 'No primary eclipse conflict'
        
    if t2egr >= 0 and t2egr <= 0.75:

        jd2egr = t2egr + jd
        
        y,m,d,h = pyasl.daycnv(jd2egr)
        
        str2 = "%2d/%2d/%4d" % (m,d,y)
        hour = dec2sex(h) 
        stre = "%2d:%2d:%.2f" % hour
        print '... observe after  UT '+stre
    else: pass

        
         
def do_obsblock(start=None,stop=None,date=[2014,6,18,0]):
    files = get_files()
    
    pdir = '/Users/jonswift/Astronomy/EBs/HIRES_RVs/'

    kics = np.loadtxt(pdir+'RV_KICs.txt',dtype='string')
    
    n = len(files)
    if stop == None:
        stop = n
    if start == None:
        start = 0

    for i in np.arange(start,stop):
        f = files[i]
        kic,path = get_kic(f)
        ind, = np.where(str(kic) == kics)
        if len(ind) == 1:
                obsblock(kic,date)

    return



def make_object(name='object',ra=None,dec=None,t0=None,period=None,M1=0.5,M2=0.4,
                inc=90.0,sma=None,RpJ=5.0,RsSun=0.5,duration=None):

        RAdeg = angcor(ra,sghms=True).d
        DECdeg = angcor(dec).d

        m1 = M1 * const.Msun
        m2 = M2 * const.Msun

        tdur = duration/24.0


        sma = ((period*24.*3600.)**2 * const.G * (m1+m2)/
               (4.0*np.pi**2))**(1.0/3.0) / const.AU
        
        object = {
                "plName": name,
                "ra": RAdeg,
                "dec": DECdeg,
                "T0": t0,
                "orbPer": period,
                "orbInc": inc,
                "SMA": sma,
                "RpJ": RpJ,
                "RsSun": RsSun,
                "Tdur": tdur  }
        
        return object



def object_etime(object,date=[2015,9,15,0],ndays=60,secondary=False,
                 observatory='Hopkins',minalt=20.0,
                 twilight='nautical',buffer=1.0,write=True):

    if observatory != None:
        otag = '_'+observatory
    else:
        latstr = '+34 08 10.0'
        lonstr = '-118 07 34.5'
        alt = 100.0
        otag = ''

    if observatory == 'Lowell':
        latstr = '+35 12 10.0'
        lonstr = '-111 39 52.0'
        alt = 2210.0
    if observatory == 'Hopkins':
        latstr = '+31 40 49.35'
        lonstr = '-110 52 44.60'
        alt = 2210.0
    if observatory == 'Caltech':
        latstr = '+34 08 10.0'
        lonstr = '-118 07 34.5'
        alt = 100.0
    if observatory == 'Keck':
        latstr = '19 49 34.9'
        lonstr =  '-155 28 30.04'
        alt = 4145.0

    lat = angcor(latstr).d
    lon = angcor(lonstr).d + 360

    d = dt.datetime(date[0],date[1],date[2],date[3])
    jd = pyasl.jdcnv(d)
    

    if secondary != False:
        tag = '_sec'
    else:
        tag = '_prim'

    fn = object['plName']+otag+tag+'_obs.dat'
    rmfile = glob.glob(fn)

    buffer /= 24.0

    if rmfile:
        os.remove(rmfile[0])
    rmfile = []
    if write == True:
        dat = pyasl.transitTimes(jd, jd+ndays, object, lon=lon, lat=lat, alt=alt, \
                                 obsOffset=buffer,minAltitude=minalt, \
                                 showTwilight=twilight,fileOutput=fn)
    else: 
        dat = pyasl.transitTimes(jd, jd+ndays, object, lon=lon, lat=lat, alt=alt, \
                                 obsOffset=buffer,minAltitude=minalt, \
                                 showTwilight=twilight)
    if dat:
        if write == True:
            pyasl.transitVisibilityPlot(dat, markTransit=True,print2file=True)
            if observatory != None:
                fold = 'transVis-'+object["plName"]+'.png'
                rmfile = glob.glob(fold)
                if rmfile:
                        os.remove(rmfile[0])
                fn = object["plName"]+tag+'.png'
                os.rename(fold,fn)
            else:
                fn = object["plName"]+otag+tag+'.png'
                rmfile = glob.glob(fn)
                if rmfile:
                    os.remove(rmfile[0])
                shutil.move(fn,path)
        else:
            pyasl.transitVisibilityPlot(dat, markTransit=True,print2file=False)
    else:
        print "No eclipses found!"

