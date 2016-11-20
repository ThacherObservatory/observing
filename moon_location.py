import ephem # See http://rhodesmill.org/pyephem/
import numpy as np
import datetime as dt 
import matplotlib.pyplot as plt

def moon_radec(year=2017,month=1,day=1,epoch=2000):
    """
    Returns the apparent geocentric Right Ascension and Declination of the 
    Moon for the given year, month and day
    """
    m = ephem.Moon()
    datestr = str(year)+'/'+str(month)+'/'+str(day)
    m.compute(datestr, epoch=str(epoch))

    return str(m.g_ra),str(m.g_dec)

