import numpy as np
import sncosmo
import cPickle as pickle
import gzip
import time
from numpy.random import normal, randint, uniform

wavelength_cuts = {'f140w.dat': [11000.,16500.], 'f160w.dat': [13500.,17000.], 'f105w.dat': [8500.,12100.], 'f814w.uvis1.dat': [6500.,10000.]}
#wavelength_cuts = {'f140w.dat': [11000.,15645.], 'f160w.dat': [13500.,15645.], 'f105w.dat': [8500.,12100.], 'f814w.uvis1.dat': [6500.,10000.]}
def filter2bandpass(filter_file):
    """Returns the sncosmo bandpass for an HST filter"""
    filter = np.loadtxt(filter_file)
    wavelength = filter[:,1]
    transmission = filter[:,2]
    index_min = np.where(wavelength == wavelength_cuts[filter_file][0])[0][0]
    index_max = np.where(wavelength ==wavelength_cuts[filter_file][1])[0][0]
    wavelength = wavelength[index_min: index_max + 1]
    transmission = transmission[index_min: index_max + 1]
    band = sncosmo.Bandpass(wavelength, transmission, name=filter_file[: -4])
    sncosmo.registry.register(band, force=True)
    return

#Only have to use once 
f105w = filter2bandpass('f105w.dat')
f140w = filter2bandpass('f140w.dat')
f160w = filter2bandpass('f160w.dat')
f814w_UVIS1 = filter2bandpass('f814w.uvis1.dat')

t0 = 0
hostr_v = 3.1 #For core collapse SNe
dust = sncosmo.CCM89Dust() #For core collapse SNe


def zero_point(filter):
    """Given a filter, returns the zero point and zero point system"""
    zpsys = 'ab'
    if filter == 'f105w':
        zp = 26.235
    elif filter == 'f140w':
        zp = 26.437
    elif filter == 'f160w':  
        zp = 25.921
    elif filter in ('f814w.uvis1', 'f814w.uvis2'):
        zp = 25.0985
    return (zp, zpsys)


def obsflux_Ia(filter, z, my_phase):
    """Given a filter, redshift z at given phase, generates the observed 
    magnitude for SNe Type Ia""" 
    #Type Ia model
    zp, zpsys = zero_point(filter)
    model_Ia = sncosmo.Model(source=sncosmo.get_source('salt2-extended', 
                                                       version='1.0'))
    
    #Type Ia: Getting parameters for salt
    alpha = 0.12
    beta = 3.
    x1 = normal(0.,1.)
    c = normal(0., 0.1)
    mabs = normal(-19.1 + alpha*x1 - beta*c, 0.15)
    model_Ia.set(z=z)
    model_Ia.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    p = {'z':z, 't0':t0, 'x1': x1, 'c': c}
    model_Ia.set(**p)

    phase = my_phase + model_Ia.source.peakphase('bessellb')
    obsflux_Ia = model_Ia.bandflux( filter, t0+(phase * (1+z)), zp, zpsys ) 
    return obsflux_Ia


def obsflux_Ibc(filter, z, my_phase): 
    """Given a filter, redshift z at given phase, generates the observed 
    magnitude for SNe Type Ibc"""
    zp, zpsys = zero_point(filter)
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper']
    obsflux_Ibc = []
    for i in model_Ibc:
        model_i = sncosmo.Model(source=sncosmo.get_source(i), effects=
                                [dust], effect_names=['host'], 
                                effect_frames=['rest'])
        mabs = normal(-17.56, 0.4)
        model_i.set(z=z)
        model_i.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        p_core_collapse = {'z':z, 't0':t0, 'hostebv': uniform(-0.1,0.65), 
                           'hostr_v': hostr_v} 
        model_i.set(**p_core_collapse)
        phase = my_phase + model_i.source.peakphase('bessellb')
        obsflux_i = model_i.bandflux( filter, t0+(phase * (1+z)), zp, zpsys ) 
        obsflux_Ibc.append(obsflux_i)
    return obsflux_Ibc[randint(len(model_Ibc))] 


def obsflux_II(filter, z, my_phase): 
    """Given a filter, redshift z at given phase, generates the 
    observed magnitude for SNe Type II"""
    zp, zpsys = zero_point(filter)
    model_II = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent-sn2l']
    obsflux_II = []
    for i in model_II:
        model_i = sncosmo.Model(source=sncosmo.get_source(i), effects=
                                [dust], effect_names=['host'], 
                                effect_frames=['rest'])
        if i == 's11-2004hx'=='nugent-sn2l': 
            mabs = normal(-17.98, 0.34)
        else:
            mabs = normal(-16.75, 0.37)
        model_i.set(z=z)
        model_i.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        p_core_collapse = {'z':z, 't0':t0, 'hostebv': uniform(-0.1,0.65), 
                           'hostr_v': hostr_v} 
        model_i.set(**p_core_collapse)  
        phase = my_phase + model_i.source.peakphase('bessellb')
        obsflux_i = model_i.bandflux( filter, t0+(phase * (1+z)), zp, zpsys ) 
        obsflux_II.append(obsflux_i)
    return obsflux_II[randint(len(model_II))] 

import math, random
random.seed(11)


def save_file(object, filename, protocol=-1):
    """Saves a compressed object to disk"""
    file = gzip.GzipFile(filename, 'wb')
    pickle.dump(object, file, protocol)
    file.close()
    return ()


def mc_file(n, filter, z):
    """Returns 3 arrays. For each SNe type, returns an array with the observed filter flux."""
    flux = np.zeros((n, 3))

    type_Ia_flux = flux[:,0]
    type_Ibc_flux = flux[:,1]
    type_II_flux = flux[:,2]
    for i in range(n):
        if i%10==0:
            print i/10, '% complete'
        my_phase = uniform(0.,5.)

        my_obsflux_Ia = obsflux_Ia(filter, z, my_phase)
        type_Ia_flux[i] = my_obsflux_Ia

        my_obsflux_Ibc= obsflux_Ibc(filter, z, my_phase)
        type_Ibc_flux[i] = my_obsflux_Ibc

        my_obsflux_II = obsflux_II(filter, z, my_phase)
        type_II_flux[i] = my_obsflux_II


    fname = '/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/files/z%0.2f' % z + '_'+ filter +'_mc.gz'
    #now we want to remove the dot
    i = fname.index('.')
    new_fname = fname[:i] + fname[i+1:]
    
    keys = ['type_Ia_flux','type_Ibc_flux', 'type_II_flux','filter']
    values = [type_Ia_flux, type_Ibc_flux, type_II_flux, filter]
    mc = dict(zip(keys, values))
    #pickle.dump( montecarlo, open( new_fname, 'wb'))
    save_file(mc, new_fname)
    return()
"""
for i in np.arange(2.00, 0.10, -0.05):
    start = time.time()
    mc_file(1000, 'f814w.uvis1', i)
    print 'It took', time.time() - start, 'seconds.'

for i in np.arange(2.00, 0.10, -0.05):
    start = time.time()
    mc_file(1000, 'f105w', i)
    print 'It took', time.time() - start, 'seconds.'

for i in np.arange(2.00, 0.10, -0.05):
    start = time.time()
    mc_file(1000, 'f140w', i)
    print 'It took', time.time() - start, 'seconds.'

for i in np.arange(2.00, 0.10, -0.05):
    start = time.time()
    mc_file(1000, 'f160w', i)
    print 'It took', time.time() - start, 'seconds.'
"""
#mc_file(100, 'f105w', 1.5)
