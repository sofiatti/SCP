import numpy as np
import sncosmo
import cPickle as pickle
import gzip
import time
from numpy.random import normal, randint, uniform
from scipy import io, interpolate, optimize

wavelength_cuts = {'f140w.dat': [11000.,16500.], 'f160w.dat': [13500.,17000.], 'f105w.dat': [8500.,12100.], 'f814w.uvis1.dat': [6500.,10000.]}

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
zero_point = {'f105w': 26.235, 'f140w': 26.437, 'f160w': 25.921,
              'f814w.uvis1': 25.0985, 'zpsys':'ab'}


def obsflux_Ia(filter, z, my_phase):
    """Given a filter, redshift z at given phase, generates the observed 
    magnitude for SNe Type Ia""" 
    #Type Ia model
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
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
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
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
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
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


def mc_file(n, filter, my_p_z, my_z=None):
    """Returns 3 arrays. For each SNe type, returns an array with the observed filter flux."""
    file = io.readsav(my_p_z)
    pz = file['p_z']
    if my_z==None:
        z = np.arange(0,5,.01)
    else:
        z = my_z
    if np.shape(pz)!=np.shape(z):
        raise ValueError("p_z array and z array are different sizes!")
    dz = z[1]-z[0]
    pz /= (dz * pz).sum()
    ecdf = np.cumsum(pz * dz)
    cdf = interpolate.interp1d(z,ecdf)
    fnew = cdf(z)
    
    def func(x, *args):
        my_cdf = args[0]
        cdf = args[1]
        return abs(my_cdf - cdf(x))
    out_z=[]
    for i in range(n):
        my_cdf = np.random.uniform(0,1)
        my_z = optimize.fmin(func, (1.5), args=(my_cdf, cdf), disp=0)
        out_z.append(my_z[0])
    
    flux = np.zeros((n, 3))

    type_Ia_flux = flux[:,0]
    type_Ibc_flux = flux[:,1]
    type_II_flux = flux[:,2]
    for i in range(n):
        if i%10==0:
            print i/10, '% complete'
        my_phase = uniform(0.,5.)

        my_obsflux_Ia = obsflux_Ia(filter, out_z[i], my_phase)
        type_Ia_flux[i] = my_obsflux_Ia

        my_obsflux_Ibc= obsflux_Ibc(filter, out_z[i], my_phase)
        type_Ibc_flux[i] = my_obsflux_Ibc

        my_obsflux_II = obsflux_II(filter, out_z[i], my_phase)
        type_II_flux[i] = my_obsflux_II


    fname = '/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff_photoz/files/' + filter + '_'+ str(n) + '_phase_0_5_mc.gz'
    
    keys = ['type_Ia_flux','type_Ibc_flux', 'type_II_flux','filter']
    values = [type_Ia_flux, type_Ibc_flux, type_II_flux, filter]
    mc = dict(zip(keys, values))
    #pickle.dump( montecarlo, open( new_fname, 'wb'))
    save_file(mc, fname)
    return()

for i in range (1000, 11000, 1000):
    start = time.time()
    mc_file(i, 'f160w', 'zprob_v0.06_09.000129')
    print 'It took', time.time() - start, 'seconds.'
