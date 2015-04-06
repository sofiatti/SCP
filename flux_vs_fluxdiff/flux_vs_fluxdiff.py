import random
import numpy as np
import sncosmo
import cPickle as pickle
import gzip
import time
import os
from numpy.random import normal, uniform, choice
from scipy import io, interpolate, optimize

my_dir = '/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/'

t0 = 0
hostr_v = 3.1
dust = sncosmo.CCM89Dust()
random.seed(11)

zero_point = {'f105w': 26.235, 'f140w': 26.437, 'f160w': 25.921,
              'uvf814w': 25.0985, 'zpsys': 'ab'}

all_model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc',
                 'nugent-hyper', 's11-2006jo', 'snana-2004fe',
                 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo',
                 'snana-sdss014475', 'snana-2006lc', 'snana-04d1la',
                 'snana-04d4jv', 'snana-2004gv', 'snana-2006ep',
                 'snana-2007y', 'snana-2004ib', 'snana-2005hm',
                 'snana-2006jo', 'snana-2007nc']


all_model_II = ['s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p',
                'snana-2004hx', 'snana-2005gi', 'snana-2006gq',
                'snana-2006kn', 'snana-2006jl', 'snana-2006iw',
                'snana-2006kv', 'snana-2006ns', 'snana-2007iz',
                'snana-2007nr', 'snana-2007nr', 'snana-2007kw',
                'snana-2007ky', 'snana-2007lj', 'snana-2007lb',
                'snana-2007ll', 'snana-2007nw', 'snana-2007ld',
                'snana-2007md', 'snana-2007lz', 'snana-2007lx',
                'snana-2007og', 'snana-2007ny', 'snana-2007nv',
                'snana-2007pg', 's11-2004hx', 'nugent-sn2l', 'nugent-sn2n',
                'snana-2006ez', 'snana-2006ix']
'''
all_model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc',
                 'nugent-hyper']
all_model_II = ['s11-2004hx', 's11-2005lc', 's11-2005gi', 's11-2006jl',
                'nugent-sn2p', 'nugent-sn2l']
'''


def all_mabs(sne_type):
    # The majority of Type II are IIp
    my_mabs = {'Ibc': normal(-17.6, scale=1),
               'II': normal(-16.80, scale=0.97),
               'IIn': normal(-18.62, scale=1.48),
               'IIl': normal(-17.98, scale=0.9)}
    return(my_mabs[sne_type])


def obsflux_Ia(filter, z, my_phase):
    """Given a filter, redshift z at given phase, generates the observed
    magnitude for SNe Type Ia"""
    alpha = 0.12
    beta = 3.
    x1 = normal(0., 1.)
    c = normal(0., 0.1)
    mabs = normal(-19.1 - alpha*x1 + beta*c, scale=0.51)
    # Type Ia model
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ia = sncosmo.Model(source=sncosmo.get_source('salt2',
                                                       version='2.4'))
    # Type Ia: Getting parameters for salt
    model_Ia.set(z=z)
    model_Ia.set_source_peakabsmag(mabs, 'bessellb', 'vega')
    p = {'z': z, 't0': t0, 'x1': x1, 'c': c}
    model_Ia.set(**p)
    phase = my_phase + model_Ia.source.peakphase('bessellb')
    obsflux_Ia = model_Ia.bandflux(filter, t0+(phase * (1+z)), zp, zpsys)
    return obsflux_Ia


def obsflux_cc(filter, z, sn_type, all_model):
    """Given a filter and redshift z, generates the observed
    magnitude for SNe Type Ibc or Type II"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    my_model = choice(all_model)
    if my_model in ['s11-2004hx', 'nugent-sn2l']:
        sn_type = 'IIl'
    elif my_model in ['nugent-sn2n', 'snana-2006ez', 'snana-20069ix']:
        sn_type = 'IIn'
    model = sncosmo.Model(source=sncosmo.get_source(my_model),
                          effects=[dust], effect_names=['host'],
                          effect_frames=['rest'])
    mabs = all_mabs(sn_type)
    model.set(z=z)
    model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    p_core_collapse = {'z': z, 't0': t0, 'hostebv': uniform(-0.1, 0.65),
                       'hostr_v': hostr_v}
    model.set(**p_core_collapse)
    my_phase = uniform(0., 5.)
    phase = my_phase + model.source.peakphase('bessellb')
    obsflux = model.bandflux(filter, t0+(phase * (1+z)), zp, zpsys)
    return obsflux


def save_file(object, filename, protocol=-1):
    """Saves a compressed object to disk"""
    file = gzip.GzipFile(filename, 'wb')
    pickle.dump(object, file, protocol)
    file.close()
    return ()


def z_from_photo_z(photo_z_file, n, my_z_array=None):
    my_p_z = io.readsav(photo_z_file)
    pz = my_p_z['p_z']
    if my_z_array is None:
        z = np.arange(0, 5, .01)
    else:
        z = my_z_array
    if np.shape(pz) != np.shape(z):
        print np.shape(pz),  np.shape(z)
        raise ValueError("p_z array and z array are different sizes")
    dz = z[1] - z[0]
    pz /= (dz * pz).sum()
    ecdf = np.cumsum(pz * dz)
    cdf = interpolate.interp1d(z, ecdf)

    def func(x, *args):
        my_cdf = args[0]
        cdf = args[1]
        return abs(my_cdf - cdf(x))
    out_z = []
    for i in range(n):
        my_cdf = np.random.uniform(0, 1)
        my_z = optimize.fmin(func, (1.5), args=(my_cdf, cdf), disp=0)
        out_z.append(my_z[0])
    out_z = np.asarray(out_z)
    return out_z


def mc_file(n, filter, min_phase, max_phase, my_z=None, photo_z_file=None,
            z_for_photo_z_file=None):
    """Returns 3 arrays. For each SNe type, returns an array with the observed
    filter flux."""

    if my_z is None:
        z_array = z_from_photo_z(photo_z_file, n, z_for_photo_z_file)
    files_dir = 'test/mc_files_%.0f_%.0f_%.0f' % (n, min_phase, max_phase)
    if not os.path.exists(files_dir):
        os.makedirs(files_dir)

    flux = np.zeros((n, 3))
    type_Ia_flux = flux[:, 0]
    type_Ibc_flux = flux[:, 1]
    type_II_flux = flux[:, 2]
    for i in range(n):
        if i % (n/100) == 0:
            print i/(n/100), '% complete'
        my_phase = uniform(min_phase, max_phase)
        my_obsflux_Ia = obsflux_Ia(filter, z_array[i], my_phase)
        type_Ia_flux[i] = my_obsflux_Ia

        my_obsflux_Ibc = obsflux_cc(filter, z_array[i], 'Ibc', all_model_Ibc)
        type_Ibc_flux[i] = my_obsflux_Ibc

        my_obsflux_II = obsflux_cc(filter, z_array[i], 'II', all_model_II)
        type_II_flux[i] = my_obsflux_II

    if my_z is None:
        new_fname = files_dir + '/' + filter + '_mc.gz'
        keys = ['type_Ia_flux', 'type_Ibc_flux', 'type_II_flux', 'z', 'filter']
        values = [type_Ia_flux, type_Ibc_flux, type_II_flux, z_array, filter]

    else:
        fname = files_dir + '/z%0.2f' % z_array[i] + '_' + filter + '_mc.gz'
        # now we want to remove the dot
        i = fname.index('.')
        new_fname = fname[:i] + fname[i+1:]

        keys = ['type_Ia_flux', 'type_Ibc_flux', 'type_II_flux', 'filter']
        values = [type_Ia_flux, type_Ibc_flux, type_II_flux, filter]
    mc = dict(zip(keys, values))
    # pickle.dump( montecarlo, open( new_fname, 'wb'))
    save_file(mc, new_fname)
    return()
'''
for i in np.arange(2.00, 0.5, -0.05):
    start = time.time()
    mc_file(1000, 'uvf814w', -5., 0., i)
    print 'It took', time.time() - start, 'seconds.'

for i in np.arange(2.00, 0.5, -0.05):
    start = time.time()
    mc_file(1000, 'f105w', -5., 0., i)
    print 'It took', time.time() - start, 'seconds.'

for i in np.arange(2.00, 0.5, -0.05):
    start = time.time()
    mc_file(1000, 'f140w', -5., 0, i)
    print 'It took', time.time() - start, 'seconds.'

for i in np.arange(2.00, 0.5, -0.05):
    start = time.time()
    mc_file(1000, 'f160w', 0., 5., i)
    print 'It took', time.time() - start, 'seconds.'
'''
start = time.time()
mc_file(10000, 'f105w', 0., 5., photo_z_file='zprob_v0.06_09.000129')
print 'It took', time.time() - start, 'seconds.'
