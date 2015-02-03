import numpy as np
import sncosmo
import matplotlib.pyplot as plt
from numpy.random import normal, randint, uniform

wavelength_cuts = {'f140w.dat': [11000.,16500.], 'f160w.dat': [13500.,17000.], 'f105w.dat': [8500.,12100.], 
                   'f814w.uvis1.dat': [6500.,10000.]}
def filter2bandpass(filter_file):
    """Returns the sncosmo bandpass for an HST filter"""
    filter = np.loadtxt(filter_file)
    wavelength = filter[:,1]
    transmission = filter[:,2]
    index_min = np.where(wavelength == wavelength_cuts[filter_file[54:]][0])[0][0]
    index_max = np.where(wavelength ==wavelength_cuts[filter_file[54:]][1])[0][0]
    wavelength = wavelength[index_min: index_max + 1]
    transmission = transmission[index_min: index_max + 1]
    band = sncosmo.Bandpass(wavelength, transmission, name=filter_file[54: -4])
    sncosmo.registry.register(band, force=True)

#Only have to use once                                                                                        
f105w = filter2bandpass('/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/f105w.dat')
f140w = filter2bandpass('/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/f140w.dat')
f160w = filter2bandpass('/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/f160w.dat')
f814w_uvis1 = filter2bandpass('/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/f814w.uvis1.dat')
t0 = 0
hostr_v = 3.1 #For core collapse SNe                                                                          
dust = sncosmo.CCM89Dust() #For core collapse SNe                                                             

zero_point = {'f105w': 26.235, 'f140w': 26.437, 'f160w': 25.921,
              'f814w.uvis1': 25.0985, 'zpsys':'ab'}

def lightcurve_Ia(filter, z, x1, c):
    """Given a filter and redshift z, generates the observed                      
    flux for SNe Type Ia"""
    #Type Ia model                                                                            
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ia = sncosmo.Model(source=sncosmo.get_source('salt2-extended',
                                                       version='1.0'))
    #Type Ia: Getting parameters for salt                                                     
    alpha = 0.12
    beta = 3.
    mabs = -19.1 + alpha*x1 - beta*c
    model_Ia.set(z=z)
    model_Ia.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    p = {'z':z, 't0':t0, 'x1': x1, 'c': c}
    model_Ia.set(**p)
    phase_array = np.linspace(model_Ia.mintime(), model_Ia.maxtime(),100)
    obsflux_Ia = model_Ia.bandflux( filter, phase_array, zp, zpsys )
    keys = ['phase_array','obsflux']
    values = [phase_array, obsflux_Ia]
    dict_Ia = dict(zip(keys,values)) 
    return (dict_Ia)

def lightcurve_Ibc(filter, z, hostebv_Ibc):
    """Given a filter, redshift z at given phase, generates the observed                      
    magnitude for SNe Type Ibc"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper']
    obsflux_Ibc = []
    phase_arrays = []
    for i in model_Ibc:
        model_i = sncosmo.Model(source=sncosmo.get_source(i), effects=
                                [dust], effect_names=['host'],
                                effect_frames=['rest'])
        mabs = -17.56
        model_i.set(z=z)
        phase_array_i = np.linspace(model_i.mintime(), model_i.maxtime(),100)
        model_i.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        p_core_collapse = {'z':z, 't0':t0, 'hostebv': hostebv_Ibc,
                           'hostr_v': hostr_v}
        model_i.set(**p_core_collapse)
        phase_arrays.append(phase_array_i)
        obsflux_i = model_i.bandflux( filter, phase_array_i, zp, zpsys )
        obsflux_Ibc.append(obsflux_i)
    keys = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper']
    values = [[obsflux_Ibc[0],phase_arrays[0]],[obsflux_Ibc[1],phase_arrays[1]],[obsflux_Ibc[2],phase_arrays[2]]
              ,[obsflux_Ibc[3],phase_arrays[3]],[obsflux_Ibc[4],phase_arrays[4]]]
    dict_Ibc = dict(zip(keys,values)) 
    return (dict_Ibc)

def lightcurve_II(filter, z, hostebv_II):
    """Given a filter and redshift z, generates the                               
    observed magnitude for SNe Type II"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_II = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent-sn2l']
    obsflux_II = []
    phase_arrays = []
    for i in model_II:
        model_i = sncosmo.Model(source=sncosmo.get_source(i), effects=
                                [dust], effect_names=['host'],
                                effect_frames=['rest'])
        if i == 's11-2004hx'=='nugent-sn2l':
            mabs = -17.98
        else:
            mabs = -16.75
        model_i.set(z=z)
        phase_array_i = np.linspace(model_i.mintime(), model_i.maxtime(),100)
        model_i.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        p_core_collapse = {'z':z, 't0':t0, 'hostebv': hostebv_II,
                           'hostr_v': hostr_v}
        model_i.set(**p_core_collapse)
        phase_arrays.append(phase_array_i)
        obsflux_i = model_i.bandflux( filter, phase_array_i, zp, zpsys )
        obsflux_II.append(obsflux_i)
    keys = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent\
-sn2l']
    values = [[obsflux_II[0],phase_arrays[0]],[obsflux_II[1],phase_arrays[1]],[obsflux_II[2],phase_arrays[2]]
              ,[obsflux_II[3],phase_arrays[3]],[obsflux_II[4],phase_arrays[4]],[obsflux_II[5],phase_arrays[5]]]
    dict_II = dict(zip(keys,values))
    return (dict_II)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

def plot_Ia(z, x1, c, data_flux, data_flux_err, phase):
    plt.figure(figsize = (3,9))
    gs1 = gridspec.GridSpec(3, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (0,1))
    ax3 = plt.subplot2grid((3,3), (0, 2))
    
    lightcurve_00 = lightcurve_Ia('f140w', z, x1, c)
    lightcurve_10 = lightcurve_Ia('f105w', z, x1, c)
    lightcurve_20 = lightcurve_Ia('f814w.uvis1', z, x1, c)
    
    ax1.plot(lightcurve_00['phase_array'], lightcurve_00['obsflux'])
    ax1.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax1.set_ylabel('Flux (counts/s)', size=14)
    ax1.set_title('F140W',size=14)
    ax1.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax2.plot(lightcurve_10['phase_array'], lightcurve_10['obsflux'])
    ax2.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax2.set_title('F105W',size=14)
    ax2.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
        
    ax3.plot(lightcurve_20['phase_array'], lightcurve_20['obsflux'])
    ax3.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax3.set_title('F814W.UVIS', size=14)
    ax3.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
     
    fig = plt.gcf()
    fig.set_size_inches(18.5,12.5)

def plot_Ibc(z, hostebv_Ibc, data_flux, data_flux_err, phase):
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper']
    colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
 (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
 (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
 (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
 (0.8, 0.7254901960784313, 0.4549019607843137),
 (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]
    
    plt.figure(figsize = (3,9))
    gs1 = gridspec.GridSpec(3, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (0,1))
    ax3 = plt.subplot2grid((3,3), (0, 2))
    
    lightcurve_00 = lightcurve_Ibc('f140w', z, hostebv_Ibc)
    lightcurve_10 = lightcurve_Ibc('f105w', z, hostebv_Ibc)    
    lightcurve_20 = lightcurve_Ibc('f814w.uvis1', z, hostebv_Ibc)
    
    for i in range(len(model_Ibc)):
        ax1.plot(lightcurve_00[model_Ibc[i]][1], lightcurve_00[model_Ibc[i]][0], color=colors[i])
        ax1.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
   # ax1.set_xlabel('x-label', fontsize=12)
    ax1.set_ylabel('Flux (counts/s)', size=14)
    ax1.set_title('F140W',size=14)
    ax1.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax2.plot(lightcurve_10[model_Ibc[i]][1], lightcurve_10[model_Ibc[i]][0], color=colors[i])
        ax2.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax2.set_title('F105W',size=14)
    ax2.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
        
    for i in range(len(model_Ibc)):
        ax3.plot(lightcurve_20[model_Ibc[i]][1], lightcurve_20[model_Ibc[i]][0], color=colors[i])
        ax3.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax3.set_title('F184.UVIS', size=14)
    ax3.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)    
    
    #plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    fig = plt.gcf()
    fig.set_size_inches(18.5,12.5)
    

def plot_II(z, hostebv_II, data_flux, data_flux_err, phase):
    model_II = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent\
-sn2l']
    colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
 (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
 (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
 (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
 (0.8, 0.7254901960784313, 0.4549019607843137),
 (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]
    
    plt.figure(figsize = (3,9))
    gs1 = gridspec.GridSpec(3, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (0,1))
    ax3 = plt.subplot2grid((3,3), (0, 2))
    
    lightcurve_00 = lightcurve_II('f140w', z, hostebv_II)
    lightcurve_10 = lightcurve_II('f105w', z, hostebv_II)
    lightcurve_20 = lightcurve_II('f814w.uvis1', z, hostebv_II)
    
    for i in range(len(model_II)):
        ax1.plot(lightcurve_00[model_II[i]][1], lightcurve_00[model_II[i]][0], color=colors[i])
        ax1.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
   # ax1.set_xlabel('x-label', fontsize=12)
    ax1.set_ylabel('Flux (counts/s)', size=14)
    ax1.set_title('F140W',size=14)
    ax1.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax2.plot(lightcurve_10[model_II[i]][1], lightcurve_10[model_II[i]][0], color=colors[i])
        ax2.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax2.set_title('F105W',size=14)
    ax2.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
        
    for i in range(len(model_II)):
        ax3.plot(lightcurve_20[model_II[i]][1], lightcurve_20[model_II[i]][0], color=colors[i])
        ax3.errorbar(phase, data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax3.set_title('F184W.UVIS', size=14)
    ax3.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    #plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    fig = plt.gcf()
    fig.set_size_inches(18.5,12.5)
    

