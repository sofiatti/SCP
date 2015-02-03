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

def lightcurve_Ia(filter, z, my_phase):
    """Given a filter, redshift z at given phase, generates the observed                      
    magnitude for SNe Type Ia"""
    #Type Ia model                                                                            
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ia = sncosmo.Model(source=sncosmo.get_source('salt2-extended',
                                                       version='1.0'))
    phase_array = np.linspace(-40,140,200)
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
    obsmag_Ia = model_Ia.bandmag(filter, 'ab', t0 + (phase * (1 + z)))
    #obsmag_Ia = model_Ia.bandmag( filter, phase_array, zp, zpsys )
    f = 'f105'
    sncosmo.plot_lc(model=my_model, bands=[f], zp=zp)
    plt.show()
    keys = ['phase_array','obsflux', 'phase']
    values = [phase_array, obsflux_Ia, phase]
    dict_Ia = dict(zip(keys,values)) 
    return (dict_Ia)

def lightcurve_Ibc(filter, z, my_phase):
    """Given a filter, redshift z at given phase, generates the observed                      
    magnitude for SNe Type Ibc"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper']
    obsflux_Ibc = []
    phase = []
    phase_array = np.linspace(-40,140,200)
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
        phase_i = my_phase + model_i.source.peakphase('bessellb')
        phase.append(phase_i)
        obsflux_i = model_i.bandflux( filter, phase_array, zp, zpsys )
        obsflux_Ibc.append(obsflux_i)
    keys = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper',
            'phase_array']
    values = [[obsflux_Ibc[0],phase[0]],[obsflux_Ibc[1],phase[1]],[obsflux_Ibc[2],phase[2]]
              ,[obsflux_Ibc[3],phase[3]],[obsflux_Ibc[4],phase[4]],phase_array]
    dict_Ibc = dict(zip(keys,values)) 
    return (dict_Ibc)

def lightcurve_II(filter, z, my_phase):
    """Given a filter, redshift z at given phase, generates the                               
    observed magnitude for SNe Type II"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_II = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent-sn2l']
    obsflux_II = []
    phase_array = np.linspace(-40,140,200)
    phase = []
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
        phase_i = my_phase + model_i.source.peakphase('bessellb')
        phase.append(phase_i)
        obsflux_i = model_i.bandflux( filter, phase_array, zp, zpsys )
        obsflux_II.append(obsflux_i)
    keys = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent\
-sn2l','phase_array']
    values = [[obsflux_II[0],phase[0]],[obsflux_II[1],phase[1]],[obsflux_II[2],phase[2]]
              ,[obsflux_II[3],phase[3]],[obsflux_II[4],phase[4]],[obsflux_II[5],phase[5]],
              phase_array]
    dict_II = dict(zip(keys,values))
    return (dict_II)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

def plot_Ia(z, data_flux, data_flux_err, phase1=-10, phase2=-5, phase3=0):
    plt.figure(figsize = (9,9))
    gs1 = gridspec.GridSpec(9, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (0,1))
    ax3 = plt.subplot2grid((3,3), (0, 2))
    ax4 = plt.subplot2grid((3,3), (1, 0))
    ax5 = plt.subplot2grid((3,3), (1, 1))
    ax6 = plt.subplot2grid((3,3), (1, 2))
    ax7 = plt.subplot2grid((3,3), (2, 0))
    ax8 = plt.subplot2grid((3,3), (2, 1))
    ax9 = plt.subplot2grid((3,3), (2, 2))
    
    lightcurve_00 = lightcurve_Ia('f140w', z, phase1)
    lightcurve_01 = lightcurve_Ia('f140w', z, phase2)
    lightcurve_02 = lightcurve_Ia('f140w', z, phase3)
    
    lightcurve_10 = lightcurve_Ia('f105w', z, phase1)
    lightcurve_11 = lightcurve_Ia('f105w', z, phase2)
    lightcurve_12 = lightcurve_Ia('f105w', z, phase3)
    
    lightcurve_20 = lightcurve_Ia('f814w.uvis1', z, phase1)
    lightcurve_21 = lightcurve_Ia('f814w.uvis1', z, phase2)
    lightcurve_22 = lightcurve_Ia('f814w.uvis1', z, phase3)
    
    ax1.plot(lightcurve_00['phase_array'], lightcurve_00['obsflux'])
    ax1.errorbar(lightcurve_00['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax1.set_ylabel('Flux for f140w (counts/s)', size=14)
    ax1.set_title('Phase = -10',size=14)
    ax1.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax2.plot(lightcurve_01['phase_array'], lightcurve_01['obsflux'])
    ax2.errorbar(lightcurve_01['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax2.set_title('Phase = -5',size=14)
    ax2.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
        
    ax3.plot(lightcurve_02['phase_array'], lightcurve_02['obsflux'])
    ax3.errorbar(lightcurve_02['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax3.set_title('Phase = 0', size=14)
    ax3.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax4.plot(lightcurve_10['phase_array'], lightcurve_10['obsflux'])
    ax4.errorbar(lightcurve_10['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax4.set_ylabel('Flux for f105w (counts/s)',size=14)
    ax4.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax4.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax5.plot(lightcurve_11['phase_array'], lightcurve_11['obsflux'])
    ax5.errorbar(lightcurve_11['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax5.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax5.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax6.plot(lightcurve_12['phase_array'], lightcurve_12['obsflux'])
    ax6.errorbar(lightcurve_12['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax6.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax6.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax7.plot(lightcurve_20['phase_array'], lightcurve_20['obsflux'])
    ax7.errorbar(lightcurve_20['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax7.set_ylabel('Flux for f814w.uvis1 (counts/s)',size=14)
    ax7.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax7.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax8.plot(lightcurve_21['phase_array'], lightcurve_21['obsflux'])
    ax8.errorbar(lightcurve_21['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax8.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax8.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    ax9.plot(lightcurve_22['phase_array'], lightcurve_22['obsflux'])
    ax9.errorbar(lightcurve_22['phase'], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax9.annotate('Type Ia',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax9.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    fig = plt.gcf()
    fig.set_size_inches(18.5,12.5)

def plot_Ibc(z, data_flux, data_flux_err, phase1=-10, phase2=-5, phase3=0):
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc','nugent-hyper']
    colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
 (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
 (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
 (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
 (0.8, 0.7254901960784313, 0.4549019607843137),
 (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]
    
    plt.figure(figsize = (9,9))
    gs1 = gridspec.GridSpec(9, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (0,1))
    ax3 = plt.subplot2grid((3,3), (0, 2))
    ax4 = plt.subplot2grid((3,3), (1, 0))
    ax5 = plt.subplot2grid((3,3), (1, 1))
    ax6 = plt.subplot2grid((3,3), (1, 2))
    ax7 = plt.subplot2grid((3,3), (2, 0))
    ax8 = plt.subplot2grid((3,3), (2, 1))
    ax9 = plt.subplot2grid((3,3), (2, 2))
    
    lightcurve_00 = lightcurve_Ibc('f140w', z, phase1)
    lightcurve_01 = lightcurve_Ibc('f140w', z, phase2)
    lightcurve_02 = lightcurve_Ibc('f140w', z, phase3)
    
    lightcurve_10 = lightcurve_Ibc('f105w', z, phase1)
    lightcurve_11 = lightcurve_Ibc('f105w', z, phase2)
    lightcurve_12 = lightcurve_Ibc('f105w', z, phase3)
    
    lightcurve_20 = lightcurve_Ibc('f814w.uvis1', z, phase1)
    lightcurve_21 = lightcurve_Ibc('f814w.uvis1', z, phase2)
    lightcurve_22 = lightcurve_Ibc('f814w.uvis1', z, phase3)
    
    for i in range(len(model_Ibc)):
        ax1.plot(lightcurve_00['phase_array'], lightcurve_00[model_Ibc[i]][0], color=colors[i])
        ax1.errorbar(lightcurve_00[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
   # ax1.set_xlabel('x-label', fontsize=12)
    ax1.set_ylabel('Flux for f140w (counts/s)', size=14)
    ax1.set_title('Phase = -10',size=14)
    ax1.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax2.plot(lightcurve_01['phase_array'], lightcurve_01[model_Ibc[i]][0], color=colors[i])
        ax2.errorbar(lightcurve_01[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax2.set_title('Phase = -5',size=14)
    ax2.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
        
    for i in range(len(model_Ibc)):
        ax3.plot(lightcurve_02['phase_array'], lightcurve_02[model_Ibc[i]][0], color=colors[i])
        ax3.errorbar(lightcurve_02[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax3.set_title('Phase = 0', size=14)
    ax3.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax4.plot(lightcurve_10['phase_array'], lightcurve_10[model_Ibc[i]][0], color=colors[i])
        ax4.errorbar(lightcurve_10[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax4.set_ylabel('Flux for f105w (counts/s)',size=14)
    ax4.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax4.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    
    for i in range(len(model_Ibc)):
        ax5.plot(lightcurve_11['phase_array'], lightcurve_11[model_Ibc[i]][0], color=colors[i])
        ax5.errorbar(lightcurve_11[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax5.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax5.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax6.plot(lightcurve_12['phase_array'], lightcurve_12[model_Ibc[i]][0], color=colors[i])
        ax6.errorbar(lightcurve_00[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax6.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax6.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax7.plot(lightcurve_20['phase_array'], lightcurve_20[model_Ibc[i]][0], color=colors[i])
        ax7.errorbar(lightcurve_20[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax7.set_ylabel('Flux for f814w.uvis1 (counts/s)',size=14)
    ax7.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax7.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax8.plot(lightcurve_21['phase_array'], lightcurve_21[model_Ibc[i]][0], color=colors[i])
        ax8.errorbar(lightcurve_21[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax8.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax8.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_Ibc)):
        ax9.plot(lightcurve_22['phase_array'], lightcurve_22[model_Ibc[i]][0], color=colors[i])
        ax9.errorbar(lightcurve_22[model_Ibc[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax9.annotate('Type Ibc',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax9.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    #plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    fig = plt.gcf()
    fig.set_size_inches(18.5,12.5)
    

def plot_II(z, data_flux, data_flux_err, phase1=-10, phase2=-5, phase3=0):
    model_II = ['s11-2004hx','s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p', 'nugent\
-sn2l']
    colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
 (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
 (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
 (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
 (0.8, 0.7254901960784313, 0.4549019607843137),
 (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]
    
    plt.figure(figsize = (9,9))
    gs1 = gridspec.GridSpec(9, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (0,1))
    ax3 = plt.subplot2grid((3,3), (0, 2))
    ax4 = plt.subplot2grid((3,3), (1, 0))
    ax5 = plt.subplot2grid((3,3), (1, 1))
    ax6 = plt.subplot2grid((3,3), (1, 2))
    ax7 = plt.subplot2grid((3,3), (2, 0))
    ax8 = plt.subplot2grid((3,3), (2, 1))
    ax9 = plt.subplot2grid((3,3), (2, 2))
    
    lightcurve_00 = lightcurve_II('f140w', z, phase1)
    lightcurve_01 = lightcurve_II('f140w', z, phase2)
    lightcurve_02 = lightcurve_II('f140w', z, phase3)
    
    lightcurve_10 = lightcurve_II('f105w', z, phase1)
    lightcurve_11 = lightcurve_II('f105w', z, phase2)
    lightcurve_12 = lightcurve_II('f105w', z, phase3)
    
    lightcurve_20 = lightcurve_II('f814w.uvis1', z, phase1)
    lightcurve_21 = lightcurve_II('f814w.uvis1', z, phase2)
    lightcurve_22 = lightcurve_II('f814w.uvis1', z, phase3)
    
    for i in range(len(model_II)):
        ax1.plot(lightcurve_00['phase_array'], lightcurve_00[model_II[i]][0], color=colors[i])
        ax1.errorbar(lightcurve_00[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
   # ax1.set_xlabel('x-label', fontsize=12)
    ax1.set_ylabel('Flux for f140w (counts/s)', size=14)
    ax1.set_title('Phase = -10',size=14)
    ax1.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax2.plot(lightcurve_01['phase_array'], lightcurve_01[model_II[i]][0], color=colors[i])
        ax2.errorbar(lightcurve_01[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax2.set_title('Phase = -5',size=14)
    ax2.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
        
    for i in range(len(model_II)):
        ax3.plot(lightcurve_02['phase_array'], lightcurve_02[model_II[i]][0], color=colors[i])
        ax3.errorbar(lightcurve_02[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax3.set_title('Phase = 0', size=14)
    ax3.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax4.plot(lightcurve_10['phase_array'], lightcurve_10[model_II[i]][0], color=colors[i])
        ax4.errorbar(lightcurve_10[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax4.set_ylabel('Flux for f105w (counts/s)',size=14)
    ax4.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax4.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    
    for i in range(len(model_II)):
        ax5.plot(lightcurve_11['phase_array'], lightcurve_11[model_II[i]][0], color=colors[i])
        ax5.errorbar(lightcurve_11[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax5.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax5.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax6.plot(lightcurve_12['phase_array'], lightcurve_12[model_II[i]][0], color=colors[i])
        ax6.errorbar(lightcurve_00[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax6.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax6.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax7.plot(lightcurve_20['phase_array'], lightcurve_20[model_II[i]][0], color=colors[i])
        ax7.errorbar(lightcurve_20[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax7.set_ylabel('Flux for f814w.uvis1 (counts/s)',size=14)
    ax7.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax7.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax8.plot(lightcurve_21['phase_array'], lightcurve_21[model_II[i]][0], color=colors[i])
        ax8.errorbar(lightcurve_21[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax8.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax8.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    for i in range(len(model_II)):
        ax9.plot(lightcurve_22['phase_array'], lightcurve_22[model_II[i]][0], color=colors[i])
        ax9.errorbar(lightcurve_22[model_II[i]][1], data_flux,xerr=0, yerr=data_flux_err, fmt='--o', c='k')
    ax9.annotate('Type II',xy=(0.75,0.85),xycoords='axes fraction', size=14)
    ax9.annotate('z=%.2f'%z,xy=(0.75,0.75),xycoords='axes fraction', size=14)
    
    #plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    fig = plt.gcf()
    fig.set_size_inches(18.5,12.5)
    

