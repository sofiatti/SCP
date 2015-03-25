import numpy as np
import sncosmo
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

t0 = 0
hostr_v = 3.1
dust = sncosmo.CCM89Dust()

zero_point = {'f105w_wfc3ir': 26.235, 'f140w_wfc3ir': 26.437, 'f160w': 25.921,
              'f814w_wfc3uvis': 25.0985, 'zpsys': 'ab'}


def filter2bandpass(filter_file):
    """Returns the sncosmo bandpass for an HST filter"""
    filter = np.loadtxt(filter_file)
    wavelength = filter[:, 0]
    transmission = filter[:, 1]
    band = sncosmo.Bandpass(wavelength, transmission, name=filter_file[: -4])
    sncosmo.registry.register(band, force=True)
    return

# Only have to use once
f105w_wfc3ir = filter2bandpass('f105w_wfc3ir.dat')
f140w_wfc3ir = filter2bandpass('f140w_wfc3ir.dat')
f814w_wfc3uvis = filter2bandpass('f814w_wfc3uvis.dat')


def lightcurve_Ia(filter, z, x1, c):
    """Given a filter and redshift z, generates the observed
    flux for SNe Type Ia"""

    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ia = sncosmo.Model(source=sncosmo.get_source('salt2',
                                                       version='2.4'))

    alpha = 0.12
    beta = 3.
    mabs = -19.1 - alpha*x1 + beta*c
    model_Ia.set(z=z)
    model_Ia.set_source_peakabsmag(mabs, 'bessellb', 'vega')
    p = {'z': z, 't0': t0, 'x1': x1, 'c': c}
    model_Ia.set(**p)
    phase_array = np.linspace(model_Ia.mintime(), model_Ia.maxtime(), 100)
    obsflux_Ia = model_Ia.bandflux(filter, phase_array, zp=zp, zpsys=zpsys)
    keys = ['phase_array', 'obsflux']
    values = [phase_array, obsflux_Ia]
    dict_Ia = dict(zip(keys, values))
    np.savetxt('test.dat', np.c_[dict_Ia['phase_array'], dict_Ia['obsflux']])
    x0 = model_Ia.get('x0')
    print x0
    return (dict_Ia)


def lightcurve_Ibc(filter, z, hostebv_Ibc):
    """Given a filter, redshift z at given phase, generates the observed
    magnitude for SNe Type Ibc"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc',
                 'nugent-hyper', 's11-2006jo', 'snana-2004fe',
                 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo',
                 'snana-sdss014475', 'snana-2006lc', 'snana-04d1la',
                 'snana-04d4jv', 'snana-2004gv', 'snana-2006ep',
                 'snana-2007y', 'snana-2004ib', 'snana-2005hm',
                 'snana-2006jo', 'snana-2007nc']
    obsflux_Ibc = []
    phase_arrays = []
    for i in model_Ibc:
        model_i = sncosmo.Model(source=sncosmo.get_source(i), effects=[dust],
                                effect_names=['host'], effect_frames=['rest'])
        mabs = -17.56
        model_i.set(z=z)
        phase_array_i = np.linspace(model_i.mintime(), model_i.maxtime(), 100)
        model_i.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        p_core_collapse = {'z': z, 't0': t0, 'hostebv': hostebv_Ibc,
                           'hostr_v': hostr_v}
        model_i.set(**p_core_collapse)
        phase_arrays.append(phase_array_i)
        obsflux_i = model_i.bandflux(filter, phase_array_i, zp, zpsys)
        obsflux_Ibc.append(obsflux_i)
    keys = model_Ibc
    values = []
    for i, item in enumerate(model_Ibc):
        values.append([obsflux_Ibc[i], phase_arrays[i]])
    dict_Ibc = dict(zip(keys, values))
    return (dict_Ibc)


def lightcurve_II(filter, z, hostebv_II):
    """Given a filter and redshift z, generates the observed magnitude for 
    SNe Type II"""
    zp = zero_point[filter]
    zpsys = zero_point['zpsys']
    model_II = ['s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p',
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
    obsflux_II = []
    phase_arrays = []
    for i in model_II:
        model_i = sncosmo.Model(source=sncosmo.get_source(i), effects=[dust],
                                effect_names=['host'],
                                effect_frames=['rest'])
        if i == 's11-2004hx' == 'nugent-sn2l':
            mabs = -17.98
        else:
            mabs = -16.75
        model_i.set(z=z)
        phase_array_i = np.linspace(model_i.mintime(), model_i.maxtime(), 100)
        model_i.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        p_core_collapse = {'z': z, 't0': t0, 'hostebv': hostebv_II,
                           'hostr_v': hostr_v}
        model_i.set(**p_core_collapse)
        phase_arrays.append(phase_array_i)
        obsflux_i = model_i.bandflux(filter, phase_array_i, zp, zpsys)
        obsflux_II.append(obsflux_i)
    keys = model_II
    values = []
    for i, item in enumerate(model_II):
        values.append([obsflux_II[i], phase_arrays[i]])
    dict_II = dict(zip(keys, values))
    return (dict_II)


def plot_Ia(z, x1, c, data_flux_filter1, data_flux_filter1_err,
            data_flux_filter2, data_flux_filter2_err,
            data_flux_filter3, data_flux_filter3_err, phase):
    plt.figure(figsize=(3, 9))
    gs1 = gridspec.GridSpec(3, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3, 3), (0, 0))
    ax2 = plt.subplot2grid((3, 3), (0, 1))
    ax3 = plt.subplot2grid((3, 3), (0, 2))

    lightcurve_00 = lightcurve_Ia('f140w', z, x1, c)
    lightcurve_10 = lightcurve_Ia('f105w', z, x1, c)
    lightcurve_20 = lightcurve_Ia('uvf814w', z, x1, c)

    ax1.plot(lightcurve_00['phase_array'], lightcurve_00['obsflux'])
    ax1.errorbar(phase, data_flux_filter1, xerr=0, yerr=data_flux_filter1_err,
                 fmt='--o', c='k')
    ax1.set_ylabel('Flux (counts/s)', size=14)
    ax1.set_title('F140W', size=14)
    ax1.annotate('Type Ia', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    ax2.plot(lightcurve_10['phase_array'], lightcurve_10['obsflux'])
    ax2.errorbar(phase, data_flux_filter2, xerr=0, yerr=data_flux_filter2_err,
                 fmt='--o', c='k')
    ax2.set_title('F105W', size=14)
    ax2.annotate('Type Ia', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    ax3.plot(lightcurve_20['phase_array'], lightcurve_20['obsflux'])
    ax3.errorbar(phase, data_flux_filter3, xerr=0, yerr=data_flux_filter3_err,
                 fmt='--o', c='k')
    ax3.set_title('F814W.UVIS', size=14)
    ax3.annotate('Type Ia', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    fig = plt.gcf()
    fig.set_size_inches(18.5, 12.5)


def plot_Ibc(z, hostebv_Ibc, data_flux_filter1, data_flux_filter1_err,
             data_flux_filter2, data_flux_filter2_err,
             data_flux_filter3, data_flux_filter3_err, phase):
    model_Ibc = ['s11-2005hl', 's11-2005hm', 's11-2006fo', 'nugent-sn1bc',
                 'nugent-hyper', 's11-2006jo', 'snana-2004fe',
                 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo',
                 'snana-sdss014475', 'snana-2006lc', 'snana-04d1la',
                 'snana-04d4jv', 'snana-2004gv', 'snana-2006ep',
                 'snana-2007y', 'snana-2004ib', 'snana-2005hm',
                 'snana-2006jo', 'snana-2007nc']
    #colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
    #          (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
    #          (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
    #          (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
    #         (0.8, 0.7254901960784313, 0.4549019607843137),
    #          (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]

    plt.figure(figsize=(3, 9))
    gs1 = gridspec.GridSpec(3, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3, 3), (0, 0))
    ax2 = plt.subplot2grid((3, 3), (0, 1))
    ax3 = plt.subplot2grid((3, 3), (0, 2))

    lightcurve_00 = lightcurve_Ibc('f140w', z, hostebv_Ibc)
    lightcurve_10 = lightcurve_Ibc('f105w', z, hostebv_Ibc)
    lightcurve_20 = lightcurve_Ibc('uvf814w', z, hostebv_Ibc)

    for i in range(len(model_Ibc)):
        ax1.plot(lightcurve_00[model_Ibc[i]][1], lightcurve_00[model_Ibc[i]][0])
        ax1.errorbar(phase, data_flux_filter1, xerr=0,
                     yerr=data_flux_filter1_err, fmt='--o', c='k')

    ax1.set_ylabel('Flux (counts/s)', size=14)
    ax1.set_title('F140W', size=14)
    ax1.annotate('Type Ibc', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    for i in range(len(model_Ibc)):
        ax2.plot(lightcurve_10[model_Ibc[i]][1], lightcurve_10[model_Ibc[i]][0])
        ax2.errorbar(phase, data_flux_filter2, xerr=0,
                     yerr=data_flux_filter2_err, fmt='--o',
                     c='k')
    ax2.set_title('F105W', size=14)
    ax2.annotate('Type Ibc', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    for i in range(len(model_Ibc)):
        ax3.plot(lightcurve_20[model_Ibc[i]][1], lightcurve_20[model_Ibc[i]][0])
        ax3.errorbar(phase, data_flux_filter3, xerr=0,
                     yerr=data_flux_filter3_err, fmt='--o',
                     c='k')
    ax3.set_title('F184.UVIS', size=14)
    ax3.annotate('Type Ibc', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    fig = plt.gcf()
    fig.set_size_inches(18.5, 12.5)


def plot_II(z, hostebv_II, data_flux_filter1, data_flux_filter1_err,
            data_flux_filter2, data_flux_filter2_err,
            data_flux_filter3, data_flux_filter3_err, phase):
    model_II = ['s11-2005lc', 's11-2005gi', 's11-2006jl', 'nugent-sn2p',
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
    #colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
    #          (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
    #          (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
    #         (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
    #          (0.8, 0.7254901960784313, 0.4549019607843137),
    #          (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]

    plt.figure(figsize=(3, 9))
    gs1 = gridspec.GridSpec(3, 9)
    gs1.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot2grid((3, 3), (0, 0))
    ax2 = plt.subplot2grid((3, 3), (0, 1))
    ax3 = plt.subplot2grid((3, 3), (0, 2))

    lightcurve_00 = lightcurve_II('f140w', z, hostebv_II)
    lightcurve_10 = lightcurve_II('f105w', z, hostebv_II)
    lightcurve_20 = lightcurve_II('uvf814w', z, hostebv_II)

    for i in range(len(model_II)):
        ax1.plot(lightcurve_00[model_II[i]][1], lightcurve_00[model_II[i]][0])
        ax1.errorbar(phase, data_flux_filter1, xerr=0,
                     yerr=data_flux_filter1_err, fmt='--o',
                     c='k')

    ax1.set_ylabel('Flux (counts/s)', size=14)
    ax1.set_title('F140W', size=14)
    ax1.annotate('Type II', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax1.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    for i in range(len(model_II)):
        ax2.plot(lightcurve_10[model_II[i]][1], lightcurve_10[model_II[i]][0])
        ax2.errorbar(phase, data_flux_filter2, xerr=0,
                     yerr=data_flux_filter2_err, fmt='--o',
                     c='k')
    ax2.set_title('F105W', size=14)
    ax2.annotate('Type II', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax2.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    for i in range(len(model_II)):
        ax3.plot(lightcurve_20[model_II[i]][1], lightcurve_20[model_II[i]][0])
        ax3.errorbar(phase, data_flux_filter3, xerr=0,
                     yerr=data_flux_filter3_err, fmt='--o',
                     c='k')
    ax3.set_title('F184W.UVIS', size=14)
    ax3.annotate('Type II', xy=(0.75, 0.85), xycoords='axes fraction', size=14)
    ax3.annotate('z=%.2f' % z, xy=(0.75, 0.75), xycoords='axes fraction',
                 size=14)

    fig = plt.gcf()
    fig.set_size_inches(18.5, 12.5)
