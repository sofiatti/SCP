import gzip
import cPickle
import mpld3
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pylab
from matplotlib import cm


def load(filename):
    """Loads a compressed object from disk"""
    file = gzip.GzipFile(filename, 'rb')
    object = cPickle.load(file)
    file.close()
    return object


def plot(filename_filter1, filename_filter2, point_flux_filter1, point_flux_filter1_err, point_flux_filter2, point_flux_filter2_err):
    """Generates a plot and transforms the figure into an HTML file from the
    montecarlo data. Adds a point to the data"""
    point_flux_diff = point_flux_filter2 - point_flux_filter1
    point_flux_diff_err = point_flux_filter2_err + point_flux_filter1_err

    dict_filter1 = load(filename_filter1)
    dict_filter2 = load(filename_filter2)
    
    type_Ia_flux_filter1 = dict_filter1['type_Ia_flux']
    type_Ibc_flux_filter1 = dict_filter1['type_Ibc_flux']
    type_II_flux_filter1 = dict_filter1['type_II_flux']
    filter1 = dict_filter1['filter']

    type_Ia_flux_filter2 = dict_filter2['type_Ia_flux']
    type_Ibc_flux_filter2 = dict_filter2['type_Ibc_flux']
    type_II_flux_filter2 = dict_filter2['type_II_flux']
    filter2 = dict_filter2['filter']

    type_Ia_flux_diff = np.subtract(type_Ia_flux_filter2,type_Ia_flux_filter1)
    type_Ibc_flux_diff = np.subtract(type_Ibc_flux_filter2,type_Ibc_flux_filter1)
    type_II_flux_diff = np.subtract(type_II_flux_filter2,type_II_flux_filter1)

    sns.set(style="white", palette="muted")
    sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})
    g = sns.JointGrid(filter2.upper() + " - " + filter1.upper() + " flux difference",
                      filter1.upper() + " flux", space=0)

    sns.kdeplot(type_Ia_flux_diff, ax=g.ax_marg_x, legend=True, label='Type Ia')
    sns.kdeplot(type_Ia_flux_filter1, ax=g.ax_marg_y, vertical=True, legend=False)
    sns.kdeplot(type_Ia_flux_diff, type_Ia_flux_filter1, cut=5, ax=g.ax_joint,
                alpha=0.8, cmap='Blues_d')

    sns.kdeplot(type_Ibc_flux_diff, ax=g.ax_marg_x, legend=True, label='Type Ib/c')
    sns.kdeplot(type_Ibc_flux_filter1, ax=g.ax_marg_y, vertical=True, legend=False)
    sns.kdeplot(type_Ibc_flux_diff, type_Ibc_flux_filter1, cut=5, ax=g.ax_joint,
                alpha=0.8, cmap='Greens_d')

    sns.kdeplot(type_II_flux_diff, ax=g.ax_marg_x, legend=True, label='Type II')
    sns.kdeplot(type_II_flux_filter1, ax=g.ax_marg_y, vertical=True, legend=False)
    sns.kdeplot(type_II_flux_diff, type_II_flux_filter1, cut=5, ax=g.ax_joint,
                alpha=0.8, cmap='Reds_d')

    g.ax_joint.errorbar(point_flux_diff, point_flux_filter1, xerr=point_flux_diff_err, yerr=point_flux_filter1_err, fmt='--o', c='k')
    g.set_axis_labels(xlabel=filter2.upper() + " - " + filter1.upper() + " flux difference (counts/s)",
                      ylabel=filter1.upper() + " flux (counts/s)")
    g.ax_joint.set_xlim(-10, 10)
    g.ax_joint.set_ylim(-1, 10)
    plt.tight_layout()
    plt.savefig(filename_filter1[:-5] + filename_filter2[5:-5]+ 'plot.png')
    mpld3.save_html(plt.gcf(), filename_filter1[:-5] + filename_filter2[5:-5] + 'plot.html')
    return ()
