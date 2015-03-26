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


def plot(filename_filter1, filename_filter2, point_flux_filter1,
         point_flux_filter1_err, point_flux_filter2, point_flux_filter2_err):
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

    type_Ia_flux_diff = np.subtract(type_Ia_flux_filter2, type_Ia_flux_filter1)
    type_Ibc_flux_diff = np.subtract(type_Ibc_flux_filter2,
                                     type_Ibc_flux_filter1)
    type_II_flux_diff = np.subtract(type_II_flux_filter2, type_II_flux_filter1)

    flux = [type_Ia_flux_filter1, type_Ibc_flux_filter1, type_II_flux_filter1]
    diff = [type_Ia_flux_diff, type_Ibc_flux_diff, type_II_flux_diff]

    for i, item in enumerate(flux):
        outlier_mask = (flux[i] < np.percentile(flux[i], 95)) & (diff[i] < np.percentile(diff[i], 95))
        flux[i] = flux[i][outlier_mask]
        diff[i] = diff[i][outlier_mask]

    sns.set(style="white", palette="muted")
    sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})
    g = sns.JointGrid(filter2.upper() + " - " + filter1.upper() + " flux difference",
                      filter1.upper() + " flux", space=0)

    sns.kdeplot(diff[0], ax=g.ax_marg_x, legend=True, label='Type Ia')
    sns.kdeplot(flux[0], ax=g.ax_marg_y, vertical=True, legend=False)
    sns.kdeplot(diff[0], flux[0], cut=4, ax=g.ax_joint,
                alpha=0.8, cmap='Blues_d')

    sns.kdeplot(diff[1], ax=g.ax_marg_x, legend=True, label='Type Ib/c')
    sns.kdeplot(flux[1], ax=g.ax_marg_y, vertical=True, legend=False)
    sns.kdeplot(diff[1], flux[1], cut=4, ax=g.ax_joint,
                alpha=0.8, cmap='Greens_d')

    sns.kdeplot(diff[2], ax=g.ax_marg_x, legend=True, label='Type II')
    sns.kdeplot(flux[2], ax=g.ax_marg_y, vertical=True, legend=False)
    sns.kdeplot(diff[2], flux[2], cut=4, ax=g.ax_joint,
                alpha=0.8, cmap='Reds_d')

    g.ax_joint.errorbar(point_flux_diff, point_flux_filter1, xerr=point_flux_diff_err, yerr=point_flux_filter1_err, fmt='--o', c='k')
    g.set_axis_labels(xlabel=filter2.upper() + " - " + filter1.upper() + " flux difference (counts/s)",
                      ylabel=filter1.upper() + " flux (counts/s)")
    plt.tight_layout()
    plt.show()
    '''
    # plt.savefig(filename_filter1[:-5] + filename_filter2[5:-5]+ 'plot.png')
    # mpld3.save_html(plt.gcf(), filename_filter1[:-5] + filename_filter2[5:-5] + 'plot.html')
    '''


def scatter_plot(filename_filter1, filename_filter2):
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

    type_Ia_flux_diff = np.subtract(type_Ia_flux_filter2, type_Ia_flux_filter1)
    type_Ibc_flux_diff = np.subtract(type_Ibc_flux_filter2,
                                     type_Ibc_flux_filter1)
    type_II_flux_diff = np.subtract(type_II_flux_filter2, type_II_flux_filter1)

    flux = [type_Ia_flux_filter1, type_Ibc_flux_filter1, type_II_flux_filter1]
    diff = [type_Ia_flux_diff, type_Ibc_flux_diff, type_II_flux_diff]

    for i, item in enumerate(flux):
        outlier_mask = (flux[i] < np.percentile(flux[i], 95)) & (diff[i] < np.percentile(diff[i], 95))
        flux[i] = flux[i][outlier_mask]
        diff[i] = diff[i][outlier_mask]

    plt.scatter(diff[0], flux[0], c='b',
                label='Type Ia')
    plt.scatter(diff[1], flux[1], c='g',
                label='Type Ib/c')
    plt.scatter(diff[2], flux[2], c='r',
                label='Type II')
    plt.title(filter2 + ' - ' + filter1 + ' vs. ' + filter1)
    plt.legend()
    plt.show()
    plt.close()
