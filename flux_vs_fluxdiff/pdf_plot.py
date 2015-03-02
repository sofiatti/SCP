import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc


def load_data(my_dir, files_dir, filter1, filter2):
    """Gets z and PDF

    Parameters
    ----------
    plot_type: str
        Name of plot type (combined or regular), optional, regular default.
    """
    z = np.loadtxt(my_dir + files_dir + 'z_' + filter1 + '_' + filter2 +
                   '.dat')

    pdf_Ia, pdf_Ibc, pdf_II = np.loadtxt(my_dir + files_dir + 'pdf_' +
                                         filter1 + '_' + filter2 + '.dat',
                                         skiprows=1)
    return (pdf_Ia, pdf_Ibc, pdf_II, z)


def plot(my_dir, files_dir, pdf_Ia, pdf_Ibc, pdf_II, z, my_title,
         my_ylabel, file_name):
    sns.set_context("talk", font_scale=1.4)
    # sns.set_context("talk", rc={"lines.linewidth": 2})
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    # for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    print np.shape(z)
    print np.shape(pdf_Ia)
    plt.plot(z, pdf_Ia, label='Type Ia')
    plt.plot(z, pdf_Ibc, label='Type Ib/c')
    plt.plot(z, pdf_II, label='Type II')

    plt.xlabel('z')
    plt.title(my_title)
    plt.ylabel(my_ylabel)
    plt.legend()
    print my_dir
    print files_dir
    print file_name
    plt.savefig(my_dir + files_dir + file_name)
    plt.close()


def make_single_plot(my_dir, files_dir, filter1, filter2):
    my_ylabel = 'PDF for (' + filter2 + ' - ' + filter1 + ') vs. ' + filter1
    indices = [i for i, x in enumerate(files_dir) if x == '_']
    phase_min = files_dir[(indices[-2] + 1):indices[-1]]
    phase_max = files_dir[(indices[-1] + 1): -1]
    n = files_dir[(indices[1] + 1):indices[2]]
    my_title = (n + ' simulated points \& Type Ia phase = (' + phase_min + ', ' +
                phase_max + ')')
    file_name = (n + '_' + phase_min + '_' + phase_max + '_' + filter1 + '_' +
                 filter2 + '_pdf_vs_z.png')
    pdf_Ia, pdf_Ibc, pdf_II, z = load_data(my_dir, files_dir, filter1,
                                           filter2)

    plot(my_dir, files_dir, pdf_Ia, pdf_Ibc, pdf_II, z, my_title, my_ylabel,
         file_name)


def make_combined_plot(my_dir, files_dir, filter1, filter2, filter3):
    my_ylabel = 'Multiplied PDFs of all filters'
    indices = [i for i, x in enumerate(files_dir) if x == '_']
    phase_min = files_dir[(indices[-2] + 1):indices[-1]]
    phase_max = files_dir[(indices[-1] + 1):-1]
    n = files_dir[(indices[1] + 1):indices[2]]
    print n
    my_title = (n + ' simulated points \& Type Ia phase = (' + phase_min + ', ' +
                phase_max + ')')
    file_name = (n + '_' + phase_min + '_' + phase_max + '_' + filter1 + '_' +
                 filter2 + '_' + filter3 + '_pdf_vs_z.png')
    pdf_Ia_I, pdf_Ibc_I, pdf_II_I, z_I = load_data(my_dir, files_dir, filter1,
                                                   filter2)
    pdf_Ia_II, pdf_Ibc_II, pdf_II_II, z_II = load_data(my_dir, files_dir,
                                                       filter1, filter3)
    pdf_Ia_III, pdf_Ibc_III, pdf_II_III, z_III = load_data(my_dir, files_dir,
                                                           filter2, filter3)
    z = z_I

    pdf_Ia = pdf_Ia_I * pdf_Ia_II * pdf_Ia_III
    pdf_Ibc = pdf_Ibc_I * pdf_Ibc_II * pdf_Ibc_III
    pdf_II = pdf_II_I * pdf_II_II * pdf_II_III

    plot(my_dir, files_dir, pdf_Ia, pdf_Ibc, pdf_II, z, my_title, my_ylabel,
         file_name)

my_dir = '/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/'
all_files_dir = ['mc_files_1000_-10_-5/', 'mc_files_1000_-5_0/',
                 'mc_files_1000_0_5/']
all_filters = ['f105w', 'f140w', 'uvf814w']

for files_dir in all_files_dir:
    make_combined_plot(my_dir, files_dir, all_filters[0], all_filters[1],
                       all_filters[2])
    make_single_plot(my_dir, files_dir, all_filters[0], all_filters[1])
    make_single_plot(my_dir, files_dir, all_filters[0], all_filters[2])
    make_single_plot(my_dir, files_dir, all_filters[1], all_filters[2])