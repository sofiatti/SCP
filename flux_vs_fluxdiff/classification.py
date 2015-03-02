import glob
import pandas as pd
import copy
import gzip
import cPickle
import numpy as np
from sklearn.ensemble import RandomForestClassifier


def load(filename):
    """Loads a compressed object from disk"""
    file = gzip.GzipFile(filename, 'rb')
    object = cPickle.load(file)
    file.close()
    return object


def make_csv(filename_filter1, filename_filter2):
    dict_filter1 = load(filename_filter1)
    dict_filter2 = load(filename_filter2)

    type_Ia_flux_filter1 = dict_filter1['type_Ia_flux']
    type_Ibc_flux_filter1 = dict_filter1['type_Ibc_flux']
    type_II_flux_filter1 = dict_filter1['type_II_flux']

    type_Ia_flux_filter2 = dict_filter2['type_Ia_flux']
    type_Ibc_flux_filter2 = dict_filter2['type_Ibc_flux']
    type_II_flux_filter2 = dict_filter2['type_II_flux']

    type_Ia_flux_diff = np.subtract(type_Ia_flux_filter2, type_Ia_flux_filter1)
    type_Ibc_flux_diff = np.subtract(type_Ibc_flux_filter2,
                                     type_Ibc_flux_filter1)
    type_II_flux_diff = np.subtract(type_II_flux_filter2,
                                    type_II_flux_filter1)

    flux = np.hstack((type_Ia_flux_filter1, type_Ibc_flux_filter1,
                     type_II_flux_filter1))
    flux_diff = np.hstack((type_Ia_flux_diff, type_Ibc_flux_diff,
                          type_II_flux_diff))

    np.savetxt('data.csv', np.transpose([flux_diff, flux]), delimiter=',',
               header='Type Ia Flux Diff,Type Ia Flux')


# Read in data generated using sncosmo and assign the correct SN Type to each
# instance


def read_data(files_dir, filename):
    indices = [i for i, x in enumerate(files_dir) if x == '_']
    n = files_dir[(indices[1] + 1):indices[2]]
    n_Ibc_start = int(n)
    n_II_start = 2 * n_Ibc_start

    all_sources = pd.read_csv(filename)
    all_sources['Type'] = 'Type Ia'
    all_sources.ix[n_Ibc_start:n_II_start, 'Type'] = 'Type Ibc'
    all_sources.ix[n_II_start:, 'Type'] = 'Type II'
    all_features = copy.copy(all_sources)
    all_label = all_sources["Type"]
    del all_features["Type"]

    X = copy.copy(all_features.values)
    Y = copy.copy(all_label.values)

    return X, Y


def obtain_proba(files_dir, flux_filter1, flux_filter2):
    my_flux_diff = flux_filter2 - flux_filter1
    my_flux = flux_filter1

    X, Y = read_data(files_dir, 'data.csv')

    Y[Y == "Type Ia"] = 0
    Y[Y == "Type Ibc"] = 1
    Y[Y == "Type II"] = 2
    Y = Y.astype(int)

    clf = RandomForestClassifier(n_estimators=200, oob_score=True)
    rfc = clf.fit(X, Y)
    proba = rfc.predict_proba([my_flux, my_flux_diff])
    return(proba)


def make_pdf_z_file(my_dir, files_dir, filter1, filter2, flux_filter1,
                    flux_filter2):
    files = []
    z = []
    pdf = []
    files.append(glob.glob(my_dir + files_dir + '*' + filter1 + '*.gz'))
    files.append(glob.glob(my_dir + files_dir + '*' + filter2 + '*.gz'))

    for a in files[0]:
        min_index = a.index('/z') + 2
        z.append(a[min_index:min_index + 3])
    my_z = np.asarray(z, dtype=float)
    my_z = np.divide(my_z, 100)

    for i in range(len(files[0])):
        make_csv(files[0][i], files[1][i])
        pdf.append(obtain_proba(files_dir, flux_filter1, flux_filter2)[0])
    np.savetxt(my_dir + files_dir + 'pdf_' + filter1 + '_' + filter2 + '.dat',
               np.transpose(pdf), header='Type Ia - Type Ib/c - Type II')
    np.savetxt(my_dir + files_dir + 'z_' + filter1 + '_' + filter2 +
               '.dat', np.transpose(my_z))

my_dir = '/Users/carolinesofiatti/projects/scp/flux_vs_fluxdiff/'
all_files_dir = ['mc_files_1000_-10_-5/', 'mc_files_1000_-5_0/', 'mc_files_1000_0_5/']

all_filters = ['f105w', 'f140w', 'uvf814w']
all_fluxes = [6.506, 4.854, 1.2]
data = dict(zip(all_filters, all_fluxes))

j = 0
for files_dir in all_files_dir:
    print files_dir
    make_pdf_z_file(my_dir, files_dir, all_filters[0], all_filters[1],
                    all_fluxes[0], all_fluxes[1])
    make_pdf_z_file(my_dir, files_dir, all_filters[0], all_filters[2],
                    all_fluxes[0], all_fluxes[2])
    make_pdf_z_file(my_dir, files_dir, all_filters[1], all_filters[2],
                    all_fluxes[1], all_fluxes[2])
    print 'done'
           