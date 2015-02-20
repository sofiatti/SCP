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


def read_data(filename):
    all_sources = pd.read_csv(filename)
    all_sources['Type'] = 'Type Ia'
    all_sources.ix[1000:2000, 'Type'] = 'Type Ibc'
    all_sources.ix[2000:, 'Type'] = 'Type II'
    all_features = copy.copy(all_sources)
    all_label = all_sources["Type"]
    del all_features["Type"]

    X = copy.copy(all_features.values)
    Y = copy.copy(all_label.values)

    return X, Y


def obtain_proba(my_flux, my_flux_diff):
    X, Y = read_data('data.csv')

    Y[Y == "Type Ia"] = 0
    Y[Y == "Type Ibc"] = 1
    Y[Y == "Type II"] = 2
    Y = Y.astype(int)

    clf = RandomForestClassifier(n_estimators=200, oob_score=True)
    rfc = clf.fit(X, Y)
    proba = rfc.predict_proba([my_flux, my_flux_diff])
    return(proba)

files = []
z = []
pdf = []
files.append(glob.glob('/Users/carolinesofiatti/projects/scp/'
             'flux_vs_fluxdiff/files_0-5/*f140w*'))
files.append(glob.glob('/Users/carolinesofiatti/projects/scp/'
             'flux_vs_fluxdiff/files_0-5/*f105w*'))

for a in files[0]:
    min_index = a.index('/z') + 2
    z.append(a[min_index:min_index + 3])
my_z = np.asarray(z, dtype=float)
my_z = np.divide(my_z, 100)

for i in range(len(files[0])):
    make_csv(files[0][i], files[1][i])
    pdf.append(obtain_proba(2, 2)[0])
np.savetxt('pdf.dat', np.transpose(pdf), header='Type Ia - Type'
           'Ib/c - Type II')
np.savetxt('z.dat', np.transpose(my_z))
