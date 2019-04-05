import PWZ
import numpy as np

bands = ['U', 'B', 'V', 'R','I','J','H','K','I1','I2','I3','I4','W1','W2','W3','W4']
plz_a = np.zeros(3*len(bands))
plz_b = np.zeros(3*len(bands))
plz_c = np.zeros(3*len(bands))
plz_ae = np.zeros(3*len(bands))
plz_be = np.zeros(3*len(bands))
plz_ce = np.zeros(3*len(bands))
plz_sig = np.zeros(3*len(bands))
plz_bands = np.zeros(3*len(bands), dtype='S2')
rel_type = np.zeros(3*len(bands), dtype='S5')

for ind, band in enumerate(bands):

    c, e = PWZ.get_PLZ(band, suppress_output=1)
    c_fun, e_fun = PWZ.get_PLZ(band, fundamentalized=1, suppress_output=1)

    inds = [3*ind, 3*ind+1, 3*ind+2]
    plz_a[inds] = [c[0], c[4], c_fun[0]]
    plz_b[inds] = [c[1], c[5], c_fun[1]]
    plz_c[inds] = [c[2], c[6], c_fun[2]]
    plz_ae[inds] = [e[0], e[3], e_fun[0]]
    plz_be[inds] = [e[1], e[4], e_fun[1]]
    plz_ce[inds] = [e[2], e[5], e_fun[2]]
    plz_sig[inds] = [c[3], c[7], c_fun[3]]
    plz_bands[inds] = [band, band, band]
    rel_type[inds] = ['FO', 'FU', 'FO+FU']

data = np.array(zip(rel_type, plz_bands, plz_a, plz_b, plz_c, plz_ae, plz_be, plz_ce,
    plz_sig), dtype=[('type', 'S5'), ('bands', 'S2'), ('a', float), ('b', float), ('c', float),
    ('ae', float), ('be', float), ('ce', float), ('sig', float)])
np.savetxt('PLZ.dat', data, fmt=['%5s', '%2s']+['%6.3f', '%6.3f', '%6.3f', '%6.3f',
    '%6.3f', '%6.3f', '%6.3f'])
