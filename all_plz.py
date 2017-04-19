import PWZ
import numpy as np

bands = ['U','B','V','R','I','J','H','K','I1','I2','I3','I4','W1','W2','W3','W4']
data = np.zeros((len(bands), 8), dtype=float)


for ind, band in enumerate(bands):
    c = PWZ.get_PLZ(band, suppress_output=1)
    data[ind,:] = c

np.savetxt('PLZ.dat',np.c_[bands, data], fmt='%s')
