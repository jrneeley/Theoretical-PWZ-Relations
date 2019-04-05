import PWZ
import numpy as np
import sys

bands = np.array(['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4',
    'W1', 'W2', 'W3'])
print_bands = np.array(['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', '3.6', '4.5', 'I3', 'I4',
    'W1', 'W2', 'W3'])
## Enter the coefficients of your reddening law here
#ext = np.array([1.348, 1.016, 0.845, 0.590, 0.289, 0.181, 0.118,
#        0.0661, 0.0507, 0.0507, 0.0507, 0.065, 0.052, 0.010])
ext = np.array([1.465, 1.332, 1.0, 0.828, 0.600, 0.288, 0.178, 0.117, 0.064,
    0.050, 0.050, 0.050, 0.065, 0.052, 0.010])
# 2 band PWZ relations
print 'Calculating 2 band PWZ relations...'

fu = open('pwz2_fu_coeff.dat', 'w')
fo = open('pwz2_fo_coeff.dat', 'w')
fofu = open('pwz2_fofu_coeff.dat', 'w')

for ind1 in np.arange(1,len(bands)):
    band1 = bands[ind1]
    pb1 = print_bands[ind1]
    ext1 = ext[ind1]

    for ind2 in np.arange(0,ind1):
        band2 = bands[ind2]
        pb2 = print_bands[ind2]
        ext2 = ext[ind2]

        rel = '{0},{1}-{0}'.format(band1, band2)
        # skip bands where where the extinction ratio is the same
        if ext[ind1] == ext[ind2]: continue
        # skip combinations that combine IRAC and WISE bands
        if (band1[0] == 'I') & (band2[0] == 'W'): continue
        if (band1[0] == 'W') & (band2[0] == 'I'): continue

        alpha = ext[ind1]/(ext[ind2] - ext[ind1])
        bands_in = [band1, band2, band1]
        ext_in = [ext1, ext2, ext1]
        c, e = PWZ.get_PWZ(bands_in, ext=ext_in, suppress_output=1)
        c_f, e_f = PWZ.get_PWZ(bands_in, ext=ext_in, fundamentalized=1, \
            suppress_output=1)
        fo.write('{:9} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n'.format(rel, \
            alpha, c[0], c[1], c[2], e[0], e[1], e[2], c[3]))
        fu.write('{:9} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n'.format(rel, \
            alpha, c[4], c[5], c[6], e[3], e[4], e[5], c[7]))
        fofu.write('{:9} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n'.format(rel, \
            alpha, c_f[0], c_f[1], c_f[2], e_f[0], e_f[1], e_f[2], c_f[3]))

        ind2 += -1

fu.close()
fo.close()
fofu.close()

# 3 band PWZ relations
print 'Calculating 3 band PWZ relations...'

fo = open('pwz3_fo_coeff.dat', 'w')
fu = open('pwz3_fu_coeff.dat', 'w')
fofu = open('pwz3_fofu_coeff.dat', 'w')

for ind2 in range(len(bands)-2):
    band2 = bands[ind2]

    for ind1 in np.arange(ind2+1, len(bands)-1):
        band1 = bands[ind1]

        for ind3 in np.arange(ind1+1, len(bands)):
            band3 = bands[ind3]

            rel = '{},{}-{}'.format(band1, band2, band3)
            #print rel
            test = np.array([ext[ind1], ext[ind2], ext[ind3]])
            # skip combinations where two filters have the same extinction
            if len(np.unique(test)) != 3: continue
            # ignore redundant band combinations
            if (band1[0] == 'I') & ((band2[0] == 'W') | (band3[0] == 'W')): continue
            if (band1[0] == 'W') & ((band2[0] == 'I') | (band3[0] == 'I')): continue

            rel = '{},{}-{}'.format(band1, band2, band3)
            bands_in = [band1, band2, band3]
            ext_in = [ext[ind1], ext[ind2], ext[ind3]]
            alpha = ext[ind1]/(ext[ind2] - ext[ind3])
            c, e = PWZ.get_PWZ(bands_in, ext=ext_in, suppress_output=1)
            c_f, e_f = PWZ.get_PWZ(bands_in, ext=ext_in, fundamentalized=1, \
                suppress_output=1)
            fo.write('{:9} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n'.format(rel, \
                alpha, c[0], c[1], c[2], e[0], e[1], e[2], c[3]))
            fu.write('{:9} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n'.format(rel, \
                alpha, c[4], c[5], c[6], e[3], e[4], e[5], c[7]))
            fofu.write('{:9} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n'.format(rel, \
                alpha, c_f[0], c_f[1], c_f[2], e_f[0], e_f[1], e_f[2], c_f[3]))

fu.close()
fo.close()
fofu.close()
