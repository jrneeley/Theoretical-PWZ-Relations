import numpy as np
import matplotlib.pyplot as mp
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg


def get_PWZ(bands, ext=0, fundamentalized=0):

    bands_all = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4',
        'W1', 'W2', 'W3']

    band1 = bands[0]
    band2 = bands[1]
    band3 = bands[2]
    col1 = bands_all.index(band1)+6
    col2 = bands_all.index(band2)+6
    col3 = bands_all.index(band3)+6
# Default reddening law Cardelli et al (1989) + Indebetouw et al. (2005)
    default_ext = [1.579, 1.348, 1.016, 0.845, 0.590, 0.289, 0.181, 0.118,
        0.0661, 0.0507, 0.0507, 0.0507, 0.065, 0.052, 0.010]
    if ext == 0 :
        ext = [default_ext[col1-6], default_ext[col2-6], default_ext[col3-6]]

    # Check for two band or three band relation
    if band1 == band3:
        band3 = band3+'.'
        # Read in theoretical mean mags for two band relation
        dtype1 = np.dtype([('Z', float), ('Y', float), ('logP', float), (band1, float),
            (band2, float), (band3, float)])
        RRab = np.loadtxt('FU-means.txt', dtype=dtype1, usecols=(0,1,5,col1,col2,col2+1), skiprows=33)
        RRc = np.loadtxt('FO-means.txt', dtype=dtype1, usecols=(0,1,5,col1,col2,col2+1), skiprows=33)
        RRab[band3] = RRab[band1]
        RRc[band3] = RRc[band1]
    else:
        # Read in theoretical mean mags for three band relation
        dtype1 = np.dtype([('Z', float), ('Y', float), ('logP', float), (band1, float),
            (band2, float), (band3, float)])
        RRab = np.loadtxt('FU-means.txt', dtype=dtype1, usecols=(0,1,5,col1,col2,col3), skiprows=33)
        RRc = np.loadtxt('FO-means.txt', dtype=dtype1, usecols=(0,1,5,col1,col2,col3), skiprows=33)

# Define reddening coefficient
    alpha = ext[0]/(ext[1] - ext[2])


    if fundamentalized == 1:
        RRc['logP'] += 0.127
        period = np.append(RRc['logP'], RRab['logP'])
        Z = np.append(RRc['Z'], RRab['Z'])
        Y = np.append(RRc['Y'], RRab['Y'])
        FeH = np.log10(Z/(1-Z-Y))+1.61 - np.log10(0.694*10**0.4+0.306)
        FeH_c = np.log10(RRc['Z']/(1-RRc['Z']-RRc['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)
        FeH_ab = np.log10(RRab['Z']/(1-RRab['Z']-RRab['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)

        m1 = np.append(RRc[band1], RRab[band1])
        m2 = np.append(RRc[band2], RRab[band2])
        m3 = np.append(RRc[band3], RRab[band3])

        W = m1 - alpha* (m2 - m3)
        W_c = RRc[band1] - alpha * (RRc[band2] - RRc[band3])
        W_ab = RRab[band1] - alpha * (RRab[band2] - RRab[band3])

        XX,YY = np.meshgrid(np.arange(-2.5, 0.5, 0.3), np.arange(-0.45, 0.25, 0.07))

    # best-fit linear plane
        A = np.c_[FeH, period, np.ones(period.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, W)    # coefficients

    # evaluate it on grid
        ZZ = C[0]*XX + C[1]*YY + C[2]
        print "\nReddening coefficient = {:0.3f}".format(alpha)
        print '\nW('+band1+','+band2+'-'+band3+') = a + b*logP + c*[Fe/H]\n'
        print 'Fundamentalized Coefficients: '
        print 'a = {:0.3f}'.format(C[2])
        print 'b = {:0.3f}'.format(C[1])
        print 'c = {:0.3f}'.format(C[0])

        fig = mp.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(FeH_c, RRc['logP'], W_c, c='b', marker='o')
        ax.scatter(FeH_ab, RRab['logP'], W_ab, c='r', marker='o')
        mp.xlabel('[Fe/H]')
        mp.ylabel('logP')
        ax.set_zlabel('W('+band1+','+band2+'-'+band3+')')
        #mp.plot(data['logP'], data['i1'], 'ro')
        mp.show()
    if fundamentalized == 0:

        FeH_c = np.log10(RRc['Z']/(1-RRc['Z']-RRc['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)
        FeH_ab = np.log10(RRab['Z']/(1-RRab['Z']-RRab['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)

        W_c = RRc[band1] - alpha* (RRc[band2] - RRc[band3])
        W_ab = RRab[band1] - alpha* (RRab[band2] - RRab[band3])

        X_c,Y_c = np.meshgrid(np.arange(-2.5, 0.5, 0.3), np.arange(-0.6, -0.2, 0.04))
        X_ab,Y_ab = np.meshgrid(np.arange(-2.5, 0.5, 0.3), np.arange(-0.45, 0.25, 0.07))

    # best-fit linear plane
        A_c = np.c_[FeH_c, RRc['logP'], np.ones(RRc.shape[0])]
        C_c,_,_,_ = scipy.linalg.lstsq(A_c, W_c)    # coefficients
        A_ab = np.c_[FeH_ab, RRab['logP'], np.ones(RRab.shape[0])]
        C_ab,_,_,_ = scipy.linalg.lstsq(A_ab, W_ab)    # coefficients


    # evaluate it on grid
        Z_c = C_c[0]*X_c + C_c[1]*Y_c + C_c[2]
        Z_ab = C_ab[0]*X_ab + C_ab[1]*Y_ab + C_ab[2]
        #print 'Reddening coefficient = ', alpha
        print "\nReddening coefficient = {:0.3f}".format(alpha)
        print '\nW('+band1+','+band2+'-'+band3+') = a + b*logP + c*[Fe/H]\n'
        print 'Fundamental Coefficients: '
        print 'a = {:0.3f}'.format(C_ab[2])
        print 'b = {:0.3f}'.format(C_ab[1])
        print 'c = {:0.3f}'.format(C_ab[0])
        print '\nFirst Overtone Coefficients: '
        print 'a = {:0.3f}'.format(C_c[2])
        print 'b = {:0.3f}'.format(C_c[1])
        print 'c = {:0.3f}'.format(C_c[0])

        fig = mp.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X_ab, Y_ab, Z_ab, rstride=1, cstride=1, alpha=0.2)
        ax.plot_surface(X_c, Y_c, Z_c, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(FeH_ab, RRab['logP'], W_ab, c='r', marker='o')
        ax.scatter(FeH_c, RRc['logP'], W_c, c='b', marker='o')
        mp.xlabel('[Fe/H]')
        mp.ylabel('logP')
        ax.set_zlabel('W('+band1+','+band2+'-'+band3+')')
        #mp.plot(data['logP'], data['i1'], 'ro')
        mp.show()

def get_PLZ(band, fundamentalized=0):
    bands_all = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4',
        'W1', 'W2', 'W3']

    col1 = bands_all.index(band)+6


    # Read in theoretical mean mags for two band relation
    dtype1 = np.dtype([('Z', float), ('Y', float), ('logP', float), (band, float)])
    RRab = np.loadtxt('FU-means.txt', dtype=dtype1, usecols=(0,1,5,col1), skiprows=33)
    RRc = np.loadtxt('FO-means.txt', dtype=dtype1, usecols=(0,1,5,col1), skiprows=33)

    if fundamentalized == 1:
        RRc['logP'] += 0.127
        period = np.append(RRc['logP'], RRab['logP'])
        Z = np.append(RRc['Z'], RRab['Z'])
        Y = np.append(RRc['Y'], RRab['Y'])
        FeH = np.log10(Z/(1-Z-Y))+1.61 - np.log10(0.694*10**0.4+0.306)
        FeH_c = np.log10(RRc['Z']/(1-RRc['Z']-RRc['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)
        FeH_ab = np.log10(RRab['Z']/(1-RRab['Z']-RRab['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)

        m1 = np.append(RRc[band], RRab[band])

        XX,YY = np.meshgrid(np.arange(-2.5, 0.5, 0.3), np.arange(-0.45, 0.25, 0.07))

    # best-fit linear plane
        A = np.c_[FeH, period, np.ones(period.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, m1)    # coefficients

    # evaluate it on grid
        ZZ = C[0]*XX + C[1]*YY + C[2]
        print '\nM('+band+') = a + b*logP + c*[Fe/H]\n'
        print 'Fundamentalized Coefficients: '
        print 'a = {:0.3f}'.format(C[2])
        print 'b = {:0.3f}'.format(C[1])
        print 'c = {:0.3f}'.format(C[0])

        fig = mp.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(FeH_c, RRc['logP'], RRc[band], c='b', marker='o')
        ax.scatter(FeH_ab, RRab['logP'], RRab[band], c='r', marker='o')
        mp.xlabel('[Fe/H]')
        mp.ylabel('logP')
        ax.set_zlabel('M('+band+')')
        #mp.plot(data['logP'], data['i1'], 'ro')
        mp.show()
    if fundamentalized == 0:

        FeH_c = np.log10(RRc['Z']/(1-RRc['Z']-RRc['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)
        FeH_ab = np.log10(RRab['Z']/(1-RRab['Z']-RRab['Y']))+1.61 - np.log10(0.694*10**0.4+0.306)

        X_c,Y_c = np.meshgrid(np.arange(-2.5, 0.5, 0.3), np.arange(-0.6, -0.2, 0.04))
        X_ab,Y_ab = np.meshgrid(np.arange(-2.5, 0.5, 0.3), np.arange(-0.45, 0.25, 0.07))

    # best-fit linear plane
        A_c = np.c_[FeH_c, RRc['logP'], np.ones(RRc.shape[0])]
        C_c,_,_,_ = scipy.linalg.lstsq(A_c, RRc[band])    # coefficients
        A_ab = np.c_[FeH_ab, RRab['logP'], np.ones(RRab.shape[0])]
        C_ab,_,_,_ = scipy.linalg.lstsq(A_ab, RRab[band])    # coefficients


    # evaluate it on grid
        Z_c = C_c[0]*X_c + C_c[1]*Y_c + C_c[2]
        Z_ab = C_ab[0]*X_ab + C_ab[1]*Y_ab + C_ab[2]
        print '\nM('+band+') = a + b*logP + c*[Fe/H]\n'
        print 'Fundamental Coefficients: '
        print 'a = {:0.3f}'.format(C_ab[2])
        print 'b = {:0.3f}'.format(C_ab[1])
        print 'c = {:0.3f}'.format(C_ab[0])
        print '\nFirst Overtone Coefficients: '
        print 'a = {:0.3f}'.format(C_c[2])
        print 'b = {:0.3f}'.format(C_c[1])
        print 'c = {:0.3f}'.format(C_c[0])

        fig = mp.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X_ab, Y_ab, Z_ab, rstride=1, cstride=1, alpha=0.2)
        ax.plot_surface(X_c, Y_c, Z_c, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(FeH_ab, RRab['logP'], RRab[band], c='r', marker='o')
        ax.scatter(FeH_c, RRc['logP'], RRc[band], c='b', marker='o')
        mp.xlabel('[Fe/H]')
        mp.ylabel('logP')
        ax.set_zlabel('M('+band+')')
        #mp.plot(data['logP'], data['i1'], 'ro')
        mp.show()
