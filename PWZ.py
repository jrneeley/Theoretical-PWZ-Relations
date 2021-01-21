import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import pandas as pd
import math


# All available bands
bands_all = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', '3.6', '4.5', '5.8', '8.0',
    'W1', 'W2', 'W3']
# Default reddening law Cardelli et al (1989) + Indebetouw et al. (2005)
cardelli = [1.579, 1.348, 1.016, 0.845, 0.590, 0.289, 0.181, 0.118,
        0.0661, 0.0507, 0.0507, 0.0507, 0.065, 0.052, 0.010]

def nCr(n, r):
    num = math.factorial(n)
    den = math.factorial(r)*math.factorial(n-r)
    return float(num)/float(den)

def get_PWZ(bands=bands_all, red_law=cardelli, fundamentalized=0,
    suppress_output=0, plot_3d=False, num_bands=2):


    # Read in theoretical mean mags
    dtype1 = np.dtype([('Z', float), ('Y', float), ('Msun', float),
        ('Teff', float), ('logL', float), ('logP', float), ('U', float),
        ('B', float), ('V', float), ('R', float), ('I', float), ('J', float),
        ('H', float), ('K', float), ('3.6', float), ('4.5', float), ('5.8', float),
        ('8.0', float), ('W1', float), ('W2', float), ('W3', float), ('W4', float)])

    RRab = np.loadtxt('FU-means.txt', dtype=dtype1, skiprows=33)
    RRc = np.loadtxt('FO-means.txt', dtype=dtype1, skiprows=33)

    # convert Z to [Fe/H]
    Z = np.append(RRc['Z'], RRab['Z'])
    Y = np.append(RRc['Y'], RRab['Y'])
    FeH = np.log10(Z/(1-Z-Y))+1.61 - 0.35

    if fundamentalized == 1:
        RRc['logP'] += 0.127
        #n_relations = math.comb(len(bands), num_bands)
        n_relations = int(nCr(len(bands), num_bands))
        dt2 = np.dtype([('band', 'U11'), ('mode', 'U4'), ('alpha', float), ('a', float),
            ('e_a', float), ('b', float), ('e_b', float), ('c', float),
            ('e_c', float), ('std', float)])
        coeff_data = np.zeros(n_relations, dtype=dt2)
    else:
        #n_relations = math.comb(len(bands), num_bands)
        n_relations = int(nCr(len(bands), num_bands))
        dt2 = np.dtype([('band', 'U11'), ('mode', 'U4'), ('alpha', float), ('a', float),
            ('e_a', float), ('b', float), ('e_b', float), ('c', float),
            ('e_c', float), ('std', float)])
        coeff_data = np.zeros(n_relations*2, dtype=dt2)

    logP = np.append(RRc['logP'], RRab['logP'])
    rrc_flag = np.repeat('FO', len(RRc['logP']))
    rrab_flag = np.repeat('FU', len(RRab['logP']))
    mode_flag = np.append(rrc_flag, rrab_flag)

    ii = 0
    if num_bands == 2:
        for ind1 in np.arange(1,len(bands)):
            band1 = bands[ind1]
            ext1 = red_law[ind1]
            mag1 = np.append(RRc[band1], RRab[band1])

            for ind2 in np.arange(0,ind1):
                band2 = bands[ind2]
                ext2 = red_law[ind2]
                mag2 = np.append(RRc[band2], RRab[band2])

                if ext1 == ext2:
                    coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band1)
                    coeff_data['mode'][ii] = 'FOFU'
                    ii += 1
                    if fundamentalized == 0:
                        coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band1)
                        coeff_data['mode'][ii] = 'FOFU'
                        ii += 1
                    continue

                alpha = ext1/(ext2 - ext1)
                if alpha < 0 :
                    coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band1)
                    coeff_data['mode'][ii] = 'FOFU'
                    ii += 1
                    if fundamentalized == 0:
                        coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band1)
                        coeff_data['mode'][ii] = 'FOFU'
                        ii += 1
                    continue
                W = mag1 - alpha*(mag2-mag1)
                #e_W = np.sqrt( ((1.+alpha)*err1)**2 + (alpha*err2)**2)

                if fundamentalized == 1:
                    A = np.c_[logP, FeH]
                    A = sm.add_constant(A)
                    res = sm.OLS(W, A).fit()

                    mag_fit = res.params[0] + res.params[1]*logP + res.params[2]*FeH
                    mag_residual = W - mag_fit
                    std = np.std(mag_residual)

                    # define the fit plane
                    XX,YY = np.meshgrid(np.arange(-0.45, 0.25, 0.07), np.arange(-3.0, 0.0, 0.3))
                    ZZ = res.params[0] + res.params[1]*XX + res.params[2]*YY

                    dt = np.dtype([('mode', 'U4'), ('logP', float), ('feh', float),
                        ('mag', float), ('mag_adj', float)])
                    data = np.zeros(len(W), dtype=dt)
                    data['mode'] = 'FOFU'
                    data['logP'] = logP
                    data['feh'] = FeH
                    data['mag'] = W
                    data['mag_adj'] = W - res.params[2]*FeH
                    df = pd.DataFrame(data)

                    coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band1)
                    coeff_data['mode'][ii] = 'FOFU'
                    coeff_data['alpha'][ii] = alpha
                    coeff_data['a'][ii] = res.params[0]
                    coeff_data['e_a'][ii] = res.bse[0]
                    coeff_data['b'][ii] = res.params[1]
                    coeff_data['e_b'][ii] = res.bse[1]
                    coeff_data['c'][ii] = res.params[2]
                    coeff_data['e_c'][ii] = res.bse[2]
                    coeff_data['std'][ii] = std
                    ii += 1

                    if suppress_output == 0:
                        print('\nW({},{}-{}) = a + b*logP + c*[Fe/H]\n'.format(band1, band2, band1))
                        print('Fundamentalized Coefficients: ')
                        print('a = {:0.3f}'.format(res.params[0]))
                        print('b = {:0.3f}'.format(res.params[1]))
                        print('c = {:0.3f}'.format(res.params[2]))
                        print('sigma = {:0.3f}'.format(std))

                        if plot_3d == False:
                            fig, ax =plt.subplots(1,2, figsize=(10,4))
                            # plot the data, color coded by metallicity
                            sns.scatterplot(data=df, x='logP', y='mag', hue='feh',
                                style='mode', ax=ax[0])
                            ax[0].invert_yaxis()
                            ax[0].legend([],[], frameon=False)
                            sns.scatterplot(data=df, x='logP', y='mag_adj', hue='feh',
                                style='mode', ax=ax[1])
                            ax[1].invert_yaxis()
                            ax[1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

                            # add in regression line
                            xx = np.linspace(-0.5, 0.2, num=20)
                            yy = res.params[0] + res.params[1]*xx
                            ax[1].plot(xx, yy, color='black')
                            plt.show()


                        else:  ### NEED TO FIX! - not currently interactive in jupyter notebook

                            fig = plt.figure()
                            ax = fig.add_subplot(111, projection='3d')
                            #ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, alpha=0.2)
                            ax.scatter(FeH, logP, W, c='black', marker='o')
                            ax.set_xlabel('[Fe/H]')
                            ax.set_ylabel('$\log$P')
                            ax.set_zlabel('M_{} mag'.format(band))
                            plt.show()

                else:
                    # Separate fits for RRc and RRab

                    logP_c = logP[mode_flag == 'FO']
                    logP_ab = logP[mode_flag == 'FU']
                    FeH_c = FeH[mode_flag == 'FO']
                    FeH_ab = FeH[mode_flag == 'FU']
                    mag_c = W[mode_flag == 'FO']
                    mag_ab = W[mode_flag == 'FU']

                    A_c = np.c_[logP_c, FeH_c]
                    A_c = sm.add_constant(A_c)
                    res_c = sm.OLS(mag_c, A_c).fit()

                    A_ab = np.c_[logP_ab, FeH_ab]
                    A_ab = sm.add_constant(A_ab)
                    res_ab = sm.OLS(mag_ab, A_ab).fit()

                    mag_fit_c = res_c.params[0] + res_c.params[1]*logP_c + res_c.params[2]*FeH_c
                    mag_residual_c = mag_c - mag_fit_c
                    std_c = np.std(mag_residual_c)

                    mag_fit_ab = res_ab.params[0] + res_ab.params[1]*logP_ab + res_ab.params[2]*FeH_ab
                    mag_residual_ab = mag_ab - mag_fit_ab
                    std_ab = np.std(mag_residual_ab)

                    # define the fit plane
                    #XX,YY = np.meshgrid(np.arange(-0.45, 0.25, 0.07), np.arange(-3.0, 0.0, 0.3))
                    #ZZ = res.params[0] + res.params[1]*XX + res.params[2]*YY

                    dt = np.dtype([('mode', 'U4'), ('logP', float), ('feh', float),
                        ('mag', float), ('mag_adj', float)])
                    data = np.zeros(len(W), dtype=dt)
                    data['mode'] = mode_flag
                    data['logP'] = logP
                    data['feh'] = FeH
                    data['mag'] = W
                    data['mag_adj'][mode_flag == 'FO'] = mag_c - res_c.params[2]*FeH_c
                    data['mag_adj'][mode_flag == 'FU'] = mag_ab - res_ab.params[2]*FeH_ab
                    df = pd.DataFrame(data)

                    coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band1)
                    coeff_data['mode'][ii] = 'FO'
                    coeff_data['alpha'][ii] = alpha
                    coeff_data['a'][ii] = res_c.params[0]
                    coeff_data['e_a'][ii] = res_c.bse[0]
                    coeff_data['b'][ii] = res_c.params[1]
                    coeff_data['e_b'][ii] = res_c.bse[1]
                    coeff_data['c'][ii] = res_c.params[2]
                    coeff_data['e_c'][ii] = res_c.bse[2]
                    coeff_data['std'][ii] = std_c
                    coeff_data['band'][ii+1] = '{},{}-{}'.format(band1, band2, band1)
                    coeff_data['mode'][ii+1] = 'FU'
                    coeff_data['alpha'][ii+1] = alpha
                    coeff_data['a'][ii+1] = res_ab.params[0]
                    coeff_data['e_a'][ii+1] = res_ab.bse[0]
                    coeff_data['b'][ii+1] = res_ab.params[1]
                    coeff_data['e_b'][ii+1] = res_ab.bse[1]
                    coeff_data['c'][ii+1] = res_ab.params[2]
                    coeff_data['e_c'][ii+1] = res_ab.bse[2]
                    coeff_data['std'][ii+1] = std_ab
                    ii += 2

                    if suppress_output == 0:
                        print('\nW({},{}-{}) = a + b*logP + c*[Fe/H]\n'.format(band1, band2, band1))
                        print('FO Coefficients: ')
                        print('a = {:.3f} +- {:.3f}'.format(res_c.params[0], res_c.bse[0]))
                        print('b = {:.3f} +- {:.3f}'.format(res_c.params[1], res_c.bse[1]))
                        print('c = {:6.3f} +- {:5.3f}'.format(res_c.params[2], res_c.bse[2]))
                        print('sigma = {:0.3f}'.format(std_c))
                        print('\nFU Coefficients: ')
                        print('a = {:.3f} +- {:.3f}'.format(res_c.params[0], res_c.bse[0]))
                        print('b = {:.3f} +- {:.3f}'.format(res_c.params[1], res_c.bse[1]))
                        print('c = {:6.3f} +- {:5.3f}'.format(res_c.params[2], res_c.bse[2]))
                        print('sigma = {:0.3f}'.format(std_ab))

                        if plot_3d == False:
                            fig, ax =plt.subplots(1,2, figsize=(11,4))
                            # plot the data, color coded by metallicity
                            sns.scatterplot(data=df, x='logP', y='mag', hue='feh',
                                style='mode', ax=ax[0])
                            ax[0].invert_yaxis()
                            ax[0].legend([],[], frameon=False)
                            ax[0].set_ylabel('W({},{}-{}) mag'.format(band1, band2, band1))

                            sns.scatterplot(data=df, x='logP', y='mag_adj', hue='feh',
                                style='mode', ax=ax[1])
                            ax[1].invert_yaxis()
                            ax[1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

                            # add in regression lines
                            xx = np.linspace(-0.5, 0.2, num=20)
                            yy = res_ab.params[0] + res_ab.params[1]*xx
                            ax[1].plot(xx, yy, color='black')

                            xx = np.linspace(-0.65, -0.3, num=20)
                            yy = res_c.params[0] + res_c.params[1]*xx
                            ax[1].plot(xx, yy, color='black')
                            ax[1].set_ylabel('W({},{}-{}) - c*[Fe/H]'.format(band1, band2, band1))

                            plt.show()
    ii = 0
    if num_bands == 3:
        for ind2 in range(len(bands)-2):
            band2 = bands[ind2]
            mag2 = np.append(RRc[band2], RRab[band2])

            for ind1 in np.arange(ind2+1, len(bands)-1):
                band1 = bands[ind1]
                mag1 = np.append(RRc[band1], RRab[band1])

                for ind3 in np.arange(ind1+1, len(bands)):
                    band3 = bands[ind3]
                    mag3 = np.append(RRc[band3], RRab[band3])

                    test = np.array([red_law[ind1], red_law[ind2], red_law[ind3]])
                    # skip combinations where two filters have the same extinction
                    if len(np.unique(test)) != 3:
                        coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band3)
                        coeff_data['mode'][ii] = 'FOFU'
                        ii += 1
                        if fundamentalized == 0:
                            coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band3)
                            coeff_data['mode'][ii] = 'FOFU'
                            ii += 1
                        continue

                    color = mag2 - mag3
                    alpha = red_law[ind1]/(red_law[ind2] - red_law[ind3])
                    if alpha < 0 :
                        coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band3)
                        coeff_data['mode'][ii] = 'FOFU'
                        ii += 1
                        if fundamentalized == 0:
                            coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band3)
                            coeff_data['mode'][ii] = 'FOFU'
                            ii += 1
                        continue
                    W = mag1 - alpha*color


                    if fundamentalized == 1:
                        A = np.c_[logP, FeH]
                        A = sm.add_constant(A)
                        res = sm.OLS(W, A).fit()

                        mag_fit = res.params[0] + res.params[1]*logP + res.params[2]*FeH
                        mag_residual = W - mag_fit
                        std = np.std(mag_residual)

                        # define the fit plane
                        XX,YY = np.meshgrid(np.arange(-0.45, 0.25, 0.07), np.arange(-3.0, 0.0, 0.3))
                        ZZ = res.params[0] + res.params[1]*XX + res.params[2]*YY

                        dt = np.dtype([('mode', 'U4'), ('logP', float), ('feh', float),
                            ('mag', float), ('mag_adj', float)])
                        data = np.zeros(len(W), dtype=dt)
                        data['mode'] = 'FOFU'
                        data['logP'] = logP
                        data['feh'] = FeH
                        data['mag'] = W
                        data['mag_adj'] = W - res.params[2]*FeH
                        df = pd.DataFrame(data)

                        coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band3)
                        coeff_data['mode'][ii] = 'FOFU'
                        coeff_data['alpha'][ii] = alpha
                        coeff_data['a'][ii] = res.params[0]
                        coeff_data['e_a'][ii] = res.bse[0]
                        coeff_data['b'][ii] = res.params[1]
                        coeff_data['e_b'][ii] = res.bse[1]
                        coeff_data['c'][ii] = res.params[2]
                        coeff_data['e_c'][ii] = res.bse[2]
                        coeff_data['std'][ii] = std
                        ii += 1

                        if suppress_output == 0:
                            print('\nW({},{}-{}) = a + b*logP + c*[Fe/H]\n'.format(band1, band2, band3))
                            print('Fundamentalized Coefficients: ')
                            print('a = {:.3f} +- {:.3f}'.format(res.params[0], res.bse[0]))
                            print('b = {:.3f} +- {:.3f}'.format(res.params[1], res.bse[1]))
                            print('c = {:6.3f} +- {:5.3f}'.format(res.params[2], res.bse[2]))
                            print('sigma = {:0.3f}'.format(std))

                            if plot_3d == False:
                                fig, ax =plt.subplots(1,2, figsize=(11,4))
                                # plot the data, color coded by metallicity
                                sns.scatterplot(data=df, x='logP', y='mag', hue='feh',
                                    style='mode', ax=ax[0])
                                ax[0].invert_yaxis()
                                ax[0].legend([],[], frameon=False)
                                ax[0].set_ylabel('W({},{}-{})'.format(band1, band2, band3))

                                sns.scatterplot(data=df, x='logP', y='mag_adj', hue='feh',
                                    style='mode', ax=ax[1])
                                ax[1].invert_yaxis()
                                ax[1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

                                # add in regression line
                                xx = np.linspace(-0.5, 0.2, num=20)
                                yy = res.params[0] + res.params[1]*xx
                                ax[1].plot(xx, yy, color='black')
                                ax[1].set_ylabel('W({},{}-{}) - c*[Fe/H]'.format(band1, band2, band3))
                                plt.show()


                            else:  ### NEED TO FIX! - not currently interactive in jupyter notebook

                                fig = plt.figure()
                                ax = fig.add_subplot(111, projection='3d')
                                #ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, alpha=0.2)
                                ax.scatter(FeH, logP, W, c='black', marker='o')
                                ax.set_xlabel('[Fe/H]')
                                ax.set_ylabel('$\log$P')
                                ax.set_zlabel('M_{} mag'.format(band))
                                plt.show()

                    else:
                        # Separate fits for RRc and RRab

                        logP_c = logP[mode_flag == 'FO']
                        logP_ab = logP[mode_flag == 'FU']
                        FeH_c = FeH[mode_flag == 'FO']
                        FeH_ab = FeH[mode_flag == 'FU']
                        mag_c = W[mode_flag == 'FO']
                        mag_ab = W[mode_flag == 'FU']

                        A_c = np.c_[logP_c, FeH_c]
                        A_c = sm.add_constant(A_c)
                        res_c = sm.OLS(mag_c, A_c).fit()

                        A_ab = np.c_[logP_ab, FeH_ab]
                        A_ab = sm.add_constant(A_ab)
                        res_ab = sm.OLS(mag_ab, A_ab).fit()

                        mag_fit_c = res_c.params[0] + res_c.params[1]*logP_c + res_c.params[2]*FeH_c
                        mag_residual_c = mag_c - mag_fit_c
                        std_c = np.std(mag_residual_c)

                        mag_fit_ab = res_ab.params[0] + res_ab.params[1]*logP_ab + res_ab.params[2]*FeH_ab
                        mag_residual_ab = mag_ab - mag_fit_ab
                        std_ab = np.std(mag_residual_ab)

                        # define the fit plane
                        #XX,YY = np.meshgrid(np.arange(-0.45, 0.25, 0.07), np.arange(-3.0, 0.0, 0.3))
                        #ZZ = res.params[0] + res.params[1]*XX + res.params[2]*YY

                        dt = np.dtype([('mode', 'U4'), ('logP', float), ('feh', float),
                            ('mag', float), ('mag_adj', float)])
                        data = np.zeros(len(W), dtype=dt)
                        data['mode'] = mode_flag
                        data['logP'] = logP
                        data['feh'] = FeH
                        data['mag'] = W
                        data['mag_adj'][mode_flag == 'FO'] = mag_c - res_c.params[2]*FeH_c
                        data['mag_adj'][mode_flag == 'FU'] = mag_ab - res_ab.params[2]*FeH_ab
                        df = pd.DataFrame(data)

                        coeff_data['band'][ii] = '{},{}-{}'.format(band1, band2, band3)
                        coeff_data['mode'][ii] = 'FO'
                        coeff_data['alpha'][ii] = alpha
                        coeff_data['a'][ii] = res_c.params[0]
                        coeff_data['e_a'][ii] = res_c.bse[0]
                        coeff_data['b'][ii] = res_c.params[1]
                        coeff_data['e_b'][ii] = res_c.bse[1]
                        coeff_data['c'][ii] = res_c.params[2]
                        coeff_data['e_c'][ii] = res_c.bse[2]
                        coeff_data['std'][ii] = std_c
                        coeff_data['band'][ii+1] = '{},{}-{}'.format(band1, band2, band3)
                        coeff_data['mode'][ii+1] = 'FU'
                        coeff_data['alpha'][ii+1] = alpha
                        coeff_data['a'][ii+1] = res_ab.params[0]
                        coeff_data['e_a'][ii+1] = res_ab.bse[0]
                        coeff_data['b'][ii+1] = res_ab.params[1]
                        coeff_data['e_b'][ii+1] = res_ab.bse[1]
                        coeff_data['c'][ii+1] = res_ab.params[2]
                        coeff_data['e_c'][ii+1] = res_ab.bse[2]
                        coeff_data['std'][ii+1] = std_ab
                        ii += 2

                        if suppress_output == 0:
                            print('\nW({},{}-{}) = a + b*logP + c*[Fe/H]\n'.format(band1, band2, band3))
                            print('FO Coefficients: ')
                            print('a = {:.3f} +- {:.3f}'.format(res_c.params[0], res_c.bse[0]))
                            print('b = {:.3f} +- {:.3f}'.format(res_c.params[1], res_c.bse[1]))
                            print('c = {:6.3f} +- {:5.3f}'.format(res_c.params[2], res_c.bse[2]))
                            print('sigma = {:0.3f}'.format(std_c))
                            print('\nFU Coefficients: ')
                            print('a = {:.3f} +- {:.3f}'.format(res_c.params[0], res_c.bse[0]))
                            print('b = {:.3f} +- {:.3f}'.format(res_c.params[1], res_c.bse[1]))
                            print('c = {:6.3f} +- {:5.3f}'.format(res_c.params[2], res_c.bse[2]))
                            print('sigma = {:0.3f}'.format(std_ab))

                            if plot_3d == False:
                                fig, ax =plt.subplots(1,2, figsize=(11,4))
                                # plot the data, color coded by metallicity
                                sns.scatterplot(data=df, x='logP', y='mag', hue='feh',
                                    style='mode', ax=ax[0])
                                ax[0].invert_yaxis()
                                ax[0].legend([],[], frameon=False)
                                ax[0].set_ylabel('W({},{}-{}) mag'.format(band1, band2, band1))

                                sns.scatterplot(data=df, x='logP', y='mag_adj', hue='feh',
                                    style='mode', ax=ax[1])
                                ax[1].invert_yaxis()
                                ax[1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

                                # add in regression lines
                                xx = np.linspace(-0.5, 0.2, num=20)
                                yy = res_ab.params[0] + res_ab.params[1]*xx
                                ax[1].plot(xx, yy, color='black')

                                xx = np.linspace(-0.65, -0.3, num=20)
                                yy = res_c.params[0] + res_c.params[1]*xx
                                ax[1].plot(xx, yy, color='black')
                                ax[1].set_ylabel('W({},{}-{}) - c*[Fe/H]'.format(band1, band2, band3))

                                plt.show()

    df_out = pd.DataFrame(coeff_data)

    return df_out

def get_PLZ(bands=bands_all, fundamentalized=0, suppress_output=0, plot_3d=False):


    # Read in theoretical mean mags
    dtype1 = np.dtype([('Z', float), ('Y', float), ('Msun', float),
        ('Teff', float), ('logL', float), ('logP', float), ('U', float),
        ('B', float), ('V', float), ('R', float), ('I', float), ('J', float),
        ('H', float), ('K', float), ('3.6', float), ('4.5', float), ('5.8', float),
        ('8.0', float), ('W1', float), ('W2', float), ('W3', float), ('W4', float)])

    RRab = np.loadtxt('FU-means.txt', dtype=dtype1, skiprows=33)
    RRc = np.loadtxt('FO-means.txt', dtype=dtype1, skiprows=33)

    # convert Z to [Fe/H]
    Z = np.append(RRc['Z'], RRab['Z'])
    Y = np.append(RRc['Y'], RRab['Y'])
    FeH = np.log10(Z/(1-Z-Y))+1.61 - 0.35

    if fundamentalized == 1:
        RRc['logP'] += 0.127
        dt2 = np.dtype([('band', 'U3'), ('mode', 'U4'), ('a', float),
            ('e_a', float), ('b', float), ('e_b', float), ('c', float),
            ('e_c', float), ('std', float)])
        coeff_data = np.zeros(len(bands), dtype=dt2)
    else:
        dt2 = np.dtype([('band', 'U3'), ('mode', 'U4'), ('a', float),
            ('e_a', float), ('b', float), ('e_b', float), ('c', float),
            ('e_c', float), ('std', float)])
        coeff_data = np.zeros(len(bands)*2, dtype=dt2)

    logP = np.append(RRc['logP'], RRab['logP'])
    rrc_flag = np.repeat('FO', len(RRc['logP']))
    rrab_flag = np.repeat('FU', len(RRab['logP']))
    mode_flag = np.append(rrc_flag, rrab_flag)



    for i in range(len(bands)):
        # do something
        band = bands[i]
        mag = np.append(RRc[band], RRab[band])



        if fundamentalized == 1:
            A = np.c_[logP, FeH]
            A = sm.add_constant(A)
            res = sm.OLS(mag, A).fit()

            mag_fit = res.params[0] + res.params[1]*logP + res.params[2]*FeH
            mag_residual = mag - mag_fit
            std = np.std(mag_residual)

            # define the fit plane
            XX,YY = np.meshgrid(np.arange(-0.45, 0.25, 0.07), np.arange(-3.0, 0.0, 0.3))
            ZZ = res.params[0] + res.params[1]*XX + res.params[2]*YY

            dt = np.dtype([('mode', 'U4'), ('logP', float), ('feh', float),
                ('mag', float), ('mag_adj', float)])
            data = np.zeros(len(mag), dtype=dt)
            data['mode'] = 'FOFU'
            data['logP'] = logP
            data['feh'] = FeH
            data['mag'] = mag
            data['mag_adj'] = mag - res.params[2]*FeH
            df = pd.DataFrame(data)

            coeff_data['band'][i] = band
            coeff_data['mode'][i] = 'FOFU'
            coeff_data['a'][i] = res.params[0]
            coeff_data['e_a'][i] = res.bse[0]
            coeff_data['b'][i] = res.params[1]
            coeff_data['e_b'][i] = res.bse[1]
            coeff_data['c'][i] = res.params[2]
            coeff_data['e_c'][i] = res.bse[2]
            coeff_data['std'][i] = std

            if suppress_output == 0:
                print('\nM('+band+') = a + b*logP + c*[Fe/H]\n')
                print('Fundamentalized Coefficients: ')
                print('a = {:.3f} +- {:.3f}'.format(res.params[0], res.bse[0]))
                print('b = {:.3f} +- {:.3f}'.format(res.params[1], res.bse[1]))
                print('c = {:6.3f} +- {:5.3f}'.format(res.params[2], res.bse[2]))
                print('sigma = {:0.3f}'.format(std))

                if plot_3d == False:
                    fig, ax =plt.subplots(1,2, figsize=(10,4))
                    # plot the data, color coded by metallicity
                    sns.scatterplot(data=df, x='logP', y='mag', hue='feh',
                        style='mode', ax=ax[0])
                    ax[0].invert_yaxis()
                    ax[0].legend([],[], frameon=False)
                    ax[0].set_ylabel('M_{}'.format(band))
                    sns.scatterplot(data=df, x='logP', y='mag_adj', hue='feh',
                        style='mode', ax=ax[1])
                    ax[1].invert_yaxis()
                    ax[1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

                    # add in regression line
                    xx = np.linspace(-0.5, 0.2, num=20)
                    yy = res.params[0] + res.params[1]*xx
                    ax[1].plot(xx, yy, color='black')
                    ax[1].set_ylabel('M_{} - c*[Fe/H]'.format(band))
                    plt.show()


                else:  ### NEED TO FIX! - not currently interactive in jupyter notebook

                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    #ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, alpha=0.2)
                    ax.scatter(FeH, logP, mag, c='black', marker='o')
                    ax.set_xlabel('[Fe/H]')
                    ax.set_ylabel('$\log$P')
                    ax.set_zlabel('M_{} mag'.format(band))
                    plt.show()

        else:
            # Separate fits for RRc and RRab

            logP_c = logP[mode_flag == 'FO']
            logP_ab = logP[mode_flag == 'FU']
            FeH_c = FeH[mode_flag == 'FO']
            FeH_ab = FeH[mode_flag == 'FU']
            mag_c = mag[mode_flag == 'FO']
            mag_ab = mag[mode_flag == 'FU']

            A_c = np.c_[logP_c, FeH_c]
            A_c = sm.add_constant(A_c)
            res_c = sm.OLS(mag_c, A_c).fit()

            A_ab = np.c_[logP_ab, FeH_ab]
            A_ab = sm.add_constant(A_ab)
            res_ab = sm.OLS(mag_ab, A_ab).fit()

            mag_fit_c = res_c.params[0] + res_c.params[1]*logP_c + res_c.params[2]*FeH_c
            mag_residual_c = mag_c - mag_fit_c
            std_c = np.std(mag_residual_c)

            mag_fit_ab = res_ab.params[0] + res_ab.params[1]*logP_ab + res_ab.params[2]*FeH_ab
            mag_residual_ab = mag_ab - mag_fit_ab
            std_ab = np.std(mag_residual_ab)

            # define the fit plane
            #XX,YY = np.meshgrid(np.arange(-0.45, 0.25, 0.07), np.arange(-3.0, 0.0, 0.3))
            #ZZ = res.params[0] + res.params[1]*XX + res.params[2]*YY

            dt = np.dtype([('mode', 'U4'), ('logP', float), ('feh', float),
                ('mag', float), ('mag_adj', float)])
            data = np.zeros(len(mag), dtype=dt)
            data['mode'] = mode_flag
            data['logP'] = logP
            data['feh'] = FeH
            data['mag'] = mag
            data['mag_adj'][mode_flag == 'FO'] = mag_c - res_c.params[2]*FeH_c
            data['mag_adj'][mode_flag == 'FU'] = mag_ab - res_ab.params[2]*FeH_ab
            df = pd.DataFrame(data)

            coeff_data['band'][2*i] = band
            coeff_data['mode'][2*i] = 'FO'
            coeff_data['a'][2*i] = res_c.params[0]
            coeff_data['e_a'][2*i] = res_c.bse[0]
            coeff_data['b'][2*i] = res_c.params[1]
            coeff_data['e_b'][2*i] = res_c.bse[1]
            coeff_data['c'][2*i] = res_c.params[2]
            coeff_data['e_c'][2*i] = res_c.bse[2]
            coeff_data['std'][2*i] = std_c
            coeff_data['band'][2*i+1] = band
            coeff_data['mode'][2*i+1] = 'FU'
            coeff_data['a'][2*i+1] = res_ab.params[0]
            coeff_data['e_a'][2*i+1] = res_ab.bse[0]
            coeff_data['b'][2*i+1] = res_ab.params[1]
            coeff_data['e_b'][2*i+1] = res_ab.bse[1]
            coeff_data['c'][2*i+1] = res_ab.params[2]
            coeff_data['e_c'][2*i+1] = res_ab.bse[2]
            coeff_data['std'][2*i+1] = std_ab

            if suppress_output == 0:
                print('\nM('+band+') = a + b*logP + c*[Fe/H]\n')
                print('FO Coefficients: ')
                print('a = {:.3f} +- {:.3f}'.format(res_c.params[0], res_c.bse[0]))
                print('b = {:.3f} +- {:.3f}'.format(res_c.params[1], res_c.bse[1]))
                print('c = {:6.3f} +- {:5.3f}'.format(res_c.params[2], res_c.bse[2]))
                print('sigma = {:0.3f}'.format(std_c))
                print('\nFU Coefficients: ')
                print('a = {:.3f} +- {:.3f}'.format(res_c.params[0], res_c.bse[0]))
                print('b = {:.3f} +- {:.3f}'.format(res_c.params[1], res_c.bse[1]))
                print('c = {:6.3f} +- {:5.3f}'.format(res_c.params[2], res_c.bse[2]))
                print('sigma = {:0.3f}'.format(std_ab))

                if plot_3d == False:
                    fig, ax =plt.subplots(1,2, figsize=(11,4))
                    # plot the data, color coded by metallicity
                    sns.scatterplot(data=df, x='logP', y='mag', hue='feh',
                        style='mode', ax=ax[0])
                    ax[0].invert_yaxis()
                    ax[0].legend([],[], frameon=False)
                    ax[0].set_ylabel('M_{} mag'.format(band))

                    sns.scatterplot(data=df, x='logP', y='mag_adj', hue='feh',
                        style='mode', ax=ax[1])
                    ax[1].invert_yaxis()
                    ax[1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

                    # add in regression lines
                    xx = np.linspace(-0.5, 0.2, num=20)
                    yy = res_ab.params[0] + res_ab.params[1]*xx
                    ax[1].plot(xx, yy, color='black')

                    xx = np.linspace(-0.65, -0.3, num=20)
                    yy = res_c.params[0] + res_c.params[1]*xx
                    ax[1].plot(xx, yy, color='black')
                    ax[1].set_ylabel('M_{} - c*[Fe/H]'.format(band))

                    plt.show()

    df_out = pd.DataFrame(coeff_data)


    return df_out
