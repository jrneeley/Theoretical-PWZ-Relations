import PWZ
import numpy as np


bands = ['B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4',
    'W1', 'W2', 'W3']
## Enter the coefficients of your reddening law here
ext = [1.348, 1.016, 0.845, 0.590, 0.289, 0.181, 0.118,
        0.0661, 0.0507, 0.0507, 0.0507, 0.065, 0.052, 0.010]

# 2 band PWZ relations
print 'Calculating 2 band PWZ relations...'
n_combinations = (len(bands)-1)*len(bands)/2

pwz2_a = np.zeros(3*n_combinations)
pwz2_b = np.zeros(3*n_combinations)
pwz2_c = np.zeros(3*n_combinations)
pwz2_ae = np.zeros(3*n_combinations)
pwz2_be = np.zeros(3*n_combinations)
pwz2_ce = np.zeros(3*n_combinations)
pwz2_sig = np.zeros(3*n_combinations)
pwz2_bands = np.zeros(3*n_combinations, dtype='S8')
rel_type = np.zeros(3*n_combinations, dtype='S5')

n_relation = 0
for ind, band in enumerate(bands):

    if ind == 0:
        continue

    band2 = ind - 1
    while band2 >= 0:

        bands_in = [bands[ind], bands[band2], bands[ind]]
        ext_in = [ext[ind], ext[band2], ext[ind]]
        if ext_in[1] == ext_in[2]:
            rel = bands[ind]+','+bands[band2]+'-'+bands[ind]
            table_ind = [n_relation, n_relation+1, n_relation+2]
            pwz2_bands[table_ind] = [rel, rel, rel]
            rel_type[table_ind] = ['FU', 'FO', 'FO+FU']
            n_relation += 3
            band2 += -1
            continue
        c, e = PWZ.get_PWZ(bands_in, ext=ext_in, suppress_output=1)
        c_fun, e_fun = PWZ.get_PWZ(bands_in, ext=ext_in, fundamentalized=1, suppress_output=1)
        table_ind = [n_relation+1, n_relation, n_relation+2]
        pwz2_a[table_ind] = [c[0], c[4], c_fun[0]]
        pwz2_b[table_ind] = [c[1], c[5], c_fun[1]]
        pwz2_c[table_ind] = [c[2], c[6], c_fun[2]]
        pwz2_ae[table_ind] = [e[0], e[3], e_fun[0]]
        pwz2_be[table_ind] = [e[1], e[4], e_fun[1]]
        pwz2_ce[table_ind] = [e[2], e[5], e_fun[2]]
        pwz2_sig[table_ind] = [c[3], c[7], c_fun[3]]
        rel = bands[ind]+','+bands[band2]+'-'+bands[ind]
        pwz2_bands[table_ind] = [rel, rel, rel]
        rel_type[table_ind] = ['FU', 'FO', 'FO+FU']
        band2 += -1
        n_relation += 3

data = np.array(zip(rel_type, pwz2_bands, pwz2_a, pwz2_b, pwz2_c, pwz2_ae, pwz2_be, pwz2_ce,
    pwz2_sig), dtype=[('type', 'S5'), ('bands', 'S8'), ('a', float), ('b', float), ('c', float),
    ('ae', float), ('be', float), ('ce', float), ('sig', float)])
np.savetxt('PWZ2.dat', data, fmt=['%5s', '%8s']+['%6.3f', '%6.3f', '%6.3f', '%6.3f',
    '%6.3f', '%6.3f', '%6.3f'])

# 3 band PWZ relations
print 'Calculating 3 band PWZ relations...'
n_combinations = 0
n_colors = len(bands)-2
while n_colors >= 1:
    n_combinations += n_colors*(n_colors+1)/2
    n_colors += -1


pwz3_a = np.zeros(3*n_combinations)
pwz3_b = np.zeros(3*n_combinations)
pwz3_c = np.zeros(3*n_combinations)
pwz3_ae = np.zeros(3*n_combinations)
pwz3_be = np.zeros(3*n_combinations)
pwz3_ce = np.zeros(3*n_combinations)
pwz3_sig = np.zeros(3*n_combinations)
pwz3_bands = np.zeros(3*n_combinations, dtype='S8')
rel_type = np.zeros(3*n_combinations, dtype='S5')

n_relation = 0
for ind in range(0,len(bands)-2):
#### NOT YET CORRECT
    band2 = ind
    for band1 in range(ind+1, len(bands)-1):
        band3 = band1 + 1
        while band3 < len(bands):
            bands_in = [bands[band1], bands[band2], bands[band3]]
            ext_in = [ext[band1], ext[band2], ext[band3]]
            if ext_in[1] == ext_in[2]:
                rel = bands[band1]+','+bands[band2]+'-'+bands[band3]
                table_ind = [n_relation, n_relation+1, n_relation+2]
                pwz3_bands[table_ind] = [rel, rel, rel]
                rel_type[table_ind] = ['FU', 'FO', 'FO+FU']
                n_relation += 3
                band3 += 1
                continue
            c, e = PWZ.get_PWZ(bands_in, ext=ext_in, suppress_output=1)
            c_fun, e_fun = PWZ.get_PWZ(bands_in, ext=ext_in, fundamentalized=1, suppress_output=1)
            rel = bands[band1]+','+bands[band2]+'-'+bands[band3]
            table_ind = [n_relation, n_relation+1, n_relation+2]
            pwz3_a[table_ind] = [c[0], c[4], c_fun[0]]
            pwz3_b[table_ind] = [c[1], c[5], c_fun[1]]
            pwz3_c[table_ind] = [c[2], c[6], c_fun[2]]
            pwz3_ae[table_ind] = [e[0], e[3], e_fun[0]]
            pwz3_be[table_ind] = [e[1], e[4], e_fun[1]]
            pwz3_ce[table_ind] = [e[2], e[5], e_fun[2]]
            pwz3_sig[table_ind] = [c[3], c[7], c_fun[3]]
            pwz3_bands[table_ind] = [rel, rel, rel]
            rel_type[table_ind] = ['FU', 'FO', 'FO+FU']
            band3 += 1
            n_relation += 3

data = np.array(zip(rel_type, pwz3_bands, pwz3_a, pwz3_b, pwz3_c, pwz3_ae, pwz3_be, pwz3_ce,
    pwz3_sig), dtype=[('type', 'S5'), ('bands', 'S8'), ('a', float), ('b', float), ('c', float),
    ('ae', float), ('be', float), ('ce', float), ('sig', float)])
np.savetxt('PWZ3.dat', data, fmt=['%5s', '%8s']+['%6.3f', '%6.3f', '%6.3f', '%6.3f',
    '%6.3f', '%6.3f', '%6.3f'])
