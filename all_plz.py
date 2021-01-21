import PWZ
import numpy as np


df_fun = PWZ.get_PLZ(fundamentalized=0, suppress_output=1)
df = PWZ.get_PLZ(fundamentalized=1, suppress_output=1)

f1 = open('plz_coeff.dat', 'w')
f2 = open('plz_coeff_fundamentalized.dat', 'w')

fmt = '%4s %4s %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %5.3f'
np.savetxt(f2, df_fun.values, fmt=fmt)
np.savetxt(f1, df.values, fmt=fmt)


f1.close()
f2.close()
