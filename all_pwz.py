import PWZ
import numpy as np



df_fun = PWZ.get_PWZ(fundamentalized=0, suppress_output=1, num_bands=2)
df = PWZ.get_PWZ(fundamentalized=1, suppress_output=1, num_bands=2)

f1 = open('pwz2_coeff.dat', 'w')
f2 = open('pwz2_coeff_fundamentalized.dat', 'w')

fmt = '%11s %4s %5.2f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %5.3f'
np.savetxt(f2, df_fun.values, fmt=fmt)
np.savetxt(f1, df.values, fmt=fmt)

f1.close()
f2.close()

df_fun = PWZ.get_PWZ(fundamentalized=0, suppress_output=1, num_bands=3)
df = PWZ.get_PWZ(fundamentalized=1, suppress_output=1, num_bands=3)

f1 = open('pwz3_coeff.dat', 'w')
f2 = open('pwz3_coeff_fundamentalized.dat', 'w')

fmt = '%11s %4s %5.2f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %5.3f'
np.savetxt(f2, df_fun.values, fmt=fmt)
np.savetxt(f1, df.values, fmt=fmt)

f1.close()
f2.close()
