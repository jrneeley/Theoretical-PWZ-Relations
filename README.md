# Theoretical-PWZ-Relations
Generate coefficients of theoretical PWZ relations


This program uses mean magnitudes of RR Lyrae stars from pulsation models to derive theoretical 
period-Wesenheit-metallicity relations. 

The Wesenheit magnitude is defined as W(m1, m2, m3) = m1 + alpha*(m2 - m3), where 
alpha is the ratio of the extinction coefficient of m1 and the color excess of bands m2 and m3,
alpha = A_m1/Av / (A_m2/Av - Am3/Av)

Enter a list of the filters you want to use in the Wesenheit magnitude, and
the relevant extinction coefficients (e.g. A_m1 / Av) for each band. If you do not enter extinction coefficients,
alpha will be calculated using the Cardelli et al. 1989 reddening law. 

Inputs: 
    bands - list of three filters to use in calculating the Wesenheit magnitude (i.e. ['U','B','V'])
    ext - extinction coefficients for your three bands (i.e [1.2,1.1,1.0])
    fundamentalized - flag to indicate whether you want separate first overtone and fundamental relations
        or fundamentalized relations. Default is 0 for separate, set to 1 for fundamentalized.

Valid filters: 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4', 'W1', 'W2', 'W3', 'W4'

UBV are Johnson
RI are Kron-Cousins
JHK are 2MASS
I1 - I4 are Spitzer IRAC 3.6, 4.5, 5.8, and 8.0 micron bands
W1 - W4 are WISE 3.4, 4.6, 12, 24 micron bands

example:

get_PWZ(['V','B',V'])
