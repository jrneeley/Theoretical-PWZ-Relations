# Theoretical-PWZ-Relations
Generate coefficients of theoretical PWZ relations

Dependencies: Numpy, matplotlib, scipy

This program uses mean magnitudes of RR Lyrae stars from pulsation models to derive theoretical 
period-Wesenheit-metallicity relations. 

The Wesenheit magnitude is defined as W(m1, m2 - m3) = m1 + alpha*(m2 - m3), where 
alpha is the ratio of the extinction coefficient of m1 and the color excess of bands m2 and m3,
alpha = A_m1/Av / (A_m2/Av - Am3/Av)

Enter a list of the filters you want to use in the Wesenheit magnitude, and
the relevant extinction coefficients (e.g. A_m1 / Av) for each band. If you do not enter extinction coefficients,
alpha will be calculated using the Cardelli et al. 1989 reddening law. 

Inputs: 
    
    bands - list of three filters to use in calculating the Wesenheit magnitude (i.e. ['U','B','V'])
            Valid filters: 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4', 'W1', 'W2', 'W3', 'W4'
            UBV: Johnson
            RI: Kron-Cousins
            JHK: 2MASS
            I1-4: Spitzer IRAC 3.6, 4.5, 5.8, and 8.0 micron bands
            W1-4: WISE 3.4, 4.6, 12, 24 micron bands

    ext - extinction coefficients for your three bands as ratios A_lambda/Av (i.e [1.2,1.1,1.0]). Default values 
        are defined according to the Cardelli et al. 1989 reddening law with Rv = 3.1, extended into the MIR bands 
        according to Indebetouw et al. 2005
    
    fundamentalized - flag to indicate whether you want separate first overtone and fundamental relations
        or fundamentalized relations. Default is 0 for separate, set to 1 for fundamentalized.


examples:

    Generate coefficients for 2 band PWZ relation W(V, B-V) with Cardelli reddening law
        > get_PWZ(['V','B',V'])

    Genereate coefficients for 3 band PWZ relation W(V, B-R) with Fitzpatrick reddening law
        > get_PWZ(['V','B','R'], [1.x, 1.x, 1.x])

    Generate coefficients of fundamentalized PWZ relation W(V, B-V) with Cardelli reddening law 
        > get_PWZ(['V','B','V'], fundamentalized=1)
