# Theoretical-PWZ-Relations
Generate coefficients of theoretical PWZ relations presented in 
Neeley et al. 2017, ApJ, 841, 84.

![Example PLZ](example_plz.png?raw=true "Example PLZ")

Dependencies: Numpy, matplotlib, statsmodels, pandas, seaborn

This program uses mean magnitudes of RR Lyrae stars from pulsation models to derive theoretical 
period-Wesenheit-metallicity relations. 

The Wesenheit magnitude is defined as W(m1, m2 - m3) = m1 + alpha*(m2 - m3), where 
alpha is the ratio of the extinction coefficient of m1 and the color excess of bands m2 and m3,
alpha = A_m1 / (A_m2 - Am3)

A few examples are shown in the jupyter notebook Examples.ipynb 

To derive the coefficients using a subset of filters, use the function get_PWZ. If you wish to derive the 
coefficients for all possible combinations of filters and save the results to a file, you may use all_pwz.py. 

To use get_PWZ():

Enter a list of all the filters you are interested in, and the relevant extinction coefficients 
(i.e. A_m1 / Av) for each of those filters. If you want to compute all possible combinations, do not 
set the keywords bands or red_law. 

Inputs: 
    
    bands - list of three filters to use in calculating the Wesenheit magnitude (i.e. ['U','B','V'])
            Valid filters: 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'I1', 'I2', 'I3', 'I4', 'W1', 'W2', 'W3', 'W4'
            UBV: Johnson
            RI: Kron-Cousins
            JHK: 2MASS
            I1-4: Spitzer IRAC 3.6, 4.5, 5.8, and 8.0 micron bands
            W1-4: WISE 3.4, 4.6, 12, 24 micron bands

    red_law - extinction coefficients for your band as ratios A_lambda/Av (i.e [1.2,1.1,1.0]). Must be same length 
        as bands. Default values are defined according to the Cardelli et al. 1989 reddening law with Rv = 3.1, 
        extended into the MIR bands according to Indebetouw et al. 2005
    
    fundamentalized - flag to indicate whether you want separate first overtone and fundamental relations
        or fundamentalized relations. Default is 0 for separate, set to 1 for fundamentalized.
        
    suppress_output - flag to either print the coefficients to the terminal and display a plot, or suppress the 
        output. Default is 0, set to 1 to suppress. 
        
    num_bands = Number of bands you want to use to define the Wesenheit magnitude. Acceptable inputs are 2 or 3, 
        default is 2. 

Outputs:
    
    df_out - A pandas array containing the coefficients a, b, c, and the standard deviation.


