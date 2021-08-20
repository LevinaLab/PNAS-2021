# Network simulation and analysis for Sukenik et al. 2021
- ABC fit of adLIF network to inter-burst interval statistics
- Code for solving to the mean-field equations

ABC_nuki_fit.py - main script to run ABC on the bursting network
ABC core is adapted from https://github.com/rcmorehead/simpleabc

BIC_simuations.py - blocking inhibtion in the fitted network (in silico
bicuculline application)
nuki/ - scripts to parse the fitting results

setup.sh - script to install NEST on ubuntu 
nestmodel/ - nestml implementation of adLIF model


Paper: https://www.pnas.org/content/118/12/e2018459118

Experimental data is available at:
https://figshare.com/projects/Neuronal_circuits_overcome_imbalance_in_excitation_and_inhibition_by_adjusting_connection_numbers/88586





