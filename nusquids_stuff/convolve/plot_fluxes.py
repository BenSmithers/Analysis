#!/usr/bin/python3.6
'''
This script plots the fluxes output by the convolve cpp script
'''

import numpy as np
import matplotlib.pyplot as plt
import os


# load the data in
data = np.loadtxt(os.path.join( os.path.dirname(__file__), 'atmosphere.txt'), dtype=float, comments='#',delimiter=' ')
n_energies = 700
n_angles = 100
assert( len(data) == (n_energies*n_angles))

# this funnny indexing is a result of the way I output the data from nuSQuIDS
# it loops through energies for each angle
energies = [data[i][0] for i in range(n_energies)]
angles = [data[n_energies*i][1] for i in range(n_angles)]

# let's fill out some flux functions
# in the data file, these data are written in a big list. But that's not a very handy format
# so I'm converting these into 2D arrays
nuE_flux = [sum([ data[energy+angle*n_energies][2] for angle in range(n_angles)]) for energy in range(n_energies)] 
nuMu_flux = [sum([ data[energy+angle*n_energies][3] for angle in range(n_angles)]) for energy in range(n_energies)] 
nuTau_flux = [sum([ data[energy+angle*n_energies][4] for angle in range(n_angles)]) for energy in range(n_energies)] 

nuEBar_flux = [sum([ data[energy+angle*n_energies][5] for angle in range(n_angles)]) for energy in range(n_energies)] 
nuMuBar_flux = [sum([ data[energy+angle*n_energies][6] for angle in range(n_angles)]) for energy in range(n_energies)] 
nuTauBar_flux = [sum([ data[energy+angle*n_energies][7] for angle in range(n_angles)]) for energy in range(n_energies)] 

cmap = plt.get_cmap('viridis')
n_colors = 6
def get_color(which):
    return( cmap( which/n_colors ) )

scale_e = 10**np.array(energies)
plt.plot( scale_e, nuE_flux, color=get_color(0), label='nuE')
plt.plot( scale_e, nuMu_flux, color=get_color(1), label='nuMu')
plt.plot( scale_e, nuTau_flux, color=get_color(2), label='nuTau')
plt.plot( scale_e, nuEBar_flux, color=get_color(3), label='nuEBar')
plt.plot( scale_e, nuMuBar_flux, color=get_color(4), label='nuMuBar')
plt.plot( scale_e, nuTauBar_flux, color=get_color(5), label='nuTauBar')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()
