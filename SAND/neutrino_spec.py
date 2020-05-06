import matplotlib.pyplot as plt
import numpy as np

"""
Plots the neutrino spectrum with errorbars for U-235 and U-238
"""

data = np.transpose(np.loadtxt("neutrino_emission", delimiter=' '))

energies = data[0]
u235_spec = data[1]
u235_error = data[1]*data[2]/100.
u238_spec = data[3]
u238_error = data[3]*data[4]/100.

plt.errorbar( energies, u235_spec, xerr=None, yerr=u235_error, capsize=5,label=r"$^{235}$U")
plt.errorbar( energies, u238_spec, xerr=None, yerr=u238_error, capsize=5,label=r"$^{238}$U")
plt.legend()
plt.title("Neutrino Emission for Uranium", size=14)
plt.xlabel(r"$\bar{\nu}_{e}$ Energy [MeV]", size=14)
plt.ylabel("dN/dE per fission [MeV$^{-1}$ boom$^{-1}$]",size=14)
plt.yscale('log')
plt.show()
