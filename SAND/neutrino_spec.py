import matplotlib.pyplot as plt
import numpy as np
plt.style.use('/home/benito/Desktop/testing/paper.mplstyle')

"""
Plots the neutrino spectrum with errorbars for U-235 and U-238
"""
pi = np.pi
power_raw = 0.98 # GW 

# define some units...
sec = 1.
mini = 60.*sec
hour = 60.*mini
MeV = 1.
eV = (1e-6)*MeV
GeV = (1e3)*MeV

joule = (6.242e18)/eV #MeV
watt = joule/sec # MeV/s
giga_watt = watt/(1e3)

# didn't actually end up using these
watt_hour = watt*hour
Gwatt_hour = watt*hour*(1e3)

# get the important units! 
# these will help convert the provided values into actual ones
reactor_out = power_raw*giga_watt
mev_per_fission = 200.*MeV
fission_rate = reactor_out/mev_per_fission # fission/sec 

# assume shape is circle
target_radius = 1.0 #meters
target_cross_section = pi*(target_radius**2)
baseline = 12.0 #meters
percent_on_target = target_cross_section / (4*pi*(baseline**2))

# load the data, get the bin widths 
data = np.transpose(np.loadtxt("neutrino_emission", delimiter=' '))
energies = data[0]
widths = np.array([0] + list(energies))
widths = widths[1:]-widths[:-1]
u235_spec = fission_rate*data[1]*percent_on_target
u235_error = u235_spec*data[2]/100.
u238_spec = fission_rate*data[3]*percent_on_target
u238_error = u238_spec*data[4]/100.

print("Total: {}".format(sum([u235_spec[i]*widths[i] for i in range(len(energies)) ])))

plt.errorbar( energies, u235_spec, xerr=None, yerr=u235_error, capsize=5,label=r"$^{235}$U")
plt.errorbar( energies, u238_spec, xerr=None, yerr=u238_error, capsize=5,label=r"$^{238}$U")
plt.legend()
plt.title("Neutrino Emission for P={:.2f}GW-th".format(power_raw), size=14)
plt.xlabel(r"$\bar{\nu}_{e}$ Energy [MeV]", size=14)
plt.ylabel("dN/dE on Target [MeV$^{-1}$ s$^{-1}$]",size=14)
plt.yscale('log')
plt.tight_layout()
plt.show()
