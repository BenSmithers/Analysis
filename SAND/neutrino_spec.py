import matplotlib.pyplot as plt
import numpy as np
plt.style.use('/home/benito/Desktop/testing/paper.mplstyle')

"""
Plots the neutrino spectrum with errorbars for U-235 and U-238
"""
pi = np.pi
power_raw = 1.21 # GW 

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

# shape is square
target_sidelen = 2.0 #meters
target_cross_section = target_sidelen**2
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

portion_235 = 0.975

net_flux = portion_235*u235_spec + (1-portion_235)*u238_spec
net_erro = portion_235*u235_error + (1-portion_235)*u235_error

#print("Total: {}".format(sum([u235_spec[i]*widths[i] for i in range(len(energies)) ])))
#print("Uncertainty: {}".format(sum([u235_error[i]*widths[i] for i in range(len(energies))])))

print("Total: {}".format(sum([net_flux[i]*widths[i] for i in range(len(energies)) ])))
print("Uncertainty: {}".format(sum([net_erro[i]*widths[i] for i in range(len(energies))])))

output = open("output_neutrinos.dat",'w')
output.write("#Total at {:.2f} GW-th with {}% U235: {}\n".format(power_raw, portion_235*100,sum([net_flux[i]*widths[i] for i in range(len(energies)) ])))
output.write("#Uncertainty: {}\n".format(sum([net_erro[i]*widths[i] for i in range(len(energies))])))
output.write("#\n")
output.write("#BinCenter[MeV]    Nus_On_Target[MeV^-1 s^-1] sigma_nus[MeV^-1 s^-1]\n")
for line in range(len(energies)):
    output.write("{} {} {}\n".format(energies[line], net_flux[line], net_flux[line]))
output.close()

plt.errorbar( energies, net_flux, xerr=None, yerr=net_erro, capsize=5)
plt.title("Neutrino Emission for P={:.2f}GW-th".format(power_raw), size=14)
plt.xlabel(r"$\bar{\nu}_{e}$ Energy [MeV]", size=14)
plt.ylabel("dN/dE on Target [MeV$^{-1}$ s$^{-1}$]",size=14)
plt.yscale('log')
plt.tight_layout()
plt.show()
