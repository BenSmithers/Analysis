import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

GeV = 1.
TeV = (1.e3)*GeV
PeV = (1.e3)*TeV

npoints = 100
energies = np.logspace(np.log10(100*GeV),np.log10( 100*PeV),npoints)

norm_scale  = (1e-18)*GeV
medianEnergy= (1e5)*GeV

def get_special_flux(energy):
    special = 2505.78207126
    norm_raw = 0.786666691303
    deltaIndex1=-0.00318415276706
    deltaIndex2=-0.000569822383113
    index1 = -2.5 + deltaIndex1
    index2 = -2.5 + deltaIndex2
    scale = 0.823108017445
    print("knee: {}".format(medianEnergy*(10**scale))) 
    indices = np.array([ (index1 if (ene > medianEnergy*(10**scale)) else index2) for ene in energy ])

    return( norm_raw*norm_scale*((energy/(medianEnergy*(10**scale)))**(indices)))

def get_broken_flux(energy):
    llh = 2507.74631697
    norm_raw    = 0.789547622204
    deltaIndex1 = 0.00384884071536
    deltaIndex2 = -0.0187919232994
    index1      = -2.5 + deltaIndex1
    index2      = -2.5 + deltaIndex2
    sigma       = 0.501858830452
    what = norm_raw*norm_scale*((energy/medianEnergy)**(index1) + (1-sigma)*( (energy/medianEnergy)**(index1 +index2)))
    return(what)

def get_reg_flux(energy):
    llh = 2505.78189606
    deltaIndex = 0.00052004866302
    index = -2.5 + deltaIndex
    norm_raw = 0.786802828312
    what = norm_raw*norm_scale*(energy/medianEnergy)**(index)
    return(what)
    
flux = get_special_flux(energies)
other_flux = get_reg_flux(energies)
sigma_flux = get_broken_flux(energies)

plt.plot(energies, flux, label="Broken")
plt.plot(energies,sigma_flux, label="Sigma Flux")
plt.plot(energies, other_flux, label="Unbroken")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Energy [GeV]",size=16)
plt.ylabel(r"Flux [cm$^{2}$ GeV s sr]$^{-1}$",size=16)
plt.tight_layout()
plt.legend()
plt.savefig("/home/benito/Dropbox/Documents/Work/Research/IceCube/presentation/WeeklyUpdate/2020July/golem.png",dpi=400)
plt.show()
