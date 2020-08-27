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
    """
    This is the one with the obv knee
    """
    special = 2505.78207126
    norm_raw = 0.784758687019
    deltaIndex1=-0.000865349953528
    deltaIndex2=-0.00663812505081
    index1 = -2.5 + deltaIndex1
    index2 = -2.5 + deltaIndex2
    scale = 0.043283585459
    print("knee: {}".format(medianEnergy*(10**scale))) 
    indices = np.array([ (index1 if (ene > medianEnergy*(10**scale)) else index2) for ene in energy ])

    return( norm_raw*norm_scale*((energy/(medianEnergy*(10**scale)))**(indices)))

def get_broken_flux(energy):
    """
    This is the sigma one
    """
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
    deltaIndex = -0.0002861886750906706
    index = -2.5 + deltaIndex
    norm_raw = 0.7863344550132751
    what = norm_raw*norm_scale*(energy/medianEnergy)**(index)
    return(what)
    
flux = get_special_flux(energies)
other_flux = get_reg_flux(energies)
sigma_flux = get_broken_flux(energies)

figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

axes[1].axvline(x=medianEnergy, lw=0.3, color=(0.7,0.7,0.7))
axes[0].plot(energies, flux, label="T2-Broken")
#axes[0].plot(energies,sigma_flux, label="T1-Broken")
axes[0].plot(energies, other_flux, label="Unbroken")
axes[0].set_xscale('log')
axes[1].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlim([min(energies),max(energies)])
axes[1].set_xlabel("Energy [GeV]",size=16)
#axes[0].set_ylabel(r"Flux [cm$^{2}$ GeV s sr]$^{-1}$",size=16)

axes[0].legend()
axes[1].plot(energies, flux/other_flux, label="T2")
#axes[1].plot(energies, sigma_flux/other_flux, label="T1")
axes[1].set_ylabel("Broken/Unbroken")
#axes[1].set_yscale('log')
#axes[1].set_yticks([0.1,1,10])
#axes[1].set_ylim([0.1,10])
axes[1].legend()
figs.tight_layout()
figs.savefig("golem.png",dpi=400)
figs.show()
