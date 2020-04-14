"""
This script will calcualte the DIS cross sections for neutrinos and convolve that with an expected flux in the ice 
"""

import nuSQuIDSpy as nsq
import numpy as np # useful for energy ranges

constants = nsq.Const()

nBins = 100
eMin = 1.*constants.GeV
eMax = (10**12)*constants.GeV

energies = np.logspace( np.log10(eMin), np.log10(eMax), nBins)

# I use these dictionaries since the keys will be helpful in making plot labels 
flavors = {'electron' : nsq.NeutrinoCrossSections_NeutrinoFlavor.electron,
           'muon' : nsq.NeutrinoCrossSections_NeutrinoFlavor.muon,
           'tau': nsq.NeutrinoCrossSections_NeutrinoFlavor.tau
           }

currents = {'CC' : nsq.NeutrinoCrossSections_Current.CC, 
            'NC' : nsq.NeutrinoCrossSections_Current.NC
            }

neut_types = {'neutrino': nsq.NeutrinoCrossSections_NeutrinoType.neutrino, 
              'antineutrino': nsq.NeutrinoCrossSections_NeutrinoType.antineutrino
              }

xs_obj = nsq.NeutrinoDISCrossSectionsFromTables()



def get_diff_flux( energy, flavor, neutrino, current, y=None, x=None):
    """
    This function is just a wrapper to the nusquids function to make sure I get my types right.... 
    """
    if not (isinstance( energy, int) or isinstance(energy, float)):
        raise TypeError("Expected {} for energy, got {}".format(float, type(energy)))
    if not isinstance( flavor, nsq.NeutrinoCrossSections_NeutrinoFlavor):
        raise TypeError("Expected {} for flavor, got {}".format( nsq.NeutrinoCrossSections_NeutrinoFlavor, type(flavor)))
    if not isinstance( current, nsq.NeutrinoCrossSections_Current):
        raise TypeError("Expected {} for current, got {}".format(nsq.NeutrinoCrossSections_Current, type(current)))
    if not isinstance( neutrino, nsq.NeutrinoCrossSections_NeutrinoType):
        raise TypeError("Expected {} for neutrino type, got {}".format(nsq.NeutrinoCrossSections_NeutrinoType, type(neutrino)))

    if (y is None) or (x is None):
        return( xs_obj.TotalCrossSection( energy, flavor, neutrino, current))
    else:
        if not (isinstance(y,float) or isinstance(y,int)):
            raise TypeError("Expected {} for y, got {}".format(float, type(y)))
        if not (isinstance(x,float) or isinstance(x,int)):
            raise TypeError("Expected {} for x, got {}".format(float, type(x)))
        return(xs_obj.

def get_total_flux( energy, flavor, neutrino, current):
    return(get_diff_flux(energy, flavor, neutrino, current))
