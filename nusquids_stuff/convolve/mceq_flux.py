from MCEq.core import config, MCEqRun
import crflux.models as crf

from bflux import HawkBPL

import time
import numpy as np
import sys

def get_angle_flux( angle , mag=0):
    """
    angle should be in degrees! 
    
    returns a dicitionary with the energies and the fluxes for each neutrino type 
    """

    if not (isinstance(angle,float) or isinstance(angle, int)):
        raise TypeError("Arg 'angle' should be {}, got {}".format(float, type(angle)))

    if not (angle<=180. and angle>=0.):
        raise ValueError("Invalid zenith angle {} deg".format(angle))

    flux = {}

    mceq = MCEqRun( 
        interaction_model = 'SIBYLL23C',
        primary_model = (HawkBPL, 0.),
#        primary_model = (crf.HillasGaisser2012, 'H3a'),
        theta_deg = angle
        )
    mceq.solve()

    flux['e_grid'] = mceq.e_grid

    flux['nue_flux'] = mceq.get_solution('nue',mag)+mceq.get_solution('pr_nue',mag)
    flux['nue_bar_flux'] = mceq.get_solution('antinue',mag)+mceq.get_solution('pr_antinue',mag)
    flux['numu_flux'] = mceq.get_solution('numu',mag)+mceq.get_solution('pr_numu',mag)
    flux['numu_bar_flux'] = mceq.get_solution('antinumu',mag)+mceq.get_solution('pr_antinumu',mag)
    flux['nutau_flux'] = mceq.get_solution('nutau',mag)+mceq.get_solution('pr_nutau',mag)
    flux['nutau_bar_flux'] = mceq.get_solution('antinutau',mag)+mceq.get_solution('pr_antinutau',mag)
    


    obj = open("temp_mceq_flux.dat".format(angle),'w')
#     obj.write("# E[Gev] Phi_nue[Gev/cm2/s/sr] Phi_numu[Gev/cm2/s/sr] Phi_nutau[Gev/cm2/s/sr] Phi_nue_bar Phi_numu_bar Phi_nutau_bar\n")
    for i in range(len(flux['e_grid'])):
        obj.write( str(flux['e_grid'][i]) )
        obj.write(' ')
        obj.write( str(flux['nue_flux'][i] ))
        obj.write(' ')
        obj.write( str(flux['numu_flux'][i] ))
        obj.write(' ')
        obj.write( str(flux['nutau_flux'][i] ))
        obj.write(' ')
        obj.write( str(flux['nue_bar_flux'][i] ))
        obj.write(' ')
        obj.write( str(flux['numu_bar_flux'][i] ))
        obj.write(' ')
        obj.write( str(flux['nutau_bar_flux'][i] ))
        obj.write('\n')
    obj.close()

get_angle_flux( float(sys.argv[1]) )

