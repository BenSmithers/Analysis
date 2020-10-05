"""
Defining my own flux for MCEq
"""

import crflux.models as crf

class HawkBPL(crf.PrimaryFlux):
    """
    Docstring
    """
    def __init__(self, ditch):
        self.name = "Hawk-BPL"
        self.params = {}

        self.e_norm = 1e5 #GeV
        
        # See https://arxiv.org/pdf/1710.00890.pdf
        # [Corsika ID] = (A, E_br, -gamma1, -gamma2)

        # the ID is calculated as (A x 100) + Z
        self.params[14] = (4.48e-2, 440.6, 2.81, 2.66) # H 
        self.params[402] = (3.31e-2, 854.6, 2.73, 2.54) # He
        self.params[1206] = (6.96e-6, 2882, 2.76, 2.55) # C
        self.params[1608] = (5.00e-6, 3843, 2.76, 2.55) # O
        self.params[2010] = (6.31e-7, 4803, 2.76, 2.55) # Ne
        self.params[2412] = (5.70e-7, 5764., 2.76, 2.55) # Mg
        self.params[2814] = (5.70e-7, 6725, 2.76, 2.55) # Si
        self.params[5426] = (2.00e-7, 13450, 2.76, 2.55) # Fe

        self.nucleus_ids = list(self.params.keys())

    def nucleus_flux(self, corsika_id, E):
        """
        Returns th eflux of nuclei corresponding to the corsika_id at energy E
        """
        corsika_id = self._find_nearby_id(corsika_id)
        
        entry = self.params[corsika_id]

        if E < entry[1]:
            return( entry[0]*pow(E/self.e_norm, entry[2]) )
        else:
            return( entry[0]*pow(entry[1]/self.e_norm, entry[2]-entry[3])*pow(E/self.e_norm, entry[3]) )
    
