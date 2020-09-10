# system tools
import os
import pickle 

# numeric tools
import numpy as np
from scipy import integrate 

# tau file verison should be updated with changes to this code! 
tau_file_version = "1.0"
tau_file_name = ".tau_branch.dat"

class TauData:
    """
    This pre-calculates some of the expectation values for the tau decay
    
    My reasoning for this system is that the integration is very costly. Since we're moving into a regime of higher dimensionality (3!), costly integration would cause tremendous slowdowns. 
    So, instead, we calculate a lot of integrals beforehand. Then, we can just interpolate between neighboring points for any arbitrary point and quickly get an accurate (enough) result 

    TODO: maybe do this with splines? 
    """
    def __init__(self):
        self._version = tau_file_version
        
        self._nodes = 1000
        self._energy_nodes = np.logspace(1, 8, self._nodes)

        self._expected_michell = np.zeros(self._nodes)
        for e_i in range(self._nodes):
            #boostedMichel( E_e, E_tau )
            vals = integrate.quad( lambda ene: BoostedMichel(ene, self._energy_nodes[e_i])*self._energy_nodes[e_i], 10, 10**8)
            self._expected_michell[e_i] = vals[0]


    def expected_electron(self, E_tau):
        """
        Returns the expectation value for the BoostedMichel function given E_tau
        Bascially tells you what the expected energy of the electron is of a decaying tau of energy E_tau
        """
        if not (isinstance(E_tau, float) or isinstance(E_tau,int)):
            raise TypeError("Expected {}, got {}".format(float,type(E_tau)))
        
        if E_tau<min(self._energy_nodes) or E_tau>max(self._energy_nodes):
            raise ValueError("E_tau of {:.2f} GeV outside of calculated bounds ({:.2f}GeV, {:.2f}GeV). Modify and update TauData, or maybe something else is bad?".format(E_tau, min(self._energy_nodes), max(self._energy_nodes)))

        upper_boundary = 1
        while E_tau>(self._energy_nodes[upper_boundary]):
            upper_boundary += 1
        lower_boundary = upper_boundary - 1

        #linear interpolation
        y2 = self._expected_michell[upper_boundary]
        y1 = self._expected_michell[lower_boundary]
        x2 = self.energy_nodes[upper_boundary]
        x1 = self.energy_nodes[lower_boundary]
        slope = (y2-y1)/(x2-x1)

        return( E_tau*slope + y2-(x2*slope) )

    @property
    def version(self):
        return(self._version)

    def __call__(self, E_tau):
        if not (isinstance(E_tau, float) or isinstance(E_tau, int)):
            raise TypeError("Expected {} for E_tau, not {}".format(float, type(E_tau)))

        return(1.0)


def args_are_floats(*args):
    """
    Simply makes sure that all the args you pass it are floats or ints 

    Returns NOTHING
    Raisese TypeError if something isn't an int or a float 
    """
    for arg in args:
        if not (isinstance(arg, float), isinstance(arg, int)):
            raise TypeError("Found an arg that's a {}, not a {}: {}".format(float, type(arg), arg))


"""
Much of this code is modified from 
https://github.com/IceCubeOpenSource/TauRunner/blob/master/python/Casino.py

So far I've only added type-checking using the above function
"""

#branching ratios of the charged tau 
RPion = 0.07856**2 
RRho = 0.43335**2 
RA1 = 0.70913**2
BrLepton = 0.18
BrPion = 0.12
BrRho = 0.26
BrA1 = 0.13
BrHad = 0.13

# Enu - outgoing tau neutrino (NC?)
# Etau - energy of tau lepton 
def BoostedMichel(E_e, E_t):
    """
    This is the differential flux rate 
    """
    args_are_floats(E_e, E_t)
    r = E_e / E_t    
    if r > 1.:
        return 0.
    return 1. / E_t * (5./3. - 3. *r**2 +4./3. * r**3)

def expected_EE(E_tau):
    """
    Calculates the expectation value of the boosted michel function for some charged tau energy

    This tells you what energy we expect of the daughter electron
    """

    exp_r = (5./3. - (3./4.) + (4./(3.*5.)))

def TauDecayToLepton(Etau, Enu, P):
    """
    """
    args_are_floats(Etau, Enu, P)
    # decays either to electron or muon
    # muon decay will be classified as tracks, so I'm throwing those away
    # therefore we suppress the spectra by the P_mu/(P_e+P_mu) ratio
    sup = 0.1739/(0.1739+0.1782)

    #z = Enu/Etau
    #g0 = (5./3.) - 3.*z**2 + (4./3.)*z**3
    #g1 = (1./3.) - 3.*z**2 + (8./3.)*z**3
    #return(g0+P*g1)

    # For a given Etau, we need the expected E_e 

def TauDecayToPion(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RPion - z)  > 0.0):
        g0 = 1./(1. - RPion)
        g1 = -(2.*z - 1. - RPion)/(1. - RPion)**2
    return(g0+P*g1)

def TauDecayToRho(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RRho - z) > 0.0):
        g0 = 1./(1. - RRho)
        g1 = -((2.*z-1.+RRho)/(1.-RRho))*((1.-2.*RRho)/(1.+2.*RRho))
    return(g0+P*g1)

def TauDecayToA1(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RA1 - z) > 0.0):
        g0 = (1./(1.-RA1))
        g1 = -((2.*z-1.+RA1)/(1.-RA1))*((1.-2.*RA1)/(1.+2.*RA1))
    return(g0 + P*g1)

def TauDecayToHadrons(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0=0.
    g1=0.
    if((0.3 - z) > 0.):
        g0 = 1./0.3
    return(g0+P*g1)

def TauDecayToAll(Etau, Enu, P):
    """
    This adds up the differential decay rates. 

    returns a (dn/dz) quantity, where dz is the width of considered E_nu/E_tau
    and n is the probability. Integrating all the (dn/dz)s is 1
    """
    args_are_floats(Etau, Enu, P)
    decay_spectra = 0
    decay_spectra+=2.0*BrLepton*TauDecayToLepton(Etau, Enu, P)
    decay_spectra+=BrPion*TauDecayToPion(Etau, Enu, P)
    decay_spectra+=BrRho*TauDecayToRho(Etau, Enu, P)
    decay_spectra+=BrA1*TauDecayToA1(Etau, Enu, P)
    decay_spectra+=BrHad*TauDecayToHadrons(Etau, Enu, P)
    return decay_spectra

Etau = 100.
zz = np.linspace(0.0,1.0,500)[1:-1]
dNTaudz = lambda z: TauDecayToAll(Etau, Etau*z, 0.)

def TauDecayToAllHadrons(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    """
    This adds up the differential decay rates 

    Enu is energy of outgoing tau neutrino
    ETau is energy of outgoing tau lepton (charged?)
    P is a polarization quantity (-1 for TauMinus)
    """
    #Etau is the energy of the tau lepton, Enu is the energy of the nu_tau after the tau decays
    decay_spectra = 0
    decay_spectra+=BrPion*TauDecayToPion(Etau, Enu, P)
    decay_spectra+=BrRho*TauDecayToRho(Etau, Enu, P)
    decay_spectra+=BrA1*TauDecayToA1(Etau, Enu, P)
    decay_spectra+=BrHad*TauDecayToHadrons(Etau, Enu, P)
    return decay_spectra/Etau

def generate_tau_file():
    """
    Generates and returns TauDecay data

    Save the file for later, quicker access
    """
    tau_file_loc = TauData()
    file_object = open(full_path, 'wb')
    pickle.dump(tau_file_loc, file_object, -1) # use the newest pickling techniques 
    file_object.close()

    return(tau_file_loc)


# try to load the file, if it's not there generate it
# if the file is old, regenerate it 
full_path = os.path.join(os.path.dirname(__file__), tau_file_name)
if os.path.exists(full_path):
    # load it, verify version number 
    file_object = open(full_path,'rb')
    tau_file = pickle.load(file_object)
    file_object.close()
    
    # this would be weird. Loaded taudata, but it's a different data type 
    # wtf? 
    if not isinstance(tau_file, TauData):
        raise TypeError("Loaded {}, not {}!".format(type(tau_file, TauData)))

    # old file, regenerate 
    if tau_file.version != tau_file_version:
        tau_file = generate_tau_file()

else: # no file, make file! 
    tau_file = generate_tau_file()


