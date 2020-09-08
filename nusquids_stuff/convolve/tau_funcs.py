import os
import pickle 

# tau file verison should be updated with changes to this code! 
tau_file_version = "1.0"
tau_file_name = ".tau_branch.dat"

class TauData:
    """
    """
    def __init__(self):
        self._version = tau_file_version

    @property
    def version(self):
        return(self._version)

    def __call__(self, E_tau):
        if not (isinstance(E_tau, float) or isinstance(E_tau, int)):
            raise TypeError("Expected {} for E_tau, not {}".format(float, type(E_tau)))

        return(1.0)


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
    
    """
    args_are_floats(E_e, E_t)
    r = E_e / E_t    
    if r > 1.:
        return 0.
    return 1. / E_t * (5./3. - 3. *r**2 +4./3. * r**3)

def TauDecayToLepton(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = (5./3.) - 3.*z**2 + (4./3.)*z**3
    g1 = (1./3.) - 3.*z**2 + (8./3.)*z**3
    return(g0+P*g1)

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
            _ tau 
    --> O < 

    Enu is energy of outgoing tau neutrino
    ETau is energy of outgoing tau lepton (charged?)
    """
    #Etau is the energy of the tau lepton, Enu is the energy of the nu_tau after the tau decays
    decay_spectra = 0
    decay_spectra+=BrPion*TauDecayToPion(Etau, Enu, P)
    decay_spectra+=BrRho*TauDecayToRho(Etau, Enu, P)
    decay_spectra+=BrA1*TauDecayToA1(Etau, Enu, P)
    decay_spectra+=BrHad*TauDecayToHadrons(Etau, Enu, P)
    return decay_spectra/Etau
