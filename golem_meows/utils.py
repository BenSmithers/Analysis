import numpy as np
import os

def set_GF(obj, config_dict):
    """
    Sets the GF object 'obj' with parameters described in the 'config_dict' dictionary 

    --> the GF object MUST already have a key for EACH key in the 'config_dict' object 
    """
    if not isinstance(config_dict, dict):
        raise TypeError("Expected {} for 'config_dict', not {}".format(dict, type(config_dict)))
    # TODO
    # need to to type-checking of 'obj'
    # ---> check if they all inherit from some GF class? 

    for key in config_dict.keys():
        if not hasattr(obj, key):
            raise KeyError("Object 'obj' has not key {}!".format(key))
        setattr(obj, key, config_dict[key])
    
def parse_point( this_point ):
    """
    This keeps parameter-space string parsing locked away in a separate function so as
        not to waste time in 
    """
    if not is_valid_point(this_point):
        raise ValueError("{} is not a valid point".format(this_point))
    broken = this_point.split("_")
    
    def round_string(number):
        return("{:.4f}".format(np.round(float(number),4)))

    dm2 = round_string(broken[1])
    th24= float(broken[3])
    s2t = "{:.4f}".format(np.sin(2*th24)**2)
    th34= round_string(broken[4])
    Uu4sq= np.sin(th24)**2
    Ut4sq= (np.cos(th24)**2)*(np.sin(float(th34))**2)
    return( Uu4sq, Ut4sq )

def is_valid_point(this_point):
    """
    Returns True if the string passed is a valid point in phase space. Otherwise returns False

    This is here just to verify that the Point is valid before the main script does any 
        heavy lifting! 
    """
    if not isinstance(this_point, str):
        return(False)
    broken = this_point.split("_")
    if not (len(broken)==7):
        return(False)
    for number in broken:
        try:
            x = float(number)
        except ValueError:
            return(False)
    return(True)
    
def check_configuration(run_options):
    """
    This should be called once the user-provided configuration is loaded from the json. 

    The point is to verify that all the filepaths are valid before doing any time-consuming
        calculations/work. It's really just a sanity check to save the user time. 
    """
    if not os.path.exists(run_options['realization']):
        raise OSError("Realization not found: {}".format(run_options['realization']))
    if not os.path.isdir(run_options['datapath']):
        raise OSError("Datapath not valid directory (does it exist?): {}".format(run_options['datapath']))
    if not os.path.isdir(run_options['fluxdir']):
        raise OSError("Fluxdir not valid directory (does it exist?): {}".format(run_options['fluxdir']))
    if not os.path.isdir(run_options['outdir']):
        raise OSError("outdir is not a valid directory (does it exist?): {}".format(run_options['outdir']))
    if not is_valid_point(run_options['point']):
        raise ValueError("{} is not a valid point in phase space.".format(run_options['point']))

def get_seed(center, width, scale, low, high):
    """
    This is a trimmed gaussian distribution. It requires five floats as parameters:
        - center: mean of the distribution
        - width:  the width of the distribution
        - scale: the width of the distribution, but not? Totally unnecessary. Not sure why it's here
        - Low:  hard cut for the low-end of the distribution
        - High: hard cut for the high-end of the distribution

    Returns: a float in the range [low, high]
    """
    seed_value = np.random.normal(center, width/scale)
    return(max(min(high, seed_value), low))

def converter():
    """
    From Brazil.py

    I have Earthly clue what this does, why it's used, or what it takes. It just seemed important. 
    I'll add more comments here once I decipher the purpose of this arcane bullshit.

    Is this just some stupid container for a few arrays? Why not just use a dictionary??
    """
    def __init__(self, golemSpit):
        self.MuExEnergy = np.asarray(golemSpit[:,0]) #?????
        self.MuExZenith = np.asarray(golemSpit[:,1])
        self.MuExAzimuth= np.asarray(golemSpit[:,2])
        self.weights    = np.asarray(golemSpit[:,-1]) #??????? 
