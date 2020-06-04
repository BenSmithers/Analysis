import numpy as np
import os

def listattr( obj ):
    """
    This function lists the attributes belonging to some object "obj."
    It ignores any attributes starting with a "_", since these are usually special
        and/or indented to be access-restricted. 

    Not sure how this will work in other versions of python 

    It returns a list of strings
    """
    # this does most of the important stuff, just without the filter 
    temp = dir(obj)

    # lambda is the anonymous function, which provides the rules to the "filter" on temp
    #   filter returns another iterable, which we cast as a list 
    temp = list(filter(lambda x:x[0]!="_", temp))
    
    return(temp)


def explicit_convert( what ):
    """
    When data is loaded from the json file, it sometimes has weird datatypes. 
    Like, you don't have strings, you have Unicode representations. 
    This function takes a some data (like a dict or otherwise) and tries to convert it.
    If it encounters a dict, it recursively converts the thing 

    @param 'what' - dict OR other

    returns: dict OR other. Depends on what was passed originally. 

    Note: this function will fail on any dicts more than 99 levels deeps. This shouldn't happen.
    If you encounter this problem... **what are you doing???**
    """
    
    if isinstance(what, dict):
        new_dict = {}
        for key in what.keys():
            new_dict[key] = explicit_convert(what[key])

        return(new_dict)
    else:
        # basically it'll either be a float, string, or bool
        as_str = str(what).lower()
        if (as_str=='true' or as_str=='false'):
            return(bool(as_str))
        
        # check if it's one of those 0 and 1 strings
        is_weird = True
        for letter in as_str:
            if not (letter=='0' or letter=='1'):
                is_weird = False
                break
        if is_weird:
            if len(as_str)>1:
                return(str(what))

        try:
            # check - does casting this as an int result in a loss of precision?
            # if not, keep this as an int
            assert(1==1.) # if this ever changes I want to know...
            as_int   = int(what)
            as_float = float(what)
            if as_int==as_float:
                return(as_int)
            else:
                return(as_float)

        except ValueError:
            return(str(what))
        except TypeError:
            return(str(what))


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
            raise KeyError("Object 'obj' has not key '{}'!".format(key))
        setattr(obj, key, config_dict[key])
    
def parse_point( this_point ):
    """
    This keeps parameter-space string parsing locked away in a separate function so as
        not to waste time in 
    """
    if not is_valid_point(this_point):
        raise ValueError("{} is not a valid point".format(this_point))
    broken = this_point.split("_") # the point is loaded in as a unicode representaion.
    #           need to do a quick conversion 
   
    bsm_params = {}
    bsm_params['index']   = int(broken[0])
    bsm_params['dm41sq']     = float(broken[1])
    bsm_params['th14'] = float(broken[2])
    bsm_params['th24'] = float(broken[3])
    bsm_params['th34'] = float(broken[4])
    bsm_params['del14'] = float(broken[5])
    bsm_params['del24'] = float(broken[6])
    
    # this function isn't really used anymore
    def round_string(number):
        return("{:.4f}".format(np.round(float(number),4)))

    return( bsm_params )

def is_valid_point(this_point):
    """
    Returns True if the string passed is a valid point in phase space. Otherwise returns False

    This is here just to verify that the Point is valid before the main script does any 
        heavy lifting! 
    """
    if not isinstance(this_point, str):
        print("Point not a string! It's a {}. Trying to convert...".format(type(this_point)))
    broken = this_point.split("_")
    if not (len(broken)==7):
        print("Too " + ("many" if len(broken)>7 else "few") + " components in point!")
        return(False)
    for number in broken:
        try:
            x = float(number)
        except ValueError:
            print("Failed to cast {} as number".format(number))
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

    Is this just some container for a few arrays? Why not just use a dictionary??
    """
    def __init__(self, golemSpit):
        self.MuExEnergy = np.asarray(golemSpit[:,0]) #?????
        self.MuExZenith = np.asarray(golemSpit[:,1])
        self.MuExAzimuth= np.asarray(golemSpit[:,2])
        self.weights    = np.asarray(golemSpit[:,-1]) #??????? 
