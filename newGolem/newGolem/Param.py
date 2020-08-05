import numpy as np

from newGolem.utils import isNumber

class Param:
    """
    Basic Fit parameter datatype
    """
    def __init__(self, **kwargs):
        """
        This accepts any number of arguments in any order. Right now it only accepts 
            name     -   a string representing the name of this parameter 
            center   -   a prior representing the center of the distribution
            width
            range
            fit
        """

        # Priors
        self.name   = ""
        self.center = 0.0 
        self.width  = 1.0
        self.range  = ( -np.inf, np.inf )

        # should this be used in the fit? 
        self.fit = False

        # try assigning these from the keywords passed to this constructor
        for key in kwargs:
            # we only care about a few of these, and we want to do some basic checks
            if key=="name":
                if not isinstance(kwargs[key], name):
                    raise TypeError("Expected {} for {}, got {}".format(str, key, type(kwargs[key])))

            elif key=="center" or key=="width":
                if not isNumber(kwargs[key]):
                    raise TypeError("{} must be {}, received {}".format(key, float, type(kwargs[key])))
                
            elif key=="range":
                if isinstance( kwargs[key], tuple):
                    if len(kwargs[key])!=2:
                        raise ValueError("{} should be length 2, received {} entries".format(key, len(kwargs[key])))
                    for entry in kwargs[key]:
                        if not isNumber(entry):
                            raise TypeError("Entries in {} should all be numbers. Found {}".format(kwargs[key], type(entry)))

            elif key=="fit":
                if not isinstance(kwargs[key], bool):
                    raise TypeError("{} should be {}, received {}".format(key, bool, type(kwargs[key])))
            else:
                continue # ignoring any arg I don't recognize

            setattr(self, key, kwargs[key])

        # make sure that the assigned values pass a very simple check
        if not self.width < abs(self.range[1] - self.range[0]):
            raise ValueError("Width of distribution is greater than width of allowed values. This must be wrong...")

        if self.center<self.range[0] or self.center>self.range[1]:
            raise ValueError("Center of distribution is outside the range of allowed values.")

        if self.range[1]<self.range[0]:
            raise ValueError("Minimum is greater than maximum...")

    
