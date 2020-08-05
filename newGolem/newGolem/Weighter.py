class Weighter:
    def __init__(self, dtype=float):
        if not isinstance(dtype, type):
            raise TypeError("Arg 'dtype' must be {}, got {}".format(type, type(dtype))) # weird...

        self._dtype = type
        self._friend = None

    @property
    def dtype(self):
        return(self._dtype)

    def __call__(self, event):
        """
        Calculate the weight of the event
        """
        return( 0.0 )

    def __add__(self, other):
        """
        Here, we combine two weighters into one super-weighter
        """
        if not isinstance(other, Weighter):
            raise TypeError("Expected {}, got {}".format(Weighter, type(other)))
        
        
class MetaWeighter(Weighter):
    def __init__(self, primary, secondary):
        if not isinstance(primary, Weighter):
            raise TypeError("Primary must be {}, got {}".format(Weighter, type(primary)))
        if not isinstance(secondary, Weighter):
            raise TypeError("Secondary must be {}, got {}".format(Weighter, type(primary)))

        # if these can be added, this should be valid
        self.dtype = type( primary.dtype() + secondary.dtype() )
        Weighter.__init__(self, self.dtype)

        self.primary = primary
        self.secondary = secondary

    def __call__(self, event):
        return( self.primary(event) + self.secondary(event) )
