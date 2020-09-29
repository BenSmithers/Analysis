from Event import Event 

class Weighter:
    def __init__(self, dtype=float):
        if not isinstance(dtype, type):
            raise TypeError("Arg 'dtype' must be {}, got {}".format(type, type(dtype))) # weird...

        self._dtype = dtype

    @property
    def dtype(self):
        return(self._dtype)

    def __call__(self, event):
        """
        Calculate the weight of the event
        """
        if not isinstance(event, Event):
            raise TypeError("Expected {}, got {}".format(Event, type(event)))

        return( self.dtype() )

    def __add__(self, other):
        """
        Here, we combine two weighters into one super-weighter
        This returns another Wighter object that evaluates the sum of the parent weighters' calculated weights 
        """
        if not isinstance(other, Weighter):
            raise TypeError("Expected {}, got {}".format(Weighter, type(other)))

        # create default event 
        ev = Event()
        dtype = type(other(ev) + self(ev))
        
        # make a little meta weighter object. It does weighting! 
        class metaWeighter(Weighter):
            def __init__(self_meta, dtype):
                Weighter.__init__(self_meta,dtype)
            def __call__(self_meta,event):
                return(self(event)+other(event))

        return(metaWeighter(dtype))

    def __mul__(self, other):
        """
        Define what happens when we multiply weighters together 

        This produces a meta-Weighter using two other weighters. This weighter weights a given event with the parent weighters, then multiplies the weights together and returns the product 
        """
        if not isinstance(other, Weighter):
            raise TypeError("Expected {}, got {}".format(Weighter, type(other)))

        ev = Event()
        dtype = type(other(ev)*self(ev))

        class metaWeighter(Weighter):
            def __init__(self_meta, dtype):
                Weighter.__init__(self_meta, dtype)
            def __call__(self_meta, event):
                return(self(event)*other(event))

        return(metaWeighter(dtype))
    
    def __div__(self, other):
        """
        Exactly the same as the multiplication, but now it's dividing 
        """
        if not isinstance(other, Weighter):
            raise TypeError("Expected {}, got {}".format(Weighter, type(other)))

        ev = Event()
        dtype = type(other(ev)/self(ev))

        class metaWeighter(Weighter):
            def __init__(self_meta, dtype):
                Weighter.__init__(self_meta, dtype)
            def __call__(self_meta, event):
                return(self(event)/other(event))

        return(metaWeighter(dtype))


