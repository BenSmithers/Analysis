from math import sqrt, log10

import os
import numpy as np
from functools import reduce
import sys

"""
This defines a few utility functions for my plotting script.

I moved this over here so the main plotter wasn't too busy
"""

flavors = ['E', 'Mu', 'Tau']
neuts = ['nu', 'nuBar']
currents = ['NC', 'CC']

def sci(number, precision=4):
    """
    Returns a string representing the number in scientific notation
    """
    if not isinstance(number, (int, float)):
        raise TypeError("Expected {}, not {}".format(float, type(number)))
    if not isinstance(precision,int):
        raise TypeError("Precision must be {}, not {}".format(int, type(precision)))

    power = int(log10(abs(number)))
    return("{0:.{1}f}".format(number/(10**power), precision)+"e{}".format( power))

def get_index( key ):
    '''
    Returns the column in the atmosphere data file for a given key
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))
    split = key.split('_')
    flavor = split[0] # E Mu Tau
    variety = split[1] # nu or nu-bar 

    flav_index = flavors.index(flavor) # 0, 1, or 2
    variety_index = neuts.index(variety) # 0 or 1
    return( 2 + int( flav_index + len(flavors)*variety_index) )

def bad_get_loc(value, edges):
    # fuck the error checking 
    if value<edges[0] or value>edges[-1]:
        return None
    else:
        scan = 0
        while not (value>=edges[scan] and value<=edges[scan+1]): # remember, the edges are sorted - so this should happen
            scan += 1
            if scan==len(edges)-1:
                raise Exception("Something bad happened with logic")
        return(scan)


def get_loc(x, domain,closest=False):
    """
    Returns the indices of the entries in domain that border 'x' 

    Raises exception if x is outside the range of domain 

    Assumes 'domain' is sorted!! 
    """
    if not isinstance(domain, (tuple,list,np.ndarray)):
        raise TypeError("'domain' has unrecognized type {}, try {}".format(type(domain), list))
    if not isinstance(x, (float,int)):
        raise TypeError("'x' should be number-like, not {}".format(type(x)))
    
    if x<domain[0] or x>domain[-1]:
        raise ValueError("x={} and is outside the domain: ({}, {})".format(sci(x), sci(domain[0]), sci(domain[-1])))

    # I think this is a binary search
    min_abs = 0
    max_abs = len(domain)-1

    lower_bin = int(abs(max_abs-min_abs)/2)
    upper_bin = lower_bin+1
    

    while not (domain[lower_bin]<=x and domain[upper_bin]>=x):
        if abs(max_abs-min_abs)<=1:
            raise Exception("Uh Oh")

        if x<domain[lower_bin]:
            max_abs = lower_bin
        if x>domain[upper_bin]:
            min_abs = upper_bin

        # now choose a new middle point for the upper and lower things
        lower_bin = min_abs + int(abs(max_abs-min_abs)/2)
        upper_bin = lower_bin + 1
    
    assert(x>=domain[lower_bin] and x<=domain[upper_bin])
    if closest:
        return( lower_bin if abs(domain[lower_bin]-x)<abs(domain[upper_bin]-x) else upper_bin )
    else:
        return(lower_bin, upper_bin)


def get_closest(x, domain, mapped):
    """
    We imagine some function maps from "domain" to "mapped"

    We have several points evaluated for this function
        domain - list-like of floats. 
        mapped - list-like of floats. Entries in domain, evaluated by the function

    The user provides a value "x," and then we interpolate the mapped value on either side of 'x' to approximate the mapped value of 'x' 
    """
    if not isinstance(domain, (tuple,list,np.ndarray)):
        raise TypeError("'domain' has unrecognized type {}, try {}".format(type(domain), list))
    if not isinstance(mapped, (tuple,list,np.ndarray)):
        raise TypeError("'mapped' has unrecognized type {}, try {}".format(type(mapped), list))
    if not isinstance(x, (float,int)):
        raise TypeError("'x' should be number-like, not {}".format(type(x)))

    if len(domain)!=len(mapped):
        raise ValueError("'len' and 'mapped' should have same length, got len(domain)={}, len(mapped)={}".format(len(domain), len(mapped)))
    
    lower_bin, upper_bin = get_loc(x, domain)
    
    # linear interp time
    x1 = domain[lower_bin]
    x2 = domain[upper_bin]
    y1 = mapped[lower_bin]
    y2 = mapped[upper_bin]

    slope = (y2-y1)/(x2-x1)
    value = (x*slope + y2 -x2*slope)

#    print("({}, {}) to ({}, {}) gave {}".format(x1,y1,x2,y2, value))
    
    return(value)

def bilinear_interp(p0, p1, p2, q11, q12, q21, q22):
    """
    Performs a bilinear interpolation on a 2D surface
    Four values are provided (the qs) relating to the values at the vertices of a square in the (x,y) domain
        p0  - point at which we want a value (len-2 tuple)
        p1  - coordinates bottom-left corner (1,1) of the square in the (x,y) domain (len-2 tuple)
        p2  - upper-right corner (2,2) of the square in the (X,y) domain (len-2 tuple)
        qs  - values at the vertices of the square (See diagram), any value supporting +/-/*
                    right now: floats, ints, np.ndarrays 


        (1,2)----(2,2)
          |        |
          |        |
        (1,1)----(2,1)

    """
    for each in [p0,p1,p2]:
        if not (isinstance(each, tuple) or isinstance(each, list) or isinstance(each, np.ndarray)):
            # I /would/ like to print out the bad value, but it might not be castable to a str?
            try:
                print("Found this: {}".format(each))
            except ValueError:
                pass
            raise TypeError("Expected {} for one of the points, got {}".format(tuple, type(each)))
        if len(each)!=2:
            raise ValueError("Points should be length {}, not {}!".format(2, len(each)))
    for val in [q11, q12, q21, q22]:
        if not isinstance(val, (float, int, np.ndarray)):
            raise TypeError("Expected {} for a value, got {}".format(float, type(val)))
    
    # check this out for the math
    # https://en.wikipedia.org/wiki/Bilinear_interpolation

    x0 = p0[0]
    x1 = p1[0]
    x2 = p2[0]
    y0 = p0[1]
    y1 = p1[1]
    y2 = p2[1]

    if not (x0>=x1 and x0<=x2):
        raise ValueError("You're doing it wrong. x0 should be between {} and {}, got {}".format(x1,x2,x0))
    if not (y0>=y1 and y0<=y2):
        raise ValueError("You're doing it wrong. y0 should be between {} and {}, got {}".format(y1,y2,y0))

    # this is some matrix multiplication. See the above link for details
    # it's not magic, it's math. Mathemagic 
    mat_mult_1 = [q11*(y2-y0) + q12*(y0-y1) , q21*(y2-y0) + q22*(y0-y1)]
    mat_mult_final = (x2-x0)*mat_mult_1[0] + (x0-x1)*mat_mult_1[1]

    return( mat_mult_final/((x2-x1)*(y2-y1)) )

class Data:
    """
    This is used as a container for the data loaded in from nuSQuIDS. 
    The main benefit of this is that the objects used by the interpolator are kept in a sterile scope, 
        and we don't have to worry about accidentally renaming an important object!  
    
    It loads it up into a convenient format for access, and provides a function for interpolating what is loaded. 
    """
    def __init__(self, filename='atmosphere.txt'):
        """
        Loads in the specified nuSQuIDS datafile. 

        Creates a "flux" dictionary for each type of neutrino and interaction. This is in units of N/s/GeV/cm2/sr
        """
        print("Extracting Data")
        data = np.loadtxt(os.path.join( os.path.dirname(__file__), 'atmosphere.txt'), dtype=float, comments='#',delimiter=' ')
        n_energies = 701
        n_angles = 100
        GeV=1e9
        if not (len(data)==n_energies*n_angles):
            raise ValueError("Datafile length error? {}!={}".format(len(data), n_energies*n_angles))

        # this funny indexing is a result of the way I output the data from nuSQuIDS
        # it loops through energies for each angle
        print("Building Flux Arrays")

        # storing the energies and the angles...
        # this was originally supposed to be agnostic to growing/shrinking energies/angles, but the code really assumes increasing.
        self._energies = [10**data[i][0] for i in range(n_energies)]
        self.growing = self._energies[1]>self._energies[0]  
        en_width = get_width(self._energies)/GeV
        self._angles = [data[n_energies*i][1] for i in range(n_angles)]
        self._ang_width = get_width(np.arccos(self._angles))
        self.ang_grow = self._angles[1]>self._angles[0]
        print(min(self._energies))
        print(type(min(self._energies)))
        print("Data spans {} GeV -> {} GeV".format(min(self._energies)/GeV, max(self._energies)/GeV))
        print("           {} rad -> {} rad in zenith".format(min(self._angles), max(self._angles)))

        print("The Energies are {} and the angles are {}".format("increasing" if self.growing else "decreasing", "increasing" if self.ang_grow else "decreasing"))
        if not (self.growing and self.ang_grow):
            raise NotImplementedError("I didn't really account for the case that the angles or energies were dereasing")

        # store the flavors, neutrinos, and currents
        self.flavors = flavors
        self.neuts = neuts
        self.currents = currents

        # let's fill out some flux functions
        # in the data file, these data are written in a big list. But that's not a very handy format
        # so I'm converting these into 2D arrays
        self._fluxes = {}
        for key in self.get_keys():
            # indexed like [energy_bin][angle_bin]
            # you may notice that NC and CC are treated as having separate fluxes, when really it'sthe same flux 
            #       this is for the most part okay since the interactions are rare enough that the fluxes are unchanged 
            self._fluxes[ key ] = [[ data[energy+angle*n_energies][get_index(key)]*2*np.pi for angle in range(n_angles)] for energy in range(n_energies)]
    
    # define a few access functions to protect the important stuff 
    # the "@property" tag makes it so these are accessed like attributes, not functions! 
    @property
    def energies(self):
        return(self._energies)
    @property
    def angles(self):
        return(self._angles)
    @property
    def ang_width(self):
        return(self._ang_width)
    @property
    def fluxes(self):
        return(self._fluxes) 

    def get_keys(self):
        """
        Returns a list of all the keys in the dictionary of fluxes 
        """
        keys = []
        for flav in self.flavors:
            for neut in self.neuts:
                for curr in self.currents:
                    key = flav+'_'+neut + '_'+curr
                    if flav=='Mu' and curr=='CC':
                        # skip tracks 
                        continue 
                    else: 
                        keys.append(key)
        return(keys)
    def get_flux(self, energy, key, use_overflow = False, angle=None):
        '''
        interpolates between entries in the flux dictionary to return the flux at arbitrary energy and angle
        Energy should be in units of eV
        Angle should be in units of radians 
        
        If an angle is provided, 
            Flux is in units of /cm2/GeV/s/sr  (incoming energy, bin width!) 
        If no angle is provided, the zenith angle is integrated over, and the
            Flux is in units of /cm2/GeV/s

        returns DOUBLE  (0.0 if beyond scope of data)
        '''
        if not (key in self._fluxes):
            raise ValueError("Bad key {}".format(key))
        if not (isinstance(energy, float) or isinstance(energy, int)):
            raise TypeError("Expected {}, not {}".format(float, type(energy)))

        if angle is not None:
            integrate = False
            if not (isinstance(angle, float) or isinstance(angle, int)):
                raise TypeError("Expected {}, not {}".format(float, type(angle)))
        else:
            integrate = True
            


        # check if it's outside the extents
        if (energy < self._energies[0] and self.growing) or (energy > self._energies[0] and not self.growing):
            # this overflow thing wasn't implemented right, so I'm disabling it for now .
            if False: # use_overflow:
                return(self._energies[0])
            else:
                return(0)
        if (energy > self._energies[-1] and self.growing) or (energy < self._energies[-1] and not self.growing):
            if False: #use_overflow:
                return(self._energies[n_energies - 1])
            else:
                return(0)

        if not integrate:
            if (angle < self._angles[0] and self.ang_grow) or (angle > self._angles[0] and not self.growing):
                return(0.)
            if (angle > self._angles[-1] and self.ang_grow) or (angle < self._angles[-1] and not self.growing):
                return(0.)
        
        lower_boundary, upper_boundary = get_loc(energy, self._energies)
        if not integrate:
            ang_lower, ang_upper = get_loc(angle, self._angles)
        
        # sanity check... 
        # essentially makes sure that the energies are monotonically increasing 
        if not ((self._energies[lower_boundary] <= energy) and (self._energies[upper_boundary] >= energy)):
            print("energy: {}".format(energy))
            print("lower bound: {}".format(self._energies[lower_boundary]))
            print("upper bound: {}".format(self._energies[upper_boundary]))
            print("indices: {}, {}".format(lower_boundary, upper_boundary))
            raise Exception()

        # if we're integrating, we get the flux at all the angles and scale it by width, otherwise it's bilinear interpolation time! 
        if integrate:
            flux_value = 0.0
            for angle_bin in range(len(self._angles)):
                # linear interpolation 
                y2 = self._fluxes[key][upper_boundary][angle_bin]
                y1 = self._fluxes[key][lower_boundary][angle_bin]
                x2 = self._energies[upper_boundary]
                x1 = self._energies[lower_boundary]
                slope = (y2-y1)/(x2-x1)

                flux_value += (energy*slope + y2 -x2*slope)*self._ang_width[angle_bin]
            return(flux_value)
        else:   
        #bilinear_interp(p0, p1, p2, q11, q12, q21, q22):
            p0 = (energy, angle)
            p1 = (self._energies[lower_boundary], self._angles[ang_lower])
            p2 = (self._energies[upper_boundary], self._angles[ang_upper])
            q11 = self._fluxes[key][lower_boundary][ang_lower]
            q21 = self._fluxes[key][upper_boundary][ang_lower]
            q12 = self._fluxes[key][lower_boundary][ang_upper]
            q22 = self._fluxes[key][upper_boundary][ang_upper]
            return(bilinear_interp(p0,p1,p2,q11,q12,q21,q22))


class IllegalArguments(ValueError):
    """
    Just using this to make it clear what the issue is! 
    """
    pass

class bhist:
    """
    It's a 1D or 2D histogram! BHist is for "Ben Hist" or "Binned Hist" depending on who's asking. 

    I made this so I could have a binned histogram that could be used for adding more stuff at arbitrary places according to some "edges" it has. The object would handle figuring out which of its bins would hold the stuff. 

    Also made with the potential to store integers, floats, or whatever can be added together and has both an additive rule and some kind of identity element correlated with the default constructor. 
        If a non-dtype entry is given, it will be explicitly cast to the dtype. 
    """
    def __init__(self,edges, dtype=float):
        """
        Arg 'edges' should be a tuple of length 1 or 2. Length 1 for 1D hist, and length 2 for 2D hist. 
        These edges represent the bin edges. 

        The type-checking could use a bit of work... Right now for 1D histograms you need to give it a length-1 list. 
        """

        if not (isinstance(edges, list) or isinstance(edges, tuple) or isinstance(edges, np.ndarray)):
            raise TypeError("Arg 'edges' must be {}, got {}".format(list, type(edges)))


        for entry in edges:
            if not (isinstance(entry, list) or isinstance(entry, tuple) or isinstance(entry, np.ndarray)):
                raise TypeError("Each entry in 'edges' should be list-like, found {}".format(type(entry)))
            if len(entry)<2:
                raise ValueError("Entries in 'edges' must be at least length 2, got {}".format(len(entry)))
        
        self._edges = np.sort(edges) # each will now be increasing. I think this is Quicksort? 
        self._dtype = dtype

        # Ostensibly you can bin strings... not sure why you would, but you could! 
        try:
            x = dtype() + dtype()
        except Exception:
            raise TypeError("It appears impossible to add {} together.".format(dtype))

        # build the function needed to register additions to the histograms.
        dims = tuple([len(self._edges[i])-1 for i in range(len(self._edges))])
        self._fill = np.zeros( shape=dims, dtype=self._dtype )
        # this function has *args, a lambda, AND arcane numpy stuff. Woo! 
        def register( amt, *args):
            """
            Tries to bin some data passed to the bhist. Arbitrarily dimensioned cause I was moving from 2D-3D and this seemed like a good opportunity 
                amt is the amount to add
                *args specifies the coordinates in our binned space 
            """
            if not len(args)==len(self._edges): 
                raise ValueError("Wrong number of args to register! Got {}, not {}".format(len(args), len(self._edges)))
            if not isinstance(amt, self._dtype):
                try:
                    amount = self._dtype(amt)
                except TypeError:
                    raise TypeError("Expected {}, got {}. Tried casting to {}, but failed.".format(self._dtype, type(amt), self._dtype))
            else:
                amount = amt

            bin_loc = tuple([get_loc( args[i], self._edges[i])[0] for i in range(len(args))]) # get the bin for each dimension

            # Verifies that nothing in the list is None-type
            if all([x is not None for x in bin_loc]):
                # itemset works like ( *bins, amount )
                try:
                    self._fill.itemset(bin_loc, self._fill.item(tuple(bin_loc))+amount)
                except TypeError:
                    print("bin_loc: {}".format(bin_loc))
                    print("amount: {}".format(amount))
                    print("previous: {}".format(self._fill.item(tuple(bin_loc))))
                    sys.exit()
                return tuple(bin_loc)
        self.register = register
 
    # some access properties. Note these aren't function calls. They are accessed like "object.centers" 
    @property
    def centers(self):
        complete = [ [0.5*(subedge[i+1]+subedge[i]) for i in range(len(subedge)-1)] for subedge in self._edges]
        return(complete[0] if len(self._edges)==1 else complete)
    @property
    def edges(self):
        complete = [[value for value in subedge] for subedge in self._edges]
        return(complete[0] if len(self._edges)==1 else complete)
    @property
    def widths(self):
        complete = [[abs(subedges[i+1]-subedges[i]) for i in range(len(subedges)-1)] for subedges in self._edges]
        return(complete[0] if len(self._edges)==1 else complete)
    @property
    def fill(self):
        return(self._fill)

def extract_from_edges(edges):
    """
    Takes list-like representing the edges of a bhist, returns the widths and centers 
    """
    if not (isinstance(edges, list) or isinstance(edges,np.ndarray) or isinstance(edges, tuple)):
        raise TypeError("Can't make bhist from {}, need {}".format(type(edges, list)))

    temp = bhist([edges])

    centers = np.array(temp.centers)
    widths = np.array(temp.widths)

    return( centers, widths )

def get_nearest_entry_to( item, array_like):
    """
    This function takes a quantity "item" and a list/array-like item "array_like"
    It returns the index of the entry in "array_like" that is closest to "item"

    Args:
        item - int
        array_like - list (or tuple, np.ndarray)
    Returns:
        index - int. =>0, <len(array_like)

    """
    # verify datatypes
    if not (isinstance(item,int) or isinstance(item, float)):
        raise TypeError("Expected number-like for arg 'item', got {}".format(type(item)))
    if not (isinstance(array_like, list) or isinstance(array_like, tuple) or isinstance(array_like,np.ndarray)):
        raise TypeError("Expected an index-able for arg 'array_like', got {}".format(type(array_like)))

    min_bin = None
    mindist = None

    # we can make no assumptions about the increasing/decreasing nature of 'array_like'
    # so we scan over it 
    for index in range(len(array_like)):
        if min_bin is None:
            min_bin = index
            mindist = abs(array_like[index] - item)
        else:
            new_distance = abs(array_like[index]-item)
            if new_distance < mindist:
                mindist = new_distance
                min_bin = index
  
    return(min_bin)

def get_width( which_list ):
    """
    Takes a list 'which_list' of floats of length N, considered the centers of some bins
    Returns a length N numpy array of floats for the widths of the bins these centers correspond to

    Arg 'which_list' must be monotonically increasing or decreasing
    """
    if not (isinstance(which_list, list) or isinstance(which_list, np.ndarray) or isinstance(which_list, tuple)):
        raise TypeError("Expected list-like, got {}".format(type(which_list)))

    # cast as a list
    if not isinstance(which_list, list):
        use = list(which_list)
    else:
        use = which_list

    n_bins = len(which_list)
    if n_bins<=1:
        raise ValueError("'which_list' should be longer than length {}".format(len(which_list)))

    increasing = which_list[1] > which_list[0]
    for i in range(n_bins-1):
        if not(increasing == (which_list[i+1] > which_list[i])):
            raise TypeError("Arg 'which_list' should be monotonically increasing or decreasing")

    n_bins = len(which_list)
    widths = [ 0. for i in range(n_bins)]
    for i in range(n_bins):
        if i==0:
            widths[i] = abs(use[1]-use[0])
        elif i==(n_bins-1):
            widths[i] = abs(use[-1]-use[-2])
        else:
            widths[i] = abs(0.5*(use[i+1]-use[i-1]))

    if not all(width>0 for width in widths):
        raise Exception("HOW {}".format(width))

    return(np.array(widths))

def get_exp_std( widths, probs, values ):
    """
    This takes a discretely binned probability density function representing a true continuous one,
    integrates it to get the expectation value and standard deviation of the distribution

    If the probability density isn't normalized, this function will normalize it

    TODO: get asymmetric error bars! 

    RETURNS; means, y+ error, y- error 
    """
    sigma = 0.5+0.341

    if not (len(widths)==len(values) and len(values)==len(probs)):
        raise ValueError("The args should have the same lengths")
   
    if not all(width>=0 for width in widths):
        raise ValueError("Found negative width, this can't be right.")

#    if not all(prob>=0 for prob in probabilities):
#       raise ValueError("Found negative probability. {}".format(min(probabilities)))

    probabilities = np.abs(probs)

    norm = sum([widths[i]*probabilities[i] for i in range(len(widths))  ])
    if norm==0.:
        return(0.,0.,0.)

    prob_density = np.array([ probabilities[i]/norm for i in range(len(probabilities))])
    if abs(1.-sum(prob_density*widths))>(1e-8):
        raise ValueError("Unable to normalize probabilities, got {}".format(sum(prob_density*widths)))
    

    mean=0.
    for item in range(len(widths)):
        # P(x)*X
        mean+= (widths[item]*prob_density[item])*values[item]

    median_bin = 0
    acc_prob = 0.
    while acc_prob<0.5:
        acc_prob+= widths[median_bin]*prob_density[median_bin]
        median_bin+=1
    if median_bin==len(prob_density):
        median_bin-=1
    median = values[median_bin]

    # get upper bound
    upper_bin = median_bin +1
    upper_b = acc_prob
    while upper_b < sigma:
        upper_b += widths[upper_bin]*prob_density[upper_bin]
        upper_bin += 1
        if upper_bin==len(prob_density):
            break

    if upper_bin==len(values):
        upper_bin-=1
    sigma_plus = values[upper_bin] - median

    lower_bin = median_bin -1
    lower_b = acc_prob
    while lower_b > (sigma - 0.5):
        lower_b -= widths[lower_bin]*prob_density[lower_bin]
        lower_bin-=1
        if lower_bin==-1:
            break

    if lower_bin==-1:
        lower_bin+=1
    sigma_minus = median - values[lower_bin]

    var=0.
    for item in range(len(widths)):
        var+= (widths[item]*prob_density[item])*(mean-values[item])**2
    sigma = sqrt(var)

    return( median, sigma_plus, sigma_minus)

