import nuSQUIDSpy as nsq
from math import sqrt

import numpy as np

"""
This defines a few utility functions for my plotting script.

I moved this over here so the main plotter wasn't too busy
"""

class IllegalArguments(ValueError):
    """
    Just using this to make it clear what the issue is! 
    """
    pass

class bhist:
    """
    It's a 1D or 2D histogram! 
    """
    def __init__(self,edges):
        """
        Arg 'edges' should be a tuple of length 1 or 2. Length 1 for 1D hist, and length 2 for 2D hist
        """

        if not (isinstance(edges, list) or isinstance(edges, tuple) or isinstance(edges, np.ndarray)):
            raise TypeError("Arg 'edges' must be {}, got {}".format(list, type(edges)))

        if not(len(edges)==1 or len(edges)==2):
            raise ValueError("Arg 'edges' should be of length 1 or 2, depending on imensionality of desired bhist")
        for entry in edges:
            if not (isinstance(entry, list) or isinstance(entry, tuple) or isinstance(entry, np.ndarray)):
                raise TypeError("Each entry in 'edges' should be list-like, found {}".format(type(entry)))
            if len(entry)<2:
                raise ValueError("Entries in 'edges' must be at least length 2, got {}".format(len(entry)))
        
        self._edges = np.sort(edges) # each will now be increasing

        # build the function needed to register additions to the hisograms 
        if len(edges)==1:
            self._fill=np.zeros(len(self._edges[0])-1)
            def register(amount, where ):
                index = self._get_loc( where, self._edges[0] )
                if index is not None:
                    self._fill[index] += amount 
            self.register = register
        
        else: # length 2
            self._fill = np.zeros((len(self._edges[0])-1, len(self._edges[1])-1))
            def register(amount, xloc, yloc ):
                xbin = self._get_loc( xloc, self._edges[0] )
                ybin = self._get_loc( yloc, self._edges[1] )
                if (xbin is not None) and (ybin is not None):
                    self._fill[xbin][ybin]+=amount
            self.register = register

    def _get_loc(self, value, edges):
        """
        Private function used by the register function. Takes a value and bin edges, both along some axis, and it returns which bin the value is in. 
        """
        if value<edges[0] or value>edges[-1]:
            return
        else:
            scan = 0
            while not (value>=edges[scan] and value<=edges[scan+1]): # remember, the edges are sorted - so this should happen
                scan += 1
                if scan==len(edges)-1:
                    raise Exception("Something bad happened with logic")
            return(scan)
            self._fill[scan]+=value 

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

def get_exp_std( widths, probabilities, values ):
    """
    This takes a discretely binned probability density function representing a true continuous one,
    integrates it to get the expectation value and standard deviation of the distribution

    If the probability density isn't normalized, this function will normalize it

    TODO: get asymmetric error bars! 

    RETURNS; means, y+ error, y- error 
    """
    sigma = 0.5+0.341

    if not (len(widths)==len(values) and len(values)==len(probabilities)):
        raise ValueError("The args should have the same lengths")
   
    if not all(width>=0 for width in widths):
        raise ValueError("Found negative width, this can't be right.")

    if not all(prob>=0 for prob in probabilities):
        raise ValueError("Found negative probability. {}".format(min(probabilities)))

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

# define a couple utility functions
def get_flavor( key ):
    '''
    take a flux dictionary key and return the nusquids flavor type
    
    The dictionary key will be like "electorn_stuff_stuff"
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))

    part = key.split('_')[0].lower()
    if part in ['e', 'eleectron']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.electron )
    elif part in ['mu', 'muon']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.muon )
    elif part in ['tau']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.tau )
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))

def get_neut( key ):
    '''
    Takes a flux dictionary key and returns the nusquids neutrino type 
    (anti neutrino or vanilla neutrino)
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))   
    part = key.split('_')[1].lower()

    if part in ['nu', 'matter','neutrino']:
        return(nsq.NeutrinoCrossSections_NeutrinoType.neutrino)
    elif part in ['nubar', 'antimatter', 'antineutrino']:
        return(nsq.NeutrinoCrossSections_NeutrinoType.antineutrino)
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))

def get_curr(key):
    '''
    Takes a flux dictionary key and returns the nusquids neutrino type 
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))   
    part = key.split('_')[2].lower()

    if part in ['neutral', 'nc']:
        return(nsq.NeutrinoCrossSections_Current.NC)
    elif part in ['charged', 'cc']:
        return(nsq.NeutrinoCrossSections_Current.CC)
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))


