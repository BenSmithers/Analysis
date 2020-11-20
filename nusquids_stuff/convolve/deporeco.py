from math import exp, sqrt, pi, log, acos
import numpy as np
import os
from utils import bhist, get_loc

"""
This script is here to approximate the uncertainties in going from "energy deposited" to "energy reconstructed"
"""

rtwo = 1./sqrt(2*pi)

use_scan = False
if use_scan:
    from imagescan import get_data
    _data_depo, _data_reco, _data_prob = get_data()
else:
    datafile = "reconstruction.txt"
    data = np.loadtxt(os.path.join(os.path.dirname(__file__), datafile), dtype=float, delimiter=",")
    data = data.transpose()
    #data[1] = data[1]*data[0] 
    # data[0] is energy
    # data[1] is sigma

from utils import get_closest

rtwolog = sqrt(2*log(2))



class DataReco:
    def __init__(self, reco_energy_edges, reco_czenith_edges, depo_energy_edges, true_czenith_edges):
        """
        Expects the energies in eV, but works in GeV
        """


        self._ereco = bhist([np.array(reco_energy_edges)*(1e-9)])
        self._edepo = bhist([np.array(depo_energy_edges)*(1e-9)])
        self._zreco = bhist([reco_czenith_edges]) # these two are in cos(Zenith)
        self._ztrue = bhist([true_czenith_edges])

        # these are now filled with the values of the probability DENSITY functions for each angle/energy combo 
        # TODO right now there is no assumed covariance ... this should be improved 
        self._energy_odds_array = np.array([[ get_odds_energy(deposited, reconstructed) for reconstructed in self.reco_energy_centers] for deposited in self.depo_energy_centers])
        self._angle_odds_array = np.array([[ get_odds_angle(true, reconstructed) for reconstructed in self.reco_czenith_centers] for true in self.true_czenith_centers]) 

        # Avoid equating two floats. Look for a sufficiently small difference! 
        max_diff = 1e-12

        # normalize these things! 
        # so for each energy deposited... the sum of (PDF*width evaluated at each reco bin) should add up to 1. 
        for depo in range(len(self._energy_odds_array)):
            self._energy_odds_array[depo] *= 1./sum(self._energy_odds_array[depo]*self.reco_energy_widths)
            assert(abs(1-sum(self._energy_odds_array[depo]*self.reco_energy_widths)) <= max_diff)

        for depo in range(len(self._angle_odds_array)):
            self._angle_odds_array[depo] *= 1./sum(self._angle_odds_array[depo]*self.reco_czenith_widths)
            assert(abs(1-sum(self._angle_odds_array[depo]*self.reco_czenith_widths)) <= max_diff)
    
    def get_energy_reco_odds(self, i_depo, i_reco ):
        return(self._energy_odds_array[i_depo][i_reco]*self.reco_energy_widths[i_reco])

    def get_czenith_reco_odds(self, i_true, i_reco, i_e_true):
        return(self._angle_odds_array[i_true][i_reco]*self.reco_czenith_widths[i_reco])
   
    @property
    def reco_energy_centers(self):
        return(self._ereco.centers)
    @property
    def reco_czenith_centers(self):
        return(self._zreco.centers)
    @property
    def depo_energy_centers(self):
        return(self._edepo.centers)
    @property
    def true_czenith_centers(self):
        return(self._ztrue.centers)

    @property
    def reco_energy_widths(self):
        return(self._ereco.widths)
    @property
    def reco_czenith_widths(self):
        return(self._zreco.widths)
    @property
    def depo_energy_widths(self):
        return(self._edepo.widths)
    @property
    def true_czenith_widths(self):
        return(self._ztrue.widths)

deg = np.pi/180.
mu = deg*6.1983e-1
fwhm = deg*(13.223 + 11.983)
sigma =  fwhm/rtwolog

def get_odds_angle(true, reconstructed):
    """
    The returns the probability density of an RA/zenith angle of 'true' being reconstructed as 'reconstructed'

    It uses the 7-year cascade analysis' evaluation of angular reconstruction from Monopod 
    """
    x = acos(reconstructed) - acos(true)
    
    value = (rtwo/sigma)*exp(-1*((x-mu)**2)/(2*sigma*sigma))
    return(value)


def get_odds_energy(deposited, reconstructed):
    """
    Takes an energy deposited and energy reconstructed.
    Loads the datafile and reads off the uncertainty from the second column. 
    The data is in %E_depo, so we scale this to be a ratio and then by that deposted energy
    """
    if not isinstance(deposited, (float,int)):
        raise Exception()
    if not isinstance(reconstructed, (float,int)):
        raise Exception()
    if use_scan:
        x_id = get_loc(deposited, _data_depo, True)
        y_id = get_loc(reconstructed, _data_reco, True)

        return(_data_prob[x_id][y_id])
    else:
        sigma = get_closest( deposited, data[0], data[1])*deposited*0.01
        
        s2 = np.log(1+ (sigma/deposited)**2)
        mu = np.log((deposited**2)/np.sqrt(deposited**2  + sigma**2))
        # now, we assume that the uncertainty follows a log normal distribution, and calculate the PDF here

        #prob = rtwo*(1./sigma)*exp(-0.5*((reconstructed - deposited)/sigma)**2)
        prob = rtwo*(1./np.sqrt(s2))*exp(-0.5*((np.log(reconstructed) - mu)**2)/s2)/reconstructed
    
    return(prob)

doplot = False
if doplot:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib import ticker

    depos_e = np.logspace(0.1, 6.9, 100)*(1e9)
    recos_e = np.logspace(0.1, 6.9, 101)*(1e9)

    ang = np.linspace(-1,1,100)

    dataobj = DataReco( recos_e, ang, depos_e, ang)
    
    depos = dataobj.depo_energy_centers
    recos = dataobj.reco_energy_centers

    odds = np.array([[dataobj.get_energy_reco_odds(depo, reco) for depo in range(len(depos_e)-1)] for reco in range(len(recos_e)-1)])
    odds = np.log10(np.ma.masked_where(odds<=1e-10, odds))
    levels = np.logspace(-3,0,10)

    print(np.min(odds))
    print(np.max(odds))

    figs, axes = plt.subplots(nrows=2, ncols=1, sharex=False, gridspec_kw={'height_ratios':[1,1]})

    ctf = axes[0].pcolormesh(depos,recos,odds, cmap='gist_yarg')
    #plt.pcolormesh(depos,recos,odds, cmap='gist_gray', locator=ticker.LogLocator())
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_xlim([10**0,10**7])
    axes[0].set_ylim([10**0,10**7])
    axes[0].set_xlabel("Deposited [GeV]")
    axes[0].set_ylabel("Reconstructed [GeV]")
    cbar = plt.colorbar(ctf) #, ticks=ticker.LogLocator())
    cbar.set_label("Probability Density")

    # slides !
    slices = 10**np.array([1,2,3,4,5])
    for cut in slices:
        


        left, right = get_loc(float(cut), depos)
        loc = left if abs(cut-depos[left])<abs(cut-depos[right]) else right
        width_factor = dataobj.depo_energy_widths[loc]
        odds = [ dataobj.get_energy_reco_odds(loc, val) for val in range(len(recos))]
        axes[1].plot(recos, odds, color=(.1,.1,.1))
    
    axes[1].set_xscale('log')
    axes[1].set_xlabel("Reconstructed [GeV]")       
    axes[1].set_ylabel("Prob. Density")
    
    figs.savefig("fig_deporeco.png",dpi=400)
    plt.show()

    ang_odds = np.array([[dataobj.get_czenith_reco_odds(true,reco) for true in range(len(ang)-1)] for reco in range(len(ang)-1)])
    ang_odds=np.log10(np.ma.masked_where(ang_odds<=1e-10,ang_odds))
    
    plt.figure(2)
    me = plt.pcolormesh(ang,ang,ang_odds,cmap='gist_yarg')
    plt.xlabel("True Zenith")
    plt.ylabel("Recon Zenith")
    cbar = plt.colorbar(me)
    cbar.set_label("Prob Density")
    plt.savefig("fig_angtruereco.png",dpi=400)
    plt.show()

