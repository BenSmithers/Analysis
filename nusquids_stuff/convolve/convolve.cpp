/*
    Ben Smithers
    
    This file will convolve fluxes with some cross sections and we'll see what we get

*/


#include <iostream>
#include <fstream>
#include <vector>
#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;

// this makes us consistent with the units we use
squids::Const un; 

double flux_function( double energy, double cos_zenith ){
    // energy should be in units of eV
    // cos_zenith should be, well, cos(zenith). 
    //       so between -1 and 1

    if (cos_zenith<-1 || cos_zenith>1){
        throw("Invalid cos_zenith received");
    }
    if (energy < 0){
        throw("I got a negative energy. Something is very bad");
    }
    
    double scale = pow(10., -18)*un.eV;
    double index = -2; //unitless 
    double knee  = 150.*un.TeV;
    
    return scale*pow(energy/knee, index );
}


int main(){
    // define some properties for our atmosphere 
    double n_nu = 3;
    double Emin = (1.e1)*un.GeV;
    double Emax = (1.e6)*un.GeV;
    double cos_zenith_min = -1.;
    double cos_zenith_max = 0.;

    int angular_bins = 40;
    int energy_bins = 100;
    bool use_earth_interactions = true;
    
    // create the atmosphere model
    
    auto zeniths = linspace(cos_zenith_min, cos_zenith_max, angular_bins);
    auto energies = logspace(Emin, Emax, energy_bins);

    std::cout << "Building NuSquids object" << std::endl;
    nuSQUIDSAtm<> nus_atm( zeniths, energies, n_nu, both, use_earth_interactions); 
    std::cout << "Done building nuSQUids object" << std::endl; 
    
    // some mixing angles
    // I don't quite understand why these have to be set manually 
    nus_atm.Set_MixingAngle(0,1,0.563942);
    nus_atm.Set_MixingAngle(0,2,0.154085);
    nus_atm.Set_MixingAngle(1,2,0.785398);
    nus_atm.Set_SquareMassDifference(1,7.65e-05);
    nus_atm.Set_SquareMassDifference(2,0.00247);

    
    nus_atm.Set_CPPhase(0, 2, 0);
    nus_atm.Set_CPPhase(1,2,-1.89); // sticking in the SK result

    // settting some zenith angle stuff 
    nus_atm.Set_rel_error(1.0e-6);
    nus_atm.Set_abs_error(1.0e-6);
    nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

    // set the initial state 
    marray<double, 4> inistate{angular_bins, energy_bins, 2, n_nu};
    std::fill( inistate.begin(), inistate.end(), 0);
    for ( int angle_bin=0; angle_bin < angular_bins; angle_bin++){
        for (int energy_bin=0; energy_bin < energy_bins; energy_bin++){
            for (int neut_type =0; neut_type<2; neut_type++){
                for (int flavor=0; flavor < n_nu; flavor++){
                    inistate[angle_bin][energy_bin][neut_type][flavor] = flux_function( energies[energy_bin], zeniths[angle_bin] );
                }
            }
        }
    }
    nus_atm.Set_initial_state( inistate, flavor);
    std::cout << "Done setting initial state" << std::endl;

    nus_atm.Set_ProgressBar(true);
    nus_atm.Set_IncludeOscillations(true);

    std::cout<<"Evolving..."<<std::endl;
    nus_atm.EvolveState();
    std::cout<<"Done. Writing!"<<std::endl;

    std::ofstream file("atmosphere.txt");
    // doing some interpolation
    int int_en = 700;
    int int_cos = 100;
    double int_min_e = log10(Emin);
    double int_max_e = log10(Emax);

    // write header
    file << "# log10(energy) cos(zenith) flux_nuE flux_nuMu flux_nuTau flux_nuEBar flux_nuMuBar flux_nuTauBar" <<std::endl;
    for(double angle=cos_zenith_min; angle<cos_zenith_max; angle+=(cos_zenith_max-cos_zenith_min)/(double)int_cos){
        for(double energy= int_min_e; energy<int_max_e; energy+=(int_max_e-int_min_e)/(double)int_en){
            // write out the angle and the energy 
            file << energy << " " << angle;
            double reg_energy = pow(10., energy);
            // write the neutrino contributions to the flux
            for( int flavor=0; flavor<n_nu; flavor++){
                file << " " << nus_atm.EvalFlavor( flavor, angle, reg_energy, 0);
            }
            // and now do it for anti-neutrinos 
            for( int flavor=0; flavor<n_nu; flavor++){
                file << " " << nus_atm.EvalFlavor( flavor, angle, reg_energy, 1);
            }
            file << std::endl;
        }
        file << std::endl;
    }

    return 0;

}
