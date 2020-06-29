/*
    Ben Smithers
    
    This file uses MCEQ to generate fluxes and nuSQUIDS to propagate them through the earth. 
    Right now the MCEQ is handled by child processes running Python... because it works and it was easy.  

*/


#include <iostream>
#include <fstream>
#include <vector>
#include <nuSQuIDS/nuSQuIDS.h>
#include <array>
#include <fstream>
#include <math.h>
using namespace nusquids;

// this makes us consistent with the units we use
squids::Const un; 

const int angular_bins = 50;
const int energy_bins = 121;
const double pi = 3.14159265359;

// takes a string, returns a pointer to a length-4 array
const int columns = 7;
std::array<double, columns> parse( std::string line ){
    std::array<double, columns> results = {0., 0., 0., 0., 0., 0., 0.};
    int index = 0;
    
    std::string working = "";
    for (const char &c: line){
        if (index>=columns){
            std::cout << "Went too far??" <<std::endl;
        }
        if (c==' '){
            results[index] = std::stod( working );
            working = "";
            index++;
        }else{
            working += c; 
        }
    }
    // last character is probably not a space! 
    results[index] = std::stod(working);
    
    return( results );

}


// interpolator for mceq fluxes
// loads in data file and can 
struct inter{
    
    // energies [GeV]
    double energies[energy_bins] = {};
    // fluxes [#/Gev/cm^2/sr]
    // these should be replaced with a map std::string --> std:array<> 
    double nue_fluxes[energy_bins] ={};
    double numu_fluxes[energy_bins] = {};
    double nutau_fluxes[energy_bins] = {};
    
    double nue_bar_fluxes[energy_bins] ={};
    double numu_bar_fluxes[energy_bins] = {};
    double nutau_bar_fluxes[energy_bins] = {};

    void load_data(std::string which){  
        /*
            This function takes a filepath and loads the data from the file into this Object

            It loads information about the Energies and all three fluxes! 
        */

        std::string holder; // holds a string while interpreting it into floats and stuff
        std::ifstream target_data;
        target_data.open( which );
        int line = 0;
        while(std::getline(target_data, holder)){
//            std::getline( target_data, holder );
            if(holder[0]=='#'){
                continue;
            }
            std::array<double,columns> values = parse( holder );
            energies[line] = values[0];
            nue_fluxes[line] = values[1];
            numu_fluxes[line] = values[2];
            nutau_fluxes[line] = values[3];
            nue_bar_fluxes[line] = values[4];
            numu_bar_fluxes[line] = values[5];
            nutau_bar_fluxes[line] = values[6];           
            //std::cout << "Interpreted: "<<energies[line]<<", "<<nue_fluxes[line]<<", "<<numu_fluxes[line]<< ", " <<nutau_fluxes[line] << std::endl;

            line++; 
        }
    }

    double get_flux(double energy, uint8_t flavor,uint8_t nutype ){
        /*
            This returns the flux of the given flavor at desired energy

            It uses the loaded fluxes and linearly interpolates between the relevant points 
            This kind of function was used so that the nusquids and mcewq binning could be different 

            Energy shoudl be in GeV

            Maybe use an enum here... 
            flavor = {0:electron, 1:muon, 2:tau}
            nutype = {0:nu, 1:nubar}
         */

        // return 0 if outside the energy range! 
        if (energy < energies[0]){
            return(0.);
        }
        if (energy > energies[energy_bins-1]){
            return(0.);
        }

        uint16_t upper_b= 1;
        while( energy > energies[upper_b] ){
            upper_b += 1;
        }
        uint16_t lower_b = upper_b - 1;
        
        double y1,y2;
        double x2 = energies[upper_b];
        double x1 = energies[lower_b];

        // access the relevant stored flux! 
        if (flavor==0 && nutype==0){
            y2 = nue_fluxes[upper_b];
            y1 = nue_fluxes[lower_b];
        }else if(flavor==1 && nutype==0){
            y2 = numu_fluxes[upper_b];
            y1 = numu_fluxes[lower_b];
        }else if(flavor==2 && nutype==0){
            y2 = nutau_fluxes[upper_b];
            y1 = nutau_fluxes[lower_b];
        }else if (flavor==0 && nutype==1){
            y2 = nue_bar_fluxes[upper_b];
            y1 = nue_bar_fluxes[lower_b];
        }else if(flavor==1 && nutype==1){
            y2 = numu_bar_fluxes[upper_b];
            y1 = numu_bar_fluxes[lower_b];
        }else if(flavor==2 && nutype==1){
            y2 = nutau_bar_fluxes[upper_b];
            y1 = nutau_bar_fluxes[lower_b];
        }else{
            throw std::invalid_argument("");
        }

        // linear interpolation 
        double flux_value = energy*((y2-y1)/(x2-x1)) + y2 -x2*((y2-y1)/(x2-x1));
        return flux_value;
    }
};


// This function isn't used anymore. It was used to dictate the overall flux before I switched to using MCEQ 
double flux_function( double energy, double cos_zenith, uint8_t falvor){
    // energy should be in units of eV
    // cos_zenith should be, well, cos(zenith). 
    //       so between -1 and 1
    // flavor: 0 e, 1 mu, 2 tau

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


// this was just used to test thatthe interpolator was working right! 
/*
int main(){
    inter Interpolator;
    Interpolator.load_data("/home/bsmithers/software_dev/Analysis/fluxes/theta_0.dat");

}*/


int main(){
    inter FluxMachine; 

    // define some properties for our atmosphere 
    long unsigned int n_nu = 3;
    double Emin = 1.*un.GeV;
    double Emax = 10*un.PeV;
    double cos_zenith_min = -0.999;
    double cos_zenith_max = 0.;

    bool use_earth_interactions = true;
    
    // create the atmosphere model
    
    auto zeniths = linspace(cos_zenith_min, cos_zenith_max, angular_bins);
    auto energies = logspace(Emin, Emax, energy_bins);
    std::cout << "start at: "<<energies[0] << std::endl;

    std::cout << "Building NuSquids object" << std::endl;
    nuSQUIDSAtm<> nus_atm( zeniths, energies, n_nu, both, use_earth_interactions); 
    std::cout << "Done building nuSQUids object" << std::endl; 
    
    //nus_atm.Set_CPPhase(1,2,-1.89); // sticking in the SK result
    nus_atm.Set_MixingAngle(0,1,0.563942);
    nus_atm.Set_MixingAngle(0,2,0.154085);
    nus_atm.Set_MixingAngle(1,2,0.785398);

    nus_atm.Set_SquareMassDifference(1,7.65e-05);
    nus_atm.Set_SquareMassDifference(2,0.00247);


    // settting some zenith angle stuff 
    nus_atm.Set_rel_error(1.0e-6);
    nus_atm.Set_abs_error(1.0e-6);
    nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

    // set the initial state 
    marray<double, 4> inistate{angular_bins, energy_bins, 2, n_nu};
    std::fill( inistate.begin(), inistate.end(), 0);
    for ( int angle_bin=0; angle_bin < angular_bins; angle_bin++){
        // load next angle file  
        // we turn this around since the nusquids zenith angle is different from the MCEQ one 
        double angle_deg = 180.-(acos(zeniths[angle_bin])*180./pi);
        if (angle_deg > 180.){
            angle_deg = 180.;
        }
        std::string command("python3 mceq_flux.py ");
        command += std::to_string(angle_deg); // append the zenith angle to the end of this as an argument to the python script 
        std::cout << " ===> Generating new flux file at "<<angle_deg<<std::endl;
        system(command.c_str()); // call the python script to generate the flux
        FluxMachine.load_data("temp_mceq_flux.dat"); // load the flux in 
        std::cout << " ===> Assigning initial state " <<std::endl;
        std::cout << std::endl;
        for (int energy_bin=0; energy_bin < energy_bins; energy_bin++){
            for (uint8_t neut_type =0; neut_type<2; neut_type++){
                for (uint8_t flavor=0; flavor < n_nu; flavor++){
                    inistate[angle_bin][energy_bin][neut_type][flavor] = FluxMachine.get_flux(energies[energy_bin]/un.GeV, flavor, neut_type);
                }
            }
        }
    }
    nus_atm.Set_initial_state( inistate, flavor);
    std::cout << " ===> Done setting initial state" << std::endl;

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
            double scale = 1.;

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
    }

    return 0;

}

