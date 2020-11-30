///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 264                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-06-06 21:05:12 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////
//#define _makeGammaPQ2_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "e_narrowResonanceCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
e_narrowResonanceCrossSection::e_narrowResonanceCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem&  bbsystem)
	:photonNucleusCrossSection(inputParametersInstance, bbsystem)
{
	_Ep         = inputParametersInstance.protonEnergy();	
	//this is in target frame
	_electronEnergy = inputParametersInstance.electronEnergy();
	//_target_beamLorentz = inputParametersInstance.targetBeamLorentzGamma();
	_target_beamLorentz = inputParametersInstance.beamLorentzGamma();
	_boost = std::acosh(inputParametersInstance.electronBeamLorentzGamma())
	  -std::acosh(inputParametersInstance.targetBeamLorentzGamma());
	_boost = _boost/2;
	// Photon energy limits in the two important frames
	_targetMaxPhotonEnergy=inputParametersInstance.targetMaxPhotonEnergy();
	_targetMinPhotonEnergy=inputParametersInstance.targetMinPhotonEnergy();
	_cmsMaxPhotonEnergy=inputParametersInstance.cmsMaxPhotonEnergy();
	_cmsMinPhotonEnergy=inputParametersInstance.cmsMinPhotonEnergy();
	//
	_VMnumEgamma = inputParametersInstance.nmbEnergyBins();
	_useFixedRange = inputParametersInstance.fixedQ2Range();
	_gammaMinQ2 = inputParametersInstance.minGammaQ2();
	_gammaMaxQ2 = inputParametersInstance.maxGammaQ2();
	_targetRadii = inputParametersInstance.targetRadius();
}


//______________________________________________________________________________
e_narrowResonanceCrossSection::~e_narrowResonanceCrossSection()
{ }


//______________________________________________________________________________
void
e_narrowResonanceCrossSection::crossSectionCalculation(const double)  // _bwnormsave (unused)
{
        // This subroutine calculates the vector meson cross section assuming
        // a narrow resonance.  For reference, see STAR Note 386.
  
        double dEgamma, minEgamma;
	double ega[3] = {0};
	double int_r,dR;
	double int_r2, dR2;
	int    iEgamma, nEgamma,beam;
	
	//Integration is done with exponential steps, in target frame
	//nEgamma = _VMnumEgamma;
	nEgamma = 1000;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy)/nEgamma;
	minEgamma = std::log(_targetMinPhotonEnergy);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	//
        printf(" gamma+nucleon threshold (CMS): %e GeV \n", _cmsMinPhotonEnergy);

        int A_1 = getbbs().electronBeam().A(); 
        int A_2 = getbbs().targetBeam().A();

	if( A_2 == 0 && A_1 >= 1 ){
	  // pA, first beam is the nucleus and is in this case the target  
	  beam = 1; 
	} else if( A_1 ==0 && A_2 >= 1){       
	  // photon energy in Target frame
	  beam = 2; 
	} else {
	  // photon energy in Target frame
	  beam = 2; 
	}
  
 	int_r=0.;
	int_r2= 0;
        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	// These might be useful in a future iteration
	int nQ2 = 1000;
        printf(" gamma+nucleon threshold (Target): %e GeV \n", _targetMinPhotonEnergy);
	for(iEgamma = 0 ; iEgamma < nEgamma; ++iEgamma){    // Integral over photon energy
	  // Target frame energies
	  ega[0] = exp(minEgamma + iEgamma*dEgamma );
	  ega[1] = exp(minEgamma + (iEgamma+1)*dEgamma );
	  ega[2] = 0.5*(ega[0]+ega[1]);

	  // Integral over Q2		
	  double dndE[3] = {0}; // Effective photon flux
	  double full_int[3] = {0}; // Full e+X --> e+X+V.M. cross section
	  //
	  for( int iEgaInt = 0 ; iEgaInt < 3; ++iEgaInt){    // Loop over the energies for the three points to integrate over Q2
	    //These are the physical limits
	    double Q2_min = std::pow(starlightConstants::mel*ega[iEgaInt],2.0)/(_electronEnergy*(_electronEnergy-ega[iEgaInt]));	    
	    double Q2_max = 4.*_electronEnergy*(_electronEnergy-ega[iEgaInt]);
	    if(_useFixedRange == true){
	      if( Q2_min < _gammaMinQ2 )
		Q2_min = _gammaMinQ2;
	      if( Q2_max > _gammaMaxQ2 )
		Q2_max = _gammaMaxQ2;
	    }
	    double lnQ2ratio = std::log(Q2_max/Q2_min)/nQ2;
	    double lnQ2_min = std::log(Q2_min);
	    //

	    for( int iQ2 = 0 ; iQ2 < nQ2; ++iQ2){     // Integral over photon virtuality
	      //
	      double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
	      double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
	      double q2_12 = (q2_2+q2_1)/2.;
	      // Integrating the effective photon flux
	      dndE[iEgaInt] +=(q2_2-q2_1)*( getcsgA_Q2_dep(q2_1)*photonFlux(ega[iEgaInt],q2_1)
					    +getcsgA_Q2_dep(q2_2)*photonFlux(ega[iEgaInt],q2_2)
					    +4.*getcsgA_Q2_dep(q2_12)*photonFlux(ega[iEgaInt],q2_12) );
	      // Integrating cross section
	      full_int[iEgaInt] += (q2_2-q2_1)*( g(ega[iEgaInt],q2_1)*getcsgA(ega[iEgaInt],q2_1,beam)
						 + g(ega[iEgaInt],q2_2)*getcsgA(ega[iEgaInt],q2_2,beam)
						 + 4.*g(ega[iEgaInt],q2_12)*getcsgA(ega[iEgaInt],q2_12,beam) );	      
	    }

	    // Finish the Q2 integration for the three end-points (Siumpsons rule)
	    dndE[iEgaInt] = dndE[iEgaInt]/6.; 
	    full_int[iEgaInt] = full_int[iEgaInt]/6.;
	  }	
	  // Finishing cross-section integral 
	  dR = full_int[0];
	  dR += full_int[1];
	  dR += 4.*full_int[2];
	  dR = dR*(ega[1]-ega[0])/6.;
	  
	  // Finishing integral over the effective photon flux
	  dR2 = dndE[0];
	  dR2 += dndE[1];
	  dR2 += 4.*dndE[2];
	  dR2 = dR2*(ega[1]-ega[0])/6.;
	  //
	  int_r = int_r+dR;
	  int_r2 = int_r2 + dR2;
	}
	cout<<endl;      
	if(_useFixedRange == true){
	  cout<<" Using fixed Q2 range "<<_gammaMinQ2 << " < Q2 < "<<_gammaMaxQ2<<endl;
	}
	printCrossSection(" Total cross section: ",int_r);
	//printCrossSection(" gamma+X --> VM+X ", int_r/int_r2);      // commented out for the mean time, not necesary in current implementation
	//
	//cout<<endl;
	setPhotonNucleusSigma(0.01*int_r);
	#ifdef _makeGammaPQ2_
	makeGammaPQ2dependence();
	#endif
}


//______________________________________________________________________________
void
e_narrowResonanceCrossSection::makeGammaPQ2dependence()
{
	// This subroutine calculates the Q2 dependence of 
        // gamma+X -> VM + X cross section for a narrow resonance
  
  /*int const nQ2bins = 19;
	double const q2Edge[nQ2bins+1] = { 0.1,1.,2.,3., 4.,5.,
					   6.,7.,8.,9.,10.,
					   11.,12.,13.,14.,15.,
					   20.,30.,40.,50.};
  */
        int const nQ2bins = 38;
	double const q2Edge[nQ2bins+1] = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
					   0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
					   1, 2, 3, 4, 5, 6, 7, 8, 9, 
					   10, 20, 30, 40, 50, 60, 70, 80, 90, 
					   100, 200 };					  
	//
	double full_x_section[nQ2bins] = {0};
	double effective_flux[nQ2bins] = {0};
	//
        ofstream gamma_flux,total_xsec;
	//
	gamma_flux.open("estarlight_gamma_flux.csv");
	total_xsec.open("estarlight_total_xsec.csv");
	//
        double dEgamma, minEgamma;
	double ega[3] = {0};
	double dR;
	double dR2;
	int    iEgamma, nEgamma,beam;

	//Integration is done with exponential steps, in target frame
	//nEgamma = _VMnumEgamma;
	nEgamma = 1000;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy)/nEgamma;
	minEgamma = std::log(_targetMinPhotonEnergy);
  
	//cout<<" Using Narrow Resonance ..."<<endl;
	
        //printf(" gamma+nucleon threshold (CMS): %e GeV \n", _cmsMinPhotonEnergy);

        int A_1 = getbbs().electronBeam().A(); 
        int A_2 = getbbs().targetBeam().A();

	if( A_2 == 0 && A_1 >= 1 ){
	  // pA, first beam is the nucleus and is in this case the target  
	  beam = 1; 
	} else if( A_1 ==0 && A_2 >= 1){       
	  // photon energy in Target frame
	  beam = 2; 
	} else {
	  // photon energy in Target frame
	  beam = 2; 
	}
        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	int nQ2 = 500;
        //printf(" gamma+nucleon threshold (Target): %e GeV \n", _targetMinPhotonEnergy);
	for( int iQ2Bin  = 0 ; iQ2Bin < nQ2bins; ++iQ2Bin){
	  for(iEgamma = 0 ; iEgamma < nEgamma; ++iEgamma){    // Integral over photon energy
	    // Target frame energies
	    ega[0] = exp(minEgamma + iEgamma*dEgamma );
	    ega[1] = exp(minEgamma + (iEgamma+1)*dEgamma );
	    ega[2] = 0.5*(ega[0]+ega[1]);
	    //
	    //
	    // Integral over Q2		
	    double dndE[3] = {0}; // Effective photon flux
	    double full_int[3] = {0}; // Full e+X --> e+X+V.M. cross section
	    //
	    for( int iEgaInt = 0 ; iEgaInt < 3; ++iEgaInt){    // Loop over the energies for the three points to integrate over Q2	   
	      //
	      double Q2_min = q2Edge[iQ2Bin];
	      double Q2_max = q2Edge[iQ2Bin+1];
	      double lnQ2ratio = std::log(Q2_max/Q2_min)/nQ2;
	      double lnQ2_min = std::log(Q2_min);
	      for( int iQ2 = 0 ; iQ2 < nQ2; ++iQ2){     // Integral over photon virtuality
		//
		double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
		double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
		double q2_12 = (q2_2+q2_1)/2.;
		// Integrating the photon flux
		dndE[iEgaInt] +=(q2_2-q2_1)*( photonFlux(ega[iEgaInt],q2_1)
					      +photonFlux(ega[iEgaInt],q2_2)
					      +4.*photonFlux(ega[iEgaInt],q2_12) );

		full_int[iEgaInt] += (q2_2-q2_1)*( g(ega[iEgaInt],q2_1)*getcsgA(ega[iEgaInt],q2_1,beam)
						   + g(ega[iEgaInt],q2_2)*getcsgA(ega[iEgaInt],q2_2,beam)
						   + 4.*g(ega[iEgaInt],q2_12)*getcsgA(ega[iEgaInt],q2_12,beam) );
	      }
	      // Finish the Q2 integration for the three end-points (Siumpsons rule)
	      dndE[iEgaInt] = dndE[iEgaInt]/6.; 
	      full_int[iEgaInt] = full_int[iEgaInt]/6.;
	    }	    
	    // Finishing cross-section integral 
	    dR = full_int[0];
	    dR += full_int[1];
	    dR += 4.*full_int[2];
	    dR = dR*(ega[1]-ega[0])/6.;
	    
	    // Finishing integral over the effective photon flux
	    dR2 = dndE[0];
	    dR2 += dndE[1];
	    dR2 += 4.*dndE[2];
	    dR2 = dR2*(ega[1]-ega[0])/6.;
	    //
	    full_x_section[iQ2Bin] += dR;
	    effective_flux[iQ2Bin] += dR2;	      
	  }
	  //cout<<gamma_x_section[iQ2Bin]/effective_flux[iQ2Bin]*1E7<<endl;
	  gamma_flux<<q2Edge[iQ2Bin]<<","<<q2Edge[iQ2Bin+1]<<","<<effective_flux[iQ2Bin]<<endl;
	  total_xsec<<q2Edge[iQ2Bin]<<","<<q2Edge[iQ2Bin+1]<<","<<full_x_section[iQ2Bin]<<endl; //No need to remove bin width
	}
	//
	//
	gamma_flux.close();
	total_xsec.close();
}

void e_narrowResonanceCrossSection::printCrossSection(const string name, const double x_section)
{
  if (0.01*x_section > 1.){
    cout<< name.c_str() <<0.01*x_section<<" barn."<<endl;
  } else if (10.*x_section > 1.){
    cout<< name.c_str() <<10.*x_section<<" mb."<<endl;
  } else if (10000.*x_section > 1.){
    cout<< name.c_str() <<10000.*x_section<<" microb."<<endl;
  } else if (10000000.*x_section > 1.){
    cout<< name.c_str() <<10000000.*x_section<<" nanob."<<endl;
  } else if (1.E10*x_section > 1.){
    cout<< name.c_str() <<1.E10*x_section<<" picob."<<endl;
  } else {
    cout<< name.c_str() <<1.E13*x_section<<" femtob."<<endl;
  }
  //cout<<endl;
}  
