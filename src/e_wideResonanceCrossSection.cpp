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
// $Author:: mlomnitz                   $: author of last commit
// $Date:: 2017-03-14 21:05:12 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////
//#define _makeGammaPQ2_

#include <iostream>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "e_wideResonanceCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
e_wideResonanceCrossSection::e_wideResonanceCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem&  bbsystem)
	: photonNucleusCrossSection(inputParametersInstance, bbsystem)//hrm
{
	_wideWmax = _wMax;
	_wideWmin = _wMin;
	_targetMaxPhotonEnergy=inputParametersInstance.targetMaxPhotonEnergy();
	_targetMinPhotonEnergy=inputParametersInstance.targetMinPhotonEnergy();
	_Ep       = inputParametersInstance.protonEnergy();
	_electronEnergy = inputParametersInstance.electronEnergy();
	_target_beamLorentz = inputParametersInstance.beamLorentzGamma();
	_VMnumEgamma = inputParametersInstance.nmbEnergyBins();
	_useFixedRange = inputParametersInstance.fixedQ2Range();
	_gammaMinQ2 = inputParametersInstance.minGammaQ2();
	_gammaMaxQ2 = inputParametersInstance.maxGammaQ2();
	_targetRadii = inputParametersInstance.targetRadius();
}


//______________________________________________________________________________
e_wideResonanceCrossSection::~e_wideResonanceCrossSection()
{

}


//______________________________________________________________________________
void
e_wideResonanceCrossSection::crossSectionCalculation(const double bwnormsave)
{
	//     This subroutine calculates the cross-section assuming a wide
	//     (Breit-Wigner) resonance.

  double W,dW, dEgamma, minEgamma;
	double ega[3] = {0};
	double int_r,dR;
	double int_r2,dR2;
	int    iW,nW,iEgamma,nEgamma,beam;

	double bwnorm = bwnormsave; //used to transfer the bwnorm from the luminosity tables
	
	// For W integration           
	nW   = 100;
	dW   = (_wideWmax-_wideWmin)/double(nW);
	// For Egamma integration
	nEgamma = 1000;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy)/nEgamma;
	minEgamma = std::log(_targetMinPhotonEnergy);
	
  
	if (getBNORM()  ==  0.){
		cout<<" Using Breit-Wigner Resonance Profile."<<endl;
	}
	else{
		cout<<" Using Breit-Wigner plus direct pi+pi- profile."<<endl;
	}
  
	cout<<" Integrating over W from "<<_wideWmin<<" to "<<_wideWmax<<endl;

        int A_1 = getbbs().electronBeam().A(); 
        int A_2 = getbbs().targetBeam().A();

	if( A_2 == 0 && A_1 >= 1 ){
          // eA, first beam is the nucleus and is in this case the target
          beam = 1;
        } else if( A_1 ==0 && A_2 >= 1){
	  // eA, second beam is the nucleus and is in this case the target
          beam = 2;
        } else {
          beam = 2;
        }

	int_r=0.;
	int_r2 = 0.;
        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Integration done using Simpson's rule
	for(iW=0;iW<=nW-1;iW++){
    
		W = _wideWmin + double(iW)*dW + 0.5*dW;
		int nQ2 = 1000;
		for(iEgamma = 0 ; iEgamma < nEgamma; ++iEgamma){    // Integral over photon energy
		  // Displaying the percentage progress
		  float ratio = float(iW)/float(nW)  + 1/float(nW)*float(iEgamma)/float(nEgamma);
		 printf("calculating cross section :%3.2f %%\r",float(ratio *100.0));
		  
		  
		  // Target frame photon energies
		  ega[0] = exp(minEgamma + iEgamma*dEgamma );
		  ega[1] = exp(minEgamma + (iEgamma+1)*dEgamma );
		  ega[2] = 0.5*(ega[0]+ega[1]);
		  // Integral over Q2				  
		  double full_int[3] = {0}; // Full e+X --> e+X+V.M. cross section
		  double dndE[3] = {0}; // Full e+X --> e+X+V.M. cross section
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
		      //			
		      // Integrating cross section
		      full_int[iEgaInt] += (q2_2-q2_1)*( g(ega[iEgaInt],q2_1)*getcsgA(ega[iEgaInt],q2_1,beam)
							 + g(ega[iEgaInt],q2_2)*getcsgA(ega[iEgaInt],q2_2,beam)
							 + 4.*g(ega[iEgaInt],q2_12)*getcsgA(ega[iEgaInt],q2_12,beam) );	      
		      // Effective flux
		      dndE[iEgaInt] +=(q2_2-q2_1)*( getcsgA_Q2_dep(q2_1)*photonFlux(ega[iEgaInt],q2_1)
						    +getcsgA_Q2_dep(q2_2)*photonFlux(ega[iEgaInt],q2_2)
						    +4.*getcsgA_Q2_dep(q2_12)*photonFlux(ega[iEgaInt],q2_12) );
		    }
		    full_int[iEgaInt] = full_int[iEgaInt]/6.;
		    dndE[iEgaInt] = dndE[iEgaInt]/6.;
		  }
		  // Finishing cross-section integral 
		  dR = full_int[0];
		  dR += full_int[1];
		  dR += 4.*full_int[2];
		  dR = dR*(ega[1]-ega[0])/6.;
		  int_r = int_r + dR*breitWigner(W,bwnorm)*dW;
		  // Finishing effective flux integral
		  	  // Finishing integral over the effective photon flux
		  dR2 = dndE[0];
		  dR2 += dndE[1];
		  dR2 += 4.*dndE[2];
		  dR2 = dR2*(ega[1]-ega[0])/6.;
		  int_r2 = int_r2 + dR2*breitWigner(W,bwnorm)*dW;
		}
	}
	cout<<endl;
	if(_useFixedRange == true){
	  cout<<" Using fixed Q2 range "<<_gammaMinQ2 << " < Q2 < "<<_gammaMaxQ2<<endl;
	}
	printCrossSection(" Total cross section: ",int_r);
	//printCrossSection(" gamma+X --> VM+X ", int_r/int_r2); 
	setPhotonNucleusSigma(0.01*int_r);
	//
#ifdef _makeGammaPQ2_
	makeGammaPQ2dependence(bwnormsave);
#endif
}

#ifdef _makeGammaPQ2_
//______________________________________________________________________________
void
e_wideResonanceCrossSection::makeGammaPQ2dependence( double bwnormsave)
{
	// This subroutine calculates the Q2 dependence of 
        // gamma+X -> VM + X cross section for a narrow resonance
  
        int const nQ2bins = 19;
	double const q2Edge[nQ2bins+1] = { 0.,1.,2.,3., 4.,5.,
					   6.,7.,8.,9.,10.,
					   11.,12.,13.,14.,15.,
					   20.,30.,40.,50.};
	double W = 0,dW,dY;
	double y1,y2,y12,ega1,ega2,ega12;
	double int_r,dR,dR2;
	double csgA1, csgA2, csgA12;
	double Eth;
	int    I,J,NW,NY,beam;

	double bwnorm = bwnormsave; //used to transfer the bwnorm from the luminosity tables
                   
	NW   = 100;
	dW   = (_wideWmax-_wideWmin)/double(NW);
  
	if (getBNORM()  ==  0.){
		cout<<" Using Breit-Wigner Resonance Profile."<<endl;
	}
	else{
		cout<<" Using Breit-Wigner plus direct pi+pi- profile."<<endl;
	}
  
	cout<<" Integrating over W from "<<_wideWmin<<" to "<<_wideWmax<<endl;

	//Lomnitz old used for XX
	Eth=0.5*(((W+protonMass)*(W+protonMass)-
	          protonMass*protonMass)/(_Ep+sqrt(_Ep*_Ep-protonMass*protonMass)));
	// Adapted for eX
	//Eth=0.5*(((W+starlightConstants::mel)*(W +starlightConstants::mel)-
	//	  starlightConstants::mel*starlightConstants::mel)/(_electronEnergy+sqrt(_electronEnergy*_electronEnergy-starlightConstants::mel*starlightConstants::mel))); 

        printf(" gamma+nucleon threshold: %e GeV \n", Eth);

        int A_1 = getbbs().electronBeam().A(); 
        int A_2 = getbbs().targetBeam().A();
  
	
        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	///
	cout<<" Lomnitz debug :: sigma_gamma_p --> VM_p "<<endl;
	cout<<" Q2+MV2 \t \t"<<" sigma_gamma_p --> VM_p (nanob)"<<endl;
	double target_cm = acosh(_target_beamLorentz);
	// another - sign from subraction in addition rule
	double exp_target_cm = exp(-target_cm);
	double int_r2;
	for( int iQ2 = 0 ; iQ2 < nQ2bins; ++iQ2){
	  int_r=0.;
	  int_r2=0.;
	  //
	  double q2_cor = getcsgA_Q2_dep( (q2Edge[iQ2+1] + q2Edge[iQ2])/2. );
	  for(I=0;I<=NW-1;I++){
	    
	    W = _wideWmin + double(I)*dW + 0.5*dW;
	    for(J=0;J<=(NY-1);J++){
	      y1  = _wideYmin + double(J)*dY;
	      y2  = _wideYmin + double(J+1)*dY;
	      y12 = 0.5*(y1+y2);  
	      double target_ega1, target_ega2, target_ega12;
	      if( A_2 == 0 && A_1 >= 1 ){
		// pA, first beam is the nucleus and is in this case the target  
		ega1  = 0.5*W*exp(-y1);
		ega2  = 0.5*W*exp(-y2);
		ega12 = 0.5*W*exp(-y12);
		beam = 1; 
	      } else if( A_1 ==0 && A_2 >= 1){
		// pA, second beam is the nucleus and is in this case the target 
		ega1  = 0.5*W*exp(y1);
		ega2  = 0.5*W*exp(y2);
		ega12 = 0.5*W*exp(y12);
		// photon energy in Target frame
		beam = 2; 
	      } else {
		ega1  = 0.5*W*exp(y1);
		ega2  = 0.5*W*exp(y2);
		ega12 = 0.5*W*exp(y12);
		// photon energy in Target frame
		beam = 2; 
	      }
	      //
	      if(ega1 < Eth || ega2 < Eth)   
		continue;
	      if(ega2 > maxPhotonEnergy() || ega1 > maxPhotonEnergy() ) 
		continue;
	      target_ega1 = ega1*exp_target_cm;
	      target_ega12 = ega12*exp_target_cm;
	      target_ega2 = ega2*exp_target_cm;
	      //cout<<"Nortmalizations "<<integrated_x_section(ega1,0,50)<<endl;
	      //		
	      csgA1=getcsgA(ega1,W,beam);
	      double full_range_1 = integrated_x_section(target_ega1);
	      //         >> Middle Point                      =====>>>
	      csgA12=getcsgA(ega12,W,beam);         
	      double full_range_12 = integrated_x_section(target_ega12);
	      //         >> Second Point                      =====>>>
	      csgA2=getcsgA(ega2,W,beam);
	      double full_range_2 = integrated_x_section(target_ega2);
	      //
	      
	      //>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
	      dR  = q2_cor*csgA1;
	      dR  = dR + 4.*q2_cor*csgA12;
	      dR  = dR + q2_cor*csgA2;
	      dR  = dR*(dY/6.)*breitWigner(W,bwnorm)*dW;
	      //
	      dR2  = full_range_1*csgA1;
	      dR2  = dR2 + 4.*full_range_12*csgA12;
	      dR2  = dR2 + full_range_2*csgA2;
	      dR2  = dR2*(dY/6.)*breitWigner(W,bwnorm)*dW;
	      //
	      int_r = int_r+dR;
	      int_r2 = int_r2 +dR2; 
	    }
	  //
	  }
	  if( iQ2 ==0 )
	    cout<<"Full range "<<int_r2*10000000<<endl;
	  cout<<(q2Edge[iQ2+1]+q2Edge[iQ2])/2.+W*W<<" ,  "<<10000000.*int_r<<endl;
	}
	
}
#endif
void e_wideResonanceCrossSection::printCrossSection(const string name, const double x_section)
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
