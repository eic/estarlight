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
// $Rev:: 45                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2011-02-27 20:52:35 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "incoherentVMCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
incoherentVMCrossSection::incoherentVMCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem&  bbsystem)
	:photonNucleusCrossSection(inputParametersInstance, bbsystem)
{
	_narrowYmax = inputParametersInstance.maxRapidity();
	_narrowYmin = -1.0*_narrowYmax;
	_narrowNumY = inputParametersInstance.nmbRapidityBins();
	_Ep         = inputParametersInstance.protonEnergy();	
}


//______________________________________________________________________________
incoherentVMCrossSection::~incoherentVMCrossSection()
{ }


//______________________________________________________________________________
void
incoherentVMCrossSection::crossSectionCalculation(const double)  // _bwnormsave (unused)
{
	// This subroutine calculates the vector meson cross section assuming
	// a narrow resonance.  For reference, see STAR Note 386.
  
	double W,dY;
	double y1,y2,y12,ega1,ega2,ega12;
	double csgA1,csgA2,csgA12,int_r,dR;
        double Wgp,csVN,csVA; 
	double Eth;
	int    J,NY,beam;
  
	NY   =  _narrowNumY;
	dY   = (_narrowYmax-_narrowYmin)/double(NY);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	W = getChannelMass();
	Eth=0.5*(((W+protonMass)*(W+protonMass)-
	          protonMass*protonMass)/(_Ep+sqrt(_Ep*_Ep-protonMass*protonMass)));
  
	// cout<<" gamma+nucleon  Threshold: "<<Eth<<endl;
        printf(" gamma+nucleon threshold: %e GeV \n", Eth);

        int A_1 = getbbs().electronBeam().A(); 
        int A_2 = getbbs().targetBeam().A();

	int_r=0.;

        // Do this first for the case when the first beam is the photon emitter 
        // Treat pA separately with defined beams 
        // The variable beam (=1,2) defines which nucleus is the target 
	for(J=0;J<=(NY-1);J++){
    
	        // This is the fdefault
		y1  = _narrowYmin + double(J)*dY;
		y2  = _narrowYmin + double(J+1)*dY;
		y12 = 0.5*(y1+y2);
    
                if( A_2 == 1 && A_1 != 1 ){
                  // pA, first beam is the nucleus and photon emitter
 		  ega1  = 0.5*W*exp(y1);
		  ega2  = 0.5*W*exp(y2);
		  ega12 = 0.5*W*exp(y12);
                  beam = 2; 
                } else if( A_1 ==1 && A_2 != 1){
                  // pA, second beam is the nucleus and photon emitter
		  ega1  = 0.5*W*exp(-y1);
		  ega2  = 0.5*W*exp(-y2);
		  ega12 = 0.5*W*exp(-y12);
                  beam = 1; 
                } else {
		  ega1  = 0.5*W*exp(y1);
		  ega2  = 0.5*W*exp(y2);
		  ega12 = 0.5*W*exp(y12);
                  beam = 2; 
                }

                // This is for checking things in the lab frame 
                // y1lab  = _narrowYmin + double(J)*dY;
                // y2lab  = _narrowYmin + double(J+1)*dY;
                // y12lab = 0.5*(y1lab+y2lab);

                // p+Pb
                // y1 = y1lab + 0.465;
                // y2 = y2lab + 0.465;
                // y12 = y12lab + 0.465; 
                // ega1  = 0.5*W*exp(y1);
                // ega2  = 0.5*W*exp(y2);
                // ega12 = 0.5*W*exp(y12);

                // Pb+p
                // y1 = y1lab - 0.465;
                // y2 = y2lab - 0.465;
                // y12 = y12lab - 0.465; 
                // ega1  = 0.5*W*exp(-y1);
                // ega2  = 0.5*W*exp(-y2);
                // ega12 = 0.5*W*exp(-y12);

    
		if(ega1 < Eth)   
			continue;
		if(ega2 > maxPhotonEnergy()) 
			continue;

		// First point 
                Wgp = sqrt(2.*ega1*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
			          +starlightConstants::protonMass*starlightConstants::protonMass);
                csVN = sigma_N(Wgp);            
                csVA = sigma_A(csVN,beam); 
                csgA1 = (csVA/csVN)*sigmagp(Wgp); 
                if( getbbs().electronBeam().A() == 1 || getbbs().targetBeam().A()==1 ){
                  csgA1 = sigmagp(Wgp);
                }

		// Middle point 
                Wgp = sqrt(2.*ega12*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
			          +starlightConstants::protonMass*starlightConstants::protonMass);
                csVN = sigma_N(Wgp);            
                csVA = sigma_A(csVN,beam); 
                csgA12 = (csVA/csVN)*sigmagp(Wgp); 
                if( getbbs().electronBeam().A() == 1 || getbbs().targetBeam().A()==1 ){
                  csgA12 = sigmagp(Wgp);
                }

		// Last point 
                Wgp = sqrt(2.*ega2*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
			          +starlightConstants::protonMass*starlightConstants::protonMass);
                csVN = sigma_N(Wgp);            
                csVA = sigma_A(csVN,beam); 
                csgA2 = (csVA/csVN)*sigmagp(Wgp); 
                if( getbbs().electronBeam().A() == 1 || getbbs().targetBeam().A()==1 ){
                  csgA2 = sigmagp(Wgp);
                }

                dR = ega1*photonFlux(ega1,beam)*csgA1;  
                dR = dR + 4*ega12*photonFlux(ega12,beam)*csgA12;
                dR = dR + ega2*photonFlux(ega2,beam)*csgA2; 
                dR = dR*(dY/6.); 

		// cout<<" y: "<<y12<<" egamma: "<<ega12<<" flux: "<<photonFlux(ega12)<<" sigma_gA: "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*dR/dY<<endl;
                // cout<<" y: "<<y12lab<<" egamma: "<<ega12<<" flux: "<<ega12*photonFlux(ega12)<<" W: "<<Wgpm<<" Wflux: "<<Wgpm*photonFlux(ega12)<<" sigma_gA (nb): "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*ega12*photonFlux(ega12)*csgA12<<endl;

		int_r = int_r+dR;

	}

        // Repeat the loop for the case when the second beam is the photon emitter. 
        // Don't repeat for pA
        if( !( (A_2 == 1 && A_1 != 1) || (A_1 == 1 && A_2 != 1) ) ){ 
	  for(J=0;J<=(NY-1);J++){
    
	        // This is the fdefault
		y1  = _narrowYmin + double(J)*dY;
		y2  = _narrowYmin + double(J+1)*dY;
		y12 = 0.5*(y1+y2);
    
                beam = 1; 
		ega1  = 0.5*W*exp(-y1);
		ega2  = 0.5*W*exp(-y2);
		ega12 = 0.5*W*exp(-y12);

                // This is for checking things in the lab frame 
                // y1lab  = _narrowYmin + double(J)*dY;
                // y2lab  = _narrowYmin + double(J+1)*dY;
                // y12lab = 0.5*(y1lab+y2lab);

                // p+Pb
                // y1 = y1lab + 0.465;
                // y2 = y2lab + 0.465;
                // y12 = y12lab + 0.465; 
                // ega1  = 0.5*W*exp(y1);
                // ega2  = 0.5*W*exp(y2);
                // ega12 = 0.5*W*exp(y12);

                // Pb+p
                // y1 = y1lab - 0.465;
                // y2 = y2lab - 0.465;
                // y12 = y12lab - 0.465; 
                // ega1  = 0.5*W*exp(-y1);
                // ega2  = 0.5*W*exp(-y2);
                // ega12 = 0.5*W*exp(-y12);

    
		if(ega2 < Eth)   
			continue;
		if(ega1 > maxPhotonEnergy()) 
			continue;

		// First point 
                Wgp = sqrt(2.*ega1*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
			          +starlightConstants::protonMass*starlightConstants::protonMass);
                csVN = sigma_N(Wgp);            
                csVA = sigma_A(csVN,beam); 
                csgA1 = (csVA/csVN)*sigmagp(Wgp); 
                if( getbbs().electronBeam().A() == 1 || getbbs().targetBeam().A()==1 ){
                  csgA1 = sigmagp(Wgp);
                }

		// Middle point 
                Wgp = sqrt(2.*ega12*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
			          +starlightConstants::protonMass*starlightConstants::protonMass);
                csVN = sigma_N(Wgp);            
                csVA = sigma_A(csVN,beam); 
                csgA12 = (csVA/csVN)*sigmagp(Wgp); 
                if( getbbs().electronBeam().A() == 1 || getbbs().targetBeam().A()==1 ){
                  csgA12 = sigmagp(Wgp);
                }

		// Last point 
                Wgp = sqrt(2.*ega2*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
			          +starlightConstants::protonMass*starlightConstants::protonMass);
                csVN = sigma_N(Wgp);            
                csVA = sigma_A(csVN,beam); 
                csgA2 = (csVA/csVN)*sigmagp(Wgp); 
                if( getbbs().electronBeam().A() == 1 || getbbs().targetBeam().A()==1 ){
                  csgA2 = sigmagp(Wgp);
                }

                dR = ega1*photonFlux(ega1,beam)*csgA1;  
                dR = dR + 4*ega12*photonFlux(ega12,beam)*csgA12;
                dR = dR + ega2*photonFlux(ega2,beam)*csgA2; 
                dR = dR*(dY/6.); 

		// cout<<" y: "<<y12<<" egamma: "<<ega12<<" flux: "<<photonFlux(ega12)<<" sigma_gA: "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*dR/dY<<endl;
                // cout<<" y: "<<y12lab<<" egamma: "<<ega12<<" flux: "<<ega12*photonFlux(ega12)<<" W: "<<Wgpm<<" Wflux: "<<Wgpm*photonFlux(ega12)<<" sigma_gA (nb): "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*ega12*photonFlux(ega12)*csgA12<<endl;

		int_r = int_r+dR;

	  }
        }

	cout<<endl;
	if (0.01*int_r > 1.){
	  cout<< " Total cross section: "<<0.01*int_r<<" barn."<<endl;
	} else if (10.*int_r > 1.){
	  cout<< " Total cross section: " <<10.*int_r<<" mb."<<endl;
        } else if (10000.*int_r > 1.){
	  cout<< " Total cross section: " <<10000.*int_r<<" microb."<<endl;
        } else if (10000000.*int_r > 1.){
	  cout<< " Total cross section: " <<10000000.*int_r<<" nanob."<<endl;
        } else if (1.E10*int_r > 1.){
	  cout<< " Total cross section: "<<1.E10*int_r<<" picob."<<endl;
        } else {
	  cout<< " Total cross section: " <<1.E13*int_r<<" femtob."<<endl;
        }
	cout<<endl;
	setPhotonNucleusSigma(0.01*int_r);
}
