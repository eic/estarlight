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
// $Rev:: 276                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-09-13 19:54:42 +0100 #$: date of last commit
//
// Description:
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "reportingUtils.h"
#include "starlightconstants.h"
#include "bessel.h"
#include "photonNucleusCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
photonNucleusCrossSection::photonNucleusCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem& bbsystem)
	: _nWbins            (inputParametersInstance.nmbWBins()          ),
	  _nYbins            (inputParametersInstance.nmbRapidityBins()   ),
	  _wMin              (inputParametersInstance.minW()              ),
	  _wMax              (inputParametersInstance.maxW()              ),
	  _yMax              (inputParametersInstance.maxRapidity()       ),
	  _beamLorentzGamma  (inputParametersInstance.beamLorentzGamma()  ),
	  _bbs               (bbsystem                                    ),
	  _protonEnergy      (inputParametersInstance.protonEnergy()      ),
	  _electronEnergy    (inputParametersInstance.electronEnergy()    ),
	  _particleType      (inputParametersInstance.prodParticleType()  ),
	  _beamBreakupMode   (inputParametersInstance.beamBreakupMode()   ),
	  _backwardsProduction(inputParametersInstance.backwardsProduction()),
      _productionMode    (inputParametersInstance.productionMode()    ),
	  _sigmaNucleus      (_bbs.targetBeam().A()          ),
	  _fixedQ2range      (inputParametersInstance.fixedQ2Range()      ),
	  _minQ2             (inputParametersInstance.minGammaQ2()        ),
	  _maxQ2             (inputParametersInstance.maxGammaQ2()        ),
	  _maxPhotonEnergy   (inputParametersInstance.cmsMaxPhotonEnergy()),
	  _cmsMinPhotonEnergy(inputParametersInstance.cmsMinPhotonEnergy()),
	  _targetRadii       (inputParametersInstance.targetRadius()      )
{
        // new options - impulse aproximation (per Joakim) and Quantum Glauber (per SK) SKQG
        _impulseSelected = inputParametersInstance.impulseVM();
	_quantumGlauber = inputParametersInstance.quantumGlauber();
	switch(_particleType) {
	case RHO:
	  _slopeParameter = 11.0;  // [(GeV/c)^{-2}]
	  _vmPhotonCoupling = 2.02;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [ GeV^{-2}]
	  _ANORM       = -2.75;
	  _BNORM       = 0.0;
	  _defaultC    = 1.0;
	  _channelMass = starlightConstants::rho0Mass; 
	  _width       = starlightConstants::rho0Width; 
	  if(_backwardsProduction){
	  	_slopeParameter = 21.8;  // [(GeV/c)^{-2}]
	  	_ANORM          = 1.;
	  }
	  break;
	case RHOZEUS:
	  _slopeParameter =11.0;
	  _vmPhotonCoupling=2.02;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75;
	  _BNORM=1.84;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::rho0Mass;
	  _width        = starlightConstants::rho0Width;
	  break;
	case FOURPRONG:
	  _slopeParameter      = 11.0;
	  _vmPhotonCoupling      = 2.02;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM       = -2.75;
	  _BNORM       = 0;  
	  _defaultC    = 11.0;
	  _channelMass  = starlightConstants::rho0PrimeMass;
	  _width        = starlightConstants::rho0PrimeWidth;
	  break;
	case OMEGA:
	case OMEGA_pipipi:
	  _slopeParameter=10.0;
	  _vmPhotonCoupling=23.13;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75;
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::OmegaMass;
	  _width        = starlightConstants::OmegaWidth;
	  if(_backwardsProduction){
	  	_slopeParameter = 21.8;  // [(GeV/c)^{-2}]
	  	_ANORM          = 1.;
	  }
	case OMEGA_pi0gamma:
	  _slopeParameter=10.0;
	  _vmPhotonCoupling=23.13;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75;
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::OmegaMass;
	  _width        = starlightConstants::OmegaWidth;
	  if(_backwardsProduction){
	  	_slopeParameter = 21.8;  // [(GeV/c)^{-2}]
	  	_ANORM          = 1.;
	  }
	  break;
	case PHI:
	  _slopeParameter=7.0;
	  _vmPhotonCoupling=13.71;
	  _vmQ2Power_c1 = 2.15;
	  _vmQ2Power_c2 = 0.0074; // [GeV^{-2}]
	  _ANORM=-2.75;
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::PhiMass;
	  _width        = starlightConstants::PhiWidth;
	  break;
	case JPSI:
	case JPSI_ee:
	  _slopeParameter=4.0;
	  _vmPhotonCoupling=10.45;
	  _vmQ2Power_c1 = 2.45;
	  _vmQ2Power_c2 = 0.00084; // [GeV^{-2}]
	  _ANORM=-2.75; 
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::JpsiMass; 
	  _width        = starlightConstants::JpsiWidth; 
	  break;
	case JPSI_mumu:
	  _slopeParameter=4.0;
	  _vmPhotonCoupling=10.45;
	  _vmQ2Power_c1 = 2.36;
	  _vmQ2Power_c2 = 0.0029; // [GeV^{-2}]
	  _ANORM=-2.75; 
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::JpsiMass; 
	  _width        = starlightConstants::JpsiWidth; 
	  break;
	case JPSI2S:
	case JPSI2S_ee:
	case JPSI2S_mumu:
	  _slopeParameter=4.3;
	  _vmPhotonCoupling=26.39;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75; 
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::Psi2SMass;
	  _width        = starlightConstants::Psi2SWidth;
	  break;
	case UPSILON:
	case UPSILON_ee:
	case UPSILON_mumu:
	  _slopeParameter=4.0;
	  _vmPhotonCoupling=125.37;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75; 
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::Upsilon1SMass;
	  _width        = starlightConstants::Upsilon1SWidth;
	  break;
	case UPSILON2S:
	case UPSILON2S_ee:
	case UPSILON2S_mumu:
	  _slopeParameter=4.0;
	  _vmPhotonCoupling=290.84;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75;
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::Upsilon2SMass;
	  _width        = starlightConstants::Upsilon2SWidth;		
	  break;
	case UPSILON3S:
	case UPSILON3S_ee:
	case UPSILON3S_mumu:
	  _slopeParameter=4.0;
	  _vmPhotonCoupling=415.10;
	  _vmQ2Power_c1 = 2.09;
	  _vmQ2Power_c2 = 0.0073; // [GeV^{-2}]
	  _ANORM=-2.75;
	  _BNORM=0.0;
	  _defaultC=1.0;
	  _channelMass  = starlightConstants::Upsilon3SMass;
	  _width        = starlightConstants::Upsilon3SWidth;
	  break;
	default:
		cout <<"No sigma constants parameterized for pid: "<<_particleType
		     <<" GammaAcrosssection"<<endl;
	}

	//Changed by Lomnitz for e case. Limit is now E_e - 100m_e
	//_maxPhotonEnergy = 12. * _beamLorentzGamma * hbarc/(_bbs.beam1().nuclearRadius()+_bbs.targetBeam().nuclearRadius());
	//_maxPhotonEnergy = _electronEnergy - 10.*starlightConstants::mel;
	/*cout<<" Lomnitz:: max energy in target frame "<< _electronEnergy - 1000.*starlightConstants::mel<<" vs electron energy "<<_electronEnergy<<endl
	    <<"           max energy in cms frame    "<<_maxPhotonEnergy<<"  vs electron energy "<<_beamLorentzGamma*starlightConstants::mel<<endl;
	    cout<<" testing original limit "<< 12. * _beamLorentzGamma * hbarc/(2.*_bbs.targetBeam().nuclearRadius())<<endl;*/
	
	  
}


//______________________________________________________________________________
photonNucleusCrossSection::~photonNucleusCrossSection()
{ }


//______________________________________________________________________________
void
photonNucleusCrossSection::crossSectionCalculation(const double)
{
	cout << "Neither narrow/wide resonance cross-section calculation.--Derived" << endl;
}

//______________________________________________________________________________
double
photonNucleusCrossSection::getcsgA(const double targetEgamma,
                                   const double Q2, 
                                   const int beam)
{
	//This function returns the cross-section for photon-nucleus interaction 
	//producing vectormesons
  
	double Av,Wgp,cs,cvma;
	double t,tmin,tmax;
	double csgA,ax,bx;
	int NGAUSS; 
	//
	double W = _channelMass; //new change, channel mass center used for the t min integration.
	  
	//     DATA FOR GAUSS INTEGRATION
	double xg[6] = {0, 0.1488743390, 0.4333953941, 0.6794095683, 0.8650633667, 0.9739065285};
	double ag[6] = {0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443};
	NGAUSS = 6;
	
	//       Note: The photon energy passed to this function is now in the target frame. The rest of the calculations are done in the
	//       CMS frame. The next lines boost the photon into the CM frame.
	double E_prime = _electronEnergy - targetEgamma;
	double cos_theta_e = 1. - Q2/(2.*_electronEnergy*E_prime);
	double theta_e = acos(cos_theta_e);
	double beam_y = acosh(_beamLorentzGamma);	
	double gamma_pt = E_prime*sin(theta_e);
	double pz_squared = targetEgamma*targetEgamma - Q2 - gamma_pt*gamma_pt;
	if( pz_squared < 0 || fabs(cos_theta_e) > 1 || 2.*targetEgamma/(Q2+W*W) < _targetRadii)
	  return 0;
	double temp_pz = sqrt(pz_squared);
	// Now boost to CM frame
	double Egamma = targetEgamma*cosh(beam_y) - temp_pz*sinh(beam_y);
	if( Egamma < _cmsMinPhotonEnergy || Egamma > _maxPhotonEnergy){
	  return 0;
	}
	//cout<<" ::: Lomnitz test in photonNucleus ::: pz^2 = "<<pz_squared << " CMS Egamma = "<<Egamma<<endl;
	//       Find gamma-proton CM energy in CMS frame in the limit Q2->0 (this is our assumption, the Q2 dependence is in the factor)
	Wgp = sqrt( 2.*(protonMass*targetEgamma)+protonMass*protonMass);
	/*Wgp = sqrt(2. * Egamma * (_protonEnergy
	                          + sqrt(_protonEnergy * _protonEnergy - protonMass * protonMass))
				  + protonMass * protonMass);*/
	
	//Used for A-A
	tmin = (W * W / (4. * Egamma * _beamLorentzGamma)) * (W * W / (4. * Egamma * _beamLorentzGamma));
  
	if ((_bbs.electronBeam().A() <= 1) && (_bbs.targetBeam().A() <= 1)){
	   // proton-proton, no scaling needed
	  csgA = getcsgA_Q2_dep(Q2)*sigmagp(Wgp);
	} else {
	   // Check if one or both beams are nuclei 
	   int A_1 = _bbs.electronBeam().A(); 
	   int A_2 = _bbs.targetBeam().A(); 
	   // coherent AA interactions
	   // Calculate V.M.+proton cross section
           // cs = sqrt(16. * pi * _vmPhotonCoupling * _slopeParameter * hbarc * hbarc * sigmagp(Wgp) / alpha); 
           cs = getcsgA_Q2_dep(Q2)*sigma_N(Wgp); //Use member function instead 
    
	   // Calculate V.M.+nucleus cross section
	   cvma = sigma_A(cs,beam); 

           // Do impulse approximation here
           if( _impulseSelected == 1){
             if( beam == 1 ){
	       cvma = A_1*cs;
	     } else if ( beam == 2 ){
               cvma = A_2*cs;
	     }   
           }	   

	   // Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
	   Av = (alpha * cvma * cvma) / (16. * pi * _vmPhotonCoupling * hbarc * hbarc);
   
	   tmax   = tmin + 0.25;
	   ax     = 0.5 * (tmax - tmin);
	   bx     = 0.5 * (tmax + tmin);
	   csgA   = 0.;
	   for (int k = 1; k < NGAUSS; ++k) { 

	       t    = ax * xg[k] + bx;
               if( A_1 <= 1 && A_2 != 1){ 
		  csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);
               }else if(A_2 <=1 && A_1 != 1){
		  csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
               }else{     
                  if( beam==1 ){
 		     csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
                  }else if(beam==2){
 		     csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);	
  		  }else{
		     cout<<"Something went wrong here, beam= "<<beam<<endl; 
                  }
               }

	       t    = ax * (-xg[k]) + bx;
               if( A_1 <= 1 && A_2 != 1){ 
			  csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);
               }else if(A_2 <=1 && A_1 != 1){
			  csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
               }else{     
                  if( beam==1 ){
 			    csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
                  }else if(beam==2){
 			    csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);	
  		  }else{
			    cout<<"Something went wrong here, beam= "<<beam<<endl; 
                  }
	       }
	   }
	   csgA = 0.5 * (tmax - tmin) * csgA;
	   csgA = Av * csgA;
	}
	return csgA;	
}


//______________________________________________________________________________
double
photonNucleusCrossSection::e_getcsgA(const double Egamma, double Q2,
                                   const double W, 
                                   const int beam)
{
	//This function returns the cross-section for photon-nucleus interaction 
	//producing vectormesons for e_starlight. Stuff will be returned in CMS frame, but
        //photon input is take in target frame
  
	double Av,Wgp,cs,cvma;
	double t,tmin,tmax;
	double csgA,ax,bx;
	int NGAUSS; 
  
	//     DATA FOR GAUSS INTEGRATION
	double xg[6] = {0, 0.1488743390, 0.4333953941, 0.6794095683, 0.8650633667, 0.9739065285};
	double ag[6] = {0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443};
	NGAUSS = 6;
  
	//       Find gamma-proton CM energy
	Wgp = sqrt( 2.*(protonMass*Egamma)
		    +protonMass*protonMass + Q2);
	
	//Used for A-A
	tmin = (W * W / (4. * Egamma * _beamLorentzGamma)) * (W * W / (4. * Egamma * _beamLorentzGamma));
  
	if ((_bbs.electronBeam().A() <= 1) && (_bbs.targetBeam().A() <= 1)){
	   // proton-proton, no scaling needed
	   csgA = sigmagp(Wgp);
	} else {
	   // coherent AA interactions
	   // Calculate V.M.+proton cross section
           // cs = sqrt(16. * pi * _vmPhotonCoupling * _slopeParameter * hbarc * hbarc * sigmagp(Wgp) / alpha); 
           cs = sigma_N(Wgp); //Use member function instead 
    
	   // Calculate V.M.+nucleus cross section
	   cvma = sigma_A(cs,beam); 

	   // Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
	   Av = (alpha * cvma * cvma) / (16. * pi * _vmPhotonCoupling * hbarc * hbarc);

           // Check if one or both beams are nuclei 
           int A_1 = _bbs.electronBeam().A(); 
           int A_2 = _bbs.targetBeam().A(); 
   
	   tmax   = tmin + 0.25;
	   ax     = 0.5 * (tmax - tmin);
	   bx     = 0.5 * (tmax + tmin);
	   csgA   = 0.;
	   for (int k = 1; k < NGAUSS; ++k) { 

	       t    = ax * xg[k] + bx;
               if( A_1 <= 1 && A_2 != 1){ 
		  csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);
               }else if(A_2 <=1 && A_1 != 1){
		  csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
               }else{     
                  if( beam==1 ){
 		     csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
                  }else if(beam==2){
 		     csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);	
  		  }else{
		     cout<<"Something went wrong here, beam= "<<beam<<endl; 
                  }
               }

	       t    = ax * (-xg[k]) + bx;
               if( A_1 <= 1 && A_2 != 1){ 
			  csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);
               }else if(A_2 <=1 && A_1 != 1){
			  csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
               }else{     
                  if( beam==1 ){
 			    csgA = csgA + ag[k] * _bbs.electronBeam().formFactor(t) * _bbs.electronBeam().formFactor(t);
                  }else if(beam==2){
 			    csgA = csgA + ag[k] * _bbs.targetBeam().formFactor(t) * _bbs.targetBeam().formFactor(t);	
  		  }else{
			    cout<<"Something went wrong here, beam= "<<beam<<endl; 
                  }
	       }
	   }
	   csgA = 0.5 * (tmax - tmin) * csgA;
	   csgA = Av * csgA;
	}
	return csgA;	
}


//______________________________________________________________________________
double
photonNucleusCrossSection::getcsgA_Q2_dep(const double Q2)
{
  double const mv2 = getChannelMass()*getChannelMass();
  double const n = vmQ2Power(Q2);
  return std::pow(mv2/(mv2+Q2),n);
}


//______________________________________________________________________________
double
photonNucleusCrossSection::photonFlux(const double Egamma, const int beam)
{
	// This routine gives the photon flux as a function of energy Egamma
	// It works for arbitrary nuclei and gamma; the first time it is
	// called, it calculates a lookup table which is used on
	// subsequent calls.
	// It returns dN_gamma/dE (dimensions 1/E), not dI/dE
	// energies are in GeV, in the lab frame
	// rewritten 4/25/2001 by SRK

        // NOTE: beam (=1,2) defines the photon TARGET

	double lEgamma,Emin,Emax;
	static double lnEmax, lnEmin, dlnE;
	double stepmult,energy,rZ;
	int nbstep,nrstep,nphistep,nstep;
	double bmin,bmax,bmult,biter,bold,integratedflux;
	double fluxelement,deltar,riter;
	double deltaphi,phiiter,dist;
	static double dide[401];
	double lnElt;
	double flux_r; 
	double Xvar;
	int Ilt;
	double RNuc=0.,RSum=0.;

	RSum=_bbs.electronBeam().nuclearRadius()+_bbs.targetBeam().nuclearRadius();
        if( beam == 1){
          rZ=double(_bbs.targetBeam().Z());
          RNuc = _bbs.electronBeam().nuclearRadius();
        } else { 
	  rZ=double(_bbs.electronBeam().Z());
          RNuc = _bbs.targetBeam().nuclearRadius();
        }

	static int  Icheck = 0;
	static int  Ibeam  = 0; 
  
	//Check first to see if pp 
	if( _bbs.electronBeam().A()==1 && _bbs.targetBeam().A()==1 ){
		int nbsteps = 400;
		double bmin = 0.5;
		double bmax = 5.0 + (5.0*_beamLorentzGamma*hbarc/Egamma);
		double dlnb = (log(bmax)-log(bmin))/(1.*nbsteps);
 
		double local_sum=0.0;

		// Impact parameter loop 
		for (int i = 0; i<=nbsteps;i++){

			double bnn0 = bmin*exp(i*dlnb);
			double bnn1 = bmin*exp((i+1)*dlnb);
			double db   = bnn1-bnn0;
        
			double ppslope = 19.0; 
			double GammaProfile = exp(-bnn0*bnn0/(2.*hbarc*hbarc*ppslope));  
			double PofB0 = 1. - (1. - GammaProfile)*(1. - GammaProfile);   
			GammaProfile = exp(-bnn1*bnn1/(2.*hbarc*hbarc*ppslope));  
			double PofB1 = 1. - (1. - GammaProfile)*(1. - GammaProfile);   

                        double loc_nofe0 = _bbs.electronBeam().photonDensity(bnn0,Egamma);
			double loc_nofe1 = _bbs.targetBeam().photonDensity(bnn1,Egamma);

			local_sum += 0.5*loc_nofe0*(1. - PofB0)*2.*starlightConstants::pi*bnn0*db; 
			local_sum += 0.5*loc_nofe1*(1. - PofB1)*2.*starlightConstants::pi*bnn1*db; 

		}
		// End Impact parameter loop 
		return local_sum;
	}

	//   first call or new beam?  - initialize - calculate photon flux
	Icheck=Icheck+1;

	// Do the numerical integration only once for symmetric systems. 
        if( Icheck > 1 && _bbs.electronBeam().A() == _bbs.targetBeam().A() && _bbs.electronBeam().Z() == _bbs.targetBeam().Z() ) goto L1000f;
        // For asymmetric systems check if we have another beam 
	if( Icheck > 1 && beam == Ibeam ) goto L1000f; 
        Ibeam = beam; 
  
  	//  Nuclear breakup is done by PofB
	//  collect number of integration steps here, in one place
  
	nbstep=1200;
	nrstep=60;
	nphistep=40;
  
	//  this last one is the number of energy steps
	nstep=100;
 
	// following previous choices, take Emin=10 keV at LHC, Emin = 1 MeV at RHIC
  	Emin=1.E-5;
	if (_beamLorentzGamma < 500) 
		Emin=1.E-3;
  
        Emax=_maxPhotonEnergy; 
	// Emax=12.*hbarc*_beamLorentzGamma/RSum;
 
	//     >> lnEmin <-> ln(Egamma) for the 0th bin
	//     >> lnEmax <-> ln(Egamma) for the last bin
  	lnEmin=log(Emin);
	lnEmax=log(Emax);
	dlnE=(lnEmax-lnEmin)/nstep; 

        printf("Calculating photon flux from Emin = %e GeV to Emax = %e GeV (CM frame) for source with Z = %3.0f \n", Emin, Emax, rZ);
	
	stepmult= exp(log(Emax/Emin)/double(nstep));
	energy=Emin;
  
	for (int j = 1; j<=nstep;j++){
		energy=energy*stepmult;
    
		//  integrate flux over 2R_A < b < 2R_A+ 6* gamma hbar/energy
		//  use exponential steps
    
		bmin=0.8*RSum; //Start slightly below 2*Radius 
		bmax=bmin + 6.*hbarc*_beamLorentzGamma/energy;
    
		bmult=exp(log(bmax/bmin)/double(nbstep));
		biter=bmin;
		integratedflux=0.;

		if( (_bbs.electronBeam().A() == 1 && _bbs.targetBeam().A() != 1) || (_bbs.targetBeam().A() == 1 && _bbs.electronBeam().A() != 1) ){
		    // This is pA 

		  if( _productionMode == PHOTONPOMERONINCOHERENT ){
		      // This pA incoherent, proton is the target

		      int nbsteps = 400;
		      double bmin = 0.7*RSum;
		      double bmax = 2.0*RSum + (8.0*_beamLorentzGamma*hbarc/energy);
		      double dlnb = (log(bmax)-log(bmin))/(1.*nbsteps);

		      double local_sum=0.0;

		      // Impact parameter loop 
		      for (int i = 0; i<=nbsteps; i++){

        		  double bnn0 = bmin*exp(i*dlnb);
			  double bnn1 = bmin*exp((i+1)*dlnb);
			  double db   = bnn1-bnn0;

                          double PofB0 = _bbs.probabilityOfBreakup(bnn0); 
                          double PofB1 = _bbs.probabilityOfBreakup(bnn1); 
      
			  double loc_nofe0 = 0.0;
			  double loc_nofe1 = 0.0; 
                          if( _bbs.electronBeam().A() == 1 ){
			    loc_nofe0 = _bbs.targetBeam().photonDensity(bnn0,energy);
			    loc_nofe1 = _bbs.targetBeam().photonDensity(bnn1,energy);
                          }
		          else if( _bbs.targetBeam().A() == 1 ){
			    loc_nofe0 = _bbs.electronBeam().photonDensity(bnn0,energy);
			    loc_nofe1 = _bbs.electronBeam().photonDensity(bnn1,energy);			    
                          }

                          // cout<<" i: "<<i<<" bnn0: "<<bnn0<<" PofB0: "<<PofB0<<" loc_nofe0: "<<loc_nofe0<<endl; 

			  local_sum += 0.5*loc_nofe0*PofB0*2.*starlightConstants::pi*bnn0*db; 
			  local_sum += 0.5*loc_nofe1*PofB1*2.*starlightConstants::pi*bnn1*db;  
		      }  // End Impact parameter loop 
	              integratedflux = local_sum; 
                    } else if ( _productionMode == PHOTONPOMERONNARROW ||  _productionMode == PHOTONPOMERONWIDE ){
                      // This is pA coherent, nucleus is the target 
                      double localbmin = 0.0;   
                      if( _bbs.electronBeam().A() == 1 ){
			localbmin = _bbs.targetBeam().nuclearRadius() + 0.7; 
		      }
                      if( _bbs.targetBeam().A() == 1 ){ 
			localbmin = _bbs.electronBeam().nuclearRadius() + 0.7; 
                      }
                      integratedflux = nepoint(energy,localbmin); 
		    }
		}else{ 
		// This is AA
		for (int jb = 1; jb<=nbstep;jb++){
		    bold=biter;
		    biter=biter*bmult;
		    // When we get to b>20R_A change methods - just take the photon flux
		    //  at the center of the nucleus.
		    if (biter > (10.*RNuc)){
		       // if there is no nuclear breakup or only hadronic breakup, which only
		       // occurs at smaller b, we can analytically integrate the flux from b~20R_A
		       // to infinity, following Jackson (2nd edition), Eq. 15.54
		       Xvar=energy*biter/(hbarc*_beamLorentzGamma);
		       // Here, there is nuclear breakup.  So, we can't use the integrated flux
		       // However, we can do a single flux calculation, at the center of the nucleus
		       // Eq. 41 of Vidovic, Greiner and Soff, Phys.Rev.C47,2308(1993), among other places
		       // this is the flux per unit area
		       fluxelement  = (rZ*rZ*alpha*energy)*
				       (bessel::dbesk1(Xvar))*(bessel::dbesk1(Xvar))/
				       ((pi*_beamLorentzGamma*hbarc)*
				       (pi*_beamLorentzGamma*hbarc));
	    
                   } else {
		       // integrate over nuclear surface. n.b. this assumes total shadowing -
		       // treat photons hitting the nucleus the same no matter where they strike
		       fluxelement=0.;
		       deltar=RNuc/double(nrstep);
		       riter=-deltar/2.;
          
		       for (int jr =1; jr<=nrstep;jr++){
			   riter=riter+deltar;
			   // use symmetry;  only integrate from 0 to pi (half circle)
			   deltaphi=pi/double(nphistep);
			   phiiter=0.;
            
			   for( int jphi=1;jphi<= nphistep;jphi++){
			       phiiter=(double(jphi)-0.5)*deltaphi;
			       // dist is the distance from the center of the emitting nucleus 
			       // to the point in question
			       dist=sqrt((biter+riter*cos(phiiter))*(biter+riter*
				     cos(phiiter))+(riter*sin(phiiter))*(riter*sin(phiiter)));
			       Xvar=energy*dist/(hbarc*_beamLorentzGamma);  
			       flux_r = (rZ*rZ*alpha*energy)*
				     (bessel::dbesk1(Xvar)*bessel::dbesk1(Xvar))/
				     ((pi*_beamLorentzGamma*hbarc)*
				     (pi*_beamLorentzGamma*hbarc));
	      
				     //  The surface  element is 2.* delta phi* r * delta r
				     //  The '2' is because the phi integral only goes from 0 to pi
				     fluxelement=fluxelement+flux_r*2.*deltaphi*riter*deltar;
				     //  end phi and r integrations
			   }//for(jphi)
		       }//for(jr)
			  //  average fluxelement over the nuclear surface
			  fluxelement=fluxelement/(pi*RNuc*RNuc);
		   }//else
			  //  multiply by volume element to get total flux in the volume element
			  fluxelement=fluxelement*2.*pi*biter*(biter-bold);
			  //  modulate by the probability of nuclear breakup as f(biter)
                          // cout<<" jb: "<<jb<<" biter: "<<biter<<" fluxelement: "<<fluxelement<<endl; 
			  if (_beamBreakupMode > 1){
			      fluxelement=fluxelement*_bbs.probabilityOfBreakup(biter);
			  }
                          // cout<<" jb: "<<jb<<" biter: "<<biter<<" fluxelement: "<<fluxelement<<endl; 
			  integratedflux=integratedflux+fluxelement;
      
		} //end loop over impact parameter 
	    }  //end of else (pp, pA, AA) 
    
	    //  In lookup table, store k*dN/dk because it changes less
	    //  so the interpolation should be better    
	    dide[j]=integratedflux*energy;                                     
	}//end loop over photon energy 
       
	//  for 2nd and subsequent calls, use lookup table immediately
  
 L1000f:
  
	lEgamma=log(Egamma);
	if (lEgamma < (lnEmin+dlnE) ||  lEgamma  > lnEmax){
		flux_r=0.0;
		cout<<"  ERROR: Egamma outside defined range. Egamma= "<<Egamma
		    <<"   "<<lnEmax<<" "<<(lnEmin+dlnE)<<endl;
	}
	else{
		//       >> Egamma between Ilt and Ilt+1
		Ilt = int((lEgamma-lnEmin)/dlnE);
		//       >> ln(Egamma) for first point 
		lnElt = lnEmin + Ilt*dlnE; 
		//       >> Interpolate
		flux_r = dide[Ilt] + ((lEgamma-lnElt)/dlnE)*(dide[Ilt+1]- dide[Ilt]);
		flux_r = flux_r/Egamma;
	}
  
	return flux_r;
}


//______________________________________________________________________________
double 
photonNucleusCrossSection::integrated_Q2_dep(double const Egamma, double const _min, double const _max)
{
  //Integration over full limits gives more accurate result
  double Q2_min =  std::pow(starlightConstants::mel*Egamma,2.0)/(_electronEnergy*(_electronEnergy-Egamma));
  double Q2_max = 4.*_electronEnergy*(_electronEnergy-Egamma);
  //double Q2_max = 2.*Egamma/_targetRadii - _wMax*_wMax;
  if( _min != 0 || _max !=0){
    if( _min > Q2_min )
      Q2_min = _min;
    if( _max < Q2_max )
      Q2_max = _max;
  }
  // Simpsons rule in using exponential step size and 10000 steps. Tested trapeve rule, linear steps and
  // nstep = 10,000 100,000 & 10,000,0000
  int nstep = 1000;
  double ln_min = std::log(Q2_min);
  double ratio = std::log(Q2_max/Q2_min)/nstep;
  double g_int = 0;
  double g_int2 = 0 ;
  double g_int3 = 0;
  //cout<<"*** Lomnitz **** Energy "<<Egamma<<" limits "<<Q2_min*1E9<<" x 1E-9 -  "<<Q2_max<<endl;
  for ( int ii = 0 ; ii< nstep; ++ii){
    double x1 =  std::exp(ln_min+(double)ii*ratio);
    double x3 =  std::exp(ln_min+(double)(ii+1)*ratio);
    double x2 =  (x3+x1)/2.;
    //cout<<"ii : "<<x1<<" "<<x2<<" "<<x3<<endl;
    g_int += (x3-x1)*( g(Egamma,x3)+g(Egamma,x1) +4.*g(Egamma,x2));
    g_int2 += (x3-x1)*( photonFlux(Egamma,x3)+photonFlux(Egamma,x1) +4.*photonFlux(Egamma,x2));
    g_int3 += (x3-x1)*( getcsgA_Q2_dep(x3)+getcsgA_Q2_dep(x1) +4.*getcsgA_Q2_dep(x2));
  }
  //return g_int2*g_int3/36.; 
  //return g_int2/6.;
  return g_int/6.;
}


//______________________________________________________________________________
double 
photonNucleusCrossSection::integrated_x_section(double const Egamma, double const _min, double const _max)
{
  //Integration over full limits gives more accurate result
  double Q2_min =  std::pow(starlightConstants::mel*Egamma,2.0)/(_electronEnergy*(_electronEnergy-Egamma));
  //double Q2_max = 2.*Egamma/_targetRadii - _wMax*_wMax;
  double Q2_max  = 4.*_electronEnergy*(_electronEnergy-Egamma);
  // Fixed ranges for plot
  if( _min != 0 || _max!=0){
    if( _min > Q2_min)
      Q2_min = _min;
    if( _max < Q2_max)
      Q2_max = _max;
  }
  // Simpsons rule in using exponential step size and 10000 steps. Tested trapeve rule, linear steps and
  // nstep = 10,000 100,000 & 10,000,0000
  int nstep = 1000;
  double ln_min = std::log(Q2_min);
  double ratio = std::log(Q2_max/Q2_min)/nstep;
  double g_int = 0;
  for ( int ii = 0 ; ii< nstep; ++ii){
    double x1 =  std::exp(ln_min+(double)ii*ratio);
    double x3 =  std::exp(ln_min+(double)(ii+1)*ratio);
    double x2 =  (x3+x1)/2.;
    //Tests from HERA https://arxiv.org/pdf/hep-ex/9601009.pdf
    //    g_int += (x3-x1)*( getcsgA_Q2_dep(x3)+getcsgA_Q2_dep(x1) +4.*getcsgA_Q2_dep(x2));
    g_int += (x3-x1)*( getcsgA_Q2_dep(x3)+getcsgA_Q2_dep(x1) +4.*getcsgA_Q2_dep(x2));
  }
  //return g_int2*g_int3/36.; 
  return g_int/6.;
}


//______________________________________________________________________________
pair< double, double >*photonNucleusCrossSection::Q2arraylimits(double const Egamma)
{
  //double Q2max = 2.*Egamma/_targetRadii - _wMax*_wMax;
  double Q2max = 4.*_electronEnergy*(_electronEnergy-Egamma);
  double Q2min= std::pow(starlightConstants::mel*Egamma,2.0)/(_electronEnergy*(_electronEnergy-Egamma));

  if( _fixedQ2range == true){
    if( Q2min < _minQ2 )
      Q2min = _minQ2;
    if( Q2max > _maxQ2 )
      Q2max = _maxQ2;
    //cout<<" Imposed limits "<<Q2min<<" - "<<Q2max<<endl;
    std::pair<double,double>* to_ret = new std::pair<double, double>(Q2min,Q2max);
    return to_ret;
  }
  int Nstep = 1000;
  //
  double ratio = std::log(Q2max/Q2min)/(double)Nstep;
  double ln_min = std::log(Q2min);
  // -- - -
  const double limit = 1E-9;
  std::vector<double>Q2_array;
  int iNstep = 0;
  double g_step = 1.;
  while( g_step>limit ){
    double Q2 = std::exp(ln_min+iNstep*ratio);
    if(Q2>Q2max) 
      break;
    g_step = g(Egamma,Q2);
    Q2_array.push_back(g_step);
    iNstep++;
  }
  if( std::exp(ln_min+iNstep*ratio) < Q2max)
    Q2max = std::exp(ln_min+iNstep*ratio);
  //cout<<Q2max<<" "<<g(Egamma,Q2max)*1E9<<endl;
  std::pair<double, double>* to_ret = new std::pair<double, double>(Q2min,Q2max);

  return to_ret;
}


//______________________________________________________________________________
double 
photonNucleusCrossSection::g(double const Egamma,
			     double const Q2)
{
  //return photonFlux(Egamma,Q2)*getcsgA_Q2_dep(Q2);
  return photonFlux(Egamma,Q2); //This could be done more elegantly in the future, but this one change should account for everything
}

//______________________________________________________________________________
double 
photonNucleusCrossSection::photonFlux(double const Egamma, 
				      double const Q2)
{
  //Need to check units later
  //double const hbar = starlightConstants::hbarc / 2.99*pow(10,14); // 6.582 * pow (10,-16) [eVs]
  //double omega = Egamma/ hbar;
  //Do we even need a lookup table for this case? This should return N(E,Q2) from dn = N(E,Q2) dE dQ2
  double const ratio = Egamma/_electronEnergy;
  double const minQ2 = std::pow( starlightConstants::mel*Egamma,2.0) / (_electronEnergy*(_electronEnergy - Egamma));
  double to_ret = alpha/(pi) *( 1- ratio + ratio*ratio/2. - (1-ratio)*( fabs(minQ2/Q2)) );
  //Comparisons:
  //  double temp = pow(2.*_electronEnergy-Egamma,2.)/(Egamma*Egamma + Q2) + 1. - 4.*starlightConstants::mel*starlightConstants::mel/Q2;
  //temp = alpha*temp*Egamma/(4.*Q2*pi*_electronEnergy*_electronEnergy);
  //  cout<<" *** Lomnitz *** Testing photon flux approximation for electron energy "<<_electronEnergy<<" gamma "<<Egamma<<" Q2 "<<Q2*1E6<<" 1E-6 "<<endl;
  //cout<<" Full expression "<<temp*1E6<<" vs. apporioximation "<<to_ret/( Egamma*fabs(Q2) )*1E6<<" --- ratio "<<temp/(to_ret/( Egamma*fabs(Q2) ))<<endl;
  return to_ret/( Egamma*fabs(Q2) );
  //return temp;
}

//______________________________________________________________________________
double
photonNucleusCrossSection::nepoint(const double Egamma,
                                   const double bmin)
{
	// Function for the spectrum of virtual photons,
	// dn/dEgamma, for a point charge q=Ze sweeping
	// past the origin with velocity gamma
	// (=1/SQRT(1-(V/c)**2)) integrated over impact
	// parameter from bmin to infinity
	// See Jackson eq15.54 Classical Electrodynamics
	// Declare Local Variables
	double beta,X,C1,bracket,nepoint_r;
  
	beta = sqrt(1.-(1./(_beamLorentzGamma*_beamLorentzGamma)));
	X = (bmin*Egamma)/(beta*_beamLorentzGamma*hbarc);
  
	bracket = -0.5*beta*beta*X*X*(bessel::dbesk1(X)*bessel::dbesk1(X)
	                              -bessel::dbesk0(X)*bessel::dbesk0(X));

	bracket = bracket+X*bessel::dbesk0(X)*bessel::dbesk1(X);
  
	// Note: NO  Z*Z!!
	C1=(2.*alpha)/pi;
  
	nepoint_r = C1*(1./beta)*(1./beta)*(1./Egamma)*bracket;
  
	return nepoint_r;
  
}


//______________________________________________________________________________
double
photonNucleusCrossSection::sigmagp(const double Wgp)
{
	// Function for the gamma-proton --> VectorMeson
	// cross section. Wgp is the gamma-proton CM energy.
	// Unit for cross section: fm**2
  
	double sigmagp_r=0.;

	// Near the threshold CM energy (between WgpMin and WgpMax), 
	// we add a linear scaling factor to bring the Xsec down to 0 
	// at the threshold CM value define by WgpMin = m_p + m_vm
	double WgpMax = 0.;
	double WgpMin = 0.;
	double thresholdScaling = 1.0;
  
	switch(_particleType)
		{ 
		case RHO:
		case RHOZEUS:
		case FOURPRONG:
			WgpMax = 1.8;
			WgpMin = 1.60; //this is the cutoff threshold for rho production. But rho has large width so it's lower
			if(Wgp<WgpMax) thresholdScaling=(Wgp-WgpMin)/(WgpMax-WgpMin);
			sigmagp_r=thresholdScaling*1.E-4*(5.0*exp(0.22*log(Wgp))+26.0*exp(-1.23*log(Wgp)));
			//This is based on the omega cross section:
			//https://arxiv.org/pdf/2107.06748.pdf
			if(_backwardsProduction)sigmagp_r=thresholdScaling*1.E-4*0.14*pow(Wgp,-2.7);
			if(_backwardsProduction)sigmagp_r=thresholdScaling*1.E-4*0.206*pow(((Wgp*Wgp-protonMass*protonMass)/(2.0*protonMass)),-2.7);
			break;
		case OMEGA:
		case OMEGA_pipipi:
			WgpMax = 1.8;
			WgpMin = 1.74; //this is the cutoff threshold for omega production: W > m_p+m_omega = 1.74
			if(Wgp<WgpMax) thresholdScaling=(Wgp-WgpMin)/(WgpMax-WgpMin);
			sigmagp_r=thresholdScaling*1.E-4*(0.55*exp(0.22*log(Wgp))+18.0*exp(-1.92*log(Wgp)));
			//if(_backwardsProduction)sigmagp_r=thresholdScaling*1.E-4*0.14*pow(Wgp,-2.7);//https://arxiv.org/pdf/2107.06748.pdf
			if(_backwardsProduction)sigmagp_r=thresholdScaling*1.E-4*0.206*pow(((Wgp*Wgp-protonMass*protonMass)/(2.0*protonMass)),-2.7);
			break;   
		case OMEGA_pi0gamma:
			WgpMax = 1.8;
			WgpMin = 1.74; //this is the cutoff threshold for omega production: W > m_p+m_omega = 1.74
			if(Wgp<WgpMax) thresholdScaling=(Wgp-WgpMin)/(WgpMax-WgpMin);
			sigmagp_r=thresholdScaling*1.E-4*(0.55*exp(0.22*log(Wgp))+18.0*exp(-1.92*log(Wgp)));
			//if(_backwardsProduction)sigmagp_r=thresholdScaling*1.E-4*0.14*pow(Wgp,-2.7);//https://arxiv.org/pdf/2107.06748.pdf
			if(_backwardsProduction)sigmagp_r=thresholdScaling*1.E-4*0.206*pow(((Wgp*Wgp-protonMass*protonMass)/(2.0*protonMass)),-2.7);
			break;                                                      
		case PHI:
			sigmagp_r=1.E-4*0.34*exp(0.22*log(Wgp));
			break;
		case JPSI:
		case JPSI_ee:
		case JPSI_mumu:
			sigmagp_r=(1.0-((_channelMass+protonMass)*(_channelMass+protonMass))/(Wgp*Wgp));
			sigmagp_r*=sigmagp_r;
			sigmagp_r*=1.E-4*0.00406*exp(0.65*log(Wgp));
			// sigmagp_r=1.E-4*0.0015*exp(0.80*log(Wgp));
			break;
		case JPSI2S:
		case JPSI2S_ee:
		case JPSI2S_mumu:
			sigmagp_r=(1.0-((_channelMass+protonMass)*(_channelMass+protonMass))/(Wgp*Wgp));
			sigmagp_r*=sigmagp_r;
			sigmagp_r*=1.E-4*0.00406*exp(0.65*log(Wgp));
			sigmagp_r*=0.166;  
			//      sigmagp_r=0.166*(1.E-4*0.0015*exp(0.80*log(Wgp)));
			break;
		case UPSILON:
		case UPSILON_ee:
		case UPSILON_mumu:
			//       >> This is W**1.7 dependence from QCD calculations
			//  sigmagp_r=1.E-10*(0.060)*exp(1.70*log(Wgp));
			sigmagp_r=(1.0-((_channelMass+protonMass)*(_channelMass+protonMass))/(Wgp*Wgp));
			sigmagp_r*=sigmagp_r;
			sigmagp_r*=1.E-10*6.4*exp(0.74*log(Wgp));
			break;
		case UPSILON2S:
		case UPSILON2S_ee:
		case UPSILON2S_mumu:
		        // sigmagp_r=1.E-10*(0.0259)*exp(1.70*log(Wgp));
		        sigmagp_r=(1.0-((_channelMass+protonMass)*(_channelMass+protonMass))/(Wgp*Wgp));
			sigmagp_r*=sigmagp_r;
			sigmagp_r*=1.E-10*2.9*exp(0.74*log(Wgp)); 
			break;
		case UPSILON3S:
		case UPSILON3S_ee:
		case UPSILON3S_mumu:
		        // sigmagp_r=1.E-10*(0.0181)*exp(1.70*log(Wgp));
		        sigmagp_r=(1.0-((_channelMass+protonMass)*(_channelMass+protonMass))/(Wgp*Wgp));
			sigmagp_r*=sigmagp_r;
			sigmagp_r*=1.E-10*2.1*exp(0.74*log(Wgp)); 
			break;
		default: cout<< "!!!  ERROR: Unidentified Vector Meson: "<< _particleType <<endl;
		}                                                                  
	return sigmagp_r;
}


//______________________________________________________________________________
double
photonNucleusCrossSection::sigma_A(const double sig_N, const int beam)
{                                                         
	// Nuclear Cross Section
	// sig_N,sigma_A in (fm**2) 

	double sum;
	double b,bmax,Pint,arg,sigma_A_r;
  
	int NGAUSS;
  
	double xg[17]=
		{.0,
		 .0483076656877383162,.144471961582796493,
		 .239287362252137075, .331868602282127650,
		 .421351276130635345, .506899908932229390,
		 .587715757240762329, .663044266930215201,
		 .732182118740289680, .794483795967942407,
		 .849367613732569970, .896321155766052124,
		 .934906075937739689, .964762255587506430,
		 .985611511545268335, .997263861849481564
		};
  
	double ag[17]=
		{.0,
		 .0965400885147278006, .0956387200792748594,
		 .0938443990808045654, .0911738786957638847,
		 .0876520930044038111, .0833119242269467552,
		 .0781938957870703065, .0723457941088485062,
		 .0658222227763618468, .0586840934785355471,
		 .0509980592623761762, .0428358980222266807,
		 .0342738629130214331, .0253920653092620595,
		 .0162743947309056706, .00701861000947009660
		};
  
	NGAUSS=16;
 
	// Check if one or both beams are nuclei 
        int A_1 = _bbs.electronBeam().A(); 
        int A_2 = _bbs.targetBeam().A(); 
        if( A_1 == 1 && A_2 == 1)cout<<" This is pp, you should not be here..."<<endl;  

	// CALCULATE P(int) FOR b=0.0 - bmax (fm)
	bmax = 25.0;
	sum  = 0.;
	for(int IB=1;IB<=NGAUSS;IB++){
    
		b = 0.5*bmax*xg[IB]+0.5*bmax;

                if( A_1 == 1 && A_2 != 1){  
                  arg=-sig_N*_bbs.targetBeam().rho0()*_bbs.targetBeam().thickness(b);
                }else if(A_2 == 1 && A_1 != 1){
                  arg=-sig_N*_bbs.electronBeam().rho0()*_bbs.electronBeam().thickness(b);
                }else{
		  // Check which beam is target 
                  if( beam == 1 ){ 
                    arg=-sig_N*_bbs.electronBeam().rho0()*_bbs.electronBeam().thickness(b);
                  }else if( beam==2 ){
                    arg=-sig_N*_bbs.targetBeam().rho0()*_bbs.targetBeam().thickness(b);
                  }else{
                    cout<<" Something went wrong here, beam= "<<beam<<endl; 
                  } 
                }
    
		Pint=1.0-exp(arg);
		// If this is a quantum Glauber calculation, use the quantum Glauber formula
		if (_quantumGlauber == 1){Pint=2.0*(1.0-exp(arg/2.0));}

		sum=sum+2.*pi*b*Pint*ag[IB];
    
		b = 0.5*bmax*(-xg[IB])+0.5*bmax;

                if( A_1 == 1 && A_2 != 1){  
                  arg=-sig_N*_bbs.targetBeam().rho0()*_bbs.targetBeam().thickness(b);
                }else if(A_2 == 1 && A_1 != 1){
                  arg=-sig_N*_bbs.electronBeam().rho0()*_bbs.electronBeam().thickness(b);
                }else{ 
		  // Check which beam is target 
                  if( beam == 1 ){ 
                    arg=-sig_N*_bbs.electronBeam().rho0()*_bbs.electronBeam().thickness(b);
                  }else if(beam==2){
                    arg=-sig_N*_bbs.targetBeam().rho0()*_bbs.targetBeam().thickness(b);
                  }else{
                    cout<<" Something went wrong here, beam= "<<beam<<endl; 
                  } 
                }

		Pint=1.0-exp(arg);
		// If this is a quantum Glauber calculation, use the quantum Glauber formula
		if (_quantumGlauber == 1){Pint=2.0*(1.0-exp(arg/2.0));}		

		sum=sum+2.*pi*b*Pint*ag[IB];

	}

	sum=0.5*bmax*sum;
  
	sigma_A_r=sum;
 
	return sigma_A_r;
}

//______________________________________________________________________________
double
photonNucleusCrossSection::sigma_N(const double Wgp)
{                                                         
        // Nucleon Cross Section in (fm**2) 
        double cs = sqrt(16. * pi * _vmPhotonCoupling * _slopeParameter * hbarc * hbarc * sigmagp(Wgp) / alpha);
        return cs;
}


//______________________________________________________________________________
double
photonNucleusCrossSection::breitWigner(const double W,
                                       const double C)
{
	// use simple fixed-width s-wave Breit-Wigner without coherent backgorund for rho'
	// (PDG '08 eq. 38.56)
	if(_particleType==FOURPRONG) {
		if (W < 4.01 * pionChargedMass)
			return 0;
		const double termA  = _channelMass * _width;
		const double termA2 = termA * termA;
		const double termB  = W * W - _channelMass * _channelMass;
		return C * _ANORM * _ANORM * termA2 / (termB * termB + termA2);
	}

	// Relativistic Breit-Wigner according to J.D. Jackson,
	// Nuovo Cimento 34, 6692 (1964), with nonresonant term. A is the strength
	// of the resonant term and b the strength of the non-resonant
	// term. C is an overall normalization.

	double ppi=0.,ppi0=0.,GammaPrim,rat;
	double aa,bb,cc;
  
	double nrbw_r;

	// width depends on energy - Jackson Eq. A.2
	// if below threshold, then return 0.  Added 5/3/2001 SRK
	// 0.5% extra added for safety margin
        // omega added here 10/26/2014 SRK
	if( _particleType==RHO ||_particleType==RHOZEUS || _particleType==OMEGA){  
		if (W < 2.01*pionChargedMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt( ((W/2.)*(W/2.)) - pionChargedMass * pionChargedMass);
		ppi0=0.358;
	}

	if( _particleType==OMEGA_pipipi){  
		if (W < 1.01*(2.0*pionChargedMass+pionNeutralMass)){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=abs((W*W - pionNeutralMass * pionNeutralMass - 2*pionChargedMass*pionChargedMass)/(2.0*W));
		ppi0=abs((_channelMass*_channelMass - pionNeutralMass * pionNeutralMass - 2*pionChargedMass*pionChargedMass)/(2.0*W));;
	}

	if( _particleType==OMEGA_pi0gamma){  
		if (W < pionNeutralMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=abs((W*W - pionNeutralMass * pionNeutralMass)/(2.0*W));
		ppi0=abs((_channelMass*_channelMass - pionNeutralMass * pionNeutralMass)/(2.0*W));
	}
  
	// handle phi-->K+K- properly
	if (_particleType  ==  PHI){
		if (W < 2.*kaonChargedMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt( ((W/2.)*(W/2.))- kaonChargedMass*kaonChargedMass);
		ppi0=sqrt( ((_channelMass/2.)*(_channelMass/2.))-kaonChargedMass*kaonChargedMass);
	}

	//handle J/Psi-->e+e- properly
	if (_particleType==JPSI || _particleType==JPSI2S){
		if(W<2.*mel){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-mel*mel);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-mel*mel);
	}
	if (_particleType==JPSI_ee){
		if(W<2.*mel){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-mel*mel);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-mel*mel);   
	}
	if (_particleType==JPSI_mumu){
		if(W<2.*muonMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-muonMass*muonMass);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-muonMass*muonMass);
	}
	if (_particleType==JPSI2S_ee){
		if(W<2.*mel){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-mel*mel);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-mel*mel);   
	}
	if (_particleType==JPSI2S_mumu){
		if(W<2.*muonMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-muonMass*muonMass);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-muonMass*muonMass);
	}

	if(_particleType==UPSILON || _particleType==UPSILON2S ||_particleType==UPSILON3S ){ 
		if (W<2.*muonMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-muonMass*muonMass);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-muonMass*muonMass);
	}
  
	if(_particleType==UPSILON_mumu || _particleType==UPSILON2S_mumu ||_particleType==UPSILON3S_mumu ){ 
		if (W<2.*muonMass){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-muonMass*muonMass);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-muonMass*muonMass);
	}
  
	if(_particleType==UPSILON_ee || _particleType==UPSILON2S_ee ||_particleType==UPSILON3S_ee ){ 
		if (W<2.*mel){
			nrbw_r=0.;
			return nrbw_r;
		}
		ppi=sqrt(((W/2.)*(W/2.))-mel*mel);
		ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-mel*mel);
	}
  
	if(ppi==0.&&ppi0==0.) 
		cout<<"Improper Gammaacrosssection::breitwigner, ppi&ppi0=0."<<endl;
  
	rat=ppi/ppi0;
	GammaPrim=_width*(_channelMass/W)*rat*rat*rat;
  
	aa=_ANORM*sqrt(GammaPrim*_channelMass*W);
	bb=W*W-_channelMass*_channelMass;
	cc=_channelMass*GammaPrim;
  
	// First real part squared 
	nrbw_r = (( (aa*bb)/(bb*bb+cc*cc) + _BNORM)*( (aa*bb)/(bb*bb+cc*cc) + _BNORM));
  
	// Then imaginary part squared 
	nrbw_r = nrbw_r + (( (aa*cc)/(bb*bb+cc*cc) )*( (aa*cc)/(bb*bb+cc*cc) ));

	//  Alternative, a simple, no-background BW, following J. Breitweg et al.
	//  Eq. 15 of Eur. Phys. J. C2, 247 (1998).  SRK 11/10/2000
	//      nrbw_r = (_ANORM*_mass*GammaPrim/(bb*bb+cc*cc))**2

  
	nrbw_r = C*nrbw_r;
  
	return nrbw_r;    
}
