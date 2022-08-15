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
// $Rev:: 1                           $: revision of last commit
// $Author:: zsweger                  $: author of last commit
// $Date:: 2022-07-29                #$: date of last commit
//
// Description:
//    Created DVCS
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "dvcs.h"
#include "e_DVCSCrossSection.h"

using namespace std;


//______________________________________________________________________________
Dvcs::Dvcs(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem):eventChannel(inputParametersInstance, bbsystem), _phaseSpaceGen(0)
{
	_VMnumega=inputParametersInstance.nmbEnergyBins();
	_VMnumQ2=inputParametersInstance.nmbGammaQ2Bins(); 
	_VMbslope=0.;//Will define in wide/narrow constructor
    _pEnergy= inputParametersInstance.protonEnergy();
	_beamNucleus = inputParametersInstance.targetBeamA();
	// electron energy in CMS frame
	_eEnergy= inputParametersInstance.electronEnergy();
	_targetMaxPhotonEnergy=inputParametersInstance.targetMaxPhotonEnergy();
	_targetMinPhotonEnergy=inputParametersInstance.targetMinPhotonEnergy();//ZS address this!
	// Now saving the photon energy limits
	_cmsMinPhotonEnergy=inputParametersInstance.cmsMinPhotonEnergy();
	_beamLorentzGamma = inputParametersInstance.beamLorentzGamma();
	_targetBeamLorentzGamma = inputParametersInstance.targetBeamLorentzGamma();
	_rap_CM=inputParametersInstance.rap_CM();
	_targetRadius = inputParametersInstance.targetRadius();
	//Turn on/off backwards production
	_backwardsProduction = inputParametersInstance.backwardsProduction();
	
}


//______________________________________________________________________________
Dvcs::~Dvcs()
{
	if (_phaseSpaceGen)
		delete _phaseSpaceGen;
	if (_dummy_pncs)
	  delete _dummy_pncs;
}


//______________________________________________________________________________
void Dvcs::momenta(double Egam,double Q2, double gamma_pz, double gamma_pt, //input conditions
				double &Y,double &E,double &px,double &py,double &pz,  //return final-state photon
				double &t_px, double &t_py, double &t_pz, double &t_E, //return pomeron
				double &e_phi,int &tcheck) //return electron (angle already known by Q2)
{
	//     This subroutine calculates momentum and energy of vector meson for electroproduction (eSTARlight)
	//     given W and photon 4-vector,   without interference.  No intereference in asymetric eX collisions
 	Q2 = 1.0*Q2;
	double phi1;
	double px1,py1;
	double xtest;
	double deltak;

	double target_px, target_py, target_pz, target_E;

	tcheck=1*tcheck;//figure out what to do with tcheck
    target_E = _beamNucleus*_pEnergy;
    target_px = 0.0;
    target_py = 0.0;
    target_pz = -_beamNucleus*sqrt(_pEnergy*_pEnergy - pow(starlightConstants::protonMass,2.) );
    phi1 = 2.*starlightConstants::pi*_randy.Rndom();
    px1 = gamma_pt*cos(phi1);
	py1 = gamma_pt*sin(phi1);
	int isbadevent = 0;
	double betax_cm = ((px1+target_px)/(Egam+target_E));
    double betay_cm = ((py1+target_py)/(Egam+target_E));
    double betaz_cm = ((gamma_pz+target_pz)/(Egam+target_E));
    transform (betax_cm,betay_cm,betaz_cm,target_E,target_px,target_py,target_pz,isbadevent);
    transform (betax_cm,betay_cm,betaz_cm,Egam, px1, py1, gamma_pz, isbadevent);
      
	e_phi = starlightConstants::pi+phi1;

    while( e_phi > 2.*starlightConstants::pi ) e_phi-= 2.*starlightConstants::pi;
	//
	if ( _bbs.targetBeam().A()!=1 ) {cout<<"DVCS is NOT yet defined for non-protons!!!"<<endl; exit(1);}
	   
L116dvcs:	    
	//Use dsig/dt= exp(-_VMbslope*t) for heavy VM
    double bslope_tdist = _VMbslope; 
	xtest = _randy.Rndom(); 
	deltak = (-1./bslope_tdist)*log(xtest);
	double u_or_t = -1.0*abs(deltak);
	//cout<<"u_or_t: "<<u_or_t<<endl;


	//rotate to put  gamma--> <--p  reaction along one axis
	double theta_pgamma = atan(sqrt(target_px*target_px + target_py*target_py)/target_pz);	
	double phi_pgamma = 0.0;
	if(target_px > 0.0) phi_pgamma = atan(target_py/target_px);
	else phi_pgamma = starlightConstants::pi - atan(target_py/abs(target_px));
	//double target_pz_rot = sqrt(target_px*target_px + target_py*target_py + target_pz*target_pz);
	//calculate scattering angle from u
	double target_ptot = sqrt(target_px*target_px + target_py*target_py + target_pz*target_pz);
	double costheta_scatter = 0.0;
	double theta_scatter = 0.0;
	double fsphoton_ptot = 0.0;
	double backwards_forwards = 1.0;

	//the variables in the following calculations are dummy variables introduced for brevity
	if(_backwardsProduction)
	{
		backwards_forwards = -1.0;
		double a = (0.5)*(pow(starlightConstants::protonMass,2.)- u_or_t -Q2);
    	double b = Egam;
    	double c = pow(starlightConstants::protonMass,2.);
    	double d = u_or_t;
    	double f = target_E;
    	double g = sqrt(target_px*target_px + target_py*target_py + target_pz*target_pz);
		costheta_scatter = (sqrt(f*f*g*g*pow((-8.0*a*a + 4.0*a*c - 4.0*a*d + 8.0*b*b*c),2) - 4.0*g*g*(4.0*a*a - 4.0*a*c + 4.0*a*d - 4.0*b*b*c + c*c - 2.0*c*d + d*d)*(4.0*a*a*f*f + b*b*(-c*c) + 2.0*b*b*c*d - 4.0*b*b*c*f*f - b*b*d*d)) - f*g*(-8.0*a*a + 4.0*a*c - 4.0*a*d + 8.0*b*b*c))/(2.0*g*g*(4.0*a*a - 4.0*a*c + 4.0*a*d - 4.0*b*b*c + c*c - 2.0*c*d + d*d));
		theta_scatter = acos(costheta_scatter);
		fsphoton_ptot = (pow(starlightConstants::protonMass,2.) - (u_or_t))/(2.0*(target_E - target_ptot*cos(theta_scatter)));
	}
	else{
		double b = target_E;
    	double c = u_or_t;
    	double d = Q2;
    	double f = Egam;
    	double g = target_ptot;
    	double h = pow(starlightConstants::protonMass,2.);
		costheta_scatter = (-1.0*sqrt(f*f*g*g*pow((-8.0*b*b*h - 2.0*c*d - 4.0*c*h + 4.0*d*h + 8.0*h*h),2) - 4.0*g*g*(4.0*b*b*h - d*d - 4.0*d*h - 4.0*h*h)*(b*b*c*c + 2.0*b*b*c*d + b*b*d*d + 4.0*b*b*f*f*h - c*c*f*f + 4.0*c*f*f*h - 4.0*f*f*h*h)) - f*g*(-8.0*b*b*h - 2.0*c*d - 4.0*c*h + 4.0*d*h + 8.0*h*h))/(2.0*g*g*(4.0*b*b*h - d*d - 4.0*d*h - 4.0*h*h));
		theta_scatter = acos(costheta_scatter);
		fsphoton_ptot = (-1.0*u_or_t - Q2)/(2.0*(Egam - target_ptot*cos(theta_scatter)));
	}
	if(costheta_scatter>1.0){goto L116dvcs;}
	//

	//final-state photon:
	double fsphoton_phi_rot = 2.*starlightConstants::pi*_randy.Rndom();
	double fsphoton_px_rot = fsphoton_ptot*sin(theta_scatter)*cos(fsphoton_phi_rot);
	double fsphoton_py_rot = fsphoton_ptot*sin(theta_scatter)*sin(fsphoton_phi_rot);
	double fsphoton_pz_rot = backwards_forwards*fsphoton_ptot*cos(theta_scatter);
	//final-state proton:
	double newion_px_rot = -1.0*fsphoton_px_rot;
	double newion_py_rot = -1.0*fsphoton_py_rot;
	double newion_pz_rot = -1.0*fsphoton_pz_rot;
	double newion_E = sqrt(newion_px_rot*newion_px_rot + newion_py_rot*newion_py_rot + newion_pz_rot*newion_pz_rot + pow(_beamNucleus*starlightConstants::protonMass,2.));
	//rotate back into normal cm frame
	double fsphoton_px = (fsphoton_px_rot*cos(theta_pgamma) + fsphoton_pz_rot*sin(theta_pgamma))*cos(phi_pgamma) - fsphoton_py_rot*sin(phi_pgamma);
	double fsphoton_py = (fsphoton_px_rot*cos(theta_pgamma) + fsphoton_pz_rot*sin(theta_pgamma))*sin(phi_pgamma) + fsphoton_py_rot*cos(phi_pgamma);
	double fsphoton_pz = fsphoton_pz_rot*cos(theta_pgamma) - fsphoton_px_rot*sin(theta_pgamma);
	double newion_px = (newion_px_rot*cos(theta_pgamma) + newion_pz_rot*sin(theta_pgamma))*cos(phi_pgamma) - newion_py_rot*sin(phi_pgamma);
	double newion_py = (newion_px_rot*cos(theta_pgamma) + newion_pz_rot*sin(theta_pgamma))*sin(phi_pgamma) + newion_py_rot*cos(phi_pgamma);
	double newion_pz = newion_pz_rot*cos(theta_pgamma) - newion_px_rot*sin(theta_pgamma);
	//
	t_px = target_px - newion_px;
	t_py = target_py - newion_py;
	t_pz = target_pz - newion_pz;
	t_E  = target_E  - newion_E;
	//
	px = fsphoton_px;
	py = fsphoton_py;
	pz = fsphoton_pz;
	E  = fsphoton_ptot;

	//transform back to e-ion CM frame
	transform (-betax_cm,-betay_cm,-betaz_cm,target_E,target_px,target_py,target_pz,isbadevent);
	transform (-betax_cm,-betay_cm,-betaz_cm,newion_E,newion_px,newion_py,newion_pz,isbadevent);
    transform (-betax_cm,-betay_cm,-betaz_cm,Egam,    px1,      py1,      gamma_pz, isbadevent);
    transform (-betax_cm,-betay_cm,-betaz_cm,E,       px,       py,       pz,       isbadevent);
    t_px = target_px-newion_px;
    t_py = target_py-newion_py;
    t_pz = target_pz-newion_pz;
    t_E  = target_E-newion_E;
	Y = 0.5*std::log( (E+fabs(pz))/(E-fabs(pz)) );
}



//______________________________________________________________________________
double Dvcs::pseudoRapidity(double px, double py, double pz)
{
	double pT = sqrt(px*px + py*py);
	double p = sqrt(pz*pz + pT*pT);
	double eta = -99.9; if((p-pz) != 0){eta = 0.5*log((p+pz)/(p-pz));}
	return eta;
}




//______________________________________________________________________________
e_Dvcs::e_Dvcs(const inputParameters& input, beamBeamSystem& bbsystem):Dvcs(input, bbsystem)
{
	cout<<"Reading in luminosity tables. e_Dvcs()"<<endl;
	e_read();
	cout<<"Creating and calculating crosssection. e_Dvcs()"<<endl;
	e_DVCSCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter(); 
}

//______________________________________________________________________________
e_Dvcs::~e_Dvcs()
{ }



//______________________________________________________________________________
void Dvcs::pickEgamq2(double &cmsEgamma, double &targetEgamma, 
				 double &Q2, double &gamma_pz, double &gamma_pt,//photon in target frame
				 double &E_prime, double &theta_e //electron
				 )
{
    double dEgamma;
	double xEgamma, xQ2, xtest, q2test; // btest;
	int  IGamma, IQ2;
	// ---------
	//	int egamma_draws = 0, cms_egamma_draws =0, q2_draws =0 ;
	// ---------
	// Following used for debugging/timing for event generation
	//std::chrono::steady_clock::time_point begin_evt = std::chrono::steady_clock::now();
	bool pick_state = false;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy);
	//cout<<_targetMaxPhotonEnergy<<endl;
	//cout<<_targetMinPhotonEnergy<<endl;
	//cout<<dEgamma<<endl;
	while( pick_state == false ){
	  
	  xEgamma = _randy.Rndom();
	  //
	  targetEgamma = std::exp(std::log(_targetMinPhotonEnergy) + xEgamma*(dEgamma));
	  IGamma = int(_VMnumega*xEgamma);
	  // Holds Q2 and integrated Q2 dependence. Array is saved in target frame
	  std::pair< double, std::vector<double> > this_energy = _g_EQ2array->operator[](IGamma);
	  double intgrated_q2 = this_energy.first;
	  // Sample single differential photon flux to obtain photon energy
	  xtest = _randy.Rndom();
	  if( xtest > intgrated_q2 ){
	    //egamma_draws+=1;
	    continue;
	  }
	  //	  btest = _randy.Rndom();
	  // Load double differential photon flux table for selected energy
	  std::vector<double> photon_flux = this_energy.second;
	  double VMQ2min = photon_flux[0];
	  double VMQ2max = photon_flux[1];
	  //
	  double ratio = std::log(VMQ2max/VMQ2min);
 	  double ln_min = std::log(VMQ2min);
	  
	  xQ2 = _randy.Rndom();
	  Q2 = std::exp(ln_min+xQ2*ratio);
	  IQ2 = int(_VMnumQ2*xQ2);	
	  // Load from look-up table. Use linear interpolation to evaluate at Q2
	  double y_1 = photon_flux[IQ2+2];
	  double y_2 = photon_flux[IQ2+3];
	  double x_1 = std::exp(ln_min+IQ2*ratio/_VMnumQ2);
	  double x_2 = std::exp(ln_min+(1+IQ2)*ratio/_VMnumQ2);
	  double m = (y_2 - y_1)/(x_2 - x_1);
	  double c = y_1-m*x_1;
	  double y = m*Q2+c;
	  q2test = _randy.Rndom();
	  if( y < q2test ){
	    continue;
	  }
	  // -- Generate electron and photon in Target frame
	  E_prime = _eEnergy - targetEgamma;
	  if(Q2>1E-6){
	  	double cos_theta_e = 1. - Q2/(2.*_eEnergy*E_prime);
	    theta_e = acos(cos_theta_e);
	  }
	  else {theta_e = sqrt(Q2/(_eEnergy*E_prime));}//updated using small angle approximation to avoid rounding to 0
	  double beam_y = acosh(_targetBeamLorentzGamma)+_rap_CM;	
	  gamma_pt = E_prime*sin(theta_e);
	  
	  double pz_squared = targetEgamma*targetEgamma + Q2 - gamma_pt*gamma_pt;
	  if( pz_squared < 0 )
	    continue;
	  double temp_pz = sqrt(pz_squared);
	  // Now boost to CMS frame to check validity
	  gamma_pz = temp_pz*cosh(beam_y) - targetEgamma*sinh(beam_y); 
	  cmsEgamma = targetEgamma*cosh(beam_y) - temp_pz*sinh(beam_y);
	  // Simple checkl, should not be needed but used for safety
	  if( cmsEgamma < _cmsMinPhotonEnergy || 2.*targetEgamma/(Q2) < _targetRadius){ //This cut is roughly RA = 0.8 fm the radius of proton and 1 eV^{-1} = 1.97 x 10 ^{-7} m
	      continue;
	  }
	  xtest = _randy.Rndom();
	  // Test against photonuclear cross section gamma+X --> VM+X
	  if( _f_WYarray[IGamma][IQ2] < xtest ){
	    continue;
	  }
	  pick_state = true;
	}
	return;
}


//______________________________________________________________________________
eXEvent Dvcs::e_produceEvent()
{
	// The new event type
	eXEvent event;
	
	int iFbadevent=0;
	int tcheck=0;

	double rapidity = 0.;
	double Q2 = 0;
	double E = 0.;
	double momx=0.,momy=0.,momz=0.;
	double targetEgamma = 0, cmsEgamma = 0 ;
	double gamma_pz = 0 , gamma_pt = 0, e_theta = 0;
	double e_E=0., e_phi=0;
	double t_px =0, t_py=0., t_pz=0, t_E;
	bool accepted = false;
	do{
	  pickEgamq2(cmsEgamma, targetEgamma, 
		   Q2, gamma_pz, gamma_pt, //photon infor in CMS frame
		   e_E, e_theta);	 //electron info in target frame  
	  //
	  momenta(cmsEgamma, Q2, gamma_pz, gamma_pt, //input
		  rapidity, E, momx, momy, momz, //final-state photon
		  t_px, t_py, t_pz, t_E, //pomeron
		  e_phi,tcheck); //electron
	  //
	  _nmbAttempts++;

	  double pt1chk = sqrt(momx*momx+momy*momy);
	  double eta1 = pseudoRapidity(momx, momy, momz);
                        

	  if(_ptCutEnabled && !_etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(!_ptCutEnabled && _etaCutEnabled){
	    if(eta1 > _etaCutMin && eta1 < _etaCutMax ){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(_ptCutEnabled && _etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax ){
	      if(eta1 > _etaCutMin && eta1 < _etaCutMax ){
	      	accepted = true;
		    _nmbAccepted++;
	      }
	    }
	  }
	  else if(!_ptCutEnabled && !_etaCutEnabled)
	    _nmbAccepted++;
	}while((_ptCutEnabled || _etaCutEnabled) && !accepted);
	if (iFbadevent==0&&tcheck==0) {

	  // - Outgoing electron - target frame
	  double e_ptot = sqrt(e_E*e_E - starlightConstants::mel*starlightConstants::mel);
	  double e_px = e_ptot*sin(e_theta)*cos(e_phi);
	  double e_py = e_ptot*sin(e_theta)*sin(e_phi);
	  double e_pz = e_ptot*cos(e_theta);
	  lorentzVector electron(e_px, e_py, e_pz, e_E);
	  event.addSourceElectron(electron);
	  // - Generated virtual photon - target frame
	  double gamma_x = gamma_pt*cos(e_phi+starlightConstants::pi);
	  double gamma_y = gamma_pt*sin(e_phi+starlightConstants::pi);
	  lorentzVector gamma(gamma_x,gamma_y,gamma_pz,cmsEgamma);
	  vector3 boostVector(0, 0, tanh(_rap_CM));
	  (gamma).Boost(boostVector);
	  event.addGamma(gamma, targetEgamma, Q2);  

	  // - Saving final-state photon
	  double E_final_photon = sqrt(momx*momx+momy*momy+momz*momz); 
	  starlightParticle particle1(momx, momy, momz, E_final_photon, starlightConstants::UNKNOWN, 22, 0);
	  event.addParticle(particle1);


	  // - Scattered target and transfered momenta at target vertex
	  double target_pz =  - _beamNucleus*sqrt(_pEnergy*_pEnergy - pow(starlightConstants::protonMass,2.) ) - t_pz;
	  
	  lorentzVector target(-t_px, -t_py, target_pz, _beamNucleus*_pEnergy - t_E);
	  double t_var = t_E*t_E - t_px*t_px - t_py*t_py - t_pz*t_pz;
	  event.addScatteredTarget(target, t_var);
	}
	return event;

}