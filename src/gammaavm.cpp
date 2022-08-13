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
// $Rev:: 277                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-09-14 10:55:55 +0100 #$: date of last commit
//
// Description:
//    Added incoherent t2-> pt2 selection.  Following pp selection scheme
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "gammaavm.h"
//#include "photonNucleusCrossSection.h"
#include "wideResonanceCrossSection.h"
#include "narrowResonanceCrossSection.h"
#include "incoherentVMCrossSection.h"
//
#include "e_narrowResonanceCrossSection.h"
#include "e_wideResonanceCrossSection.h"

using namespace std;


//______________________________________________________________________________
Gammaavectormeson::Gammaavectormeson(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem):eventChannel(inputParametersInstance, bbsystem), _phaseSpaceGen(0)
{
	_VMNPT=inputParametersInstance.nmbPtBinsInterference();
	_VMWmax=inputParametersInstance.maxW();
	_VMWmin=inputParametersInstance.minW();
	//_VMYmax=inputParametersInstance.maxRapidity();
	//_VMYmin=-1.*_VMYmax;
	_VMnumw=inputParametersInstance.nmbWBins();
	//_VMnumy=inputParametersInstance.nmbRapidityBins();
	_VMnumega=inputParametersInstance.nmbEnergyBins();
	_VMnumQ2=inputParametersInstance.nmbGammaQ2Bins(); 
	_VMgamma_em=inputParametersInstance.beamLorentzGamma();
	_VMinterferencemode=inputParametersInstance.interferenceEnabled();
	_VMbslope=0.;//Will define in wide/narrow constructor
        _bslopeDef=inputParametersInstance.bslopeDefinition();
	_bslopeVal=inputParametersInstance.bslopeValue();
	_pEnergy= inputParametersInstance.protonEnergy();
	_beamNucleus = inputParametersInstance.targetBeamA();
	// electron energy in CMS frame
	_eEnergy= inputParametersInstance.electronEnergy();
	_VMpidtest=inputParametersInstance.prodParticleType();
	_VMptmax=inputParametersInstance.maxPtInterference();
	_VMdpt=inputParametersInstance.ptBinWidthInterference();
        _ProductionMode=inputParametersInstance.productionMode();
	_targetMaxPhotonEnergy=inputParametersInstance.targetMaxPhotonEnergy();
	_targetMinPhotonEnergy=inputParametersInstance.targetMinPhotonEnergy();
	// Now saving the photon energy limits
	_cmsMaxPhotonEnergy=inputParametersInstance.cmsMaxPhotonEnergy();
	_cmsMinPhotonEnergy=inputParametersInstance.cmsMinPhotonEnergy();
	_beamLorentzGamma = inputParametersInstance.beamLorentzGamma();
	_targetBeamLorentzGamma = inputParametersInstance.targetBeamLorentzGamma();
	_rap_CM=inputParametersInstance.rap_CM();
	_targetRadius = inputParametersInstance.targetRadius();
	//Turn on/off backwards production
	_backwardsProduction = inputParametersInstance.backwardsProduction();
	
        N0 = 0; N1 = 0; N2 = 0; 
	if (_VMpidtest == starlightConstants::FOURPRONG){
	  // create n-body phase-spage generator
	  _phaseSpaceGen = new nBodyPhaseSpaceGen(_randy);
	}
	if(_ProductionMode == 12 || _ProductionMode == 13) // Narrow and wide coherent photon-pommeron in eSTARlight
	  _dummy_pncs = new photonNucleusCrossSection(inputParametersInstance, bbsystem);
	//
	const double r_04_00 = 0.674;
	const double cos_delta = 0.925;
	for( int ii = 0; ii < 100; ++ii){ //epsilon 0-1
	  double epsilon = 0.01*ii;
	  const double R = (1./epsilon)*r_04_00/(1.-r_04_00);
	  for(int jj = 0; jj < 200; ++jj){ //psi 0 - 2pi
	    double psi = jj*starlightConstants::pi/100.;
	    double max_bin;
	    for( int kk = 0; kk < 200; ++kk){ //temp
	      double theta = kk*starlightConstants::pi/100.;
	      // Fin max
	      double this_test = std::pow(std::sin(theta),2.)*(1+epsilon*cos(2.*psi)) + 2.*epsilon*R*std::pow(std::cos(theta),2.)
		+std::sqrt(2.*epsilon*(1+epsilon))*cos_delta*std::sin(2.*theta)*std::cos(psi);
	      if(this_test >  max_bin)
		max_bin = this_test;
	    }
	    _angular_max[ii][jj] = max_bin;
	  }
	}
}


//______________________________________________________________________________
Gammaavectormeson::~Gammaavectormeson()
{
	if (_phaseSpaceGen)
		delete _phaseSpaceGen;
	if (_dummy_pncs)
	  delete _dummy_pncs;
}


//______________________________________________________________________________
void Gammaavectormeson::pickwy(double &W, double &Y)
{
        double dW, dY, xw,xy,xtest;
	int  IW,IY;
  
	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	dY = (_VMYmax-_VMYmin)/double(_VMnumy);
  
 L201pwy:

	xw = _randy.Rndom();
	W = _VMWmin + xw*(_VMWmax-_VMWmin);

	if (W < 2 * starlightConstants::pionChargedMass)
		goto L201pwy;
  
	IW = int((W-_VMWmin)/dW);
	xy = _randy.Rndom();
	Y = _VMYmin + xy*(_VMYmax-_VMYmin);
	IY = int((Y-_VMYmin)/dY); 
	xtest = _randy.Rndom();

	if( xtest > _Farray[IW][IY] )
		goto L201pwy;

        N0++; 
	// Determine the target nucleus 
	// For pA this is given, for all other cases use the relative probabilities in _Farray1 and _Farray2 
	_TargetBeam=2; // Always true for eX
}         



//______________________________________________________________________________                                               
void Gammaavectormeson::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,
                                     double  W,
                                     double  px0, double  py0, double  pz0,
				     double  spin_element,
                                     double& px1, double& py1, double& pz1,
                                     double& px2, double& py2, double& pz2,
                                     int&    iFbadevent)
{
	// This routine decays a particle into two particles of mass mdec,
	// taking spin into account
	double pmag;
	double phi,theta,Ecm;
	double betax,betay,betaz;
	double mdec1=0.0,mdec2=0.0;
	double E1=0.0,E2=0.0;

	//    set the mass of the daughter particles
	mdec1=getDaughterMass(ipid);

	if(_VMpidtest == starlightConstants::OMEGA_pi0gamma){
		mdec2=starlightConstants::pionNeutralMass;
		if(W < mdec2){
			cout<<" ERROR: W="<<W<<endl;
			iFbadevent = 1;
			return;
		}
		//  calculate the magnitude of the momenta
		pmag = (W*W - mdec2*mdec2)/(2.0*W);
	}else{
		mdec2=mdec1;
		if(W < 2*mdec1){
			cout<<" ERROR: W="<<W<<endl;
			iFbadevent = 1;
			return;
		}
		//  calculate the magnitude of the momenta
		pmag = sqrt(W*W/4. - mdec2*mdec2);
	}
	  
	//  pick an orientation, based on the spin
	//  phi has a flat distribution in 2*pi
	phi = _randy.Rndom()*2.*starlightConstants::pi;
                                                                                                                
	//  find theta, the angle between one of the outgoing particles and
	//  the beamline, in the frame of the two photons
    //cout<<"spin element: "<<spin_element<<endl;
	theta=getTheta(ipid, spin_element);
 
	//  compute unboosted momenta
	px1 = sin(theta)*cos(phi)*pmag;
	py1 = sin(theta)*sin(phi)*pmag;
	pz1 = cos(theta)*pmag;
	px2 = -px1;
	py2 = -py1;
	pz2 = -pz1;

	Ecm = sqrt(W*W+px0*px0+py0*py0+pz0*pz0);
	E1 = sqrt(mdec1*mdec1+px1*px1+py1*py1+pz1*pz1);
	E2 = sqrt(mdec2*mdec2+px2*px2+py2*py2+pz2*pz2);

	//cout<<"Ecm: "<<Ecm<<endl;
	//cout<<"W: "<<W<<endl;

	betax = -(px0/Ecm);
	betay = -(py0/Ecm);
	betaz = -(pz0/Ecm);

	transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
	transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);

	if(iFbadevent == 1)
	   return;

}



// Code for returning a point in the Dalitz plot corresponding to a 3 body decay
// Draws uniformly from the plot
void Gammaavectormeson::dalitzPlot
(double M, double m1, double m2, double m3, double &s12, double &s23)
{

//M is mass of parent particle
//m_i is mass of decay products
double S = pow(M,2);
bool dalitzFound = false;

while (!dalitzFound){

	//generates random numbers between 0 and .5, correspnding to (x,y) on Dalitz plot
	double x =  _randy.Rndom()/2;
	double y =  _randy.Rndom()/2;

	
	if (x==0) continue;

	double lamb_1 = pow(x,2)+pow(S,2)+pow(m3,4)-2*(x*S+x*pow(m3,2)+pow(m3,2)*S);
	double lamb_2 = pow(x,2)+pow(m1,4)+pow(m2,4)-2*(x*pow(m1,2)+x*pow(m2,2)+pow(m2,2)*pow(m1,2));

	double y_max = pow(m1,2)+pow(m3,2)+(1/(2*x)*((S-x-pow(m1,2))*(x-pow(m3,2)+pow(m2,2))+sqrt(lamb_1*lamb_2)));
	double y_min = pow(m1,2)+pow(m3,2)+(1/(2*x)*((S-x-pow(m1,2))*(x-pow(m3,2)+pow(m2,2))-sqrt(lamb_1*lamb_2)));

		//The +/-.01 terms are to prevent this from chooisng points that lie on the boundary plot surface, 
		//which can lead to unphysical results 
	if (x>pow(m1+m1,2)+.01 && x<pow(M-m3,2)-.01){
		if (y>y_min+.01 && y<y_max-.01)
		{
			s12 = x; 
			s23 = y;

			dalitzFound = true;
		}

		else continue;
	}

else continue;

}


return;
}




//Simulates omega -> Pi+Pi-P0 decay
void Gammaavectormeson::threePionDecay
(                                    starlightConstants::particleTypeEnum &ipid,
				     double  px0, double  py0, double  pz0,
                                     double& px1, double& py1, double& pz1,
                                     double& px2, double& py2, double& pz2,
                                     double& px3, double& py3, double& pz3,
                                     int&    iFbadevent,
                                     double spin_element)
{


ipid = starlightConstants::PION;

//Mass of omega particle 
double M_omega = starlightConstants::OmegaMass;

//Masses of Pi_0,Pi_plus/Pi_minus
double m_pi0 = starlightConstants::pionNeutralMass;
double m_pip = starlightConstants::pionChargedMass;
double m_pim = m_pip;

//Dalitz plot variables 
double s12 = 0;
double s23 = 0;


dalitzPlot(M_omega,m_pip,m_pim,m_pi0,s12,s23);


double E0 = M_omega;


double E1 = (E0*E0+m_pip*m_pip-s23)/(2*E0);
double E2 = ((s12+s23)-(m_pip*m_pip+m_pi0*m_pi0))/(2*E0);
double E3 = (E0*E0+m_pi0*m_pi0-s12)/(2*E0);


//Computes a random althimuthal angle for the decay in the x-axis
double az = _randy.Rndom()*2.*starlightConstants::pi;


double p1mag = sqrt(E1*E1-m_pip*m_pip);
double p2mag = sqrt(E2*E2-m_pim*m_pim);
double p3mag = sqrt(E3*E3-m_pi0*m_pi0);



//Calculates the angle between particle 1 and 2 in the frame in which particle 3 has angle 0
double a12 = acos(-(s12-pow((E1+E2),2)+p1mag*p1mag+p2mag*p2mag)/(2*p1mag*p2mag));


//Calculates the angle between particle 2 and 3 in the frame in which particle 3 has angle 0
double a23 = acos(-(s23-pow((E3+E2),2)+p3mag*p3mag+p2mag*p2mag)/(2*p2mag*p3mag));



// if (isnan(a12)) 
// { 	cout<<"BAD ANGLE"<<endl;
// 	cout<<-(s12-pow((E1+E2),2)+p1mag*p1mag+p2mag*p2mag)/(2*p1mag*p2mag)<<endl;
// 	if(-(s12-pow((E1+E2),2)+p1mag*p1mag+p2mag*p2mag)/(2*p1mag*p2mag)>1){a12 = 0;}
// 	else if (-(s12-pow((E1+E2),2)+p1mag*p1mag+p2mag*p2mag)/(2*p1mag*p2mag)<-1){a12 = starlightConstants::pi;}
// 	else{a12 = -starlightConstants::pi/2;}


// }


// if (isnan(a23)) 
// { 	cout<<"BAD ANGLE"<<endl;
// 	cout<<-(s23-pow((E3+E2),2)+p3mag*p3mag+p2mag*p2mag)/(2*p2mag*p3mag)<<endl;
// 	if(-(s23-pow((E3+E2),2)+p3mag*p3mag+p2mag*p2mag)/(2*p2mag*p3mag)>1){a23 = 0;}
// 	else if (-(s23-pow((E3+E2),2)+p3mag*p3mag+p2mag*p2mag)/(2*p2mag*p3mag)<-1){a23 = starlightConstants::pi;}
// 	else{a23 = starlightConstants::pi/2;}

// }



//Calculates the angle of particle 1 and 2 in this frame
double unrotated_phi2 = -(starlightConstants::pi-a23);
double unrotated_phi1 = a12+unrotated_phi2;


//Calulates the angle particle 3 makes with the z-axis
double theta = getTheta(starlightConstants::PION, spin_element);


//Make theta such that the projected particle lies on the negative x axis
if(theta>=0 && theta<=starlightConstants::pi){theta = -(starlightConstants::pi-theta);}
if (theta>starlightConstants::pi && theta<=2*starlightConstants::pi){theta = -(theta-starlightConstants::pi);}


px3 =  -p3mag*sin(az)*sin(theta);
py3 =   p3mag*sin(theta)*cos(az);
pz3 =  -p3mag*cos(theta);


//We have to relate the azimuthal angle we rotated Pi0 by to the angle needed to rotate the charged pions in order to ensure momentum conseveration in the x direction

//This quanitiy is the sine of the angle above 
double x_mom_angle = -px3/(p1mag*(-cos(unrotated_phi1)*sin(theta)+sin(unrotated_phi1)*cos(theta))+p2mag*(-cos(unrotated_phi2)*sin(theta)+sin(unrotated_phi2)*cos(theta)));

//Cosine of the angle needed to rotate the charged pions in order to ensure momentum conseveration
//in the y direction
//This will of course be the same as above, modulo possible phase shifts that led to some of y momenta being off by a neagtive sign 
double y_mom_factor = -py3/((p1mag*(-cos(unrotated_phi1)*sin(theta)+sin(unrotated_phi1)*cos(theta))+p2mag*(-cos(unrotated_phi2)*sin(theta)+sin(unrotated_phi2)*cos(theta))));



//Computed unboosted momenta
px1 = p1mag*x_mom_angle*(-cos(unrotated_phi1)*sin(theta)+sin(unrotated_phi1)*cos(theta));
py1 = p1mag*(-cos(unrotated_phi1)*sin(theta)+sin(unrotated_phi1)*cos(theta))*y_mom_factor;
pz1 = p1mag*(cos(theta)*cos(unrotated_phi1)+(sin(theta)*sin(unrotated_phi1)));
//py1 = sqrt(pow(p1mag,2)-pow(px1,2)-pow(pz1,2));


px2 =  p2mag*x_mom_angle*(-cos(unrotated_phi2)*sin(theta)+sin(unrotated_phi2)*cos(theta));
py2 =  p2mag*(-cos(unrotated_phi2)*sin(theta)+sin(unrotated_phi2)*cos(theta))*y_mom_factor;
pz2 =  p2mag*(cos(theta)*cos(unrotated_phi2)+(sin(theta)*sin(unrotated_phi2)));



//  cout<<"Conservation of total momentum check"<<endl;
//  cout<<"Net mom_x = "<<(px1+px2+px3)<<endl;
//  cout<<"Net mom_y = "<<(py1+py2+py3)<<endl;
//  cout<<"Net mom_z = "<<(pz1+pz2+pz3)<<endl;

//  cout<<"Conservation of total energy check"<<endl;
//  cout<<"Net energy = " <<(M_omega-(E1+E2+E3))<<endl;
// cout<<"###############"<<endl;



double  Ecm = sqrt(M_omega*M_omega+px0*px0+py0*py0+pz0*pz0);
double	betax = -(px0/Ecm);
double	betay = -(py0/Ecm);
double	betaz = -(pz0/Ecm);




//Transforms boosted momenta and energy
transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);
transform (betax,betay,betaz,E3,px3,py3,pz3,iFbadevent);

// cout<<"Energy conseveration after transformation"<<endl;
// cout<<"Pi+:"<<E1-sqrt(pow(starlightConstants::pionChargedMass,2)+px1*px1+py1*py1+pz1*pz1)<<endl; 
// cout<<"Pi-:"<<E2-sqrt(pow(starlightConstants::pionChargedMass,2)+px2*px2+py2*py2+pz2*pz2)<<endl; 
// cout<<"Pi0:"<<E3-sqrt(pow(starlightConstants::pionNeutralMass,2)+px3*px3+py3*py3+pz3*pz3)<<endl; 

//  cout<<"####"<<endl;


//Decays the pi0 into two photons
// double e_g1, px_g1, py_g1, pz_g1, e_g2, px_g2, py_g2, pz_g2;
// pi0Decay(px3,py3,pz3,e_g1,px_g1,py_g1,pz_g1,e_g2,px_g2,py_g2,pz_g2,iFbadevent);

	if(iFbadevent == 1)
	   return;
	
}




//______________________________________________________________________________                                               
// decays a particle into four particles with isotropic angular distribution
bool Gammaavectormeson::fourBodyDecay
(starlightConstants::particleTypeEnum& ipid,
 const double                  E,           // E (unused)
 const double                  W,          // mass of produced particle
 const double*                 p,          // momentum of produced particle; expected to have size 3
 lorentzVector*                decayVecs,  // array of Lorentz vectors of daughter particles; expected to have size 4
 int&                          iFbadevent)
{
	const double parentMass = W;

	// set the mass of the daughter particles
	const double daughterMass = getDaughterMass(ipid);
	if (parentMass < 4 * daughterMass){
		cout << " ERROR: W=" << parentMass << " GeV too small" << endl;
		iFbadevent = 1;
		return false;
	}

	// construct parent four-vector
	const double        parentEnergy = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
	                                        + parentMass * parentMass);
	const lorentzVector parentVec(p[0], p[1], p[2], parentEnergy);

	// setup n-body phase-space generator
	assert(_phaseSpaceGen);
	static bool firstCall = true;
	if (firstCall) {
		const double m[4] = {daughterMass, daughterMass, daughterMass, daughterMass};
		_phaseSpaceGen->setDecay(4, m);
		// estimate maximum phase-space weight
		_phaseSpaceGen->setMaxWeight(1.01 * _phaseSpaceGen->estimateMaxWeight(_VMWmax));
		firstCall = false;
	}

	// generate phase-space event
	if (!_phaseSpaceGen->generateDecayAccepted(parentVec))
		return false;

	// set Lorentzvectors of decay daughters
	for (unsigned int i = 0; i < 4; ++i)
		decayVecs[i] = _phaseSpaceGen->daughter(i);
	return true;
}

//______________________________________________________________________________                                               
void Gammaavectormeson::pi0Decay(double& px_pi0, double& py_pi0, double& pz_pi0,
                                 double& e_g1, double& px_g1, double& py_g1, double& pz_g1,
                                 double& e_g2, double& px_g2, double& py_g2, double& pz_g2,
                                 int&    iFbadevent)
{
	double pmag;
	double phi,theta,Ecm;
	double betax,betay,betaz;
	double E1=0.0,E2=0.0;

	// This routine decays a pi0 into two isotropically produced photons
	double m_pi0=starlightConstants::pionNeutralMass;
	
	//  calculate the magnitude of the momenta
	//pmag = sqrt(m_pi0*m_pi0/4.);
	  
	pmag = m_pi0/2;  
	//  pick an orientation, based on the spin
	//  phi has a flat distribution in 2*pi
	phi = _randy.Rndom()*2.*starlightConstants::pi;
                                                                                                                
	//  find theta, the angle between one of the outgoing photons and
	//  the beamline, in the frame of the pi0
	theta=getTheta(starlightConstants::PION, 1.0/3.0);
 
	//  compute unboosted momenta
	px_g1 = sin(theta)*cos(phi)*pmag;
	py_g1 = sin(theta)*sin(phi)*pmag;
	pz_g1 = cos(theta)*pmag;
	px_g2 = -px_g1;
	py_g2 = -py_g1;
	pz_g2 = -pz_g1;

	Ecm = sqrt(m_pi0*m_pi0+px_pi0*px_pi0+py_pi0*py_pi0+pz_pi0*pz_pi0);
	E1 = sqrt(px_g1*px_g1+py_g1*py_g1+pz_g1*pz_g1);
	E2 = sqrt(px_g2*px_g2+py_g2*py_g2+pz_g2*pz_g2);

	betax = -(px_pi0/Ecm);
	betay = -(py_pi0/Ecm);
	betaz = -(pz_pi0/Ecm);


	transform (betax,betay,betaz,E1,px_g1,py_g1,pz_g1,iFbadevent);
	transform (betax,betay,betaz,E2,px_g2,py_g2,pz_g2,iFbadevent);

	e_g1=E1;
	e_g2=E2;

	if(iFbadevent == 1)
	   return;
}


//______________________________________________________________________________
double Gammaavectormeson::getDaughterMass(starlightConstants::particleTypeEnum &ipid)
{
	//This will return the daughter particles mass, and the final particles outputed id...
	double mdec=0.;
  
	switch(_VMpidtest){
	case starlightConstants::RHO:
	case starlightConstants::RHOZEUS:
	case starlightConstants::FOURPRONG:
	case starlightConstants::OMEGA:
		mdec = starlightConstants::pionChargedMass;
		ipid = starlightConstants::PION;
		break;
	case starlightConstants::OMEGA_pi0gamma:
		mdec = 0.0;
		ipid = starlightConstants::PHOTON;
		break;
	case starlightConstants::PHI:
		mdec = starlightConstants::kaonChargedMass;
		ipid = starlightConstants::KAONCHARGE;
		break;
	case starlightConstants::JPSI:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break; 
	case starlightConstants::JPSI2S_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI2S_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break; 

	case starlightConstants::JPSI2S:
	case starlightConstants::UPSILON:
	case starlightConstants::UPSILON2S:
	case starlightConstants::UPSILON3S:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break;
	case starlightConstants::UPSILON_ee:
	case starlightConstants::UPSILON2S_ee:
	case starlightConstants::UPSILON3S_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break;
	case starlightConstants::UPSILON_mumu:
	case starlightConstants::UPSILON2S_mumu:
	case starlightConstants::UPSILON3S_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;   
		break;
	default: cout<<"No daughtermass defined, gammaavectormeson::getdaughtermass"<<endl;
	}
  
	return mdec;
}



//______________________________________________________________________________
double Gammaavectormeson::getTheta(starlightConstants::particleTypeEnum ipid, double r_04_00)
{
	//This depends on the decay angular distribution
	//Valid for rho, phi, omega.
	double theta=0.;
	double xtest=1.;
	double dndtheta=0.;

	while(xtest > dndtheta){
    
	  theta = starlightConstants::pi*_randy.Rndom();
	  xtest = _randy.Rndom();
	  //  Follow distribution for helicity +/-1 with finite Q2
	  //  Physics Letters B 449, 328; The European Physical Journal C - Particles and Fields 13, 37; 
	  //  The European Physical Journal C - Particles and Fields 6, 603
  	

	  switch(ipid){
	    
	  case starlightConstants::MUON:
	  case starlightConstants::ELECTRON:
	    //primarily for upsilon/j/psi.  VM->ee/mumu
	    dndtheta = sin(theta)*(1 + r_04_00+( 1-3.*r_04_00 )*cos(theta)*cos(theta));
	    break;
	    
	  case starlightConstants::PION:
	  	dndtheta=  sin(theta)*(1 - r_04_00+( 3.*r_04_00-1 )*cos(theta)*cos(theta));
	  	break;

	  case starlightConstants::KAONCHARGE:
	    //rhos etc
	    dndtheta=  sin(theta)*(1 - r_04_00+( 3.*r_04_00-1 )*cos(theta)*cos(theta));
	    break;
	    
	  default: if(!_backwardsProduction) 
	  {
	  	cout<<"No proper theta dependence defined, check gammaavectormeson::gettheta"<<endl;
			exit( 3 ); //User needs to handle this case
			} 

	  }//end of switch

	  // Assume unpolarized vector-mesons for now
	  if(_backwardsProduction){
	  	dndtheta = sin(theta);
	  }
	}

	return theta;

}


//______________________________________________________________________________
double Gammaavectormeson::getSpinMatrixElement(double W, double Q2, double epsilon)
{
  double m2 = 0.6*(W*W);
  double R = starlightConstants::pi/2.*m2/Q2 - std::pow(m2,3./2.)/sqrt(Q2)/(Q2+m2) - m2/Q2*atan(sqrt(m2/Q2));
  //////////needed for problem with double precision
  double R1 = starlightConstants::pi/2.*m2/Q2;
  if(R1>1.0E10)R=0.0;
  ////////// 
  R = (Q2+m2)*(Q2+m2)*R*R/m2;
  return epsilon*R/(1+epsilon*R);
}


//______________________________________________________________________________
double Gammaavectormeson::getWidth()
{
	return _width;
}


//______________________________________________________________________________
double Gammaavectormeson::getMass()
{
	return _mass;
}


//______________________________________________________________________________
double Gammaavectormeson::getSpin()
{
	return 1.0; //VM spins are the same
}

//______________________________________________________________________________
void Gammaavectormeson::momenta(double W,double Egam,double Q2, double gamma_pz, double gamma_pt, //input conditions
				double &Y,double &E,double &px,double &py,double &pz,  //return vm
				double &t_px, double &t_py, double &t_pz, double &t_E, //return pomeron
				double &e_phi,int &tcheck) //return electron (angle already known by Q2)
{
	//     This subroutine calculates momentum and energy of vector meson for electroproduction (eSTARlight)
	//     given W and photon 4-vector,   without interference.  No intereference in asymetric eX collisions
 
	double Epom,Pom_pz,tmin,pt2,phi1,phi2;
	double px1,py1;
	double xt,xtest,ytest;
	double t2;

	double target_px, target_py, target_pz, target_E;

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
	// Pomeron pz is != than its energy in eSTARlight, in order to conserve energy/momentum of scattered
	// target
    //Pom_pz = 0.5*(W*W-Q2)/(Egam + gamma_pz);
    Epom = 0.5*(W*W+Q2)/(Egam + target_E);

    L522vm:
	while( e_phi > 2.*starlightConstants::pi ) e_phi-= 2.*starlightConstants::pi;
	//
	if ( ( _bbs.targetBeam().A()==1 ) 
          || (_ProductionMode == 4) ) {
	    if( (_VMpidtest == starlightConstants::RHO) || (_VMpidtest == starlightConstants::RHOZEUS) || (_VMpidtest == starlightConstants::OMEGA)){
	      // Use dipole form factor for light VM
	    L613vm:
	      xtest = 2.0*_randy.Rndom();
              double ttest = xtest*xtest; 
              ytest = _randy.Rndom();
              double t0 = 1./2.23; 
              double yprob = xtest*_bbs.electronBeam().dipoleFormFactor(ttest,t0)*_bbs.electronBeam().dipoleFormFactor(ttest,t0); 
              if( ytest > yprob ) goto L613vm; 
              t2 = ttest; 
              pt2 = xtest;              
	    }else{
		//Use dsig/dt= exp(-_VMbslope*t) for heavy VM
                double bslope_tdist = _VMbslope; 
		double Wgammap = 0.0; 
                switch(_bslopeDef){
		  case 0:
		    //This is the default, as before
		    bslope_tdist = _VMbslope;
		    break;
		  case 1:
		    //User defined value of bslope. BSLOPE_VALUE default is 4.0 if not set. 
                    bslope_tdist = _bslopeVal;
		    if( N0 <= 1 )cout<<" ATTENTION: Using user defined value of bslope = "<<_bslopeVal<<endl;
                    break; 
		  case 2:
                    //This is Wgammap dependence of b from H1 (Eur. Phys. J. C 46 (2006) 585)
		    Wgammap = sqrt(4.*Egam*_pEnergy); 
		    bslope_tdist = 4.63 + 4.*0.164*log(Wgammap/90.0);
		    if( N0 <= 1 )cout<<" ATTENTION: Using energy dependent value of bslope!"<<endl; 
		    break;
		  default:
		    cout<<" Undefined setting for BSLOPE_DEFINITION "<<endl;
		}

	        xtest = _randy.Rndom(); 
		// t2 = (-1./_VMbslope)*log(xtest);
		t2 = (-1./bslope_tdist)*log(xtest);
		pt2 = sqrt(1.*t2);
	    }
	} else {
	    // >> Check tmin
	    tmin = ((Epom/_VMgamma_em)*(Epom/_VMgamma_em));

	    if(tmin > 0.5){
		cout<<" WARNING: tmin= "<<tmin<<endl;
                cout<< " Y = "<<Y<<" W = "<<W<<" Epom = "<<Epom<<" gamma = "<<_VMgamma_em<<endl; 
		cout<<" Will pick a new W,Y "<<endl;
		tcheck = 1;
		return;
	    }
 L663vm:
	    xt = _randy.Rndom(); 
            if( _bbs.targetBeam().A() != 1){ 
	      pt2 = 8.*xt*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius();
	    } 
	    else{
	      std::cout<<"Can't find the electron for eX"<<std::endl;
	    }  

	    xtest = _randy.Rndom();
	    t2 = tmin + pt2*pt2;

	    double comp=0.0; 
            if( _bbs.targetBeam().A() != 1){ 
	      comp = _bbs.targetBeam().formFactor(t2)*_bbs.targetBeam().formFactor(t2)*pt2;
	    }
	    else 
	      std::cout<<"Can't find the electron for eX"<<std::endl;
            if( xtest > comp ) goto L663vm;
       		
	}
	phi2 = 2.*starlightConstants::pi*_randy.Rndom();

	//
	t_px = pt2*cos(phi2);
	t_py = pt2*sin(phi2);
	// Used to return the pomeron pz to generator
	// Compute scattered target kinematics p_(initial ion) = p_(pomeron) + p_(final ion)
	double newion_E = target_E-Epom;
	double newion_px = target_px - t_px;
	double newion_py = target_py - t_py;
	double newion_pz_squared = newion_E*newion_E - newion_px*newion_px - newion_py*newion_py - pow(_beamNucleus*starlightConstants::protonMass,2.);
	if(newion_pz_squared<0)goto L522vm;
	double newion_pz = -sqrt( newion_pz_squared );
	Pom_pz = target_pz - newion_pz;
	t_pz = Pom_pz;
	// Compute vector sum Pt = Pt1 + Pt2 to find pt for the vector meson
	px = px1 + t_px;
	py = py1 + t_py;

	t_E = Epom;
	// Finally V.M. energy, pz and rapidity from photon + pommeron.
	E = Egam + t_E;
	pz = gamma_pz + t_pz;
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
	//	cout << "pz_final = " << pz << endl;


	// Backwards Production... updated by Zachary Sweger on 9/21/2021
	if (_backwardsProduction) {
	  double is_proton_E  = target_E;
	  double is_proton_px = target_px;
	  double is_proton_py = target_py;
	  double is_proton_pz = target_pz;
	  double is_gamma_E  = Egam;
	  double is_gamma_px = px1;
	  double is_gamma_py = py1;
	  double is_gamma_pz = gamma_pz;
	  double is_tot_E  = is_proton_E  + is_gamma_E;
	  double is_tot_px = is_proton_px + is_gamma_px;
	  double is_tot_py = is_proton_py + is_gamma_py;
	  double is_tot_pz = is_proton_pz + is_gamma_pz;

	  double fs_proton_E  = newion_E;
	  double fs_proton_px = newion_px;
	  double fs_proton_py = newion_py;
	  double fs_proton_pz = newion_pz;
	  
	  double fs_vm_E  = E;
	  double fs_vm_px = px;
	  double fs_vm_py = py;
	  double fs_vm_pz = pz;

      int isbadevent = 0;
	  double betax_cm = (is_tot_px/is_tot_E);
      double betay_cm = (is_tot_py/is_tot_E);
      double betaz_cm = (is_tot_pz/is_tot_E);
      transform (betax_cm,betay_cm,betaz_cm,is_proton_E,is_proton_px,is_proton_py,is_proton_pz,isbadevent);
      transform (betax_cm,betay_cm,betaz_cm,is_gamma_E, is_gamma_px, is_gamma_py, is_gamma_pz, isbadevent);
      transform (betax_cm,betay_cm,betaz_cm,fs_proton_E,fs_proton_px,fs_proton_py,fs_proton_pz,isbadevent);
      transform (betax_cm,betay_cm,betaz_cm,fs_vm_E,    fs_vm_px,    fs_vm_py,    fs_vm_pz,    isbadevent);

      double u_fs_proton_px = -1.0*fs_proton_px;
	  double u_fs_proton_py = -1.0*fs_proton_py;
	  double u_fs_proton_pz = -1.0*fs_proton_pz;
	  double u_fs_proton_E  = fs_proton_E;
	  
	  double u_fs_vm_px = -1.0*fs_vm_px;
	  double u_fs_vm_py = -1.0*fs_vm_py;
	  double u_fs_vm_pz = -1.0*fs_vm_pz;
	  double u_fs_vm_E  = fs_vm_E;

      betax_cm = -1.0*betax_cm;
	  betay_cm = -1.0*betay_cm;
	  betaz_cm = -1.0*betaz_cm;
      transform (betax_cm,betay_cm,betaz_cm,is_proton_E,  is_proton_px,  is_proton_py,  is_proton_pz,  isbadevent);
      transform (betax_cm,betay_cm,betaz_cm,is_gamma_E,   is_gamma_px,   is_gamma_py,   is_gamma_pz,   isbadevent);
      transform (betax_cm,betay_cm,betaz_cm,u_fs_proton_E,u_fs_proton_px,u_fs_proton_py,u_fs_proton_pz,isbadevent);
      transform (betax_cm,betay_cm,betaz_cm,u_fs_vm_E,    u_fs_vm_px,    u_fs_vm_py,    u_fs_vm_pz,    isbadevent);

      px=u_fs_vm_px;
      py=u_fs_vm_py;
      pz=u_fs_vm_pz;
      E =u_fs_vm_E;

      t_px=-u_fs_proton_px;
      t_py=-u_fs_proton_py;
      t_pz=is_proton_pz-u_fs_proton_pz;
      t_E =is_proton_E-u_fs_proton_E;
      //cout<<"Energy Violation: "<<(E+u_fs_proton_E)-(is_tot_E)<<endl;
      //cout<<"Pz Violation: "<<(pz+u_fs_proton_pz)-(is_tot_pz)<<endl;
      //cout<<"Px Violation: "<<(px+u_fs_proton_px)-(is_tot_px)<<endl;
      //cout<<"Py Violation: "<<(py+u_fs_proton_py)-(is_tot_py)<<endl;
      //cout<<"proton mass: "<< sqrt(u_fs_proton_E*u_fs_proton_E-u_fs_proton_pz*u_fs_proton_pz-u_fs_proton_py*u_fs_proton_py-u_fs_proton_px*u_fs_proton_px)<<endl;
      //cout << "VM MASS = " << sqrt(E*E-px*px - py*py-pz*pz) << endl;

	  //Calculate the rapidity
	  Y = 0.5*std::log( (E+fabs(pz))/(E-fabs(pz)) );
	  if(Y!=Y){
	  	//FIXME the following should be upated if not using omegas
	  	if(_VMpidtest==starlightConstants::OMEGA_pi0gamma || _VMpidtest==starlightConstants::OMEGA) E=sqrt(starlightConstants::OmegaMass*starlightConstants::OmegaMass+px*px+py*py+pz*pz); 
		else if(_VMpidtest==starlightConstants::RHO) E=sqrt(starlightConstants::rho0Mass*starlightConstants::rho0Mass+px*px+py*py+pz*pz);
		//cout<<"Energy Violation: "<<(E+u_fs_proton_E)-(is_tot_E)<<endl;
	  }
	}
	
}


//______________________________________________________________________________
double Gammaavectormeson::pTgamma(double E)
{
    // returns on random draw from pp(E) distribution
    double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
    double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
    int satisfy =0;
        
    ereds = (E/_VMgamma_em)*(E/_VMgamma_em);
    //sqrt(3)*E/gamma_em is p_t where the distribution is a maximum
    Cm = sqrt(3.)*E/_VMgamma_em;
    // If E is very small, the drawing of a pT below is extremely slow. 
    // ==> Set pT = sqrt(3.)*E/_VMgamma_em for very small E. 
    // Should have no observable consequences (JN, SRK 11-Sep-2014)
    if( E < 0.0005 )return Cm; 
 
    //the amplitude of the p_t spectrum at the maximum

    if( _bbs.electronBeam().A()==1 && _bbs.targetBeam().A() != 1){ 
      if( _ProductionMode == 2 || _ProductionMode ==3 ){
         singleformfactorCm=_bbs.electronBeam().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.targetBeam().formFactor(Cm*Cm+ereds);
      }  
    } else if( _bbs.targetBeam().A()==1 && _bbs.electronBeam().A() != 1 ){
      if( _ProductionMode == 2 || _ProductionMode ==3){
         singleformfactorCm=_bbs.targetBeam().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.electronBeam().formFactor(Cm*Cm+ereds);
      }  
    } else if (_TargetBeam == 1) {
      singleformfactorCm=_bbs.targetBeam().formFactor(Cm*Cm+ereds);
    } else {
      singleformfactorCm=_bbs.electronBeam().formFactor(Cm*Cm+ereds);
    }

    Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
    //pick a test value pp, and find the amplitude there
    x = _randy.Rndom();

    if( _bbs.electronBeam().A()==1 && _bbs.targetBeam().A() != 1){ 
      if( _ProductionMode == 2 || _ProductionMode ==3){
        pp = x*5.*starlightConstants::hbarc/_bbs.electronBeam().nuclearRadius(); 
        singleformfactorpp1=_bbs.electronBeam().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius(); 
        singleformfactorpp1=_bbs.targetBeam().formFactor(pp*pp+ereds);
      }  
    } else if( _bbs.targetBeam().A()==1 && _bbs.electronBeam().A() != 1 ){
      if( _ProductionMode == 2 || _ProductionMode ==3){
        pp = x*5.*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius(); 
        singleformfactorpp1=_bbs.targetBeam().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.electronBeam().nuclearRadius(); 
        singleformfactorpp1=_bbs.electronBeam().formFactor(pp*pp+ereds);
      }  
    } else if (_TargetBeam == 1) {
        pp = x*5.*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius(); 
        singleformfactorpp1=_bbs.electronBeam().formFactor(pp*pp+ereds);
    } else {
        pp = x*5.*starlightConstants::hbarc/_bbs.electronBeam().nuclearRadius(); 
        singleformfactorpp1=_bbs.electronBeam().formFactor(pp*pp+ereds);
    }

    test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

    while(satisfy==0){
	u = _randy.Rndom();
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =_randy.Rndom();
            if( _bbs.electronBeam().A()==1 && _bbs.targetBeam().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pp = x*5.*starlightConstants::hbarc/_bbs.electronBeam().nuclearRadius(); 
                singleformfactorpp2=_bbs.electronBeam().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius(); 
                singleformfactorpp2=_bbs.targetBeam().formFactor(pp*pp+ereds);
              }  
            } else if( _bbs.targetBeam().A()==1 && _bbs.electronBeam().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pp = x*5.*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius(); 
                singleformfactorpp2=_bbs.targetBeam().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.electronBeam().nuclearRadius(); 
                singleformfactorpp2=_bbs.electronBeam().formFactor(pp*pp+ereds);
              }  
            } else if (_TargetBeam == 1) {
              pp = x*5.*starlightConstants::hbarc/_bbs.targetBeam().nuclearRadius(); 
              singleformfactorpp2=_bbs.electronBeam().formFactor(pp*pp+ereds);
            } else {
              pp = x*5.*starlightConstants::hbarc/_bbs.electronBeam().nuclearRadius(); 
              singleformfactorpp2=_bbs.electronBeam().formFactor(pp*pp+ereds);
            }
	    test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
	}
    }

    return pp;
}


//______________________________________________________________________________
void Gammaavectormeson::vmpt(double W,double Y,double &E,double &px,double &py, double &pz,
                             int&) // tcheck (unused)
{
	//    This function calculates momentum and energy of vector meson
	//    given W and Y, including interference.
	//    It gets the pt distribution from a lookup table.
	double dY=0.,yleft=0.,yfract=0.,xpt=0.,pt1=0.,ptfract=0.,pt=0.,pt2=0.,theta=0.;
	int IY=0,j=0;
  
	dY  = (_VMYmax-_VMYmin)/double(_VMnumy);
  
	//  Y is already fixed; choose a pt
	//  Follow the approach in pickwy
	//  in  _fptarray(IY,pt) IY=1 corresponds to Y=0, IY=numy/2 corresponds to +y
 	// Changed,  now works -y to +y.
	IY=int((Y-_VMYmin)/dY);
	if (IY > (_VMnumy)-1){
        	IY=(_VMnumy)-1;
	}

	yleft=(Y-_VMYmin)-(IY)*dY;

	yfract=yleft*dY;
  
	xpt=_randy.Rndom();
	for(j=0;j<_VMNPT+1;j++){
		if (xpt < _fptarray[IY][j]) goto L60;
	}
 L60:
  
	//  now do linear interpolation - start with extremes
  	if (j == 0){
		pt1=xpt/_fptarray[IY][j]*_VMdpt/2.;
		goto L80;
	}
	if (j == _VMNPT){
		pt1=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY][j])/(1.-_fptarray[IY][j]);
		goto L80;
	}
  
	//  we're in the middle
  	ptfract=(xpt-_fptarray[IY][j])/(_fptarray[IY][j+1]-_fptarray[IY][j]);
	pt1=(j+1)*_VMdpt+ptfract*_VMdpt;
  
	//  at an extreme in y?
	if (IY == (_VMnumy/2)-1){
		pt=pt1;
		goto L120;
	}
 L80:

	//  interpolate in y repeat for next fractional y bin      
	for(j=0;j<_VMNPT+1;j++){
		if (xpt < _fptarray[IY+1][j]) goto L90;
	}
 L90:
  
	//  now do linear interpolation - start with extremes
	if (j == 0){
		pt2=xpt/_fptarray[IY+1][j]*_VMdpt/2.;
		goto L100;
	}
	if (j == _VMNPT){
		pt2=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY+1][j])/(1.-_fptarray[IY+1][j]);
		goto L100;
	}
  
	//  we're in the middle
	ptfract=(xpt-_fptarray[IY+1][j])/(_fptarray[IY+1][j+1]-_fptarray[IY+1][j]);
	pt2=(j+1)*_VMdpt+ptfract*_VMdpt;
 L100:

	//  now interpolate in y  
	pt=yfract*pt2+(1-yfract)*pt1;
 L120:

	//  we have a pt 
	theta=2.*starlightConstants::pi*_randy.Rndom();
	px=pt*cos(theta);
	py=pt*sin(theta);

	E  = sqrt(W*W+pt*pt)*cosh(Y);
	pz = sqrt(W*W+pt*pt)*sinh(Y);
	//      randomly choose to make pz negative 50% of the time
	if(_randy.Rndom()>=0.5) pz = -pz;
}



//______________________________________________________________________________
double Gammaavectormeson::pseudoRapidity(double px, double py, double pz)
{
	double pT = sqrt(px*px + py*py);
	double p = sqrt(pz*pz + pT*pT);
	double eta = -99.9; if((p-pz) != 0){eta = 0.5*log((p+pz)/(p-pz));}
	return eta;
}

//______________________________________________________________________________
Gammaanarrowvm::Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaanarrowvm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaanarrowvm()"<<endl;
	narrowResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaanarrowvm::~Gammaanarrowvm()
{ }


//______________________________________________________________________________
Gammaaincoherentvm::Gammaaincoherentvm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
        cout<<"Reading in luminosity tables. Gammaainkoherentvm()"<<endl;
        read();
        cout<<"Creating and calculating crosssection. Gammaaincoherentvm()"<<endl;
        incoherentVMCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
        _VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaaincoherentvm::~Gammaaincoherentvm()
{ }


//______________________________________________________________________________
Gammaawidevm::Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaawidevm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaawidevm()"<<endl;
	wideResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter();
}


//______________________________________________________________________________
Gammaawidevm::~Gammaawidevm()
{ }


//______________________________________________________________________________
e_Gammaanarrowvm::e_Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. e_Gammaanarrowvm()"<<endl;
	e_read();
	cout<<"Creating and calculating crosssection. e_Gammaanarrowvm()"<<endl;
	e_narrowResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
e_Gammaanarrowvm::~e_Gammaanarrowvm()
{ }


//______________________________________________________________________________
e_Gammaawidevm::e_Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. e_Gammaawidevm()"<<endl;
	e_read();
	cout<<"Creating and calculating crosssection. e_Gammaawidevm()"<<endl;
	e_wideResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter();
}


//______________________________________________________________________________
e_Gammaawidevm::~e_Gammaawidevm()
{ }


//______________________________________________________________________________
void Gammaavectormeson::pickwEgamq2(double &W, double &cmsEgamma, double &targetEgamma, 
				 double &Q2, double &gamma_pz, double &gamma_pt,//photon in target frame
				 double &E_prime, double &theta_e //electron
				 )
{
        double dW, dEgamma;
	double xw,xEgamma, xQ2, xtest, q2test; // btest;
	int  IW,IGamma, IQ2;
	// ---------
	//	int egamma_draws = 0, cms_egamma_draws =0, q2_draws =0 ;
	// ---------
	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	// Following used for debugging/timing for event generation
	//std::chrono::steady_clock::time_point begin_evt = std::chrono::steady_clock::now();
	bool pick_state = false;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy);
	while( pick_state == false ){

	  xw = _randy.Rndom();
	  W = _VMWmin + xw*(_VMWmax-_VMWmin);
	  double w_test = _randy.Rndom();
	  if (W < 2 * starlightConstants::pionChargedMass)
	    continue;
	  IW = int((W-_VMWmin)/dW);
	  if (_BWarray[IW] < w_test)
	    continue;
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
	  N0++; 
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
	    //q2_draws++;
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
	  if( cmsEgamma < _cmsMinPhotonEnergy || 2.*targetEgamma/(Q2+W*W) < _targetRadius){ //This cut is roughly RA = 0.8 fm the radius of proton and 1 eV^{-1} = 1.97 x 10 ^{-7} m
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
eXEvent Gammaavectormeson::e_produceEvent()
{
	// The new event type
	eXEvent event;

	int iFbadevent=0;
	int tcheck=0;
	starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
        starlightConstants::particleTypeEnum vmpid = starlightConstants::UNKNOWN; 
	// at present 4 prong decay is not implemented
	double comenergy = 0.;
	double rapidity = 0.;
	double Q2 = 0;
	double E = 0.;
	double momx=0.,momy=0.,momz=0.;
	double targetEgamma = 0, cmsEgamma = 0 ;
	double gamma_pz = 0 , gamma_pt = 0, e_theta = 0;
	double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
	double px3=0.,py3=0.,pz3=0.;
	double pt1chk = 0., pt2chk = 0., pt3chk = 0.;
	double eta1 = 0., eta2 = 0., eta3 = 0.;
	double e_E=0., e_phi=0;
	double t_px =0, t_py=0., t_pz=0, t_E;
	bool accepted = false;
	bool threeBody = false;
	do{

	  pickwEgamq2(comenergy,cmsEgamma, targetEgamma, 
		   Q2, gamma_pz, gamma_pt, //photon infor in CMS frame
		   e_E, e_theta);	 //electron info in target frame  
	  //

	
	  momenta(comenergy,cmsEgamma, Q2, gamma_pz, gamma_pt, //input
		  rapidity, E, momx, momy, momz, //VM
		  t_px, t_py, t_pz, t_E, //pomeron
		  e_phi,tcheck); //electron
	  //
	  // inelasticity: used for angular distributions
	  double col_y = 1. - (e_E/_eEnergy)*std::pow(std::cos(e_theta/2.),2.);
	  double col_polarization = (1 - col_y)/(1-col_y+col_y*col_y/2.);
	  double spin_element = getSpinMatrixElement(comenergy, Q2, col_polarization);
	  _nmbAttempts++;
	  
	  vmpid = ipid; 

	  if (_VMpidtest == starlightConstants::OMEGA_pipipi){
	  	threePionDecay(ipid,momx,momy,momz,px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,iFbadevent,spin_element);

	  	double pt1chk = sqrt(px1*px1+py1*py1);
	  	double pt2chk = sqrt(px2*px2+py2*py2);
	  	double pt3chk = sqrt(px3*px3+py3*py3+pz3*pz3);
	  	double eta1 = pseudoRapidity(px1, py1, pz1);
	  	double eta2 = pseudoRapidity(px2, py2, pz2);
	  	double eta3 = pseudoRapidity(px3, py3, pz3);

	  	threeBody = true;

	  }

	  else{
	  	// Two body dedcay in eSTARlight includes the angular corrections due to finite virtuality
	  
	  twoBodyDecay(ipid,comenergy,momx,momy,momz,spin_element,
		       px1,py1,pz1,px2,py2,pz2,iFbadevent);
	  double pt1chk = sqrt(px1*px1+py1*py1);
	  double pt2chk = sqrt(px2*px2+py2*py2);
	  double eta1 = pseudoRapidity(px1, py1, pz1);
	  double eta2 = pseudoRapidity(px2, py2, pz2);
}


                        
if (!threeBody){
	  if(_ptCutEnabled && !_etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(!_ptCutEnabled && _etaCutEnabled){
	    if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(_ptCutEnabled && _etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
	      if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
		accepted = true;
		_nmbAccepted++;
	      }
	    }
	  }
	  else if(!_ptCutEnabled && !_etaCutEnabled)
	    _nmbAccepted++;
}

else{

if(_ptCutEnabled && !_etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax && pt3chk > _ptCutMin && pt3chk < _ptCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(!_ptCutEnabled && _etaCutEnabled){
	    if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax && eta3 > _etaCutMin && eta3 < _etaCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(_ptCutEnabled && _etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax && pt3chk > _ptCutMin && pt3chk < _ptCutMax){
	      if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax && eta3 > _etaCutMin && eta3 < _etaCutMax){
		accepted = true;
		_nmbAccepted++;
	      }
	    }
	  }
	  else if(!_ptCutEnabled && !_etaCutEnabled)
	    _nmbAccepted++;

}


	}while((_ptCutEnabled || _etaCutEnabled) && !accepted);
	

	if (iFbadevent==0&&tcheck==0) {
	  int q1=0,q2=0;
	  int ipid1,ipid2=0;
	  
	  double xtest = _randy.Rndom(); 
	  if (xtest<0.5)
	    {
	      q1=1;
	      q2=-1;
	    }
	  else {
	    q1=-1;
	    q2=1;
	  }
	  
	  
	  if ( ipid == 11 || ipid == 13 ){
	    ipid1 = -q1*ipid;
	    ipid2 = -q2*ipid;
	  } else if (ipid == 22){ // omega --> pi0 + gamma
	  	ipid1 = 22;
	  	ipid2 = 111;
	  } else {
	    ipid1 = q1*ipid;
	    ipid2 = q2*ipid;
	  }


	  // - Outgoing electron - target frame
	  double e_ptot = sqrt(e_E*e_E - starlightConstants::mel*starlightConstants::mel);
	  double e_px = e_ptot*sin(e_theta)*cos(e_phi);
	  double e_py = e_ptot*sin(e_theta)*sin(e_phi);
	  double e_pz = e_ptot*cos(e_theta);
	  lorentzVector electron(e_px, e_py, e_pz, e_E);
	  event.addSourceElectron(electron);
	  // - Generated photon - target frame
	  double gamma_x = gamma_pt*cos(e_phi+starlightConstants::pi);
	  double gamma_y = gamma_pt*sin(e_phi+starlightConstants::pi);
	  lorentzVector gamma(gamma_x,gamma_y,gamma_pz,cmsEgamma);
	  vector3 boostVector(0, 0, tanh(_rap_CM));
	  (gamma).Boost(boostVector);
	  event.addGamma(gamma, targetEgamma, Q2);   
	  // - Saving V.M. daughters
	  double md = 0; 

	  if (!threeBody){
	  double md = getDaughterMass(vmpid);
	  bool isOmegaNeutralDecay = (vmpid == starlightConstants::PHOTON);
	  if(isOmegaNeutralDecay){
	  	//double mpi0 = starlightConstants::pionNeutralMass;
	  	double e_gamma1,px_gamma1,py_gamma1,pz_gamma1,e_gamma2,px_gamma2,py_gamma2,pz_gamma2;
	    double Ed1 = sqrt(px1*px1+py1*py1+pz1*pz1); 
	  	starlightParticle particle1(px1, py1, pz1, Ed1, starlightConstants::UNKNOWN, ipid1, q1);
	  	event.addParticle(particle1);
	  	//

	  	pi0Decay(px2,py2,pz2,e_gamma1,px_gamma1,py_gamma1,pz_gamma1,e_gamma2,px_gamma2,py_gamma2,pz_gamma2,iFbadevent);
	  	starlightParticle particle2(px_gamma1, py_gamma1, pz_gamma1, e_gamma1, starlightConstants::UNKNOWN, ipid1, q2);
	  	starlightParticle particle3(px_gamma2, py_gamma2, pz_gamma2, e_gamma2, starlightConstants::UNKNOWN, ipid1, q2);
	  	event.addParticle(particle2);
	  	event.addParticle(particle3);
	  	//double invmass = sqrt(mpi0*mpi0 + 2.0*(Ed1*Ed2 - (px1*px2+py1*py2+pz1*pz2)));
	  } 


	 
	  else {
	  	double Ed1 = sqrt(md*md+px1*px1+py1*py1+pz1*pz1); 
	    starlightParticle particle1(px1, py1, pz1, Ed1, starlightConstants::UNKNOWN, ipid1, q1);
	    event.addParticle(particle1);
	    //
	    double Ed2 = sqrt(md*md+px2*px2+py2*py2+pz2*pz2); 
	    starlightParticle particle2(px2, py2, pz2, Ed2, starlightConstants::UNKNOWN, ipid2, q2);
	    event.addParticle(particle2);
	  }
	}

	if (threeBody){

	    double e_gamma1,px_gamma1,py_gamma1,pz_gamma1,e_gamma2,px_gamma2,py_gamma2,pz_gamma2;
	    double Ed1 = sqrt(pow(starlightConstants::pionChargedMass,2)+px1*px1+py1*py1+pz1*pz1); 
	    double Ed2 = sqrt(pow(starlightConstants::pionChargedMass,2)+px2*px2+py2*py2+pz2*pz2);
	    double Ed3 = sqrt(pow(starlightConstants::pionNeutralMass,2)+px3*px3+py3*py3+pz3*pz3);

	    double E0 = sqrt(pow(starlightConstants::OmegaMass,2)+momx*momx+momy*momy+momz*momz);


	  	starlightParticle particle1(px1, py1, pz1, Ed1, starlightConstants::UNKNOWN, 211, q1);
	  	event.addParticle(particle1);

	  	starlightParticle particle2(px2, py2, pz2, Ed2, starlightConstants::UNKNOWN, 211, q2);
	  	event.addParticle(particle2);



	  	//

	  	pi0Decay(px3,py3,pz3,e_gamma1,px_gamma1,py_gamma1,pz_gamma1,e_gamma2,px_gamma2,py_gamma2,pz_gamma2,iFbadevent);
	  	starlightParticle particle3(px_gamma1, py_gamma1, pz_gamma1, e_gamma1, starlightConstants::UNKNOWN, ipid1, 22);
	  	starlightParticle particle4(px_gamma2, py_gamma2, pz_gamma2, e_gamma2, starlightConstants::UNKNOWN, ipid1, 22);


	  	//cout<<"Looking at particles after boosted into beam frame"<<endl;
	  	// cout<<"Conservation of Omega->3Pi momentum check"<<endl;
	  	// cout<<"Net mom_x = "<<momx-(px1+px2+px3)<<endl;
	  	// cout<<"Net mom_y = "<<momy-(py1+py2+py3)<<endl;
	  	// cout<<"Net mom_z = "<<momz-(pz1+pz2+pz3)<<endl;


	  	// cout<<"Conservation of total energy check"<<endl;
	  	//cout<<"Net energy = " <<sqrt(pow(starlightConstants::OmegaMass,2)+momx*momx+momy*momy+momz*momz)-(Ed1+Ed2+sqrt(pow(starlightConstants::pionNeutralMass,2)+px3*px3+py3*py3+pz3*pz3))<<endl;

	  	// cout<<"########################"<<endl;


	 //  	cout<<"Now looking at particles after boosted into beam frame"<<endl;
		// cout<<"Conservation of total momentum check"<<endl;
		// cout<<"Net mom_x = "<<momx-(px1+px2+px_gamma1+px_gamma2)<<endl;
		// cout<<"Net mom_y = "<<momy-(py1+py2+py_gamma1+py_gamma2)<<endl;
		// cout<<"Net mom_z = "<<momz-(pz1+pz2+pz_gamma1+pz_gamma2)<<endl;


		// cout<<"Conservation of total energy check"<<endl;
		// cout<<"Energy conseveration from Omega->3Pi decay: " <<E0-(Ed1+Ed2+Ed3)<<endl;

		// cout<<"Energy conseveration from Pi0 decay:"<<Ed3-(e_gamma1+e_gamma2)<<endl;
		
		// cout<<"Total net energy = " <<E0-(Ed1+Ed2+e_gamma1+e_gamma2)<<endl;
		// cout<<"########################"<<endl;



	  	event.addParticle(particle3);
	  	event.addParticle(particle4);

	}

	  // - Scattered target and transfered momenta at target vertex
	  double target_pz =  - _beamNucleus*sqrt(_pEnergy*_pEnergy - pow(starlightConstants::protonMass,2.) ) - t_pz;
	  //Sign of t_px in following equation changed to fix sign error and conserve p_z.  Change made by Spencer Klein based on a bug report from Ya-Ping Xie.  Nov. 14, 2019
	  lorentzVector target(-t_px, -t_py, target_pz, _beamNucleus*_pEnergy - t_E);
	  double t_var = t_E*t_E - t_px*t_px - t_py*t_py - t_pz*t_pz;
	  event.addScatteredTarget(target, t_var);
	}
	return event;

}
string Gammaavectormeson::gammaTableParse(int ii, int jj)
{
  ostringstream tag1, tag2;
  tag1<<ii;
  tag2<<jj;
  string to_ret = tag1.str()+","+tag2.str();
  return to_ret;
}
