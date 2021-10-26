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
//
///////////////////////////////////////////////////////////////////////////


#ifndef GAMMAAVM_H
#define GAMMAAVM_H


#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "eventchannel.h"
#include "upcevent.h"
#include "eXevent.h"
#include "nBodyPhaseSpaceGen.h"
//Now here for eSTARlight
#include "photonNucleusCrossSection.h"

class Gammaavectormeson : public eventChannel
{
  
 public:
  Gammaavectormeson(const inputParameters& ipnut, beamBeamSystem& bbsystem);
  virtual ~Gammaavectormeson();
  
  eXEvent e_produceEvent();

  void pickwy(double &W, double &Y);
  void pickwEgamq2(double &W, double &cmsEgamma, double &targetEgamma, 
		double &Q2, double &gamma_pz, double & gamma_pt, 
		double &E_prime, double &cos_theta_e);
  void momenta(double W,double Y,double &E,double &px,double &py,double &pz,int &tcheck);
  void momenta(double W,double Egam,double Q2, double gamma_pz, double gamma_pt,
	       double &rapidity, double &E,double &px,double &py,double &pz,
	       double &t_px, double &t_py, double &t_pz, double &t_E,
	       double &e_phi,int &tcheck);
  double pTgamma(double E); 
  void vmpt(double W,double Y,double &E,double &px,double &py, double &pz,int &tcheck);
  void twoBodyDecay(starlightConstants::particleTypeEnum &ipid,double W,double px0,double py0,double pz0,double spin_element,
		    double &px1,double &py1,double&pz1,double &px2,double &py2,double &pz2,int &iFbadevent);
  bool fourBodyDecay(starlightConstants::particleTypeEnum& ipid, const double E, const double W, const double* p, lorentzVector* decayMoms, int& iFbadevent);
  void pi0Decay(double& px_pi0, double& py_pi0, double& pz_pi0,double& e_g1, double& px_g1, double& py_g1, double& pz_g1,double& e_g2, double& px_g2, double& py_g2, double& pz_g2,int& iFbadevent);
  double getMass();
  double getWidth();
  virtual double getTheta(starlightConstants::particleTypeEnum ipid, double r_04_00);
  double getSpinMatrixElement(double W, double Q2, double epsilon);
  double getSpin();
  double _VMbslope;
  virtual double getDaughterMass(starlightConstants::particleTypeEnum &ipid);                
  double pseudoRapidity(double px, double py, double pz);
  
 private:
  std::string gammaTableParse(int ii, int jj);
  starlightConstants::particleTypeEnum _VMpidtest;
  int _VMnumw;
  int _VMnumy;
  int _VMnumega;
  int _VMnumQ2;
  int _VMinterferencemode;
  int _ProductionMode;
  double _targetMaxPhotonEnergy;
  double _targetMinPhotonEnergy;
  double _cmsMaxPhotonEnergy;
  double _cmsMinPhotonEnergy;
  double _beamLorentzGamma;
  double _beam2LorentzGamma;
  double _rap_CM;
  double _targetRadius;
  int _TargetBeam; 
  int N0;
  int N1;
  int N2; 
  double _VMgamma_em;
  double _VMNPT;
  double _VMWmax;
  double _VMWmin;
  double _VMYmax;
  double _VMYmin;
  double _mass;
  double _width;
  double _VMptmax;
  double _VMdpt;
  int    _bslopeDef;
  double _bslopeVal;
  double _pEnergy;
  int _beamNucleus;
  double _eEnergy;
  nBodyPhaseSpaceGen* _phaseSpaceGen;
  // eSTARlight
  photonNucleusCrossSection* _dummy_pncs;
  double _angular_max[100][200];
};

class Gammaanarrowvm : public Gammaavectormeson
{
 public:
  Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~Gammaanarrowvm();
};

class Gammaawidevm : public Gammaavectormeson
{  
 public:
  Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~Gammaawidevm();
};

class Gammaaincoherentvm : public Gammaavectormeson
{  
 public:
  Gammaaincoherentvm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~Gammaaincoherentvm();
};
class e_Gammaanarrowvm : public Gammaavectormeson
{
 public:
  e_Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~e_Gammaanarrowvm();
};

class e_Gammaawidevm : public Gammaavectormeson
{  
 public:
  e_Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~e_Gammaawidevm();
};

#endif  // GAMMAAVM_H
