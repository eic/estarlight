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


#ifndef PI0_H
#define PI0_H


#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "eventchannel.h"
#include "eXevent.h"
#include "nBodyPhaseSpaceGen.h"
//Now here for eSTARlight
#include "photonNucleusCrossSection.h"

class Pi0 : public eventChannel
{
  
 public:
  Pi0(const inputParameters& ipnut, beamBeamSystem& bbsystem);
  virtual ~Pi0();
  
  eXEvent e_produceEvent();

  void pickEgamq2(double &cmsEgamma, double &targetEgamma, 
		double &Q2, double &gamma_pz, double & gamma_pt, 
		double &E_prime, double &cos_theta_e);
  void momenta(double Egam,double Q2, double gamma_pz, double gamma_pt,
	       double &rapidity, double &E,double &px,double &py,double &pz,
	       double &t_px, double &t_py, double &t_pz, double &t_E,
	       double &e_phi,int &tcheck);
  void twoPhotonDecay(double  px0, double  py0, double  pz0,
                           double& px1, double& py1, double& pz1,
                           double& px2, double& py2, double& pz2,
                           int&    iFbadevent);
  double _VMbslope;              
  double pseudoRapidity(double px, double py, double pz);
  double UfromCosTheta(double cosTheta, double W, double Q2);
  double cosThetaFromU(double u, double W, double Q2);
  double TfromCosTheta(double cosTheta, double W, double Q2);
  double cosThetaFromT(double t, double W, double Q2);

  
 private:
  int _VMnumega;
  int _VMnumQ2;
  double _targetMaxPhotonEnergy;
  double _targetMinPhotonEnergy;
  double _cmsMinPhotonEnergy;
  double _beamLorentzGamma;
  double _targetBeamLorentzGamma;
  double _rap_CM;
  double _targetRadius;
  double _pEnergy;
  int _beamNucleus;
  double _eEnergy;
  nBodyPhaseSpaceGen* _phaseSpaceGen;
  // eSTARlight
  photonNucleusCrossSection* _dummy_pncs;
  bool _backwardsProduction;
};

class e_Pi0 : public Pi0
{
 public:
  e_Pi0(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~e_Pi0();
};



#endif  // GAMMAPI0_H
