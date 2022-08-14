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


#ifndef DVCS_H
#define DVCS_H


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

class Dvcs : public eventChannel
{
  
 public:
  Dvcs(const inputParameters& ipnut, beamBeamSystem& bbsystem);
  virtual ~Dvcs();
  
  eXEvent e_produceEvent();

  void pickEgamq2(double &cmsEgamma, double &targetEgamma, 
		double &Q2, double &gamma_pz, double & gamma_pt, 
		double &E_prime, double &cos_theta_e);
  void momenta(double Egam,double Q2, double gamma_pz, double gamma_pt,
	       double &rapidity, double &E,double &px,double &py,double &pz,
	       double &t_px, double &t_py, double &t_pz, double &t_E,
	       double &e_phi,int &tcheck);
  double _VMbslope;              
  double pseudoRapidity(double px, double py, double pz);
  
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

class e_Dvcs : public Dvcs
{
 public:
  e_Dvcs(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~e_Dvcs();
};



#endif  // GAMMAAVM_H
