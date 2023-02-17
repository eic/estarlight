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
// $Rev:: 259                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-04-19 01:58:25 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef PHOTONNUCLEUSCROSSSECTION_H
#define PHOTONNUCLEUSCROSSSECTION_H


#include "starlightconstants.h"
#include "beambeamsystem.h"
#include "inputParameters.h"

class photonNucleusCrossSection{

public:

	photonNucleusCrossSection(const inputParameters& input, const beamBeamSystem&  bbsystem);
	~photonNucleusCrossSection();
  
	double         slopeParameter  () const { return _slopeParameter;   }  ///< returns slope of t-distribution [(GeV/c)^{-2}]
	double         getChannelMass  () const { return _channelMass;      }  ///< returns mass of the produced system [GeV/c^2]
	double         getBNORM        () const { return _BNORM;            }
	beamBeamSystem getbbs          () const { return _bbs;              }  ///< returns beamBeamSystem
	double         vmPhotonCoupling() const { return _vmPhotonCoupling; }  ///< vectormeson-photon coupling constant f_v / 4 pi (cf. Eq. 10 in KN PRC 60 (1999) 014903)
	double         vmQ2Power       (double Q2) const { return _vmQ2Power_c1+ _vmQ2Power_c2*(_channelMass*_channelMass + Q2);        }
	double         getDefaultC     () const { return _defaultC;         }
	double         maxPhotonEnergy () const { return _maxPhotonEnergy;  }  ///< returns max photon energy in lab frame [GeV] (for vectormesons only)

	void crossSectionCalculation(const double bwnormsave);
	double backwardsPropagationOmegaCrossSection(const double targetEgamma);

	// Use the wide or narrow constructor to calculate sigma
	// wide/narrow will inherit this.
	double getcsgA(const double Egamma,
	               const double Q2,
                       const int beam);
	double e_getcsgA(const double Egamma, double Q2,
	               const double W,
                       const int beam);
	// Midification to csg due to virtuality
	double getcsgA_Q2_dep(const double Q2);
	double photonFlux(const double Egamma,
                       const int beam);
	// --- Added for electron 
	double photonFlux(const double Egamma,
			  const double Q2);
	double integrated_Q2_dep(const double Egamma, const double _min = 0 , const double _max = 0);
	double integrated_x_section(const double Egamma, const double _min = 0 , const double _max = 0);
	std::pair<double,double>* Q2arraylimits(double const Egamma);
	double g(double const Egamma, double const Q2);
	// ---
	double sigmagp(const double Wgp);
	double sigma_A(const double sig_N, 
                       const int beam);
        double sigma_N(const double Wgp);
	double breitWigner(const double W,
	                   const double C);
	double nepoint(const double Egamma,
	               const double bmin);

	double getPhotonNucleusSigma () const {return _photonNucleusSigma;}
	void   setPhotonNucleusSigma (double sigma) {_photonNucleusSigma = sigma;}
	
protected:
	const unsigned int _nWbins;
	const unsigned int _nYbins;
	
	const double _wMin;
	const double _wMax;
	const double _yMax;

	const double _beamLorentzGamma;

	double _photonNucleusSigma; 

	int    _printDef; 
        int    _impulseSelected;
	int    _quantumGlauber;  // from input parameter; 1 for Quantum Glauber, 0 for classical Glauber
	
private:

	beamBeamSystem _bbs;
  
	// copied from inputParameters
	double                               _protonEnergy;
	double                               _electronEnergy;
	starlightConstants::particleTypeEnum _particleType;
	int                                  _beamBreakupMode;     ///< breakup mode for beam particles
	bool _backwardsProduction;
        int                                  _productionMode; 
	int                                  _sigmaNucleus; 

	// locally defined parameters
	double _slopeParameter;    ///< slope of t-distribution [(GeV/c)^{-2}]
	double _vmPhotonCoupling;  ///< vectormeson-photon coupling constant f_v / 4 pi (cf. Eq. 10 in KN PRC 60 (1999) 014903)
	double _vmQ2Power_c1;  ///< VM power law Q2 correction (Mv2/(Q2+Mv2))^n; with n = c1 + c2*(Q2+Mv2)
	double _vmQ2Power_c2;  ///< VM power law Q2 correction (Mv2/(Q2+Mv2))^
	double _ANORM;
	double _BNORM;
	double _defaultC;
	double _width;            ///< width of the produced system  [GeV/c^2]
	double _channelMass;      ///< mass of the produced system  [GeV/c^2]
	double _fixedQ2range;     ///< loads state of Q2 range
	double _minQ2;            ///< min Q2 range in case it is fixed
	double _maxQ2;            ///< max Q3 range in case it is fixed
	double _maxPhotonEnergy;  ///< max photon energy in lab frame [GeV] (for vectormesons only)
	double _cmsMinPhotonEnergy;
	double _targetRadii;
	double _maxW_GP;		  ///< max W_GP energy
	double _minW_GP; 		  ///< min W_GP energy
	
};


#endif  // PHOTONNUCLEUSCROSSSECTION_H
