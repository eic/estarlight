/////////////////////////////////////////////////////////////////////////// 
// 
// Copyright 2010 
// 
// This file is part of starlight. 
// 
// starlight is free software: you can redistribute it and/or modify 
// it under the terms of the GNU General Public License as published by 
// the Free Software Foundation, either version 3 of the License, or 
// (at your option) any later version. 
// 
// starlight is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
// GNU General Public License for more details. 
// 
// You should have received a copy of the GNU General Public License 
// along with starlight. If not, see <http://www.gnu.org/licenses/>. 
// 
/////////////////////////////////////////////////////////////////////////// 
// 
// File and Version Information: 
// $Rev:: 276 $: revision of last commit // $Author:: jnystra#$: author of last commit 
// $Date:: 2016-09-13 19:54:42 +0100 #$: date of last commit // 
// Description: 
// 
// 
// 
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "reportingUtils.h"
#include "starlightconstants.h"
#include "inputParameters.h"
#include "inputParser.h"
#include "starlightconfig.h"
#include <cmath>
#include <cstring>
#include "randomgenerator.h"

using namespace std;
using namespace starlightConstants;

parameterlist parameterbase::_parameters;

#define REQUIRED true
#define NOT_REQUIRED false

//______________________________________________________________________________
inputParameters::inputParameters()
        : _baseFileName          ("baseFileName","slight"),
	  _targetBeamZ                ("TARGET_BEAM_Z",0),
	  _targetBeamA                ("TARGET_BEAM_A",0),
	  _electronBeamLorentzGamma     ("ELECTRON_BEAM_GAMMA",0),
	  _targetBeamLorentzGamma     ("TARGET_BEAM_GAMMA",0),
	  _maxW                  ("W_MAX",0),
	  _minW                  ("W_MIN",0),
	  _nmbWBins              ("W_N_BINS",0),
	  _maxRapidity           ("RAP_MAX",9,NOT_REQUIRED),
	  _nmbRapidityBins       ("RAP_N_BINS",100,NOT_REQUIRED),
	  _nmbEnergyBins         ("EGA_N_BINS",0),
	  _ptCutEnabled          ("CUT_PT",false, NOT_REQUIRED),
	  _ptCutMin              ("PT_MIN",0, NOT_REQUIRED),
	  _ptCutMax              ("PT_MAX",0, NOT_REQUIRED),
	  _etaCutEnabled         ("CUT_ETA",false, NOT_REQUIRED),
	  _etaCutMin             ("ETA_MIN",0, NOT_REQUIRED),
	  _etaCutMax             ("ETA_MAX",0, NOT_REQUIRED),
	  _productionMode        ("PROD_MODE",0),
	  _nmbEventsTot          ("N_EVENTS",0),
	  _prodParticleId        ("PROD_PID",0),
	  _randomSeed            ("RND_SEED",0),
	  _beamBreakupMode       ("BREAKUP_MODE",0),
	  _interferenceEnabled   ("INTERFERENCE",false,NOT_REQUIRED),
	  _interferenceStrength  ("IF_STRENGTH",0,NOT_REQUIRED),
	  _maxPtInterference     ("INT_PT_MAX",0,NOT_REQUIRED),
	  _nmbPtBinsInterference ("INT_PT_N_BINS",0),
	  _ptBinWidthInterference("INT_PT_WIDTH",0),
	  _protonEnergy          ("PROTON_ENERGY",0, NOT_REQUIRED),
	  _electronEnergy        ("ELECTRON_ENERGY",0, NOT_REQUIRED),
	  _minGammaEnergy	 ("MIN_GAMMA_ENERGY",6.0, NOT_REQUIRED),
	  _maxGammaEnergy	 ("MAX_GAMMA_ENERGY",600000.0, NOT_REQUIRED),
	  _minGammaQ2            ("MIN_GAMMA_Q2",0,NOT_REQUIRED),
	  _maxGammaQ2            ("MAX_GAMMA_Q2",0,NOT_REQUIRED),
	  _nmbGammaQ2Bins       ("INT_GAMMA_Q2_BINS",400,NOT_REQUIRED),
	  _pythiaParams          ("PYTHIA_PARAMS","", NOT_REQUIRED),
      _outputFormat           ("OUTPUT_FORMAT",0, NOT_REQUIRED),
	  _backwardsProduction    ("BACKWARDS_PRODUCTION",false, NOT_REQUIRED),

	  _xsecCalcMethod	 ("XSEC_METHOD",0, NOT_REQUIRED),
          _axionMass             ("AXION_MASS",50, NOT_REQUIRED),  // AXION HACK
          _bslopeDefinition      ("BSLOPE_DEFINITION",0, NOT_REQUIRED),
	  _bslopeValue           ("BSLOPE_VALUE",4.0,NOT_REQUIRED),
	  _impulseVM             ("SELECT_IMPULSE_VM",0,NOT_REQUIRED),
	  _quantumGlauber        ("QUANTUM_GLAUBER",0,NOT_REQUIRED)
{
  // All parameters must be initialised in initialisation list! 
  // If not: error: 'parameter<T, validate>::parameter() [with T = unsigned int, bool validate = true]' is private
  // or similar
	
        _ip.addParameter(_baseFileName);

	_ip.addParameter(_targetBeamZ);
	_ip.addParameter(_targetBeamA);

	_ip.addParameter(_targetBeamLorentzGamma);
	_ip.addParameter(_electronBeamLorentzGamma);
	
	_ip.addParameter(_maxW);
	_ip.addParameter(_minW);
	
	_ip.addParameter(_nmbWBins);

	_ip.addParameter(_maxRapidity);
	_ip.addParameter(_nmbRapidityBins);
	//
	_ip.addParameter(_nmbEnergyBins);

	_ip.addParameter(_ptCutEnabled);
	_ip.addParameter(_ptCutMin);
	_ip.addParameter(_ptCutMax);
	
	_ip.addParameter(_etaCutEnabled);
	_ip.addParameter(_etaCutMax);
	_ip.addParameter(_etaCutMin);
	
	_ip.addParameter(_productionMode);
	_ip.addParameter(_nmbEventsTot);
	_ip.addParameter(_prodParticleId);
	_ip.addParameter(_randomSeed);
	_ip.addParameter(_beamBreakupMode);
	_ip.addParameter(_interferenceEnabled);
	_ip.addParameter(_interferenceStrength);
	_ip.addParameter(_maxPtInterference);
	_ip.addParameter(_nmbPtBinsInterference);
	_ip.addParameter(_minGammaEnergy);
	_ip.addParameter(_maxGammaEnergy);
	_ip.addParameter(_minGammaQ2);
	_ip.addParameter(_maxGammaQ2);
	_ip.addParameter(_nmbGammaQ2Bins);
	_ip.addParameter(_pythiaParams);
	_ip.addParameter(_outputFormat);
	_ip.addParameter(_backwardsProduction);
	_ip.addParameter(_xsecCalcMethod);
        _ip.addParameter(_axionMass);     // AXION HACK
        _ip.addParameter(_bslopeDefinition); 
        _ip.addParameter(_bslopeValue); 
	_ip.addParameter(_impulseVM);
	_ip.addParameter(_quantumGlauber);
}


//______________________________________________________________________________
inputParameters::~inputParameters()
{ }


//______________________________________________________________________________
bool
inputParameters::configureFromFile(const std::string &_configFileName)
{

	int nParameters = _ip.parseFile(_configFileName);
	
	if(nParameters == -1) 
	{
		printWarn << "could not open file '" << _configFileName << "'" << endl;
	  return false;
	}
	
	
	if(_ip.validateParameters(cerr))
		printInfo << "successfully read input parameters from '" << _configFileName << "'" << endl;
	else {
 		printWarn << "problems reading input parameters from '" << _configFileName << "'" << endl
 		          << *this;
 		return false;
 	}
 	return true;
}
 bool inputParameters::init()
 {
	
 	// Calculate beam gamma in CMS frame
 	double rap1 = acosh(electronBeamLorentzGamma());
	double rap2 = -acosh(targetBeamLorentzGamma());

	double _electronEnergy_lab = (electronBeamLorentzGamma())*(starlightConstants::mel);
	double _ionEnergy_lab = targetBeamLorentzGamma()*(targetBeamA())*protonMass;
	double _ionPz_lab = -1.0*sqrt(_ionEnergy_lab*_ionEnergy_lab - (targetBeamA()*targetBeamA())*protonMass*protonMass);
	double _electronPz_lab = 1.0*sqrt(_electronEnergy_lab*_electronEnergy_lab - starlightConstants::mel*(starlightConstants::mel));
	double _totalEnergy_lab = 1.0*_electronEnergy_lab + _ionEnergy_lab;
	double _totalPz_lab = _ionPz_lab + _electronPz_lab;
	_rap_CM = (1.0/2.0)*log((_totalEnergy_lab + _totalPz_lab)/(_totalEnergy_lab - _totalPz_lab));

	_beamLorentzGamma = cosh(_rap_CM-rap2);
	_targetLorentzGamma = cosh(rap1-rap2);

	if( targetBeamA() == 1) //proton case 0.87 fm = 4.4 GeV^{-1}
	  _targetR = 4.4;
	else
	  _targetR = 6.1 * pow(targetBeamA(), 1./3.);
	  //_targetR = 4.4/2.;

	_fixedQ2Range = false;
	std::cout << "Rapidity electron beam: " << rap1 << ", rapidity target beam: " << rap2 << ", rapidity CMS system: " << _rap_CM << ", beam gamma in CMS: " << _beamLorentzGamma<< std::endl;
	std::cout << "Rapidity beam 1 in beam 2 frame: " << rap1-rap2 << ", beam 1 gamma in beam 2 frame: " << _targetLorentzGamma<< std::endl;
	_ptBinWidthInterference = maxPtInterference() / nmbPtBinsInterference();
	_protonEnergy           = _beamLorentzGamma * protonMass;
	// Electron energy in the target frame of reference
	_electronEnergy         = _targetLorentzGamma * starlightConstants::mel;

	if( (targetBeamZ()==1) && (targetBeamA()==3) ){
		printWarn << "tritium is not currently supported" << endl;
		return false;}
	// check that rho production uses wide resonance option
	if(_prodParticleId.value()==113 && _productionMode.value()==2){
		printWarn << endl<< "For rho meson production, you should choose the wide resonance option (production mode = 3)" << endl;
		return false;}

	// define interaction type
	switch (productionMode()) {
	case 1:
		_interactionType = PHOTONPHOTON;
		break;
	case 2:
		_interactionType = PHOTONPOMERONNARROW;
		break;
	case 3:
		_interactionType = PHOTONPOMERONWIDE;
		break;
        case 4:
                _interactionType = PHOTONPOMERONINCOHERENT;
                break;
	case 5:
		_interactionType = PHOTONUCLEARSINGLE;
		break;
	case 6:
		_interactionType = PHOTONUCLEARDOUBLE;
		break;
	case 7:
		_interactionType = PHOTONUCLEARSINGLEPA;
		break;
	case 8:
		_interactionType = PHOTONUCLEARSINGLEPAPY;
		break;
// 	case 9:
// 		_interactionType = PHOTONPHOTONKNIEHL;
// 		break;
//         case 10:
//                 _interactionType = PHOTONPHOTONKNIEHLMODIFIED;
//                 break;
	case 12:
		_interactionType = E_PHOTONPOMERONNARROW;
		break;
	case 13:
		_interactionType = E_PHOTONPOMERONWIDE;
		break;
	  
	default:
		printWarn << "unknown production mode '" << _productionMode << "'" << endl;
		return false;
	}
	if( (_productionMode.value() ==4) && (_interferenceEnabled.value())) {
		printWarn << " cannot enable interference for incoherent production " << endl;
		return false;
	}
  
	double mass        = 0;
	double width       = 0;
	double defaultMinW = 0;  // default for _minW, unless it is defined later [GeV/c^2]
	double defaultMaxW = 0;  // default for _maxW, unless it is defined later [GeV/c^2]
	switch (prodParticleId()) {

// 	case 24:  // W+W- pair
// 		_particleType = W;
// 		_decayType    = WW;
// 		defaultMinW   = 2 * muonMass;
// 		break;
	case 115:  // a_2(1320)
		_particleType = A2;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::a2Mass;
		width         = starlightConstants::a2Width;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW   = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0; 
		break;
	case 221:  // eta
		_particleType = ETA;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::etaMass;
		width         = starlightConstants::etaWidth;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW   = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0; 
		break;
	case 225:  // f_2(1270)
		_particleType = F2;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::f2Mass;
		width         = starlightConstants::f2Width;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = starlightConstants::f2BrPiPi; 
		break;
	case 331:  // eta'(958)
		_particleType = ETAPRIME;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::etaPrimeMass;
		width         = starlightConstants::etaPrimeWidth;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0; 
		break;
	case 335:  // f_2'(1525)
		_particleType = F2PRIME;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::f2PrimeMass;
		width         = starlightConstants::f2PrimeWidth;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = starlightConstants::f2PrimeBrKK; 
		break;
	case 441:  // eta_c(1s)
		_particleType = ETAC;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::etaCMass;
		width         = starlightConstants::etaCWidth;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0; 
		break;
	case 9010221:  // f_0(980), was orginally called 10221? updated to standard number
		_particleType = F0;
		_decayType    = SINGLEMESON;
		mass          = starlightConstants::f0Mass;
		width         = starlightConstants::f0Width;
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
               _inputBranchingRatio = starlightConstants::f0BrPiPi; 
		break;
	case 33:  // Z"/Z0  This is the rho^0 rho^0 final state SRK
		_particleType = ZOVERZ03;
		_decayType    = SINGLEMESON;
                defaultMinW   = 4*pionChargedMass;
                defaultMaxW         = 1.6; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0; 
		break;
        case 88:  // axion// AXION HACK, till break statement
                _particleType = AXION;
                _decayType    = SINGLEMESON;
                mass          = _axionMass.value();
                width         = 1/(64*starlightConstants::pi)*mass*mass*mass/(1000*1000);//Fix Lambda=1000 GeV, rescaling is trivial.
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                break;    // AXION HACK, end
// 	case 25: // Higgs
// 		_particleType = HIGGS;
// 		_decayType    = SINGLEMESON;
// 		break;
	case 113:  // rho(770)
		_particleType = RHO;
		_decayType    = WIDEVMDEFAULT;
		mass          = starlightConstants::rho0Mass;
		width         = starlightConstants::rho0Width;
		defaultMinW   = 2 * pionChargedMass;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::rho0BrPiPi; 
		break;
	case 913:  // rho(770) with direct pi+pi- decay, interference given by ZEUS data
		_particleType = RHOZEUS;
		_decayType    = WIDEVMDEFAULT;
		mass          = starlightConstants::rho0Mass;
		width         = starlightConstants::rho0Width;
		defaultMinW   = 2 * pionChargedMass;
		defaultMaxW         = mass + 5 * width;  // use the same 1.5GeV max mass as ZEUS
		_inputBranchingRatio = starlightConstants::rho0BrPiPi; 
		break;
	case 999:  // pi+pi-pi+pi- phase space decay
		_particleType = FOURPRONG;
		_decayType    = WIDEVMDEFAULT;
		mass          = starlightConstants::rho0PrimeMass;
		width         = starlightConstants::rho0PrimeWidth;		
		defaultMinW   = 4 * pionChargedMass;
        defaultMaxW   = 10.0;	    
		//defaultMaxW     = sqrt(electronBeamLorentzGamma()*targetBeamLorentzGamma())*2*(starlightConstants::hbarc)/(pow(float(targetBeamA()),1./6.)); // JES 6.17.2015 to avoid problems with no default
    	_inputBranchingRatio = 1.0; 
		break;

	case 223:  // omega(782)
		_particleType = OMEGA;
		_decayType    = NARROWVMDEFAULT;  
		mass          = starlightConstants::OmegaMass;
		width         = starlightConstants::OmegaWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::OmegaBrPiPi; 
		break;
	case 223022:  // omega(782) -> pi0 + gamma
		_particleType = OMEGA_pi0gamma;
		_decayType    = NARROWVMDEFAULT;  
		mass          = starlightConstants::OmegaMass;
		width         = starlightConstants::OmegaWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::OmegaBrPi0Gamma; 
		break;
	case 333:  // phi(1020)
		_particleType = PHI;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::PhiMass;
		width         = starlightConstants::PhiWidth;
		defaultMinW   = 2 * kaonChargedMass;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::PhiBrKK; 
		break;
	case 443:  // J/psi
		_particleType = JPSI;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::JpsiMass;
		width         = starlightConstants::JpsiWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = (starlightConstants::JpsiBree + starlightConstants::JpsiBrmumu)/2.; 
		break;
   	case 443011:  // J/psi
		_particleType = JPSI_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::JpsiMass;
		width         = starlightConstants::JpsiWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::JpsiBree; 
		break;
	case 443013:  // J/psi
	        cout<<"In inputParameters setting J/psi mass!"<<endl; 
		_particleType = JPSI_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::JpsiMass;
		width         = starlightConstants::JpsiWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::JpsiBrmumu; 
		break;
	case 444:  // psi(2S) 
		_particleType = JPSI2S;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Psi2SMass;
		width         = starlightConstants::Psi2SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = (starlightConstants::Psi2SBree + starlightConstants::Psi2SBrmumu)/2.; 
		break;
	case 444011: // psi(2S)
		_particleType = JPSI2S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Psi2SMass;
		width         = starlightConstants::Psi2SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Psi2SBree; 
		break;
	case 444013:  // psi(2S)
		_particleType = JPSI2S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Psi2SMass;
		width         = starlightConstants::Psi2SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Psi2SBrmumu; 
		break;
	case 553:  // Upsilon(1S) 
		_particleType = UPSILON;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon1SMass;
		width         = starlightConstants::Upsilon1SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = (starlightConstants::Upsilon1SBree + starlightConstants::Upsilon1SBrmumu)/2.; 
		break;
	case 553011:  // Upsilon
		_particleType = UPSILON_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon1SMass;
		width         = starlightConstants::Upsilon1SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Upsilon1SBree; 
		break;
	case 553013:  // Upsilon
		_particleType = UPSILON_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon1SMass;
		width         = starlightConstants::Upsilon1SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Upsilon1SBrmumu; 
		break;
	case 554:  // Upsilon(2S)
		_particleType = UPSILON2S;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon2SMass;
		width         = starlightConstants::Upsilon2SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = (starlightConstants::Upsilon2SBree + starlightConstants::Upsilon2SBrmumu)/2.; 
		break;
	case 554011:  // Upsilon(2S)
		_particleType = UPSILON2S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon2SMass;
		width         = starlightConstants::Upsilon2SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Upsilon2SBree; 
		break;
	case 554013:  // Upsilon(2S)
		_particleType = UPSILON2S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon2SMass;
		width         = starlightConstants::Upsilon2SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Upsilon2SBrmumu; 
		break;
	case 555:  // Upsilon(3S)
		_particleType = UPSILON3S;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon3SMass;
		width         = starlightConstants::Upsilon3SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = (starlightConstants::Upsilon3SBree + starlightConstants::Upsilon3SBrmumu)/2.; 
		break;
	case 555011:  // Upsilon(3S)
		_particleType = UPSILON3S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon3SMass;
		width         = starlightConstants::Upsilon3SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Upsilon3SBree; 
		break;
	case 555013:  // Upsilon(3S)
		_particleType = UPSILON3S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = starlightConstants::Upsilon3SMass;
		width         = starlightConstants::Upsilon3SWidth;
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = starlightConstants::Upsilon3SBrmumu; 
		break;
	default:
		printWarn << "unknown particle ID " << _prodParticleId << endl;
		return false;
	}  // _prodParticleId

	if (_minW.value() == -1)
		_minW = defaultMinW;
	if (_maxW.value() == -1)
		_maxW = defaultMaxW;
	if ( _maxW.value() <= _minW.value() ) {
		printWarn << "maxW must be greater than minW" << endl;
		return false;
	}

	// Sanity check on Q2 range in case it is specified by user
	if( _minGammaQ2.value() != 0 || _maxGammaQ2.value() != 0){
	  if( _minGammaQ2.value() <0 || _maxGammaQ2.value() <=_minGammaQ2.value() )
	    printWarn << "Input values for min and max Q2 are inconsistent: "<<_minGammaQ2.value()<<","<<_maxGammaQ2.value()
		      <<". Continuing with default "<<endl;
	  else
	    _fixedQ2Range = true;
	}
	//Photon energy limits in C.M.S frame, used for some basic safety chekcs
	_cmsMinPhotonEnergy  = 0.5*(((mass+protonMass)*(mass+protonMass)-
				     protonMass*protonMass)/(_protonEnergy.value()+sqrt(_protonEnergy.value()*_protonEnergy.value()-protonMass*protonMass)));
	_cmsMaxPhotonEnergy        = 0.5*mass*exp(9);
	// Photon limits in target frame: this is where the photon flux is well described
	_targetMaxPhotonEnergy        = (_targetLorentzGamma - 10. ) *starlightConstants::mel;
	_targetMinPhotonEnergy = mass;
	printInfo << "using the following " << *this;
	
	return true;
}


//______________________________________________________________________________
ostream&
inputParameters::print(ostream& out) const
{
	out << "starlight parameters:" << endl
	    << "    base file name  ...................... '"  << _baseFileName.value() << "'" << endl
	    << "    target beam atomic number .............. " << _targetBeamZ.value() << endl
	    << "    target beam atomic mass number ......... " << _targetBeamA.value() << endl
	    << "    Lorentz gamma of beams in CM frame ..... " << _beamLorentzGamma << endl
	    << "    mass W of produced hadronic system ..... " << _minW.value() << " < W < " << _maxW.value() << " GeV/c^2" << endl
	    << "    # of W bins ............................ " << _nmbWBins.value() << endl
	    << "    maximum absolute value for rapidity .... " << _maxRapidity.value() << endl
	    << "    # of rapidity bins ..................... " << _nmbRapidityBins.value() << endl
            << "    # of Egamma bins ....................... " << _nmbEnergyBins.value() << endl
	    << "    cut in pT............................... " << yesNo(_ptCutEnabled.value()) << endl;
    if (_ptCutEnabled.value()) {
	out << "        minumum pT.......................... " << _ptCutMin.value() << " GeV/c" << endl
	    << "        maximum pT.......................... " << _ptCutMax.value() << " GeV/c" << endl;}
	out << "    cut in eta.............................. " << yesNo(_etaCutEnabled.value()) << endl;
    if (_etaCutEnabled.value()) {
	out << "        minumum eta......................... " << _etaCutMin.value() << endl
	    << "        maximum eta......................... " << _etaCutMax.value() << endl;}
        out << "    production mode ........................ " << _productionMode.value() << endl
	    << "    number of events to generate ........... " << _nmbEventsTot.value() << endl
	    << "    PDG ID of produced particle ............ " << _prodParticleId.value() << endl
	    << "    seed for random generator .............. " << _randomSeed.value() << endl
	    << "    breakup mode for beam particles ........ " << _beamBreakupMode.value() << endl
	    << "    interference enabled ................... " << yesNo(_interferenceEnabled.value()) << endl;
    if (_interferenceEnabled.value()) {
	out << "    interference strength .................. " << _interferenceStrength.value() << endl
	    << "    maximum p_T for interference calc. ..... " << _maxPtInterference.value() << " GeV/c" << endl
	    << "    # of p_T bins for interference calc. ... " << _nmbPtBinsInterference.value() << endl;}
    if (_productionMode.value()!=1) {
		if (_productionMode.value()==4) {
	 	  out  << "    coherent scattering off nucleus ........ no" << endl;}
		else {
	 	  out  << "    coherent scattering off nucleus ........ yes" << endl;
		}
    }
    if (_fixedQ2Range == true ){
      out << "    fixed photon Q2 range .................. " << _minGammaQ2.value() << " < Q2 < "
	  <<_maxGammaQ2.value() << " GeV/c^2 "<<endl;
    }
    out<< " Q2_BINS"       <<nmbGammaQ2Bins        () <<endl;
    out     <<"    Quantum Glauber parameter...............  "<<_quantumGlauber.value()<<endl;
    out     <<"    Impulse VM parameter....................  "<<_impulseVM.value()<<endl;
    return out;
}


//______________________________________________________________________________
ostream&
inputParameters::write(ostream& out) const
{
  
        out << "baseFileName"  << baseFileName         () <<endl
	    << "TARGET_BEAM_Z" << targetBeamZ          () <<endl
	    << "TARGET_BEAM_A" << targetBeamA          () <<endl
	    << "BEAM_GAMMA"    << beamLorentzGamma     () <<endl
	    << "W_MAX"         << maxW                 () <<endl
	    << "W_MIN"         << minW                 () <<endl
	    << "W_N_BINS"      << nmbWBins             () <<endl
	    << "RAP_MAX"       << maxRapidity          () <<endl
	    << "RAP_N_BINS"    << nmbRapidityBins      () <<endl
	    << "CUT_PT"        << ptCutEnabled         () <<endl
	    << "PT_MIN"        << ptCutMin             () <<endl
	    << "PT_MAX"        << ptCutMax             () <<endl
	    << "CUT_ETA"       << etaCutEnabled        () <<endl
	    << "ETA_MIN"       << etaCutMin            () <<endl
	    << "ETA_MAX"       << etaCutMax            () <<endl
	    << "PROD_MODE"     << productionMode       () <<endl
	    << "N_EVENTS"      << nmbEvents            () <<endl
	    << "PROD_PID"      << prodParticleId       () <<endl
	    << "RND_SEED"      << randomSeed           () <<endl
	    << "BREAKUP_MODE"  << beamBreakupMode      () <<endl
	    << "INTERFERENCE"  << interferenceEnabled  () <<endl
	    << "IF_STRENGTH"   << interferenceStrength () <<endl
	    << "INT_PT_MAX"    << maxPtInterference    () <<endl
	    << "INT_PT_N_BINS" << nmbPtBinsInterference() <<endl;
	out<< "USING FIXED RANGE"<<_fixedQ2Range<<endl;
	if( _fixedQ2Range == true ){
	  out << "Q2_MIN"      << minGammaQ2           () <<endl
	      << "Q2_MAX"      << maxGammaQ2           () <<endl;
	}
	out<< " Q2_BINS"       <<nmbGammaQ2Bins        () <<endl;
	return out;
}

std::string
inputParameters::parameterValueKey() const 
{
  return parameterbase::_parameters.validationKey();
}
