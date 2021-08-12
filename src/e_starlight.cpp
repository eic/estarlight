///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2017
//
//    This file is part of estarlight.
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
// $Rev:: 263                         $: revision of last commit
// $Author:: mlomnitz                  $: author of last commit
// $Date:: 02/28/2017 #$: date of last commit
//
// Description:
//
// Main interface for eSTARlight code inherited from STARlight
//
//
///////////////////////////////////////////////////////////////////////////
 

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "starlightconfig.h"

#ifdef ENABLE_PYTHIA
#include "PythiaStarlight.h"
#endif

#ifdef ENABLE_PYTHIA6
#include "starlightpythia.h"
#endif

#include "reportingUtils.h"
#include "inputParameters.h"
#include "eventchannel.h"
#include "gammagammasingle.h"
#include "gammaavm.h"
#include "twophotonluminosity.h"
#include "gammaaluminosity.h"
#include "gammaeluminosity.h"
#include "incoherentPhotonNucleusLuminosity.h"
#include "e_starlight.h"


using namespace std;
using namespace starlightConstants;


e_starlight::e_starlight() :
		_beamSystem            (0),
		_eventChannel          (0),
		_nmbEventsPerFile      (100),
		_isInitialised         (false),
		_inputParameters       (0)
{ }


e_starlight::~e_starlight()
{ }


bool
e_starlight::init()
{
  // ??? 
        if(Starlight_VERSION_MAJOR == 9999)
	{
	  cout  << endl << "#########################################" << endl
	     	<< " Initialising Starlight version: trunk..." << endl
	     	<< "#########################################" << endl << endl;
	}
	else
	{
	  cout  << endl << "#########################################" << endl
	     	<< " Initialising Starlight version: " << Starlight_VERSION_MAJOR << "."
	     	<< Starlight_VERSION_MINOR << "." << Starlight_VERSION_MINOR_MINOR << "..." << endl
	        << "#########################################" << endl << endl;
	}

	_nmbEventsPerFile    = _inputParameters->nmbEvents();  // for now we write only one file...

	_beamSystem = new beamBeamSystem(*_inputParameters);
	
	// This sets the precision used by cout for printing
	cout.setf(ios_base::fixed,ios_base::floatfield);
	cout.precision(3);

        std::string _baseFileName;
        _baseFileName = _inputParameters->baseFileName();
       _lumLookUpTableFileName = _baseFileName + ".txt";

	const bool lumTableIsValid = luminosityTableIsValid();

	// Do some sanity checks of the input parameters here.
	if( _inputParameters->targetBeamA() == 0 ){
	  //	  printErr << endl << "One of the two beams must be an electron in eSTARlight" << endl;
	  return false;
	}

	bool createChannel = true;
	switch (_inputParameters->interactionType())	{
	  // Photon-photon scattering is not implemented in eSTARlight as of yet
	  /*
	    case PHOTONPHOTON:
		if (!lumTableIsValid) {
			printInfo << "creating luminosity table for photon-photon channel" << endl;
			twoPhotonLuminosity(*_inputParameters, _beamSystem->beam1(), _beamSystem->beam2());
		}
		break;		
	  */
	case E_PHOTONPOMERONNARROW:  // narrow and wide resonances use
	case E_PHOTONPOMERONWIDE:    // the same luminosity function
		if (!lumTableIsValid) {
			printInfo << "creating luminosity table for coherent photon-Pomeron channel" <<endl;
			photonElectronLuminosity lum(*_inputParameters, *_beamSystem);
		}
		break;
		// Incoherent electro-production is still pending 
		/*        case PHOTONPOMERONINCOHERENT:  // narrow and wide resonances use
                if (!lumTableIsValid) {
                        printInfo << "creating luminosity table for incoherent photon-Pomeron channel" << endl;
                        incoherentPhotonNucleusLuminosity lum(*_inputParameters, *_beamSystem);
			}
			break;*/
	default:
		{
			printWarn << "unknown interaction type '" << _inputParameters->interactionType() << "'."
			          << " cannot initialize eSTARlight." << endl;
			return false;
		}
	}
	
	if(createChannel)
	{
	  if (!createEventChannel())
		  return false;
	}

	_isInitialised = true;
	return true;
}


eXEvent
e_starlight::produceEvent()
{
	if (!_isInitialised) {
	  //		printErr << "trying to generate event but Starlight is not initialised. aborting." << endl;
		exit(-1);
	}
	return _eventChannel->e_produceEvent();
}


bool
e_starlight::luminosityTableIsValid() const
{
	printInfo << "using random seed = " << _inputParameters->randomSeed() << endl;

	ifstream lumLookUpTableFile(_lumLookUpTableFileName.c_str());
	lumLookUpTableFile.precision(15);
	if ((!lumLookUpTableFile) || (!lumLookUpTableFile.good())) {
		return false;
	}

	unsigned int targetBeamZ, targetBeamA;
	double       beamLorentzGamma = 0;
	double       maxW = 0, minW = 0;
	unsigned int nmbWBins;
	double       maxRapidity = 0;
	unsigned int nmbRapidityBins;
	int          productionMode, beamBreakupMode;
	//Remove interference options for eX collisions (no logner symmetric)
	//	bool         interferenceEnabled = false;
	//	double       interferenceStrength = 0;
	bool         coherentProduction = false;
	double       incoherentFactor = 0, deuteronSlopePar = 0, maxPtInterference = 0;
	int          nmbPtBinsInterference;
	std::string  validationKey;
	if (!(lumLookUpTableFile
	      >> validationKey
	      >> targetBeamZ >> targetBeamA
	      >> beamLorentzGamma
	      >> maxW >> minW >> nmbWBins
	      >> maxRapidity >> nmbRapidityBins
	      >> productionMode
	      >> beamBreakupMode
	      >> maxPtInterference
	      >> nmbPtBinsInterference
	      >> coherentProduction >> incoherentFactor
	      >> deuteronSlopePar))
	  // cannot read parameters from lookup table file
		return false;
        	
	std::string validationKeyEnd;
	while(!lumLookUpTableFile.eof())
	{
	  lumLookUpTableFile >> validationKeyEnd; 
	}
	lumLookUpTableFile.close();
	return (validationKey == _inputParameters->parameterValueKey() && validationKeyEnd == validationKey);
	return true;
}


bool
e_starlight::createEventChannel()
{
	switch (_inputParameters->prodParticleType()) {
	case ELECTRON:
	case MUON:
	case TAUON:
        case TAUONDECAY:
	  printWarn << "cannot construct Gammagammaleptonpair event channel." << endl;
	  return false;
	
	case A2:        // jetset/pythia
	case ETA:       // jetset/pythia
	case ETAPRIME:  // jetset/pythia
	case ETAC:      // jetset/pythia
	case F0:        // jetset/pythia
		{

		  //#ifdef ENABLE_PYTHIA
			// PythiaOutput = true;
		  // 		        cout<<"Pythia is enabled!"<<endl;
// 			return true;
//#else
//			printWarn << "Starlight is not compiled against Pythia8; "
//			          << "jetset event channel cannot be used." << endl;
// 			return false;
//#endif
		}
	case F2:
	case F2PRIME:
	case ZOVERZ03:
    case AXION: // AXION HACK
		{
		  //  #ifdef ENABLE_PYTHIA
	 	        cout<<" This is f2, f2prim, rho^0 rho^0, or axion "<<endl; 
			_eventChannel= new Gammagammasingle(*_inputParameters, *_beamSystem);
			if (_eventChannel)
				return true;
			else {
				printWarn << "cannot construct Gammagammasingle event channel." << endl;
				return false;
			}
		}
	case RHO:
	case RHOZEUS:
	case FOURPRONG:
	case OMEGA:  
	case OMEGA_pi0gamma:
	case PHI:
	case JPSI:
	case JPSI2S:
	case JPSI2S_ee:
	case JPSI2S_mumu:
	case JPSI_ee:
	case JPSI_mumu:
	case UPSILON:
	case UPSILON_ee:
	case UPSILON_mumu:
	case UPSILON2S:
	case UPSILON2S_ee:
	case UPSILON2S_mumu:
	case UPSILON3S:
	case UPSILON3S_ee:
	case UPSILON3S_mumu:
		{
			if (_inputParameters->interactionType() == E_PHOTONPOMERONNARROW) {
				_eventChannel = new e_Gammaanarrowvm(*_inputParameters, *_beamSystem);
				if (_eventChannel)
					return true;
				else {
					printWarn << "cannot construct Gammaanarrowvm event channel." << endl;
					return false;
				}
			}

			if (_inputParameters->interactionType() == E_PHOTONPOMERONWIDE) {
				_eventChannel = new e_Gammaawidevm(*_inputParameters, *_beamSystem);
				if (_eventChannel)
					return true;
				else {
					printWarn << "cannot construct Gammaawidevm event channel." << endl;
					return false;
				}
			}


			printWarn << "interaction type '" << _inputParameters->interactionType() << "' "
			          << "cannot be used with particle type '" << _inputParameters->prodParticleType() << "'. "
			          << "cannot create event channel." << endl;
			return false;
		}
	default:
		{
			printWarn << "unknown event channel '" << _inputParameters->prodParticleType() << "'."
			          << " cannot create event channel." << endl;
			return false;
		}
	}
}
