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
// Main interfqce for eSTARlight code inherited from STARlight
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef eSTARLIGHT_H
#define eSTARLIGHT_H


#include <string>

#include "eventchannel.h"


class inputParameters;
class beam;
class beamBeamSystem;


class e_starlight {

public:
      
	e_starlight();
	~e_starlight();
      
	bool init();

	eXEvent produceEvent();
      
        std::string   baseFileName  () const { return _baseFileName;		  }
	unsigned long nmbAttempts   () const { return _eventChannel->nmbAttempts(); }
	unsigned long nmbAccepted   () const { return _eventChannel->nmbAccepted(); } 
	double getTotalCrossSection () const { return _eventChannel->getTotalChannelCrossSection(); } 
	void setInputParameters(inputParameters* inputParams) { _inputParameters = inputParams; }   

private:
      
	bool luminosityTableIsValid() const;

	bool createEventChannel();
      
	beamBeamSystem*    _beamSystem;
	eventChannel*      _eventChannel;
	unsigned int       _nmbEventsPerFile;
	std::string        _baseFileName;
	std::string        _lumLookUpTableFileName;
	bool               _isInitialised;
	inputParameters*   _inputParameters;
};


#endif  // e_STARLIGHT_H
