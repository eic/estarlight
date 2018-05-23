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
// $Rev:: 211                         $: revision of last commit
// $Author:: mlomnitz                   $: author of last commit
// $Date:: 2017-03-14 03:05:09 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef e_WIDERESONANCECROSSSECTION_H
#define e_WIDERESONANCECROSSSECTION_H


#include "photonNucleusCrossSection.h"
#include "inputParameters.h"

class e_wideResonanceCrossSection : public photonNucleusCrossSection {

public:

	e_wideResonanceCrossSection(const inputParameters& input, const beamBeamSystem& bbsystem);
	~e_wideResonanceCrossSection();

	void crossSectionCalculation(const double bwnormsave);
	void makeGammaPQ2dependence(const double bwnormsave);
	void printCrossSection(const std::string name, const double x_section);
private:

	double _wideWmax;
	double _wideWmin;
	double _targetMaxPhotonEnergy;
	double _targetMinPhotonEnergy;
	double _Ep;
	double _electronEnergy;
	double _target_beamLorentz;
	double _VMnumEgamma;
	double _useFixedRange;
	double _gammaMinQ2;
	double _gammaMaxQ2;
	double _targetRadii;

};


#endif  // WIDERESONANCECROSSSECTION_H
