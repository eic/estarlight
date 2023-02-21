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
// $Rev:: 0                            $: revision of last commit
// $Author:: zsweger                   $: author of last commit
// $Date:: 2022-07-29                 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef e_PI0CROSSSECTION_H
#define e_PI0CROSSSECTION_H


#include "photonNucleusCrossSection.h"


class e_Pi0CrossSection : public photonNucleusCrossSection {

public:

	e_Pi0CrossSection(const inputParameters& input, const beamBeamSystem&  bbsystem);
	~e_Pi0CrossSection();

	void crossSectionCalculation(const double bwnormsave);
	void makeGammaPQ2dependence();
	void printCrossSection(const std::string name, const double x_section);

private:
	
	double _Ep;
	double _electronEnergy;
	double _target_beamLorentz;
	double _boost;
	//
	double _targetMaxPhotonEnergy;
	double _targetMinPhotonEnergy;
	double _cmsMaxPhotonEnergy;
	double _cmsMinPhotonEnergy;
	//
	double _VMnumEgamma;
	double _useFixedRange;
	double _gammaMinQ2;
	double _gammaMaxQ2;
	double _targetRadii;
	bool _backwardsProduction;
};


#endif  // PI0CROSSSECTION_H
