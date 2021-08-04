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
// $Rev:: 264                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-06-06 21:05:12 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "inputParameters.h"
#include "beambeamsystem.h"
#include "beam.h"
#include "starlightconstants.h"
#include "nucleus.h"
#include "bessel.h"
#include "gammaeluminosity.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
photonElectronLuminosity::photonElectronLuminosity(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem)
  : photonNucleusCrossSection(inputParametersInstance, bbsystem)
  ,_protonEnergy(inputParametersInstance.protonEnergy())
  ,_electronEnergy(inputParametersInstance.electronEnergy())
  ,_beamLorentzGamma(inputParametersInstance.beamLorentzGamma())
  ,_baseFileName(inputParametersInstance.baseFileName())
  ,_maxW(inputParametersInstance.maxW())
  ,_minW(inputParametersInstance.minW())
  ,_nmbWBins(inputParametersInstance.nmbWBins())
  ,_maxRapidity(inputParametersInstance.maxRapidity())
  ,_nmbRapidityBins(inputParametersInstance.nmbRapidityBins())
  ,_nEBins(inputParametersInstance.nmbEnergyBins())
  ,_minGammaQ2(inputParametersInstance.minGammaQ2())
  ,_maxGammaQ2(inputParametersInstance.maxGammaQ2())
  ,_nmbGammaQ2Bins(inputParametersInstance.nmbGammaQ2Bins()) 
  ,_cmsMaxPhotonEnergy(inputParametersInstance.cmsMaxPhotonEnergy())
  ,_cmsMinPhotonEnergy(inputParametersInstance.cmsMinPhotonEnergy())
  ,_targetMaxPhotonEnergy(inputParametersInstance.targetMaxPhotonEnergy())
  ,_targetMinPhotonEnergy(inputParametersInstance.targetMinPhotonEnergy())
  ,_productionMode(inputParametersInstance.productionMode())
  ,_beamBreakupMode(inputParametersInstance.beamBreakupMode())
{
  cout <<"Creating Luminosity Tables."<<endl;
  photonNucleusDifferentialLuminosity();
  cout <<"Luminosity Tables created."<<endl;
}


//______________________________________________________________________________
photonElectronLuminosity::~photonElectronLuminosity()
{ }


//______________________________________________________________________________
void photonElectronLuminosity::photonNucleusDifferentialLuminosity()
{
  double W,dW, dY;
  double Egamma,Y;
  //
  double testint;
  // 
  double g_E;
  double csgA;
  double C;  
  int beam; 
  //
  //int nQ2steps =100;
  
  std::string wyFileName;
  wyFileName = _baseFileName +".txt";

  ofstream wylumfile;
  wylumfile.precision(15);
  
  std::string EQ2FileName;
  EQ2FileName = "e_"+_baseFileName +".txt";

  ofstream EQ2lumfile;
  EQ2lumfile.precision(15);

  double  bwnorm;

  dW = (_wMax-_wMin)/_nWbins;
  dY  = (_yMax-(-1.0)*_yMax)/_nYbins;

  // Look-up table storing double and single photon flux differentials
  EQ2lumfile.open(EQ2FileName.c_str());
  EQ2lumfile << _nmbGammaQ2Bins <<endl;
  // Look-up table storing photo-nuclear cross section
  // Write the values of W used in the calculation to slight.txt.  
  wylumfile.open(wyFileName.c_str());
  wylumfile << getbbs().beam1().Z() <<endl;
  wylumfile << getbbs().beam1().A() <<endl;
  wylumfile << getbbs().beam2().Z() <<endl;
  wylumfile << getbbs().beam2().A() <<endl;
  wylumfile << _beamLorentzGamma <<endl;
  wylumfile << _maxW <<endl;
  wylumfile << _minW <<endl;
  wylumfile << _nmbWBins <<endl;
  wylumfile << _maxRapidity <<endl;
  wylumfile << _nmbRapidityBins <<endl;
  wylumfile << _productionMode <<endl;
  wylumfile << _beamBreakupMode <<endl;
  wylumfile << starlightConstants::deuteronSlopePar <<endl;
  
  // Normalize the Breit-Wigner Distribution and write values of W to slight.txt
  testint=0.0;
  // Grabbing default value for C in the breit-wigner calculation
  C=getDefaultC();
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {
    W = _wMin + double(i)*dW + 0.5*dW;
    testint = testint + breitWigner(W,C)*dW;
    wylumfile << W << endl;
  }
  bwnorm = 1./testint;
  
  // Write the values of Y used in the calculation to slight.txt.
  for(unsigned int i = 0; i <= _nYbins - 1; ++i) {
    Y = -1.0*_yMax + double(i)*dY + 0.5*dY;
    wylumfile << Y << endl;
  }
    
  int A_1 = getbbs().beam1().A(); 
  int A_2 = getbbs().beam2().A();
  if( A_2 == 0 && A_1 != 0 ){
    beam = 1;
  } else if( A_1 ==0 && A_2 != 0){
    beam = 2;
  } else{
    beam = 2;
  }
  // Exponential steps for both tables: Target frame for photon generation and CMS frame for final state generation
  double dE_target = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy)/_nEBins;
  // - - - - Calculations in Target frame
  for(unsigned int i = 0; i < _nWbins; ++i) {

    W = _wMin + double(i)*dW + 0.5*dW;
    wylumfile << breitWigner(W,bwnorm) <<endl;
  }
  for(int j = 0; j < _nEBins ; ++j) { 
    // Used for calculation of g(Egamma) which must be done in the ion (target) rest frame
    Egamma = std::exp( std::log(_targetMinPhotonEnergy)+j*dE_target); 
    g_E = 0;
    // Photon energy limits are determnined in target frame - multiply by Egamma to speed up event generation
    std::pair< double, double >* this_energy = Q2arraylimits(Egamma);
    double Q2min = this_energy->first;
    double Q2max = this_energy->second;
    if( Q2min > Q2max)
      continue;
    //Accounts for the phase space factor 
    g_E = Egamma*integrated_Q2_dep(Egamma, Q2min, Q2max);
    EQ2lumfile << g_E<<endl;
    // Write out Q2 range used for generation
    EQ2lumfile << Q2min << endl;
    EQ2lumfile << Q2max << endl;
    for( unsigned int iQ2 =0 ;iQ2<_nmbGammaQ2Bins; ++iQ2){
      double Q2 = std::exp( std::log(Q2min)+iQ2*std::log(Q2max/Q2min)/_nmbGammaQ2Bins );
      EQ2lumfile<< g(Egamma,Q2) <<endl;
      csgA=getcsgA(Egamma,Q2,beam);
      wylumfile << csgA << endl;
    }
  }
  
  EQ2lumfile.close();

  wylumfile << bwnorm << endl;
  wylumfile.close();
   
}
string photonElectronLuminosity::gammaTableParse(int ii, int jj)
{
  ostringstream tag1, tag2;
  tag1<<ii;
  tag2<<jj;
  string to_ret = tag1.str()+","+tag2.str();
  return to_ret;
}

