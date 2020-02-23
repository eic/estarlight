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
// $Rev:: 265                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2016-06-06 21:37:26 +0100 #$: date of last commit
//
// Description:
//    Added 18->19 for reading in the luminosity table
//    Incoherent factor added to table --Joey
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "readinluminosity.h"
#include "starlightconstants.h"
#include "inputParameters.h"

using namespace std;


//______________________________________________________________________________
readLuminosity::readLuminosity(const inputParameters& inputParametersInstance)
  : _Warray(0), _Yarray(0), _Farray(0), _Farray1(0), _Farray2(0),
    _f_WYarray(0), _g_Earray(0)/*, _g_EQ2array(0)*/
  , _ReadInputNPT(inputParametersInstance.nmbPtBinsInterference())
  , _ReadInputnumy(inputParametersInstance.nmbRapidityBins())
  , _ReadInputnumw(inputParametersInstance.nmbWBins())
  , _ReadInputnumega(inputParametersInstance.nmbEnergyBins())
  , _ReadInputnumQ2(inputParametersInstance.nmbGammaQ2Bins())
  , _ReadInputgg_or_gP(inputParametersInstance.productionMode())
  , _ReadInputinterferencemode(inputParametersInstance.interferenceEnabled())
  , _baseFileName(inputParametersInstance.baseFileName())
{

}


//______________________________________________________________________________
readLuminosity::~readLuminosity()
{ 
  if(_Warray) delete [] _Warray;
  if(_Yarray) delete [] _Yarray;
  if(_Farray) delete [] _Farray;
  if(_Farray1) delete [] _Farray1;
  if(_Farray2) delete [] _Farray2; 
  // For eSTARLIGHT
  if(_f_WYarray) delete [] _f_WYarray;
  if(_g_Earray) delete [] _g_Earray;
  if(_g_EQ2array) delete _g_EQ2array;
}


//______________________________________________________________________________
void readLuminosity::read()
{
  
  if(!_Warray) _Warray = new double[_ReadInputnumw];
  if(!_Yarray) _Yarray = new double[_ReadInputnumy];
  if(!_Farray) 
  {
    _Farray = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _Farray[i] = new double[_ReadInputnumy];
    }
  }
  if(!_Farray1) 
  {
    _Farray1 = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _Farray1[i] = new double[_ReadInputnumy];
    }
  }
  if(!_Farray2) 
  {
    _Farray2 = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _Farray2[i] = new double[_ReadInputnumy];
    }
  }

  double dummy[17]; //number of lines used to read in input parameters saved to lookup table[slight.txt].


  std::string wyFileName;
  wyFileName = _baseFileName +".txt";
  
//  cout << "wyFileName being read in" << wyFileName << endl;

  double fpart =0.;
  double fptsum=0.;
  ifstream wylumfile;

  _f_max=0.0;
  _f_max1=0.0;
  _f_max2=0.0;

  wylumfile.open(wyFileName.c_str());

  for(int i=0;i < 17;i++){ 
    wylumfile >> dummy[i];
  }
  int A_1 = dummy[1];
  int A_2 = dummy[3];

  for(int i=0;i<_ReadInputnumw;i++){
    wylumfile >> _Warray[i];
  }
  for(int i=0;i<_ReadInputnumy;i++){
    wylumfile >> _Yarray[i];
  }

  if( (_ReadInputgg_or_gP == 1) || (A_2 == 1 && A_1 != 1) || (A_1 ==1 && A_2 != 1) ){ 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _Farray[i][j];
        if( _Farray[i][j] > _f_max ) _f_max=_Farray[i][j];
      }
    }
    //Normalize farray 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        _Farray[i][j]=_Farray[i][j]/_f_max;
      }
    }
  } else {
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _Farray1[i][j];
      }
    }
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _Farray2[i][j];
        if( _Farray1[i][j] + _Farray2[i][j] > _f_max ) _f_max=(_Farray1[i][j] + _Farray2[i][j]);
      }
    }
    //Normalize farray, farray1, farray2 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        _Farray1[i][j]=_Farray1[i][j]/_f_max;
        _Farray2[i][j]=_Farray2[i][j]/_f_max;
        _Farray[i][j]=_Farray1[i][j]+_Farray2[i][j];
      }
    }
  }
  wylumfile >> _bwnormsave;

  if (_ReadInputgg_or_gP != 1 && _ReadInputinterferencemode != 0) {
        // only numy/2 y bins here, from 0 (not -ymax) to ymax
        double **finterm  = new double*[starlightLimits::MAXWBINS];
        for (int i = 0; i < starlightLimits::MAXWBINS; i++) finterm[i] = new double[starlightLimits::MAXYBINS];
        for (int i=0;i<_ReadInputnumy;i++) {
            //fmax=0;
            //we want to convert _fptarray to an integral array where fpt(i,j) is near 0, and fpt(j,NPT) ~1. This will facilitate a simple table loookup
            fptsum=0.;
            for (int j=0;j<_ReadInputNPT;j++) {
                wylumfile >> fpart;
                finterm[i][j] = fpart;
                _fptarray[i][j]=0.;
                fptsum=fptsum+fpart;
            }
            //convert array to integral
            _fptarray[i][0]=finterm[i][0]/fptsum;
            for (int j=1;j<_ReadInputNPT;j++) {
                for (int k=0;k<j;k++) {
                    _fptarray[i][j]=_fptarray[i][j]+finterm[i][k];
                }
                _fptarray[i][j]=_fptarray[i][j]/fptsum;
            }
        }
        delete [] finterm;

    }
  wylumfile.close();
  return;
}


//______________________________________________________________________________
void readLuminosity::e_read()
{
  if(!_Warray) _Warray = new double[_ReadInputnumw];
  if(!_Yarray) _Yarray = new double[_ReadInputnumy];
  if(!_BWarray) _BWarray = new double[_ReadInputnumw];
  if(!_f_WYmax) 
  {
    _f_WYarray = new double*[_ReadInputnumega];
    for(int i = 0; i < _ReadInputnumega; i++)
    {
      _f_WYarray[i] = new double[_ReadInputnumQ2];
    }
  }
  if(!_g_Earray) 
  {
    _g_Earray = new double*[_ReadInputnumega];
    for(int i = 0; i < _ReadInputnumega; i++)
    {
      _g_Earray[i] = new double[_ReadInputnumQ2];
    }
  }

  double dummy[13]; //number of lines used to read in input parameters saved to lookup table[slight.txt].


  std::string wyFileName;
  wyFileName = _baseFileName +".txt";
  
//  cout << "wyFileName being read in" << wyFileName << endl;

  ifstream wylumfile;

  _f_WYmax=0.0;
  _g_Emax=0.0;

  wylumfile.open(wyFileName.c_str());

  for(int i=0;i < 13;i++){ 
    wylumfile >> dummy[i];
  }
  int A_1 = dummy[1];
  int A_2 = dummy[3];

  for(int i=0;i<_ReadInputnumw;i++){
    wylumfile >> _Warray[i];
  }
  //No longer needed, can be removed later
  for(int i=0;i<_ReadInputnumy;i++){
    wylumfile >> _Yarray[i];
  }

  //New table saving the Breit-Wigner function for form factor
  double bw_max = 0 ;
  for(int i=0;i<_ReadInputnumw;i++){
    wylumfile >> _BWarray[i];
    if( _BWarray[i] > bw_max )
      bw_max = _BWarray[i];
  }
  for(int i=0;i<_ReadInputnumw;i++){
    _BWarray[i] = _BWarray[i]/bw_max;
  }


  if( (A_2 == 0 && A_1 >= 1) || (A_1 ==0 && A_2 >= 1) ){ 
    for(int i=0;i<_ReadInputnumega;i++){
      for(int j=0;j<_ReadInputnumQ2;j++){
        wylumfile >> _f_WYarray[i][j];
        if( _f_WYarray[i][j] > _f_WYmax ) _f_WYmax=_f_WYarray[i][j];
	//
	//wylumfile >> _g_Earray[i][j];
	//if( _g_Earray[i][j] > _g_Emax ) _g_Emax = _g_Earray[i][j];
      }
    }
    //Normalize f_WY array, g does not need to be normalized, it is used for normalization
    for(int i=0;i<_ReadInputnumega;i++){
      for(int j=0;j<_ReadInputnumQ2;j++){
        _f_WYarray[i][j] = _f_WYarray[i][j]/( _f_WYmax );
	//_g_Earray[i][j] = _g_Earray[i][j]/_g_Emax;
      }
    }
  }
  wylumfile >> _bwnormsave;

  wylumfile.close();
  cout<<"Done reading wylumi file"<<endl;
  std::string EQ2FileName;
  EQ2FileName = "e_"+_baseFileName+".txt";
  ifstream EQlumfile;
  
  EQlumfile.open(EQ2FileName.c_str());
  int n_steps;
  EQlumfile >> n_steps;
  double integrated_max = 0 ;
  //  _g_EQ2array = new map<string,std::vector<double> >();
  _g_EQ2array = new vector< std::pair< double, vector<double> > > ();
  while( !EQlumfile.eof() ){
    _g_EQ2max = 0 ;
    double integral;
    std::vector<double> p;
    for( int iQ2 = 0 ; iQ2 < _ReadInputnumQ2+3 ; ++iQ2){      
      if(iQ2 == 0 ){
	EQlumfile >> integral;
	if( integral > integrated_max)
	  integrated_max = integral;
      }
      else{
	double temp;
	EQlumfile >> temp;
	p.push_back(temp);
	//cout<<p[iQ2-1]<<endl;
      }
      if( iQ2 > 2 && p[iQ2-1] > _g_EQ2max){	
	_g_EQ2max = p[iQ2-1];
      }
    }
    for( unsigned int iQ2=2; iQ2 < p.size(); ++iQ2)
      p[iQ2] = p[iQ2]/_g_EQ2max;
    
    _g_EQ2array->push_back(std::pair< double, std::vector<double> >(integral,p));
  }
  for(  std::vector< std::pair<double, std::vector<double> > >::iterator it =_g_EQ2array->begin() ; it != _g_EQ2array->end(); ++it){
    //cout<<it->first<<endl;
    it->first = it->first/integrated_max;
    //cout<<"Now"<<it->first<<endl;
  }

  //normalize
  /*for(  std::map<string,std::vector<double> >::iterator it =_g_EQ2array->begin() ; it != _g_EQ2array->end(); ++it){
    std::string key = it->first;
    if(key=="")
      continue;
    for(unsigned int iQ2=2; iQ2 < it->second.size(); ++iQ2)
      it->second[iQ2] = it->second[iQ2]/_g_EQ2max;
      }*/
  EQlumfile.close();
  return;  
  
}

