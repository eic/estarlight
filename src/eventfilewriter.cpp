
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
// $Rev:: 228                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2016-01-18 17:10:17 +0000 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////

#ifdef HEPMC3_ON
#include "hepmc3writer.h"
#endif

#include <iostream>
#include <exception>
#include <cstdlib>

#include "eventfilewriter.h"
#include "starlightparticlecodes.h"

using namespace std;

//______________________________________________________________________________
eventFileWriter::eventFileWriter()
: fileWriter()
,_writeFullPythia(false)
,_writeFullHepMC3(false)
{ }

//______________________________________________________________________________
eventFileWriter::eventFileWriter(std::string filename)
: fileWriter(filename)
{ }

//______________________________________________________________________________
int eventFileWriter::writeInit(inputParameters &_p)
{
  _fileStream<<"CONFIG_OPT: "<<_p.productionMode()<<" "<<_p.prodParticleId()<<" "<<_p.nmbEvents()
	     <<" "<<_p.quantumGlauber()<<" "<<_p.impulseVM()<<" "<<_p.randomSeed()<<std::endl;
  _fileStream<<"BEAM_1: "<<_p.beam1Z()<<" "<<_p.beam1A()<<" "<<_p.beam1LorentzGamma()<<std::endl;
  _fileStream<<"BEAM_2: "<<_p.beam2Z()<<" "<<_p.beam2A()<<" "<<_p.beam2LorentzGamma()<<std::endl;
  _fileStream<<"PHOTON: "<<_p.nmbEnergyBins()<<" "<<_p.fixedQ2Range()<<" "<<_p.nmbGammaQ2Bins()
	     <<" "<<_p.minGammaQ2()<<" "<<_p.maxGammaQ2()<<std::endl;
  
  if ( _writeFullHepMC3 )
    {
#ifdef HEPMC3_ON
      _hepmc3writer = new hepMC3Writer();
      _hepmc3writer->initWriter(_p);
#else
      std::cout << "***HepMC3 file format not written --- eStarlight not compiled with HepMC3 Enabled***" << std::endl;
      _writeFullHepMC3 = false;
#endif
    }
  
  return 0;
}

//______________________________________________________________________________
int eventFileWriter::writeEvent(eXEvent &event, int eventnumber)
{
   
    int numberoftracks = 0;
    if(_writeFullPythia)
    {
        numberoftracks = event.getParticles()->size();
    }
    else
    {
        for(unsigned int i = 0; i<event.getParticles()->size(); ++i)
        {
            if(event.getParticles()->at(i).getStatus() >= 0) numberoftracks++;
        }
    }
    

    // sometimes we don't have tracks due to wrongly picked W , check it
    if(numberoftracks){
      eventnumber++;
      
      _fileStream << "EVENT: " << eventnumber << " " << numberoftracks << " " << 1 << std::endl;

      _fileStream <<"VERTEX: "<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<" "<<1<<" "<<0<<" "<<0<<" "<<numberoftracks<<std::endl;

      for( uint igam = 0 ; igam < event.getGammaEnergies()->size(); ++igam){
	_fileStream <<"GAMMA: "<<event.getGammaEnergies()->at(igam)<<" "<<event.getGammaMasses()->at(igam)<<std::endl;
	lorentzVector gam = event.getGamma()->at(igam);
	// Write the photon 4-vector out to file. Might be used in later iterations, so I keep it here
	//_fileStream <<"GAMMA VECTOR: "<<gam.GetPx()<<" "<<gam.GetPy()<<" "<<gam.GetPz()<<" "<<gam.GetE()<<" "<<-temp<<std::endl;
      }
      for( uint itarget = 0 ; itarget < event.getTarget()->size(); ++itarget){
	lorentzVector target = event.getTarget()->at(itarget);
	_fileStream <<"t: "<<event.getVertext()->at(itarget)<<std::endl;
	_fileStream <<"TARGET: "<<target.GetPx()<<" "<<target.GetPy()<<" "<<target.GetPz()<<" "<<target.GetE()<<std::endl;
      }
      for( uint iel = 0 ; iel < event.getSources()->size(); ++iel){
	lorentzVector el = event.getSources()->at(iel); 
	_fileStream <<"SOURCE: "<<el.GetPx()<<" "<<el.GetPy()<<" "<<el.GetPz()<<" "<<el.GetE()<<std::endl;
      }
	     
      int ipart = 0;
      std::vector<starlightParticle>::const_iterator part = (event.getParticles())->begin();
      
      for (part = event.getParticles()->begin(); part != event.getParticles()->end(); part++, ipart++)
	{
          if(!_writeFullPythia) 
          {
              if((*part).getStatus() < 0) continue;
          }
	  _fileStream << "TRACK: " << " "<< starlightParticleCodes::jetsetToGeant((*part).getPdgCode()) <<" "<< (*part).GetPx() << " " << (*part).GetPy()
		      << " "<< (*part).GetPz() << " " << eventnumber << " " << ipart << " " << 0 << " "
		      << (*part).getPdgCode();
		      
	  if(_writeFullPythia)
	  {
	    lorentzVector vtx = (*part).getVertex();
	    _fileStream << " " << vtx.GetPx() << " " << vtx.GetPy() << " " << vtx.GetPz() << " " << vtx.GetE();
	    _fileStream << " " << (*part).getFirstParent() << " " << (*part).getLastParent() << " " << (*part).getFirstDaughter() << " " << (*part).getLastDaughter() << " " << (*part).getStatus();
	  }
		      
	  _fileStream <<std::endl;
	}
    }

#ifdef HEPMC3_ON
    if ( _writeFullHepMC3 ){
      _hepmc3writer->writeEvent(event,eventnumber);
    }
#endif
    
    return 0;
}

int eventFileWriter::close()
{
  
#ifdef HEPMC3_ON
    if ( _writeFullHepMC3 ){
      _hepmc3writer->close();
    }
#endif
  
    try
    {
        _fileStream.close();
    }
    catch (const ios::failure & error)
    {
        cerr << "I/O exception: " << error.what() << endl;
        return EXIT_FAILURE;
    }
    return 0;
}
