
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
  _fileStream<<"ELECTRON_BEAM: "<<" "<<_p.electronBeamLorentzGamma()<<std::endl;
  _fileStream<<"TARGET_BEAM: "<<_p.targetBeamZ()<<" "<<_p.targetBeamA()<<" "<<_p.targetBeamLorentzGamma()<<std::endl;
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
int eventFileWriter::writeInitLUND(inputParameters &_p)
{
  _fileStream<<"PYTHIA EVENT FILE"<<std::endl;
  _fileStream<<"============================================"<<std::endl;
  _fileStream<<"I, ievent, genevent, subprocess, nucleon,                 targetparton, xtargparton, beamparton, xbeamparton,               thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu, leptonphi,   s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R, sigma_rad,       SigRadCor, EBrems, photonflux, t-diff, nrTracks"<<std::endl;
  _fileStream<<"============================================"<<std::endl;
  _fileStream<<"I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)"<<std::endl;
  _fileStream<<"============================================"<<std::endl;

  //Proton and Electron Mass
  double p_mass = 0.93827;
  double e_mass = 0.000511;

  //Calculate Four Vector
  double _electronBeam_E  = e_mass*_p.electronBeamLorentzGamma();
  double _targetBeam_E  = p_mass*_p.targetBeamLorentzGamma();
  double _electronBeam_pz = -1.0*sqrt(_electronBeam_E*_electronBeam_E - e_mass*e_mass);
  double _targetBeam_pz = sqrt(_targetBeam_E*_targetBeam_E - p_mass*p_mass);

  //Define Four Vector
  _electronBeam_four_vector_={0,0,_electronBeam_pz,_electronBeam_E};
  _targetBeam_four_vector_ ={0,0,_targetBeam_pz,_targetBeam_E};
  _electronBeam_pdg_id_ = 11;
  int nuclear_pid = _p.targetBeamA()*10 + _p.targetBeamZ()*10000 + 1000000000;// Form 100ZZZAAAl; l=0
  _targetBeam_pdg_id_ = ( _p.targetBeamA() > 1 ) ? nuclear_pid : 2212; //Everything is a proton right now
  
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


//______________________________________________________________________________
int eventFileWriter::writeEventLUND(eXEvent &event, int eventnumber)
{
    _fileStream << "0         "<<eventnumber<<"         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0"<<std::endl;
    _fileStream<<"============================================"<<std::endl;

    int particleIndex=1;
    uint lineNumberOfFirstElectronDaughter=3;
    uint lineNumberOfLastElectronDaughter=4;
    uint lineNumberOfFirstIonDaughter=6;
    uint lineNumberOfLastIonDaughter=lineNumberOfFirstIonDaughter+event.getParticles()->size();
    uint numberOfElectrons = event.getSources()->size();
    uint numberOfIons = event.getTarget()->size();
    uint numberOfGammas = event.getGammaEnergies()->size();

    auto sp0 = std::setw(8);//column widths
    auto sp1 = std::setw(14);

    double electronmass = sqrt(_electronBeam_four_vector_[3]*_electronBeam_four_vector_[3] - _electronBeam_four_vector_[2]*_electronBeam_four_vector_[2] - _electronBeam_four_vector_[1]*_electronBeam_four_vector_[1] - _electronBeam_four_vector_[0]*_electronBeam_four_vector_[0]);
    _fileStream <<particleIndex << sp0 <<21 << sp0 <<_electronBeam_pdg_id_<<"       "<<0<<"       "<<lineNumberOfFirstElectronDaughter<<"       "<<lineNumberOfLastElectronDaughter<<"  " << std::setprecision(6) << sp1 <<_electronBeam_four_vector_[0]<<"  " << sp1 <<_electronBeam_four_vector_[1]<<"  " << sp1 <<_electronBeam_four_vector_[2]<<"  " <<sp1<<_electronBeam_four_vector_[3]<<"   "<<sp1<<electronmass<<"       0       0       0"<<std::endl;
    particleIndex++;

    double ionmass = sqrt(_targetBeam_four_vector_[3]*_targetBeam_four_vector_[3] - _targetBeam_four_vector_[2]*_targetBeam_four_vector_[2] - _targetBeam_four_vector_[1]*_targetBeam_four_vector_[1] - _targetBeam_four_vector_[0]*_targetBeam_four_vector_[0]);
    _fileStream <<particleIndex <<sp0<<21 << sp0 <<_targetBeam_pdg_id_<<"       "<<0<<"       "<<lineNumberOfFirstIonDaughter<<"       "<<lineNumberOfLastIonDaughter<<"  " << sp1 <<_targetBeam_four_vector_[0]<<"  " <<sp1<<_targetBeam_four_vector_[1]<<"  " <<sp1<<_targetBeam_four_vector_[2]<<"  " <<sp1<<_targetBeam_four_vector_[3]<<"   "<<sp1<<ionmass<<"       0       0       0"<<std::endl;
    particleIndex++;

    //FOR ALEX
    for( uint iel = 0 ; iel < numberOfElectrons; ++iel){
      lorentzVector el = event.getSources()->at(iel); 
      _fileStream <<particleIndex<<sp0<<21 <<sp0<<11<<"       "<<1<<"       "<<0<<"       "<<0<<"  " <<sp1<<el.GetPx()<<"  " <<sp1<<el.GetPy()<<"  " <<sp1<<el.GetPz()<<"  " <<sp1<<el.GetE()<<"   "<<sp1<<electronmass<<"       0       0       0"<<std::endl;
      particleIndex++;
    }

    for( uint igam = 0 ; igam < numberOfGammas; ++igam){
      lorentzVector gam = event.getGamma()->at(igam); 
      _fileStream <<particleIndex <<sp0<<21 <<sp0<<22<<"       "<<1<<"       "<<0<<"       "<<0<<"  " <<sp1<<gam.GetPx()<<"  " <<sp1<<gam.GetPy()<<"  " <<sp1<<gam.GetPz()<<"  " <<sp1<<event.getGammaEnergies()->at(igam)<<"   "<<sp1<<event.getGammaMasses()->at(igam)<<"       0       0       0"<<std::endl;
      particleIndex++;
    }



    for( uint iel = 0 ; iel < numberOfElectrons; ++iel){
      lorentzVector el = event.getSources()->at(iel); 
      _fileStream <<particleIndex<<sp0<<1 <<sp0<<11<<"       "<<3<<"       "<<0<<"       "<<0<<"  " <<sp1<<el.GetPx()<<"  " <<sp1<<el.GetPy()<<"  " <<sp1<<el.GetPz()<<"  " <<sp1<<el.GetE()<<"   "<<sp1<<electronmass<<"       0       0       0"<<std::endl;
      particleIndex++;
    }

    //FIXME This is only set up for protons right now
    for( uint itarget = 0 ; itarget < numberOfIons; ++itarget){
      lorentzVector target = event.getTarget()->at(itarget); 
      _fileStream <<particleIndex <<sp0<<1 <<sp0<<2212<<"       "<<2<<"       "<<lineNumberOfFirstIonDaughter<<"       "<<lineNumberOfLastIonDaughter<<"  " <<sp1<<target.GetPx()<<"  " <<sp1<<target.GetPy()<<"  " <<sp1<<target.GetPz()<<"  " <<sp1<<target.GetE()<<"   "<<sp1<<ionmass<<"       0       0       0"<<std::endl;
      particleIndex++;
    }

       
    int ipart = 0;
    std::vector<starlightParticle>::const_iterator part = (event.getParticles())->begin();
    
    for (part = event.getParticles()->begin(); part != event.getParticles()->end(); part++, ipart++)
    {
      _fileStream <<particleIndex <<sp0<<1 <<sp0<<(*part).getPdgCode()<<"       "<<2<<"       "<<0<<"       "<<0<<"  " <<sp1<<(*part).GetPx()<<"  " <<sp1<<(*part).GetPy()<<"  " <<sp1<<(*part).GetPz()<<"  " <<sp1<<(*part).GetE()<<"   "<<sp1<<(*part).M()<<"       0       0       0"<<std::endl;      
      particleIndex++;
    }
    

    _fileStream<<"=============== Event finished ==============="<<std::endl;
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
