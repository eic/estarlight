#include "HepMC3/GenVertex.h"
#include "HepMC3/GenVertex_fwd.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "inputParameters.h"
#include "eXevent.h"
#include "hepmc3writer.h"

using HepMC3::FourVector;
using HepMC3::Print;

hepMC3Writer::hepMC3Writer()
{
}

int hepMC3Writer::initWriter(const inputParameters &param)
{

  std::string hepmc3_filename = param.baseFileName() + ".hepmc";
  
  /** need to pass the parameters to HepMC3 **/
  param.beam1Z(); //just do something for now 
  _hepmc3_output = new HepMC3::WriterAscii(hepmc3_filename);
  return 0;
}

int hepMC3Writer::writeEvent(const eXEvent &event, int eventnumber)
{

  /** Make HepMC3 Event and Vertex ... Currently only 1 Vertex per Event  **/ 
  HepMC3::GenEvent hepmc3_evt(HepMC3::Units::GEV , HepMC3::Units::MM);
  HepMC3::GenVertexPtr hepmc3_vertex_to_write = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));

  /** Get starlight particle vector **/
  const std::vector<starlightParticle> * particle_vector = event.getParticles();
  lorentzVector gamma_lorentzVec = (*event.getGamma())[0];
  FourVector hepmc3_gamma_four_vector = FourVector(gamma_lorentzVec.GetPx(),
						   gamma_lorentzVec.GetPy(),
						   gamma_lorentzVec.GetPz(),
						   gamma_lorentzVec.GetE());
  HepMC3::GenParticlePtr hepmc3_gamma = std::make_shared<HepMC3::GenParticle>( hepmc3_gamma_four_vector, 22, 0 ); // status currently set to 0... need to double check
  hepmc3_vertex_to_write->add_particle_in( hepmc3_gamma );
  
  /** Takes e_starlight events and converts to hepmc3 format **/
  for ( std::vector<starlightParticle>::const_iterator particle_iter = (*particle_vector).begin();
	particle_iter != (*particle_vector).end();
	++particle_iter)
    {
      int hepmc3_pid = (*particle_iter).getPdgCode();
      /** pass to HepMC3 FourVector **/
      FourVector hepmc3_four_vector = FourVector( (*particle_iter).GetPx(),
						  (*particle_iter).GetPy(),
						  (*particle_iter).GetPz(),
						  (*particle_iter).GetE());
      
      HepMC3::GenParticlePtr hepmc3_particle = std::make_shared<HepMC3::GenParticle>( hepmc3_four_vector, hepmc3_pid, 0 ); // status currently set to 0... need to double check
      hepmc3_vertex_to_write->add_particle_out( hepmc3_particle );
    }

  lorentzVector electron_lorentzVec = (*event.getSources())[0];
  FourVector hepmc3_electron_four_vector = FourVector(electron_lorentzVec.GetPx(),
						   electron_lorentzVec.GetPy(),
						   electron_lorentzVec.GetPz(),
						   electron_lorentzVec.GetE());

  HepMC3::GenParticlePtr hepmc3_electron = std::make_shared<HepMC3::GenParticle>( hepmc3_electron_four_vector, 11, 0 ); // status currently set to 0... need to double check

  // add code to output ion  SRK August 8 2021 
 lorentzVector ion_lorentzVec = (*event.getTarget())[0];
  FourVector hepmc3_ion_four_vector = FourVector(ion_lorentzVec.GetPx(),
						   ion_lorentzVec.GetPy(),
						   ion_lorentzVec.GetPz(),
						   ion_lorentzVec.GetE());
  
 //  PDG code for protons is 2212; not clear what to do for ions, but use proton ID for now  SRK August 8, 2021
HepMC3::GenParticlePtr hepmc3_ion = std::make_shared<HepMC3::GenParticle>( hepmc3_ion_four_vector, 2212, 0 ); // status currently set to 0... need to double check  

  hepmc3_vertex_to_write->add_particle_out( hepmc3_ion );
  
  /** add vertex to HepMC3 Event**/
  hepmc3_evt.add_vertex( hepmc3_vertex_to_write );
  _hepmc3_output->write_event(hepmc3_evt);
  //  Print::listing(hepmc3_evt);
  //  Print::content(hepmc3_evt);
  hepmc3_evt.clear();
  
  return eventnumber;
}

