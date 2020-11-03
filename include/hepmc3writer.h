
#ifndef HEPMC3WRITER_H
#define HEPMC3WRITER_H

#include "inputParameters.h"
#include "eXevent.h"
#include  "HepMC3/WriterAscii.h"
#include "HepMC3/FourVector.h"

using HepMC3::WriterAscii;
using HepMC3::FourVector;

class hepMC3Writer
{

 public:
  /** Default Constructor **/
  hepMC3Writer();
  ~hepMC3Writer();
  
  /** Init HepMC3 Writer with input parameters **/
  int initWriter(const inputParameters & param);
  
  /** Write an eX event to file  **/
  int writeEvent(const eXEvent &event, int eventnumber);

  /* closes file */
  void close(){ _hepmc3_output->close(); };
  
 private:

  /** Init HepMC3 Beam Four Vectors **/
  int initBeamHepMC3(const inputParameters &param);
  
  WriterAscii * _hepmc3_output;
  FourVector hepmc3_beam1_four_vector_;
  FourVector hepmc3_beam2_four_vector_;
  int beam1_pdg_id_;
  int beam2_pdg_id_;
  
};

#endif // HEPMC3WRITER_H
