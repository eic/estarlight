
#ifndef HEPMC3WRITER_H
#define HEPMC3WRITER_H

#include "inputParameters.h"
#include "eXevent.h"
#include  "HepMC3/WriterAscii.h"

using HepMC3::WriterAscii;

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

  WriterAscii * _hepmc3_output;

};

#endif // HEPMC3WRITER_H
