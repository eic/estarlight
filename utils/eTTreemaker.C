// this very simple macro takes as input an Ascii eSTARlight output file, slight.out, and creates a standard TTree named ntuple.root
// S. Klein, June, 2018
// Fixed Bjorken-x calc Jan. 2020

#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

void eTTreemaker(void)
{
using namespace std;
  using namespace std;
  double me=0.000511;
  double mfinal;
// create output TTree here
  
  double fspx,fspy,fspz,fse,fsrap,fspt,fsmass,p1x,p1y,p1z,p1e,p1prap;
  double p2x,p2y,p2z,p2e,p2prap,p1pt,p2pt;

// new variables for eSTARlight
  double kphoton,kperp,qsquared,xbj,ttransfer;

// ttransfer will be imlemented later

  TNtuple *esNTuple = new TNtuple("esNT","eslightNtuple","fspt:fspz:fsrap:fsmass:p1pt:p1z:p1prap:p2pt:p2z:p2prap:kg:qsq:xbj");
// This is the electron, photon, final state (particle combination), particle 1 and particle 2
// These are pseudorapidities, except for the final state, where it is rapidity
// For the photon there is the photon energy, Q^2 the bjorken x of the target (calculated with x = Q^2/(2km_p) and t, the momentum transfer from the target

   cout <<" Opening slight.out"<< endl;

   ifstream inFile;
   inFile.open("slight.out");
   cout << "slight.out open"<<endl;
   int countLines=0;

  
 // Define variables for event loop
   int ntrk,nvtx;
   int i1,i2,i3,i4;  // junk variables for fscanf lines

	string line;
	stringstream lineStream;
  
 // event card
	string label;
	int eventNmb,nmbTracks, pcode;
	int nev,ntr,stopv,pdgpid1,pdgpid2;   

 //  go through input cards, which include the gammas of the electron and ion.
// We expect a CONFIG_OPT, BEAM_1 (electron), BEAM_2 (ion) and PHOTON
	double gamma_e,gamma_ion,junk1;


	//lineStream=stringstream(line);
	//assert(lineStream>>label);
	//assert(label == "VERTEX:");

	
	if (!getline(inFile,line)) {cout <<" Error reading CONFIG_OPT line"<<endl;}
	countLines++;
	lineStream=stringstream(line);
	assert(lineStream>>label);// CONFIG_OPT

	if (!getline(inFile,line)) {cout <<" Error reading BEAM_1 line"<<endl;}
	countLines++;
	lineStream=stringstream(line);
	assert(lineStream>>label>>i1>>i2>>gamma_e);  // BEAM_1 is the electron
	cout<<"Electron Lorentz boost is "<<gamma_e<<endl;
	
	if (!getline(inFile,line)) {cout <<" Error reading BEAM_2 line"<<endl;}
	countLines++;
	lineStream=stringstream(line);
       	assert(lineStream>>label>>i1>>i2>>gamma_ion);  // BEAM_2 is the ion
	cout<<"Ion Lorentz boost is "<<gamma_ion<<endl;

	if (!getline(inFile,line)) {cout <<" Error reading PHOTON line"<<endl;}
	countLines++;
	lineStream=stringstream(line); // PHOTON
	
  
// begin event loop here.  eSTARlight events have the format:
// Event card, Vertex card, photon card, source card, track 1 card, track 2 card

   while (inFile.good())
      {
 
	if (!getline(inFile,line)) {break;}
	countLines++;
	lineStream=stringstream(line);
	assert(lineStream>>label>>eventNmb>>nmbTracks);
	if (!(label =="EVENT:")) continue;
	if (eventNmb < 5) {cout <<"Reached Event "<<eventNmb<<endl;}
	
 // vertex card should be followed by gamma card, t card and source card.
	if (!getline(inFile,line)) {break;}
    	countLines++;
	//lineStream.str(line);
	//lineStream.clear();
	lineStream=stringstream(line);
	assert(lineStream>>label);
	assert(label == "VERTEX:");
	
	if (!getline(inFile,line)) {break;}
	lineStream=stringstream(line);
	assert(lineStream>>label>>kphoton>>qsquared);
	assert(label == "GAMMA:");

	if (!getline(inFile,line)) {break;}
	lineStream=stringstream(line);
	assert(lineStream>>label>>ttransfer);
	assert(label == "t:");

	if (!getline(inFile,line)) {break;}	
	lineStream=stringstream(line);
	assert(lineStream>>label);
	assert(label == "TARGET:");

       	if (!getline(inFile,line)) {break;}
	lineStream=stringstream(line);
	assert(lineStream>>label);
	assert(label == "SOURCE:");
	
// two track cards
	if (!getline(inFile,line)) {break;}
    	countLines++;
	//lineStream.str(line);
	lineStream=stringstream(line);
	assert(lineStream>>label>>pcode>>p1x>>p1y>>p1z>>i1>>i2>>i3>>pdgpid1);
	assert(label == "TRACK:");
	
	
	if (!getline(inFile,line)) {break;}
    	countLines++;
	lineStream=stringstream(line);
	//lineStream.str(line);
	assert(lineStream>>label>>pcode>>p2x>>p2y>>p2z>>i1>>i2>>i3>>pdgpid2);
	assert(label == "TRACK:");

	// get the final state masses  should be particle anti-particle, so pdgcodes should be opposite
	if (pdgpid1 != -pdgpid2)
	  { cout<<"Error pdgpid codes don't match"<<pdgpid1<<" "<<pdgpid2<<endl;
	    exit(-1);
	  }
	double pdgabs=abs(pdgpid1);
	mfinal=0.;
	if (pdgabs == 211) mfinal=0.139;
	if (pdgabs == 11) mfinal=0.000511;
	if (pdgabs == 13) mfinal=0.105;
	if (pdgabs == 321) mfinal=0.494;
	if (mfinal ==0)
	  {cout<<" Error final mass=0.  pdgabs="<<pdgabs<<endl;
	    exit(-1);
	  }
// Now do needed kinematics computations

	double p1=sqrt(p1x*p1x+p1y*p1y+p1z*p1z);
	double e1=sqrt(p1*p1+mfinal*mfinal);
        p1prap=0.5*log((p1+p1z)/(p1-p1z));
	p1pt=sqrt(p1x*p1x+p1y*p1y);

        double p2=sqrt(p2x*p2x+p2y*p2y+p2z*p2z);
	double e2=sqrt(p2*p2+mfinal*mfinal);
        p2prap=0.5*log((p2+p2z)/(p2-p2z));
	p2pt=sqrt(p2x*p2x+p2y*p2y);

 // final state
        fspx=p1x+p2x;
        fspy=p1y+p2y;
        fspz=p1z+p2z;
	double fse=e1+e2;
	// need to determine energy from pdgpid
        fspt=sqrt(fspx*fspx+fspy*fspy);
        fsmass=sqrt(fse*fse-(fspx*fspx+fspy*fspy+fspz*fspz));
        fsrap= 0.5*log((fse+fspz)/(fse-fspz));
	// done - now book NTuple

	// to calculate X_bjorken, need photon energy in lab frame
	double kphotontargetframe=kphoton*gamma_ion;
	//	double xbjorken=qsquared/(2.*kphotontargetframe*gamma_ion);
	double xbjorken=fsmass/(2.*gamma_ion*0.939)*exp(-fsrap);
	
	esNTuple->Fill(fspt,fspz,fsrap,fsmass,p1pt,p1z,p1prap,p2pt,p2z,p2prap,kphoton,qsquared,xbjorken);

	if (eventNmb < 5)
	   {cout<<fspx<<fspy<<fspz<<fse<<fsrap<<fspt<<fsmass<<fsrap<<endl;
	     cout<<p1x<<p1y<<p1z<<p1prap<<" "<<p1<<"  pt1="<<p1pt<<endl;
	     cout<<p2x<<p2y<<p2z<<p2prap<<" "<<p2<<"  pt2="<<p2pt<<endl;
	   }
	
      } //Done with event loop.
		   
// write out NTuple

   TFile *NTfile = new TFile("ntuple.root","recreate");
 esNTuple->Write();
 NTfile->Close();
}
		   
    
