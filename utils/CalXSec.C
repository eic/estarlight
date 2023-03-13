#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

#include <iostream>
#include <fstream>
#include <cmath>


//This function returns the photon flux within Q2min<Q2<Q2max and k_low<photon energy<k_high, (unit is GeV), photon energy is in proton rest framework 
//The differential photon flux is written as the function of y and Q2, see reference P. Fleischmann, PhD thesis, DESY-THESIS-2004-013
//s is the center of mass energy for proton and electron
//MJpsi is the vector meson mass (default set is J/Psi), also unit is GeV/c^2
double calPhotonFluxY(double Q2min, double Q2max, double k_low, double k_high, double s, double MJpsi=3.0969)
{
  double me = 0.000511; // electron mass
  double m_I = 0.93827; // proton mass

  double dndE = 0;
  double alpha = 1./137.036, pi = TMath::Pi();

  double minQ2 = Q2min;
  double maxQ2 = Q2max;
  //Q2 and y nsteps
  int const nQ2 = 100;
  int const ny = 100;
  double lnQ2ratio = std::log(maxQ2/minQ2)/(1.*nQ2);
  double lnQ2_min = std::log(minQ2);

  for (int iQ2=0;iQ2<nQ2;iQ2++)
  {
    double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
    double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
    double mQ2 = (q2_2+q2_1)/2.; //mean
    double dQ2 = q2_2-q2_1;

    double yH = (k_high*k_high + mQ2-m_I)/(s-m_I);
    double yL = (k_low*k_low+ mQ2-m_I)/(s-m_I);
    double lnyratio = std::log(yH/yL)/(1.*ny); //in log 
    double lny_min = std::log(yL);

    for (int iy=0;iy<ny;iy++){
      double y2_1 = exp( lny_min + iy*lnyratio);	    
      double y2_2 = exp( lny_min + (iy+1)*lnyratio);
      double dy = y2_2-y2_1;
      double my = (y2_2+y2_1)/2.; //mean
      //kenimatic limit
      double min_Q2 = me*me*my*my/(1-my);
      if (mQ2<min_Q2) {continue; cout<<"continue" <<endl;}

      dndE+=dQ2*dy*( alpha/2./ pi/mQ2 ) * ( (1+(1-my)*(1-my))/my - 2*me*me*my/mQ2 );

      if((iQ2%10)==0&&(iy==99)) {cout << "min_Q2: " << min_Q2 << endl;
      	cout << "yL, lny_min, y2_1: "<<yL<<", " <<lny_min << ", "<< y2_1 << endl;
      	cout <<"dQ2, dy, mQ2, my, dndE: "<< dQ2 <<", " << dy <<", "<<mQ2 << ", "<< my << ", " <<dndE << endl;}
    }

  }
  return dndE;
}




// //This function returns the photon flux within Q2min<Q2<Q2max and k_low<photon energy<k_high, (unit is GeV), photon energy is in proton rest framework 
// //The differential photon flux is written as the function of photon energy k and Q2, see reference arXiv:1803.06420 
// //Ee is the electron beam energy at proton rest framework
// double calPhotonFluxEg(double Q2min, double Q2max, double k_low, double k_high, double Ee)
// {
//   double me = 0.000511;  //electron mass

//   double dndE = 0;
//   double alpha = 1./137.036, pi = TMath::Pi();
  
//   //Q2 and photon energy nsteps
//   int const nk=100;
//   int const nQ2 = 100;
//   double kH = k_high;
//   double kL = k_low;
//   double lnkratio = std::log(kH/kL)/(1.*nk); //in log 
//   double lnk_min = std::log(kL);
//   // double test=0;
//   for (int ik=0;ik<nk;ik++){
//     double k2_1 = exp( lnk_min + ik*lnkratio);	    
//     double k2_2 = exp( lnk_min + (ik+1)*lnkratio);
//     double dk = k2_2-k2_1;
//     double mk = (k2_2+k2_1)/2.; //mean

//     double min_Q2 = me*me*mk*mk/(Ee*(Ee-mk));
//     double minQ2 = min_Q2>Q2min? min_Q2:Q2min;
//     double maxQ2 = Q2max; // note: the maximum Q2 at kenimatic limit is 4*Ee*(Ee-k)
//     double lnQ2ratio = std::log(maxQ2/minQ2)/(1.*nQ2);
//     double lnQ2_min = std::log(minQ2);

//     for (int iQ2=0;iQ2<nQ2;iQ2++)
//     {
//       double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
//       double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
//       double mQ2 = (q2_2+q2_1)/2.; //mean
//       double dQ2 = q2_2-q2_1;

//       // Integrating the effective photon flux
//       dndE+=dQ2*dk*( alpha/ pi/mk/mQ2 ) * ( 1 - mk/Ee + mk*mk/2./Ee/Ee - (1-mk/Ee)*fabs(min_Q2/mQ2));
//     }
//   }
//   return dndE;
// }


void CalXSec(){

	// ******** Start measuring time ******** //
	clock_t start, end, cpu_time;
	start = clock();

	// need a header file
  	setTDRStyle();

  	TFile *ntuple = TFile::Open("ntuple_slightJpsi_27p5x920_1000000_Ee_Q2_2to5.root", "read");
  	TNtuple *esNtuple = (TNtuple*)gDirectory -> Get("esNT");

	// ******** Define variables in the tree ******** //
    float_t p1pt;
    float_t p1eta;
    float_t p2pt;
    float_t p2eta;
    float_t Jpsipt;
    float_t Jpsiy;
    float_t W;
    float_t PhotonE;
    float_t Q2;
    float_t Ep;
    float_t Ee;
  

    double me=0.000511; //in GeV
    double mp=0.939;


	esNtuple -> SetBranchAddress("p1pt", &p1pt);
	esNtuple -> SetBranchAddress("p1prap", &p1eta);
	esNtuple -> SetBranchAddress("p2pt", &p2pt);
	esNtuple -> SetBranchAddress("p2prap", &p2eta);
	esNtuple -> SetBranchAddress("fspt", &Jpsipt);
	esNtuple -> SetBranchAddress("fsrap", &Jpsiy);
	esNtuple -> SetBranchAddress("W", &W);
	esNtuple -> SetBranchAddress("kg", &PhotonE);
	esNtuple -> SetBranchAddress("qsq", &Q2);
	esNtuple -> SetBranchAddress("Ep", &Ep);
	esNtuple -> SetBranchAddress("Ee", &Ee);


	// ******** Print out the total number of entries ******** //	
	Long64_t totEntries = esNtuple->GetEntries();
	cout << "Total Entries: " << totEntries << endl;

	TCanvas *WCan = new TCanvas("WCan", "WCan", 1500, 1000);
	TH1D *WHist = new TH1D("WHist", "Photon Energy(W);W(GeV);counts", 20, 0, 200);
	esNtuple -> Draw("W>>WHist");
	int Nbins = WHist -> GetNbinsX();

	TH1D *photonFluxHist = new TH1D("photonFluxHist", "Photon Flux;W(GeV);#Phi_{Photon}(1/GeV)", 20, 0, 200);


	TH1D *XSecHist2to5 = new TH1D("XSecHist2to5", "Cross section;W(GeV);#sigma#gamma p #rightarrow J/#psi(nb)", 20, 0, 200);
	TH1D *XSecHist5to10 = new TH1D("XSecHist5to10", "Cross section;W(GeV);#sigma#gamma p to J/#psi(nb)", 20, 0, 200);

	//crosssection from estarlight when 10x100GeV
	// double pXsec2to5 = 6.072; //pb
	// double pXsec5to10 = 4.913; //pb

	

	// //crosssection from estarlight when 18x275GeV
	// double pXsec2to5 = 12.540; //pb
	// double pXsec5to10 = 4.913; //pb

	// //crosssection from estarlight when 27.5x920GeV
	double pXsec2to5 = 23.788*pow(10.,-3); //nb



	int content = 0, count;
	double Q2min=0, Q2max=0, k_low=0, k_high=0;
	// ******** Start the event loop - read onia tree and save values for bit 14 in Ntuples ******** //
	for(int iBin = 1; iBin<=5/*Nbins*/ ; iBin++){
		count = 0;
		if((WHist->GetBinContent(iBin))==0) continue;

		for(int iEvent=0; iEvent<(totEntries); iEvent++){

			// if(iEvent > 5){
			// 	break;
			// }

			esNtuple -> GetEntry(iEvent);
			// cout << "Q2: " << Q2 << endl;
			// cout << "PhotonE: " << PhotonE << endl;

			if(W > (WHist->GetBinLowEdge(iBin)) && W < (WHist->GetBinLowEdge(iBin+1)) ){

				if(count == 0){
					Q2min = Q2;
					Q2max = Q2;
					k_low = PhotonE;
					k_high = PhotonE;
				}

				if(Q2 < Q2min) Q2min=Q2;
				else if(Q2 > Q2max) Q2max=Q2;

				if(PhotonE < k_low) k_low = PhotonE;
				else if(PhotonE > k_high) k_high = PhotonE;
				
				// cout <<"Q2 min, max: " << Q2min << ", " << Q2max << endl;
				// cout <<"k min, max: " << k_low << ", " << k_high << endl;
				// cout << "*****************" << endl;
				count++;
			}
		}
		

		double s = pow(Ep+Ee,2)-pow(sqrt(pow(Ep,2)-pow(mp,2))-sqrt(pow(Ee,2)-pow(me,2)),2);
		cout << "Ep: " << Ep << endl;
		cout << "Ee: " << Ee << endl;
		cout << "s: " <<s << endl;

		// double PhotonFlux = calPhotonFluxEg(Q2min, Q2max, k_low, k_high, Ee);
		double PhotonFlux = calPhotonFluxY(Q2min, Q2max, k_low, k_high, s, 3.0969);

		photonFluxHist -> SetBinContent(iBin, PhotonFlux);
		XSecHist2to5 -> SetBinContent(iBin, (WHist->GetBinContent(iBin))*1/PhotonFlux*pXsec2to5);
	
		cout << "Q2 min, max: " << Q2min << ", " << Q2max << endl;
		cout << "k min, max: " << k_low << ", " << k_high << endl;

		cout << "Photon Flux: " << PhotonFlux << endl;
		cout << "****************************" << endl;

	}

	TCanvas *PhotonFluxCan = new TCanvas("PhotonFluxCan", "PhotonFluxCan", 1500, 1000);
	PhotonFluxCan -> SetLogy();
	photonFluxHist -> Draw();


	TCanvas *XSecCan = new TCanvas("XSecCan", "XSecCan", 1500, 1000);
	XSecCan -> SetLogy();
	XSecHist2to5 -> Draw("p");
	XSecHist2to5 -> SetMarkerStyle(20);
	XSecHist2to5 -> SetMarkerColor(kRed);
	XSecHist2to5 -> SetMarkerSize(3);

	// XSecHist5to10 -> Draw("same");
  	
  	// calPhotonFluxEg()

	

  	// ******** End measuring time ******** //
	end = clock();
	cpu_time = (double)(end - start)/CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time/60. << "minutes" << endl;
	

	return;

}
