//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 17 14:39:53 2014 by ROOT version 5.34/14
// from TTree starlightTree/starlightTree
// found on file: starlight.root
//////////////////////////////////////////////////////////

#ifndef AnalyzeTree_h
#define AnalyzeTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalyzeTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *parent;
   TClonesArray    *daughters;

   // List of branches
   TBranch        *b_parent;   //!
   TBranch        *b_daughters;   //!

   AnalyzeTree(TString fileName,TString treeName = "esNT");
   virtual ~AnalyzeTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //Data members to be filled
   float fspt_, fspz_, fsrap_, fsmass_, p1pt_, p1z_, p1prap_, p2pt_, p2z_, p2prap_, kg_, qsq_, xbj_;
   
};

#endif

#ifdef AnalyzeTree_cxx
AnalyzeTree::AnalyzeTree(TString fileName,TString treeName) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TFile * infile = new TFile(fileName,"read");
  TTree * tree = (TTree*) infile->Get( treeName );
  Init(tree);
}

AnalyzeTree::~AnalyzeTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalyzeTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   parent = 0;
   daughters = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fspt",&fspt_);
   fChain->SetBranchAddress("fsrap",&fsrap_);
   fChain->SetBranchAddress("fspz",&fspz_);
   fChain->SetBranchAddress("fsmass",&fsmass_);
   fChain->SetBranchAddress("p1pt",&p1pt_);
   fChain->SetBranchAddress("p1z",&p1z_);
   fChain->SetBranchAddress("p1prap",&p1prap_);
   fChain->SetBranchAddress("p2pt",&p2pt_);
   fChain->SetBranchAddress("p2z",&p2z_);
   fChain->SetBranchAddress("p2prap",&p2prap_);
   fChain->SetBranchAddress("kg",&kg_);
   fChain->SetBranchAddress("qsq",&qsq_);
   fChain->SetBranchAddress("xbj",&xbj_);
      
   Notify();
}

Bool_t AnalyzeTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyzeTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyzeTree_cxx
