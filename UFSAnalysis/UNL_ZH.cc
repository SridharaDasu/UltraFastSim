#include "UNL_ZH.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <string>
#include <math.h>
#include <algorithm>
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"
#include "UltraFastSim.h"
#include "LinkDef.h"

// Utilities
bool sort_TLorentzVector_by_pt(TLorentzVector i,TLorentzVector j) { return (i.Pt() > j.Pt()); }

// Constructor opens the file with the list of root files to chain
// and opens the output histogram file
UNL_ZH::UNL_ZH(TChain *t)
  : outFile_(0), chain_(0), nEvents_(0) 
{
  // Main object to read in
  ufs_ = new UltraFastSim();

  if (t) nEvents_ = t->GetEntries();
}
  
UNL_ZH::~UNL_ZH() {
  outFile_->Write();
  outFile_->Close();
  delete ufs_;
}

bool UNL_ZH::SetOutputFile(std::string outFile) {
  // Open output file with histograms
  outFile_ = new TFile(outFile.c_str(),"RECREATE");
  return true;
}


void UNL_ZH::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   chain_ = tree;

   nEvents_ = chain_->GetEntries();

   // Setup the branches to read in
   b_ufs_ = chain_->GetBranch("UltraFastSim");
   b_ufs_->SetAddress(&ufs_);
}

Bool_t UNL_ZH::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  
  return true;
}

Int_t UNL_ZH::GetEntry(Long64_t entry, Int_t getall) { 
  int localentry=-1;
  if (entry >=0 ) localentry = chain_->GetTree()->GetEntry(entry,getall);
  return localentry;
}
 
void UNL_ZH::Begin(TTree * t)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   outFile_->cd();
   outTree_ = new TTree("zh","Simple ZH analysis from UNL");
   outTree_->Branch("sel",&preselU_,"preselU/O:preselS:selU_tight:selS_tight:selU_med:selS_med:selU_loose:selS_loose");
   outTree_->Branch("btag",&nU_tight_,"nU_tight/I:nS_tight:nU_med:nS_med:nU_loose:nS_loose");
   outTree_->Branch("pt",&ptU_tight_,"ptU_tight/D:ptS_tight:ptU_med:ptS_med:ptU_loose:ptS_loose");
   outTree_->Branch("mH",&mHU_tight_,"mHU_tight/D:mHS_tight:mHU_med:mHS_med:mHU_loose:mHS_loose");
   outTree_->Branch("muon",&mu_nU_,"nU/I:nS/I");
   outTree_->Branch("muonpT", &muonptU_, "muonpTU[2]/D:muonpTS[2]");
   outTree_->Branch("ptZ", &ptZU_, "ptZU/D:ptZS");
   outTree_->Branch("mZ", &mZU_, "mZU/D:mZS");
   outTree_->Branch("dPhiHZ", &dPhiHZU_tight_, "dPhiHZU_tight/D:dPhiHZS_tight:dPhiHZU_med:dPhiHZS_med:dPhiHZU_loose:dPhiHZS_loose");
   outTree_->Branch("ptB",&ptBU_tight_,"ptBU_tight[2]/D:ptBS_tight[2]:ptBU_med[2]:ptBS_med[2]:ptBU_loose[2]:ptBS_loose[2]");
   outTree_->Branch("flavB",&flavU_tight_,"flavU_tight[2]/I:flavS_tight[2]:flavU_med[2]:flavS_med[2]:flavU_loose[2]:flavS_loose[2]");
   outTree_->SetWeight(1.0/nEvents_);
}

void UNL_ZH::ResetOutputTree() {
  preselU_=preselS_=selU_tight_=selS_tight_=selU_med_=selS_med_=selU_loose_=selS_loose_=false;
  nU_tight_=nS_tight_=nU_med_=nS_med_=nU_loose_=nS_loose_=-1;
  ptU_tight_=ptS_tight_=ptU_med_=ptS_med_=ptU_loose_=ptS_loose_=-1;
  mHU_tight_=mHS_tight_=mHU_med_=mHS_med_=mHU_loose_=mHS_loose_=-1;
  mu_nU_=mu_nS_=-1;
  for (int i=0; i!=2; i++){
    muonptU_[i] = -1.;
    muonptS_[i] = -1.;
  }
  ptZU_=ptZS_=mZU_=mZS_=-1.;
  dPhiHZU_tight_=dPhiHZS_tight_=dPhiHZU_med_=dPhiHZS_med_=dPhiHZU_loose_=dPhiHZS_loose_=-999.;
  for (int i=0; i<2; i++) {
    ptBU_tight_[i]=ptBS_tight_[i]=ptBU_med_[i]=ptBS_med_[i]=ptBU_loose_[i]=ptBS_loose_[i]=-1;
    flavU_tight_[i]=flavS_tight_[i]=flavU_med_[i]=flavS_med_[i]=flavU_loose_[i]=flavS_loose_[i]=0;
  }
}

void UNL_ZH::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void UNL_ZH::SortByDecreasingPt() {
  std::sort(muonsUpgrade_.begin(), muonsUpgrade_.end(), sort_TLorentzVector_by_pt);
  std::sort(muonsStdGeom_.begin(), muonsStdGeom_.end(), sort_TLorentzVector_by_pt);
  std::sort(btagsLooseUpgrade_.begin(), btagsLooseUpgrade_.end(), sort_TLorentzVector_by_pt);
  std::sort(btagsLooseStdGeom_.begin(), btagsLooseStdGeom_.end(), sort_TLorentzVector_by_pt);
  std::sort(btagsMediumUpgrade_.begin(), btagsMediumUpgrade_.end(), sort_TLorentzVector_by_pt);
  std::sort(btagsMediumStdGeom_.begin(), btagsMediumStdGeom_.end(), sort_TLorentzVector_by_pt);
  std::sort(btagsTightUpgrade_.begin(), btagsTightUpgrade_.end(), sort_TLorentzVector_by_pt);
  std::sort(btagsTightStdGeom_.begin(), btagsTightStdGeom_.end(), sort_TLorentzVector_by_pt);
  std::sort(bquarks_.begin(), bquarks_.end(), sort_TLorentzVector_by_pt);
  std::sort(cquarks_.begin(), cquarks_.end(), sort_TLorentzVector_by_pt);
}

Bool_t UNL_ZH::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either UNL_ZH::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  int nb = GetEntry(entry);

  // Get the TLorentzVector from the TParticle for muons
  muonsUpgrade_.clear(); muonsStdGeom_.clear();
  for (int i=0; i<ufs_->muonList().size(); ++i) {
    TLorentzVector m;
    ufs_->muonList()[i].Momentum(m);
    muonsUpgrade_.push_back(m);
  }
  for (int i=0; i<ufs_->muonListStdGeom().size(); ++i) {
    TLorentzVector m;
    ufs_->muonListStdGeom()[i].Momentum(m);
    muonsStdGeom_.push_back(m);
  }

  btagsLooseUpgrade_ = ufs_->bJetListLoose();
  btagsLooseStdGeom_ = ufs_->bJetListLooseStdGeom();
  btagsMediumUpgrade_ = ufs_->bJetListMedium();
  btagsMediumStdGeom_ = ufs_->bJetListMediumStdGeom();
  btagsTightUpgrade_ = ufs_->bJetListTight();
  btagsTightStdGeom_ = ufs_->bJetListTightStdGeom();

  bquarks_.clear(); cquarks_.clear();
  for (int i=0; i<ufs_->bQuarkList().size(); ++i) {
    TLorentzVector b;
    ufs_->bQuarkList()[i].Momentum(b);
    bquarks_.push_back(b);
  }
  for (int i=0; i<ufs_->cQuarkList().size(); ++i) {
    TLorentzVector c;
    ufs_->cQuarkList()[i].Momentum(c);
    cquarks_.push_back(c);
  }

  SortByDecreasingPt();

  ResetOutputTree();

  if (! Preselection()) return true;
  
//  std::cout << "Event " << entry << " muonsUpgrade " << muonsUpgrade_[0].Pt() << " " << muonsUpgrade_[1].Pt()
//	    << " muonsStdGeom " << muonsStdGeom_[0].Pt() << " " << muonsStdGeom_[1].Pt()
//	    << " btagsMediumUpgrade " << btagsMediumUpgrade_.size()
//	    << " btagsMediumStdGeom " << btagsMediumStdGeom_.size()
//	    << std::endl;
//
  Selection();
  outTree_->Fill();

  return kTRUE;
}

void UNL_ZH::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void UNL_ZH::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}


bool UNL_ZH::Preselection() {
  if (muonsUpgrade_.size() > 1) {
    if (muonsUpgrade_[0].Pt() > 20.0 && muonsUpgrade_[1].Pt() > 20.0){
      TLorentzVector mumu = muonsUpgrade_[0] + muonsUpgrade_[1];
      muonptU_[0] = muonsUpgrade_[0].Pt();
      muonptU_[1] = muonsUpgrade_[1].Pt();
      ptZU_ = mumu.Pt();
      mZU_ = mumu.M();
      mu_nU_ = muonsUpgrade_.size();
      preselU_ = true;
    }
  }
  if (muonsStdGeom_.size() > 1) {
    if (muonsStdGeom_[0].Pt() > 20.0 && muonsStdGeom_[1].Pt() > 20.0){
      TLorentzVector mumu = muonsStdGeom_[0] + muonsStdGeom_[1];
      muonptS_[0] = muonsStdGeom_[0].Pt();
      muonptS_[1] = muonsStdGeom_[1].Pt();
      ptZS_ = mumu.Pt();
      mZS_ = mumu.M();
      mu_nS_ = muonsStdGeom_.size();
      preselS_ = true;
    }
  }
  
  return preselU_ || preselS_;
}

bool UNL_ZH::Selection() {
  bool selected=false;
  
  nU_tight_ = btagsTightUpgrade_.size();
  nS_tight_ = btagsTightStdGeom_.size();
  nU_med_ = btagsMediumUpgrade_.size();
  nS_med_ = btagsMediumStdGeom_.size();
  nU_loose_ = btagsLooseUpgrade_.size();
  nS_loose_ = btagsLooseStdGeom_.size();

  if (nU_tight_ > 1) {
    TLorentzVector bb = btagsTightUpgrade_[0] + btagsTightUpgrade_[1];
    TLorentzVector mumu = muonsUpgrade_[0] + muonsUpgrade_[1];
    ptU_tight_ = bb.Pt();
    mHU_tight_ = bb.M();
    dPhiHZU_tight_ = bb.DeltaPhi(mumu);
    ptBU_tight_[0]=btagsTightUpgrade_[0].Pt();  ptBU_tight_[1]=btagsTightUpgrade_[1].Pt(); 
    AssignFlavor(flavU_tight_, 2, bquarks_, cquarks_, btagsTightUpgrade_);
    if (btagsTightUpgrade_[0].Pt() > 20.0 && btagsTightUpgrade_[1].Pt() > 20.0 ) {
      if (bb.Pt() > 100 ) {
     	selU_tight_ = true;
      }
    }
  }
  if (nS_tight_ > 1) {
    TLorentzVector bb = btagsTightStdGeom_[0] + btagsTightStdGeom_[1];
    TLorentzVector mumu = muonsStdGeom_[0] + muonsStdGeom_[1];
    ptS_tight_ = bb.Pt();
    mHS_tight_ = bb.M();
    dPhiHZS_tight_ = bb.DeltaPhi(mumu);
    ptBS_tight_[0]=btagsTightStdGeom_[0].Pt();  ptBS_tight_[1]=btagsTightStdGeom_[1].Pt(); 
    AssignFlavor(flavS_tight_, 2, bquarks_, cquarks_, btagsTightStdGeom_);
    if (btagsTightStdGeom_[0].Pt() > 20.0 && btagsTightStdGeom_[1].Pt() > 20.0 ) {
      if (bb.Pt() > 100 ) {
	selS_tight_ = true;
      }
    }
  }
  if (nU_med_ > 1) {
    TLorentzVector bb = btagsMediumUpgrade_[0] + btagsMediumUpgrade_[1];
    TLorentzVector mumu = muonsUpgrade_[0] + muonsUpgrade_[1];
    ptU_med_ = bb.Pt();
    mHU_med_ = bb.M();
    dPhiHZU_med_ = bb.DeltaPhi(mumu);
    ptBU_med_[0]=btagsMediumUpgrade_[0].Pt();  ptBU_med_[1]=btagsMediumUpgrade_[1].Pt(); 
    AssignFlavor(flavU_med_, 2, bquarks_, cquarks_, btagsMediumUpgrade_);
    if (btagsMediumUpgrade_[0].Pt() > 20.0 && btagsMediumUpgrade_[1].Pt() > 20.0 ) {
      if (bb.Pt() > 100 ) {
	selU_med_ = true;
      }
    }
  }
  if (nS_med_ > 1) {
    TLorentzVector bb = btagsMediumStdGeom_[0] + btagsMediumStdGeom_[1];
    TLorentzVector mumu = muonsStdGeom_[0] + muonsStdGeom_[1];
    ptS_med_ = bb.Pt();
    mHS_med_ = bb.M();
    dPhiHZS_med_ = bb.DeltaPhi(mumu);
    ptBS_med_[0]=btagsMediumStdGeom_[0].Pt();  ptBS_med_[1]=btagsMediumStdGeom_[1].Pt(); 
    AssignFlavor(flavS_med_, 2, bquarks_, cquarks_, btagsMediumStdGeom_);
    if (btagsMediumStdGeom_[0].Pt() > 20.0 && btagsMediumStdGeom_[1].Pt() > 20.0 ) {
      if (bb.Pt() > 100 ) {
	selS_med_ = true;
      }
    }
  }
  
  if (nU_loose_ > 1) {
    TLorentzVector bb = btagsLooseUpgrade_[0] + btagsLooseUpgrade_[1];
    TLorentzVector mumu = muonsUpgrade_[0] + muonsUpgrade_[1];
    ptU_loose_ = bb.Pt();
    mHU_loose_ = bb.M();
    dPhiHZU_loose_ = bb.DeltaPhi(mumu);
    ptBU_loose_[0]=btagsLooseUpgrade_[0].Pt();  ptBU_loose_[1]=btagsLooseUpgrade_[1].Pt(); 
    AssignFlavor(flavU_loose_, 2, bquarks_, cquarks_, btagsLooseUpgrade_);
    if (btagsLooseUpgrade_[0].Pt() > 20.0 && btagsLooseUpgrade_[1].Pt() > 20.0 ) {
      if (bb.Pt() > 100 ) {
	selU_loose_ = true;
      }
    }
  }
  if (nS_loose_ > 1) {
    TLorentzVector bb = btagsLooseStdGeom_[0] + btagsLooseStdGeom_[1];
    TLorentzVector mumu = muonsStdGeom_[0] + muonsStdGeom_[1];
    ptS_loose_ = bb.Pt();
    mHS_loose_ = bb.M();
    dPhiHZS_loose_ = bb.DeltaPhi(mumu);
    ptBS_loose_[0]=btagsLooseStdGeom_[0].Pt();  ptBS_loose_[1]=btagsLooseStdGeom_[1].Pt(); 
    AssignFlavor(flavS_loose_, 2, bquarks_, cquarks_, btagsLooseStdGeom_);
    if (btagsLooseStdGeom_[0].Pt() > 20.0 && btagsLooseStdGeom_[1].Pt() > 20.0 ) {
      if (bb.Pt() > 100 ) {
	selS_loose_ = true;
      }
    }
  }
  
  return selected;
}

void UNL_ZH::AssignFlavor(int *flav, int n, 
			  const std::vector<TLorentzVector>& b, const std::vector<TLorentzVector>& c,
			  const std::vector<TLorentzVector>& tag) {
  for (int i=0; i<n; ++i) {
    bool found=false;
    for (int ib=0; ib<b.size() && !found; ++ib) {
      double dR = b[ib].DeltaR(tag[i]);
      if (dR < 0.3) {
	flav[i] = 5;
	found = true;
      }
    }
    for (int ic=0; ic<c.size() && !found; ++ic) {
      double dR = c[ic].DeltaR(tag[i]);
      if (dR < 0.3) {
	flav[i] = 4;
	found = true;
      }
    }
    if (!found) flav[i]=1;
  }
}


			  
