#ifndef _UNL_ZH_h_
#define _UNL_ZH_h_
//
// Simple ZH analysis that runs over Sridhara's UltraFastSim files
//
#include <string>
#include <TSelector.h>
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TChain.h"
#include "TFile.h"
class UltraFastSim;

class UNL_ZH : public TSelector {
 public:
  UNL_ZH(TChain *t=0);
  ~UNL_ZH();

   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   bool SetOutputFile(std::string outFile);
  
 private:
  // Preselection
  bool Preselection();
  // Basic analysis cuts
  bool Selection();
  void ResetOutputTree();
  void SortByDecreasingPt();
  void AssignFlavor(int *flav, int n, 
		    const std::vector<TLorentzVector>& b, const std::vector<TLorentzVector>& c,
		    const std::vector<TLorentzVector>& tag);

  // Data members
  TFile *outFile_;
  TTree *chain_;
  int nEvents_;
  std::vector<TLorentzVector> muonsUpgrade_, muonsStdGeom_;
  std::vector<TLorentzVector> btagsLooseUpgrade_, btagsLooseStdGeom_, btagsMediumUpgrade_, btagsMediumStdGeom_, btagsTightUpgrade_, btagsTightStdGeom_;
  std::vector<TLorentzVector> bquarks_, cquarks_;
  UltraFastSim *ufs_;
  TBranch *b_ufs_;

  // Output tree
  TTree *outTree_;
  bool preselU_, preselS_, selU_tight_, selS_tight_, selU_med_, selS_med_, selU_loose_, selS_loose_;
  int nU_tight_, nS_tight_, nU_med_, nS_med_, nU_loose_, nS_loose_;
  double ptU_tight_, ptS_tight_, ptU_med_, ptS_med_, ptU_loose_, ptS_loose_;
  double mHU_tight_, mHS_tight_, mHU_med_, mHS_med_, mHU_loose_, mHS_loose_;
  int mu_nU_, mu_nS_;
  double muonptU_[2], muonptS_[2];
  double ptZU_, ptZS_, mZU_, mZS_;
  double dPhiHZU_tight_, dPhiHZS_tight_, dPhiHZU_med_, dPhiHZS_med_, dPhiHZU_loose_, dPhiHZS_loose_;
  double ptBU_tight_[2], ptBS_tight_[2], ptBU_med_[2], ptBS_med_[2], ptBU_loose_[2], ptBS_loose_[2];
  int flavU_tight_[2], flavS_tight_[2], flavU_med_[2], flavS_med_[2], flavU_loose_[2], flavS_loose_[2];
  


  ClassDef(UNL_ZH,0);
};
#endif
