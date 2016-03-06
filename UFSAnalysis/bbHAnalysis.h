#ifndef bbHAnalysis_H
#define bbHAnalysis_H
#include "TH1D.h"

namespace Pythia8 {
  class Pythia;
}

class UltraFastSim;

class TFile;
class TTree;

class bbHAnalysis {

public:
  bbHAnalysis(TFile *outFile, Pythia8::Pythia *pythia, UltraFastSim *ufs, bool verbosity);

  virtual ~bbHAnalysis();

  bool run();
  bool end();

private:

  bbHAnalysis();

  TFile* outFile;

  //ml added 

  double bJetPtCut; 
  double muonE;
  double muonPx;
  double muonPy;
  double muonPz;
  double taumE;
  double taumPx;
  double taumPy;
  double taumPz;
  
  double hadE;
  double hadPx;
  double hadPy;
  double hadPz;
  double tauhE;
  double tauhPx;
  double tauhPy;
  double tauhPz;
 
  int nbJetCut;
  int nbJetCutStdGeom;
  int nbasicCut;
  int nbasicCutStdGeom; 
  
  
  double tau1E;
  double tau1Px;
  double tau1Py;
  double tau1Pz;
  double tau2E;
  double tau2Px;
  double tau2Py;
  double tau2Pz;
 
  double rTau1E;
  double rTau1Px;
  double rTau1Py;
  double rTau1Pz;
  double rTau2E;
  double rTau2Px;
  double rTau2Py;
  double rTau2Pz;
  
  double hadronPt;// for hadron loop

  double bJet1PtStdGeom;
  double bJet1EtaStdGeom;
  double bJet1PhiStdGeom;
  double bJet2PtStdGeom;
  double bJet2EtaStdGeom;
  double bJet2PhiStdGeom;

  Pythia8::Pythia *pythia;
  UltraFastSim *ufs;
  bool verbose_;
  int iEvent;

  TTree *tree;

  int nElecs;
  int nMuons;
  int nTaus;
  int ncQrks;
  int nbQrks;
  int nJets;
  int nBJets;
  int nBJetsStdGeom;
  int nRTaus;
  double muonPt;
  double muonEta;
  double muonPhi;
  double tauPt;
  double tauEta;
  double tauPhi;
  double taumPt;
  double taumEta;
  double taumPhi;
  double tauhPt;
  double tauhEta;
  double tauhPhi;
  double rTauPt;
  double rTauEta;
  double rTauPhi;
  double bQrk1Pt;
  double bQrk1Eta;
  double bQrk1Phi;
  double bQrk2Pt;
  double bQrk2Eta;
  double bQrk2Phi;
  double bJet1Pt;
  double bJet1Eta;
  double bJet1Phi;
  double bJet2Pt;
  double bJet2Eta;
  double bJet2Phi;
  double gVisMass;
  double rVisMass;
  double ET;
  double HT;
  double MET;
  double MHT;


};

#endif
