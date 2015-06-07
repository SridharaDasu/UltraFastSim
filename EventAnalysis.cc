#include "EventAnalysis.h"

using namespace std;
#include <math.h>
#include <string>
#include <map>
#include <iostream>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"

#include "UltraFastSim.h"

EventAnalysis::EventAnalysis(string histFileName) {
  file = new TFile(histFileName.c_str(), "recreate");
  // Define histograms
  histograms["NVisibleParticles"] = new TH1F("NVisibleParticles", "Number of Visible Particles", 50, 0., 1000.);
  histograms["RecoPhotonET"] = new TH1F("RecoPhotonET", "Transverse energy of reconstructed photons (GeV)", 100, 0., 100.);
  histograms["MuPairMass@PreSelection"] = new TH1F("MuPairMass@PreSelection", "MuPair Mass @ PreSelection (GeV)", 100, 0., 200.);
  histograms["MuPairMass@FullSelection"] = new TH1F("MuPairMass@FullSelection", "MuPair Mass @ FullSelection (GeV)", 100, 0., 200.);
}

void EventAnalysis::processEvent(UltraFastSim *ufs) {
  
  // Process the generator level event

  int nVisibleParticles = 
    ufs->photonList().size() + 
    ufs->electronList().size() + 
    ufs->muonList().size() + 
    ufs->chargedHadronList().size() + 
    ufs->neutralHadronList().size();
  histograms["NVisibleParticles"]->Fill(nVisibleParticles);

  // Process the reconstructed event
  for(unsigned int i = 0; i < ufs->photonList().size(); i++) {
    double recoPhotonET = ufs->photonList()[i].Pt();
    histograms["RecoPhotonET"]->Fill(recoPhotonET);
  }

  vector<TParticle> muons;
  TLorentzVector MuPair;
  muons=ufs->muonList();
  if(muons.size()>1){
    MuPair=TLorentzVector(muons[0].Px(), muons[0].Py(), muons[0].Pz(), muons[0].Energy()) +
      TLorentzVector(muons[1].Px(), muons[1].Py(), muons[1].Pz(), muons[1].Energy());
    double MuPairMass=MuPair.M();
    double MuPairPt=MuPair.Pt();
    double MuPairPhi=MuPair.Phi();
    histograms["MuPairMass@PreSelection"]->Fill(MuPairMass);
    if(muons[0].Pt()>20. && muons[1].Pt()>20.){  
      histograms["MuPairMass@FullSelection"]->Fill(MuPairMass);
    }//muon pt>20 GeV
  }//muon.size>1
  
}

void EventAnalysis::finalize() {
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events          : "<<histograms["NVisibleParticles"]->GetEntries() << endl;
  cout<<" Lepton pair preselection : "<<histograms["MuPairMass@PreSelection"]->GetEntries() << endl;
  cout<<" Lepton pair    selection : "<<histograms["MuPairMass@FullSelection"]->GetEntries() << endl;
  file->cd();
  file->Write();
  file->Close();
}
