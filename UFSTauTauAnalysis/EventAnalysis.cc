#include "EventAnalysis.h"

using namespace std;
#include <math.h>
#include <string>
#include <map>

#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"

#include "UltraFastSim.h"

EventAnalysis::EventAnalysis(string histFileName) {
  file = new TFile(histFileName.c_str(), "recreate");
  // Define histograms
  histograms["NVisibleParticles"] = new TH1F("NVisibleParticles", "Number of Visible Particles", 100, 0., 1000.);
  histograms["muTauMass"] = new TH1F("muTauMass", "Visible Invariant Mass of (mu,tau) (GeV)", 100, 0., 200.);
  histograms["muTauMass_lowMET"] = new TH1F("muTauMass_lowMET", "Visible Invariant Mass of (mu,tau) (GeV) for MET < 50 GeV", 100, 0., 200.);
  histograms["muTauMass_highPT"] = new TH1F("muTauMass_highPT", "Visible Invariant Mass of (mu,tau) (GeV) for PT_muTau > 100 GeV", 100, 0., 200.);
  histograms["muTauMass_allCuts"] = new TH1F("muMuMass_allCuts", "Visible Invariant Mass of (mu,tau) (GeV) after all cuts", 100, 0., 200.);
}

void EventAnalysis::processEvent(UltraFastSim *ufs) {
  
  // Visible particle count

  int nVisibleParticles = 
    ufs->photonList().size() + 
    ufs->electronList().size() + 
    ufs->muonList().size() + 
    ufs->chargedHadronList().size() + 
    ufs->neutralHadronList().size();
  histograms["NVisibleParticles"]->Fill(nVisibleParticles);

  // muTau invariant mass

  if(ufs->visTauList().size() > 0 && ufs->muonList().size() > 0) {
    if(ufs->visTauList()[0].Pt() > 30 && ufs->visTauList()[0].Eta() < 2.5 && 
       ufs->muonList()[0].Pt() > 20 && ufs->muonList()[0].Eta() < 2.5) {
      int charge = 
	(ufs->visTauList()[0].GetPdgCode() * ufs->muonList()[0].GetPdgCode()) /
	abs(ufs->visTauList()[0].GetPdgCode() * ufs->muonList()[0].GetPdgCode());
      if(charge < 0) {
	TLorentzVector muTau = 
	  TLorentzVector(ufs->muonList()[0].Px(), ufs->muonList()[0].Py(), ufs->muonList()[0].Pz(), ufs->muonList()[0].Energy()) +
	  TLorentzVector(ufs->visTauList()[0].Px(), ufs->visTauList()[0].Py(), ufs->visTauList()[0].Pz(), ufs->visTauList()[0].Energy());
	double muTauMass = muTau.M();
	histograms["muTauMass"]->Fill(muTauMass);
	if(muTau.Pt() > 100) {
	  histograms["muTauMass_highPT"]->Fill(muTauMass);
	}
	if(ufs->getMET().Pt() < 50.) {
	  histograms["muTauMass_lowMET"]->Fill(muTauMass);
	}
	if(muTau.Pt() > 100 && ufs->getMET().Pt() < 50.) {
	  histograms["muTauMass_allCuts"]->Fill(muTauMass);
	}
	double a = 
	  (ufs->getMET().Px()*ufs->visTauList()[0].Py()-ufs->getMET().Py()*ufs->visTauList()[0].Px())/
	  (ufs->muonList()[0].Px()*ufs->visTauList()[0].Py()-ufs->muonList().Py()*ufs->visTauList()[0].Px());
	double b = 
	  (ufs->getMET().Px()*ufs->muonList()[0].Py()-ufs->getMET().Py()*ufs->muList()[0].Px())/
	  (ufs->visTauList()[0].Px()*ufs->muonList()[0].Py()-ufs->visTauList().Py()*ufs->muonList()[0].Px());
	if(a >= 0 && b >= 0) {
	}
	else if(a < 0 && b < 0) a = -a; b = -b;
	else {
	  cout << "collinear approximation failed; a = " << a << " b = " << b << endl;
	  a = 0; b = 0;
	}
	TLorentzVector tauTau =
      }
    }
  }

}

void EventAnalysis::finalize() {
  file->cd();
  file->Write();
  file->Close();
}
