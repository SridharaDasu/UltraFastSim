#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <math.h>
using namespace std;
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "EventData.h"
#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath

double deltaR(double eta1, double eta2, double phi1, double phi2) 
{
  double dphi = fabs(phi1-phi2);
  if(dphi > M_PI) 
    {
      dphi = (2*M_PI - dphi);
    }
  double deta = fabs(eta1-eta2);
  return sqrt(dphi*dphi + deta*deta);
}

vector<string> cutLevel;

vector<TH1F*> hMuonPt;  // For top two Muons
vector<TH1F*> hMuonPhi; // For top two Muons
vector<TH1F*> hMuonEta; // For top two Muons
vector<TH1F*> hElectronPt;  // For top two Electrons
vector<TH1F*> hElectronPhi; // For top two Electrons
vector<TH1F*> hElectronEta; // For top two Electrons
vector<TH1F*> hZPt;
vector<TH1F*> hZPhi;
vector<TH1F*> hZEta;
vector<TH1F*> hZMass;
vector<TH1F*> het;
vector<TH1F*> hmet;
vector<TH1F*> hht;
vector<TH1F*> hmht;
vector<TH1F*> hnJets;
vector<TH1F*> hleadJetET;
vector<TH1F*> hnextJetET;
vector<TH1F*> hSumPxRatio;
vector<TH1F*> hSumPyRatio;
vector<TH1F*> hMETRatio;

vector<TParticle> bQs;
vector<TParticle> cQs;
vector<TParticle> muons;
vector<TParticle> electrons;
vector<TLorentzVector> jets;

TLorentzVector Z;
double ET;
TLorentzVector MET;
double HT;
TLorentzVector MHT;

void bookHists() {

  cutLevel.push_back("Zmm");
  cutLevel.push_back("Zmm+ZPt");
  cutLevel.push_back("Zmm+ZPt+JetVeto");
  cutLevel.push_back("Zmm+ZPt+JetVeto+PxPyMETRatio");
  cutLevel.push_back("Zee");
  cutLevel.push_back("Zee+ZPt");
  cutLevel.push_back("Zee+ZPt+JetVeto+PxPyMETRatio");

  string hName;
  TH1F* h;
  for(unsigned int c = 0; c < cutLevel.size(); c++) {
    hName = "MuonPt-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hMuonPt.push_back(h);
    hName = "MuonEta-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hMuonEta.push_back(h);
    hName = "MuonPhi-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hMuonPhi.push_back(h);	
    hName = "ElectronPt-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hElectronPt.push_back(h);
    hName = "ElectronEta-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hElectronEta.push_back(h);
    hName = "ElectronPhi-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hElectronPhi.push_back(h);	
    hName = "ZPt-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hZPt.push_back(h);
    hName = "ZEta-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hZEta.push_back(h);
    hName = "ZPhi-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hZPhi.push_back(h);	
    hName = "ZMass-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 200.);
    hZMass.push_back(h);
    hName = "et" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 800.);
    het.push_back(h);	
    hName = "met" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hmet.push_back(h);	
    hName = "ht" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 800.);
    hht.push_back(h);	
    hName = "mht" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hmht.push_back(h);	
    hName = "nJets" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 40.);
    hnJets.push_back(h);
    hName = "leadJetET" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hleadJetET.push_back(h);	
    hName = "nextJetET" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hnextJetET.push_back(h);	
    hName = "SumPxRatio-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 3.);
    hSumPxRatio.push_back(h);
    hName = "SumPyRatio-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 3.);
    hSumPyRatio.push_back(h);
    hName = "METRatio-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 3.);
    hMETRatio.push_back(h);
  }

}

void fillHists(int c) {
  if(muons.size() > 0) {
    hMuonPt[c]->Fill(muons[0].Pt());
    hMuonEta[c]->Fill(muons[0].Eta());
    hMuonPhi[c]->Fill(muons[0].Phi());
    if(muons.size() > 1) {
      hMuonPt[c]->Fill(muons[1].Pt());
      hMuonEta[c]->Fill(muons[1].Eta());
      hMuonPhi[c]->Fill(muons[1].Phi());
    }
  }
  if(electrons.size() > 0) {
    hElectronPt[c]->Fill(electrons[0].Pt());
    hElectronPt[c]->Fill(electrons[1].Pt());
    hElectronEta[c]->Fill(electrons[0].Eta());
    if(electrons.size() > 1) {
      hElectronEta[c]->Fill(electrons[1].Eta());
      hElectronPhi[c]->Fill(electrons[0].Phi());
      hElectronPhi[c]->Fill(electrons[1].Phi());
    }
  }
  hZPt[c]->Fill(Z.Pt());
  hZEta[c]->Fill(Z.Eta());
  hZPhi[c]->Fill(Z.Phi());
  hZMass[c]->Fill(Z.M());
  het[c]->Fill(ET);
  hmet[c]->Fill(MET.Pt());
  hht[c]->Fill(HT);
  hmht[c]->Fill(MHT.Pt());
  unsigned int nJets = 0;
  for(nJets = 0; nJets < jets.size(); nJets++) {
    if(jets[nJets].Pt() < 30.) {
      break;
    }
  }
  hnJets[c]->Fill(nJets);
  if(jets.size() > 0) hleadJetET[c]->Fill(jets[0].Pt());
  if(jets.size() > 1) hnextJetET[c]->Fill(jets[1].Pt());
  hSumPxRatio[c]->Fill(-MET.Px()/Z.Px());
  hSumPyRatio[c]->Fill(-MET.Py()/Z.Py());
  hMETRatio[c]->Fill(MET.Pt()/Z.Pt());
}

int main(int argc, char **argv){

  int AnalyzedEvents=0;
  int MuonPairPreSelection=0;
  int MuonPairSelection=0;
  int ElectronPairPreSelection=0;
  int ElectronPairSelection=0;
  int ZMassWindow=0;
  int ZptSelection=0;
  int JetVetoSelection=0;
  int PxPyMETRatioSelection=0;

  //================================================
  //event loop

  ifstream ConfigFile;
  ConfigFile.open(argv[1]);
  if(!ConfigFile){cerr<<"Unable to read the config file";exit(1);}
  string OutputRootFile;
  ConfigFile>>OutputRootFile;

  //file should be created before tree to avoid memory resident trees
  TFile fout(OutputRootFile.c_str(),"recreate");
  bookHists();

  //Just to know the total statistics of the root files
  //This is needed for the event normalization
  string filename;
  while(ConfigFile>>filename && filename!="EOF"){
    TFile *file=TFile::Open(filename.c_str());
    TTree *rootTree = (TTree*)file->Get("EventData");
    int nev = int(rootTree->GetEntries());
    cerr<<"Reading : "<<filename<<endl;
    cerr<<"Number of entries is : "<<nev<<endl;
    TBranch *branch = rootTree->GetBranch("EventData");
    EventData *ufs=new EventData();
    branch->SetAddress(&ufs);
    for(int i = 0; i < nev; i++){
      AnalyzedEvents++;
      if(!(AnalyzedEvents%1000))cerr<<"event number "<<AnalyzedEvents<<endl;
      rootTree->GetEvent(i);
      bQs=ufs->bQuarkList();
      cQs=ufs->cQuarkList();
      muons=ufs->muonList();
      electrons=ufs->electronList();
      jets=ufs->jetList();
      ET=ufs->getET();
      MET=ufs->getMET();
      HT=ufs->getHT();
      MHT=ufs->getMHT();
      if(muons.size()>1){
	MuonPairPreSelection++;
	if(muons[0].Pt()>20. && muons[1].Pt()>20.){	  
	  MuonPairSelection++;	  
	  Z=TLorentzVector(muons[0].Px(), muons[0].Py(), muons[0].Pz(), muons[0].Energy()) +
	    TLorentzVector(muons[1].Px(), muons[1].Py(), muons[1].Pz(), muons[1].Energy());
	  double Zinvmass=Z.M();
	  double Zpt=Z.Pt();
	  fillHists(0);
	  // Take events with a ~real Z with significant PT
	  if(Zinvmass>70. && Zinvmass<110.){
	    ZMassWindow++;
	    if(Zpt>100.){
	      fillHists(1);
	      ZptSelection++;
	      // Make sure that the top two jets are soft -- consistent with UE or PU "Jet Veto"
	      if((jets.size() == 0) ||
		 (jets.size() > 0 && jets[0].Pt() < 80.) ||
		 (jets.size() > 1 && jets[1].Pt() < 50.)) {
		fillHists(2);
		JetVetoSelection++;
		// Make sure that MET is balancing the visible Z both in magnitude and direction
		if((-MET.Px()/Z.Px())>0.7 && (-MET.Py()/Z.Py())>0.7 && (MET.Pt()/Z.Pt())>0.7 &&
		   (-MET.Px()/Z.Px())<1.3 && (-MET.Py()/Z.Py())<1.3 && (MET.Pt()/Z.Pt())<1.3) {
		  fillHists(3);
		  PxPyMETRatioSelection++;
		}
	      }//Jet Veto
	    }//Z pt cut
	  }//Z mass window
	}//muon pt>20 GeV
      }//muon.size>1
      if(electrons.size()>1){
	ElectronPairPreSelection++;
	if(electrons[0].Pt()>20. && electrons[1].Pt()>20.){	  
	  ElectronPairSelection++;	  
	  Z=TLorentzVector(electrons[0].Px(), electrons[0].Py(), electrons[0].Pz(), electrons[0].Energy()) +
	    TLorentzVector(electrons[1].Px(), electrons[1].Py(), electrons[1].Pz(), electrons[1].Energy());
	  double Zinvmass=Z.M();
	  double Zpt=Z.Pt();
	  fillHists(4);
	  // Take events with a ~real Z with significant PT
	  if(Zinvmass>70. && Zinvmass<110.){
	    ZMassWindow++;
	    if(Zpt>100.){
	      fillHists(5);
	      ZptSelection++;
	      // Make sure that the top two jets are soft -- consistent with UE or PU "Jet Veto"
	      if((jets.size() == 0) ||
		 (jets.size() > 0 && jets[0].Pt() < 80.) ||
		 (jets.size() > 1 && jets[1].Pt() < 50.)) {
		fillHists(2);
		JetVetoSelection++;
		// Make sure that MET is balancing the visible Z both in magnitude and direction
		if((-MET.Px()/Z.Px())>0.7 && (-MET.Py()/Z.Py())>0.7 && (MET.Pt()/Z.Pt())>0.7 &&
		   (-MET.Px()/Z.Px())<1.3 && (-MET.Py()/Z.Py())<1.3 && (MET.Pt()/Z.Pt())<1.3) {
		  fillHists(3);
		  PxPyMETRatioSelection++;
		}
	      }//Jet Veto
	    }//Z pt cut
	  }//Z mass window
	}//muon pt>20 GeV
      }//muon.size>1
    }//for loop on events
    delete rootTree;
    delete ufs;
    delete file;
  }//for loop on files
  
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events            : "<<AnalyzedEvents<<endl;
  cout<<" Muon pair preselection     : "<<MuonPairPreSelection<<endl;
  cout<<" Muon pair    selection     : "<<MuonPairSelection<<endl;
  cout<<" Electron pair preselection : "<<ElectronPairPreSelection<<endl;
  cout<<" Electron pair    selection : "<<ElectronPairSelection<<endl;
  cout<<" Z Invariant Mass Window    : "<<ZMassWindow<<endl;
  cout<<" Z pt                       : "<<ZptSelection<<endl;
  cout<<" Jet Veto                   : "<<JetVetoSelection<<endl;
  cout<<" PxPyMETRatio selection     : "<<PxPyMETRatioSelection<<endl;
  cout<<" ======================================== " <<endl;

  fout.cd();
  fout.Write();
  fout.Close();

  return 0;
}
