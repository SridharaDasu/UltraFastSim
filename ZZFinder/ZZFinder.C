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
vector<TH1F*> hZ1Pt;
vector<TH1F*> hZ1Phi;
vector<TH1F*> hZ1Eta;
vector<TH1F*> hZ1Mass;
vector<TH1F*> hZ2Pt;
vector<TH1F*> hZ2Phi;
vector<TH1F*> hZ2Eta;
vector<TH1F*> hZ2Mass;
vector<TH1F*> hZZ4LMass;
vector<TH1F*> het;
vector<TH1F*> hmet;
vector<TH1F*> hht;
vector<TH1F*> hmht;
vector<TH1F*> hnJets;
vector<TH1F*> hleadJetET;
vector<TH1F*> hnextJetET;

vector<TParticle> genTaus;
vector<TParticle> muons;
vector<TParticle> electrons;
vector<TLorentzVector> jets;

TLorentzVector ZMM;
TLorentzVector ZEE;

TLorentzVector Z1;
TLorentzVector Z2;

TLorentzVector ZZ4L;

double ET;
TLorentzVector MET;
double HT;
TLorentzVector MHT;

void bookHists() {

  cutLevel.push_back("mmee");
  cutLevel.push_back("mmee+PTCuts");
  cutLevel.push_back("mmee+PTCuts+Z1");
  cutLevel.push_back("mmee+PTCuts+Z1+tau");
  cutLevel.push_back("mmmm");
  cutLevel.push_back("mmmm+PTCuts");
  cutLevel.push_back("mmmm+PTCuts+Z1");
  cutLevel.push_back("mmmm+PTCuts+Z1+tau");
  cutLevel.push_back("eeee");
  cutLevel.push_back("eeee+PTCuts");
  cutLevel.push_back("eeee+PTCuts+Z1");
  cutLevel.push_back("eeee+PTCuts+Z1+tau");

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

    hName = "Z1Pt-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hZ1Pt.push_back(h);
    hName = "Z1Eta-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hZ1Eta.push_back(h);
    hName = "Z1Phi-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hZ1Phi.push_back(h);	
    hName = "Z1Mass-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 200.);
    hZ1Mass.push_back(h);

    hName = "Z2Pt-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hZ2Pt.push_back(h);
    hName = "Z2Eta-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hZ2Eta.push_back(h);
    hName = "Z2Phi-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hZ2Phi.push_back(h);	
    hName = "Z2Mass-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 200.);
    hZ2Mass.push_back(h);

    hName = "ZZ4LMass-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 80, 0., 400.);
    hZZ4LMass.push_back(h);

    hName = "et-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 800.);
    het.push_back(h);	
    hName = "met-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hmet.push_back(h);	
    hName = "ht-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 800.);
    hht.push_back(h);	
    hName = "mht-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hmht.push_back(h);	
    hName = "nJets-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 40.);
    hnJets.push_back(h);
    hName = "leadJetET-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hleadJetET.push_back(h);	
    hName = "nextJetET-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 400.);
    hnextJetET.push_back(h);	
  }

}

void fillHists(int c) {

  for(unsigned int i = 0; i < muons.size(); i++) {
    hMuonPt[c]->Fill(muons[i].Pt());
    hMuonEta[c]->Fill(muons[i].Eta());
    hMuonPhi[c]->Fill(muons[i].Phi());
  }

  for(unsigned int i = 0; i < electrons.size(); i++) {
    hElectronPt[c]->Fill(electrons[i].Pt());
    hElectronEta[c]->Fill(electrons[i].Eta());
    hElectronPhi[c]->Fill(electrons[i].Phi());
  }

  hZ1Pt[c]->Fill(Z1.Pt());
  hZ1Eta[c]->Fill(Z1.Eta());
  hZ1Phi[c]->Fill(Z1.Phi());
  hZ1Mass[c]->Fill(Z1.M());

  hZ2Pt[c]->Fill(Z2.Pt());
  hZ2Eta[c]->Fill(Z2.Eta());
  hZ2Phi[c]->Fill(Z2.Phi());
  hZ2Mass[c]->Fill(Z2.M());

  hZZ4LMass[c]->Fill(ZZ4L.M());

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

}

int main(int argc, char **argv){

  int AnalyzedEvents=0;
  int MMMMPreSelection=0;
  int MMMMSelection=0;
  int MMEEPreSelection=0;
  int MMEESelection=0;
  int EEEEPreSelection=0;
  int EEEESelection=0;
  int Z1MassWindow=0;
  int WithGenTaus=0;

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
      genTaus=ufs->genTauList();
      muons=ufs->muonList();
      electrons=ufs->electronList();
      jets=ufs->jetList();
      ET=ufs->getET();
      MET=ufs->getMET();
      HT=ufs->getHT();
      MHT=ufs->getMHT();
      //MMEE
      if(muons.size()>1 && electrons.size()>1){
	MMEEPreSelection++;
	ZMM=TLorentzVector(muons[0].Px(), muons[0].Py(), muons[0].Pz(), muons[0].Energy()) +
	  TLorentzVector(muons[1].Px(), muons[1].Py(), muons[1].Pz(), muons[1].Energy());
	ZEE=TLorentzVector(electrons[0].Px(), electrons[0].Py(), electrons[0].Pz(), electrons[0].Energy()) +
	  TLorentzVector(electrons[1].Px(), electrons[1].Py(), electrons[1].Pz(), electrons[1].Energy());
	double ZMMinvmass=ZMM.M();
	double ZEEinvmass=ZEE.M();
	if(ZMMinvmass > ZEEinvmass) {
	  Z1 = ZMM;
	  Z2 = ZEE;
	}
	else {
	  Z1 = ZEE;
	  Z2 = ZMM;
	}
	ZZ4L=ZMM+ZEE;
	fillHists(0);
	if((muons[0].Pt()>20. && muons[1].Pt()>10. && electrons[0].Pt()>7.  && electrons[1].Pt()>7. ) ||
	   (muons[0].Pt()>5.  && muons[1].Pt()>5.  && electrons[0].Pt()>20. && electrons[1].Pt()>10.) ||
	   (muons[0].Pt()>20. && muons[1].Pt()>5.  && electrons[0].Pt()>10. && electrons[1].Pt()>7. ) ||
	   (muons[0].Pt()>10. && muons[1].Pt()>5.  && electrons[0].Pt()>20. && electrons[1].Pt()>7.)) {
	  MMEESelection++;
	  fillHists(1);
	  // Take events with a ~real Z
	  if((ZMMinvmass>40. && ZMMinvmass<110.) || (ZEEinvmass>40. && ZEEinvmass<110.)) {
	    Z1MassWindow++;
	    fillHists(2);
	    // Increment separately those with taus
	    if(genTaus.size() > 0) {
	      if(genTaus[0].Pt() > 20) {
		WithGenTaus++;
		fillHists(3);
	      }
	    }
	  }//Z mass window
	}//Lepton Pt selection
      }//Four lepton
      //MMMM
      if(muons.size()>3){
	bool found = false;
	for(int j = 0; j < 3; j++) {
	  for(int k = j + 1; k < 4; k++) {
	    if((muons[j].GetPdgCode() * muons[k].GetPdgCode()) < 0) {
	      ZMM=TLorentzVector(muons[j].Px(), muons[j].Py(), muons[j].Pz(), muons[j].Energy()) +
		TLorentzVector(muons[k].Px(), muons[k].Py(), muons[k].Pz(), muons[k].Energy());
	      if(!found || (fabs(ZMM.M()-91.) < fabs(Z1.M()-91.))) {
		found = true;
		Z1 = ZMM;
		int l, m;
		if(j == 0 && k == 1) {l = 2; m = 3;}
		else if(j == 0 && k == 2) {l = 1; m = 3;}
		else if(j == 0 && k == 3) {l = 1; m = 2;}
		else if(j == 1 && k == 2) {l = 0; m = 3;}
		else if(j == 1 && k == 3) {l = 0; m = 2;}
		else if(j == 2 && k == 3) {l = 0; m = 1;}
		Z2 = TLorentzVector(muons[l].Px(), muons[l].Py(), muons[l].Pz(), muons[l].Energy()) +
		  TLorentzVector(muons[m].Px(), muons[m].Py(), muons[m].Pz(), muons[m].Energy());
	      }
	    }
	  }
	  if(found) {
	    ZZ4L = Z1 + Z2;
	    MMMMPreSelection++;
	    fillHists(4);
	    if(muons[0].Pt()>20. && muons[1].Pt()>10. && muons[2].Pt()>5. && muons[3].Pt()>5.){
	      MMMMSelection++;
	      fillHists(5);
	      double Z1invmass=Z1.M();
	      if((Z1invmass>40. && Z1invmass<110.)) {
		Z1MassWindow++;
		fillHists(6);
		// Increment separately those with taus
		if(genTaus.size() > 0) {
		  if(genTaus[0].Pt() > 20) {
		    WithGenTaus++;
		    fillHists(7);
		  }
		}
	      }
	    }
	  }
	}
      }//Four Muon
      //EEEE
      if(electrons.size()>3){
	bool found = false;
	for(int j = 0; j < 3; j++) {
	  for(int k = j + 1; k < 4; k++) {
	    if((electrons[j].GetPdgCode() * electrons[k].GetPdgCode()) < 0) {
	      ZMM=TLorentzVector(electrons[j].Px(), electrons[j].Py(), electrons[j].Pz(), electrons[j].Energy()) +
		TLorentzVector(electrons[k].Px(), electrons[k].Py(), electrons[k].Pz(), electrons[k].Energy());
	      if(!found || (fabs(ZMM.M()-91.) < fabs(Z1.M()-91.))) {
		found = true;
		Z1 = ZMM;
		int l, m;
		if(j == 0 && k == 1) {l = 2; m = 3;}
		else if(j == 0 && k == 2) {l = 1; m = 3;}
		else if(j == 0 && k == 3) {l = 1; m = 2;}
		else if(j == 1 && k == 2) {l = 0; m = 3;}
		else if(j == 1 && k == 3) {l = 0; m = 2;}
		else if(j == 2 && k == 3) {l = 0; m = 1;}
		Z2 = TLorentzVector(electrons[l].Px(), electrons[l].Py(), electrons[l].Pz(), electrons[l].Energy()) +
		  TLorentzVector(electrons[m].Px(), electrons[m].Py(), electrons[m].Pz(), electrons[m].Energy());
	      }
	    }
	  }
	  if(found) {
	    ZZ4L = Z1 + Z2;
	    EEEEPreSelection++;
	    fillHists(8);
	    if(electrons[0].Pt()>20. && electrons[1].Pt()>10. && electrons[2].Pt()>5. && electrons[3].Pt()>5.){
	      EEEESelection++;
	      fillHists(9);
	      double Z1invmass=Z1.M();
	      if((Z1invmass>40. && Z1invmass<110.)) {
		Z1MassWindow++;
		fillHists(10);
		// Increment separately those with taus
		if(genTaus.size() > 0) {
		  if(genTaus[0].Pt() > 20) {
		    WithGenTaus++;
		    fillHists(11);
		  }
		}
	      }
	    }
	  }
	}
      }//Four Electron
    }//for loop on events
    delete rootTree;
    delete ufs;
    delete file;
  }//for loop on files
  
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events            : "<<AnalyzedEvents<<endl;
  cout<<" MMMM Preselection          : "<<MMMMPreSelection<<endl;
  cout<<" MMMM    selection          : "<<MMMMSelection<<endl;
  cout<<" MMEE Preselection          : "<<MMEEPreSelection<<endl;
  cout<<" MMEE    selection          : "<<MMEESelection<<endl;
  cout<<" EEEE Preselection          : "<<EEEEPreSelection<<endl;
  cout<<" EEEE    selection          : "<<EEEESelection<<endl;
  cout<<" Z Invariant Mass Window    : "<<Z1MassWindow<<endl;
  cout<<" With Generator Taus        : "<<WithGenTaus<<endl;
  cout<<" ======================================== " <<endl;

  fout.cd();
  fout.Write();
  fout.Close();

  return 0;
}
