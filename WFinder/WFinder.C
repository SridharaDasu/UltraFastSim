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

vector<string> cutLevel;

vector<TH1F*> hMuonPt; 
vector<TH1F*> hMuonPhi;
vector<TH1F*> hMuonEta;
vector<TH1F*> hElectronPt; 
vector<TH1F*> hElectronPhi;
vector<TH1F*> hElectronEta;
vector<TH1F*> hmet;
vector<TH1F*> hcosdeltaphi;
vector<TH1F*> hWPt;
vector<TH1F*> hWPhi;
vector<TH1F*> hWMass;

vector<TParticle> bQs;
vector<TParticle> cQs;
vector<TParticle> muons;
vector<TParticle> electrons;
vector<TLorentzVector> jets;

TLorentzVector W;
TLorentzVector MET;

double cosdeltaphi;
double mT;

void bookHists() {

  cutLevel.push_back("-m");
  cutLevel.push_back("-Wmn");
  cutLevel.push_back("-e");
  cutLevel.push_back("-Wen");

  string hName;
  TH1F* h;
  for(unsigned int c = 0; c < cutLevel.size(); c++) {
    hName = "MuonPt" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hMuonPt.push_back(h);
    hName = "MuonEta" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hMuonEta.push_back(h);
    hName = "MuonPhi" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, 0, 7.0);
    hMuonPhi.push_back(h);
    hName = "ElectronPt" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hElectronPt.push_back(h);
    hName = "ElectronEta" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hElectronEta.push_back(h);
    hName = "ElectronPhi" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, 0, 7.0);
    hElectronPhi.push_back(h);
    hName = "met" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 200.);
    hmet.push_back(h);
    hName = "cosdeltaphi" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, -1.0, 1.0);
    hcosdeltaphi.push_back(h);
    hName = "WMass" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 200.);
    hWMass.push_back(h);
    hName = "WPt" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hWPt.push_back(h);
    hName = "WPhi" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hWPhi.push_back(h);	
  }

}

void fillHists(int c) {
  if(muons.size() > 0) {
    hMuonPt[c]->Fill(muons[0].Pt());
    hMuonEta[c]->Fill(muons[0].Eta());
    hMuonPhi[c]->Fill(muons[0].Phi());
  }
  if(electrons.size() > 0) {
    hElectronPt[c]->Fill(electrons[0].Pt());
    hElectronEta[c]->Fill(electrons[0].Eta());
    hElectronPhi[c]->Fill(electrons[0].Phi());
  }
  hmet[c]->Fill(MET.Pt());
  hcosdeltaphi[c]->Fill(cosdeltaphi);
  hWPt[c]->Fill(W.Pt());
  hWPhi[c]->Fill(W.Phi());
  hWMass[c]->Fill(mT);
}

int main(int argc, char **argv){

  int AnalyzedEvents=0;
  int MuonPreSelection=0;
  int MuonSelection=0;
  int ElectronPreSelection=0;
  int ElectronSelection=0;
  int WMassWindow=0;

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
      MET=ufs->getMET();
      if(muons.size()>=1){
	MuonPreSelection++;
	TLorentzVector MuonMomentum;
	muons[0].Momentum(MuonMomentum);
	W=MET+MuonMomentum;
	if(muons[0].Pt()>20.){	  
	  MuonSelection++;
	  cosdeltaphi= (MET.Px()*muons[0].Px()+MET.Py()*muons[0].Py())/(MET.Pt()*muons[0].Pt()); // dot product
	  mT=sqrt(2*muons[0].Pt()*MET.Pt()*(1-cosdeltaphi)); // calculate transverse mass
	  double Winvmass= mT;
	  fillHists(0);
	  if(Winvmass>30.){
	    fillHists(1);
	    WMassWindow++;
	  }//W mass window
	}//muon pt>20 GeV
      }//muons.size>1
      if(electrons.size()>=1){
	ElectronPreSelection++;
	TLorentzVector ElectronMomentum;
	electrons[0].Momentum(ElectronMomentum);
	W=MET+ElectronMomentum;
	if(electrons[0].Pt()>20.){	  
	  ElectronSelection++;
	  cosdeltaphi= (MET.Px()*electrons[0].Px()+MET.Py()*electrons[0].Py())/(MET.Pt()*electrons[0].Pt()); // dot product
	  mT=sqrt(2*electrons[0].Pt()*MET.Pt()*(1-cosdeltaphi)); // calculate transverse mass
	  double Winvmass= mT;
	  fillHists(2);
	  if(Winvmass>30.){
	    fillHists(3);
	    WMassWindow++;
	    // Make sure there are two jets above 30 GeV
	    // We take the top two jets only, which results in requiring that other event activity is lower
	  }//W mass window
	}//electron pt>20 GeV
      }//electrons.size>1
    }
    delete ufs;
    delete rootTree;
    delete file;
  }//for loop on events
  
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events            : "<<AnalyzedEvents<<endl;
  cout<<" Muon preselection          : "<<MuonPreSelection<<endl;
  cout<<" Muon selection             : "<<MuonSelection<<endl;
  cout<<" Electron preselection      : "<<ElectronPreSelection<<endl;
  cout<<" Electron selection         : "<<ElectronSelection<<endl;
  cout<<" W Invariant Mass Window    : "<<WMassWindow<<endl;
  cout<<" ======================================== " <<endl;

  fout.cd();
  fout.Write();
  fout.Close();

  return 0;
}
