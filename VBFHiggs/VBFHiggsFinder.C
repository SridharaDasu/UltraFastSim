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

vector<TH1F*> het;
vector<TH1F*> hmet;
vector<TH1F*> hht;
vector<TH1F*> hmht;
vector<TH1F*> hnJets;
vector<TH1F*> hleadJetET;
vector<TH1F*> hnextJetET;
vector<TH1F*> hleadJetEta;
vector<TH1F*> hnextJetEta;
vector<TH1F*> hdeltaEta;
vector<TH1F*> hdiJetMass;
vector<TH1F*> hdiJetPt;
vector<TH1F*> hdiJetEta;
vector<TH1F*> hdiJetPhi;

vector<TLorentzVector> jets;

TLorentzVector diJet;
double ET;
TLorentzVector MET;
double HT;
TLorentzVector MHT;

TH1F* etHist;
TH1F* metHist;

void bookHists() {

  cutLevel.push_back("diJet");
  cutLevel.push_back("diJet+deltaEta");
  cutLevel.push_back("diJet+deltaEta+diJetMass");
  cutLevel.push_back("diJet+deltaEta+diJetMass+JetVeto");

  string hName;
  TH1F* h;

  etHist = new TH1F("ET", "ET", 40, 0., 800.);
  metHist = new TH1F("MET", "MET", 40, 0., 800.);

  for(unsigned int c = 0; c < cutLevel.size(); c++) {
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
    hName = "leadJetEta" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hleadJetEta.push_back(h);	
    hName = "nextJetEta" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hnextJetEta.push_back(h);	
    hName = "deltaEta" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, 0, 7.0);
    hdeltaEta.push_back(h);
    hName = "diJetPt-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
    hdiJetPt.push_back(h);
    hName = "diJetEta-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hdiJetEta.push_back(h);
    hName = "diJetPhi-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
    hdiJetPhi.push_back(h);	
    hName = "diJetMass-" + cutLevel[c];
    h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 2000.);
    hdiJetMass.push_back(h);
  }

}

void fillHists(int c) {
  double METValue = MET.Pt();
  double MHTValue = MHT.Pt();
  double leadJetPT = jets[0].Pt();
  double nextJetPT = jets[1].Pt();
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
  hleadJetET[c]->Fill(jets[0].Pt());
  hnextJetET[c]->Fill(jets[1].Pt());
  hleadJetEta[c]->Fill(jets[0].Eta());
  hnextJetEta[c]->Fill(jets[1].Eta());
  hdeltaEta[c]->Fill(fabs(jets[0].Eta()-jets[1].Eta()));
  hdiJetPt[c]->Fill(diJet.Pt());
  hdiJetEta[c]->Fill(diJet.Eta());
  hdiJetPhi[c]->Fill(diJet.Phi());
  hdiJetMass[c]->Fill(diJet.M());
}

int main(int argc, char **argv){

  int AnalyzedEvents=0;
  int DiJetPreSelection=0;
  int DiJetSelection=0;
  int DeltaEtaWindow=0;
  int DiJetMassSelection=0;
  int JetVetoSelection=0;

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
      jets=ufs->jetList();
      ET=ufs->getET();
      etHist->Fill(ET);
      MET=ufs->getMET();
      metHist->Fill(MET.Pt());
      HT=ufs->getHT();
      MHT=ufs->getMHT();
      if(jets.size() > 1 && MET.Pt() > 100.){
	DiJetPreSelection++;  // Assumed trigger Level Cuts
	if(jets[0].Pt() > 30. && jets[1].Pt() > 30.) {
	  DiJetSelection++;
	  fillHists(0);
	  double jetEtaMin = jets[0].Eta();
	  double jetEtaMax = jets[1].Eta();
	  if(jetEtaMin > jetEtaMax) {
	    jetEtaMin = jets[1].Eta();
	    jetEtaMax = jets[0].Eta();
	  }
	  double deltaEta = jetEtaMax - jetEtaMin;
	  diJet=TLorentzVector(jets[0].Px(), jets[0].Py(), jets[0].Pz(), jets[0].Energy()) +
	    TLorentzVector(jets[1].Px(), jets[1].Py(), jets[1].Pz(), jets[1].Energy());
	  double diJetInvmass=diJet.M();
	  double diJetPt=diJet.Pt();
	  // Take events with a ~real Z with significant PT
	  if(deltaEta > 4.) {
	    DeltaEtaWindow++;
	    fillHists(1);
	    if(diJetInvmass > 1000.) {
	      DiJetMassSelection++;
	      fillHists(2);
	      // Make sure that there is no jet in between the first two jets in eta
	      bool jetVeto = false;
	      for(unsigned int i = 2; i < jets.size(); i++) {
		if((jets[i].Pt() > 30) && (jets[i].Eta() > jetEtaMin && jets[i].Eta() < jetEtaMax)) {
		  jetVeto = true;
		  break;
		}
	      }
	      if(!jetVeto) {
		JetVetoSelection++;
		fillHists(3);
	      }//Jet Veto
	    }//diJet Mass Window
	  }//Delta Eta
	}//Tag jets 1 and 2 ET
      }//jets.size>1 && MET > 100
    }//for loop on events
    delete rootTree;
    delete ufs;
    delete file;
  }//for loop on files
  
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events            : "<<AnalyzedEvents<<endl;
  cout<<" DiJet preselection         : "<<DiJetPreSelection<<endl;
  cout<<" DiJet    selection         : "<<DiJetSelection<<endl;
  cout<<" DeltaEta Window            : "<<DeltaEtaWindow<<endl;
  cout<<" DiJet Invariant Mass Window: "<<DiJetMassSelection<<endl;
  cout<<" Jet Veto                   : "<<JetVetoSelection<<endl;
  cout<<" ======================================== " <<endl;

  fout.cd();
  fout.Write();
  fout.Close();

  return 0;
}
