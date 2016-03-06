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
#include "UltraFastSim.h"
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
vector<string> bJetOP;
vector<string> geom;

vector< vector< vector<TH1F*> > > hMuonPt;  // For top two Muons
vector< vector< vector<TH1F*> > > hMuonPhi; // For top two Muons
vector< vector< vector<TH1F*> > > hMuonEta; // For top two Muons
vector< vector< vector<TH1F*> > > hZPt;
vector< vector< vector<TH1F*> > > hZPhi;
vector< vector< vector<TH1F*> > > hZEta;
vector< vector< vector<TH1F*> > > hZMass;
vector< vector< vector<TH1F*> > > hJet1Pt;
vector< vector< vector<TH1F*> > > hJet2Pt;
vector< vector< vector<TH1F*> > > hJet3Pt;
vector< vector< vector<TH1F*> > > hJetPhi;  // For top two jets
vector< vector< vector<TH1F*> > > hJetEta;  // For top two jets
vector< vector< vector<TH1F*> > > hHPt;
vector< vector< vector<TH1F*> > > hHPhi;
vector< vector< vector<TH1F*> > > hHEta;
vector< vector< vector<TH1F*> > > hHMass;

vector<TParticle> bQs;
vector<TParticle> cQs;
vector<TParticle> muons;
vector<TLorentzVector> jets;
vector<TLorentzVector> bjets[2][3];
TLorentzVector Z;
TLorentzVector H;

void bookHists() {

  cutLevel.push_back("Zmm");
  cutLevel.push_back("Zmm+ZPt");
  cutLevel.push_back("Zmm+ZPt-2Jets-JPt");
  cutLevel.push_back("Zmm+ZPt-2Jets-JPt-bTags");
  cutLevel.push_back("Zmm+ZPt-2Jets-JPt-bTags-HPt");
  cutLevel.push_back("Zmm+ZPt-2Jets-JPt-bTags-HPt-DPhi");
  cutLevel.push_back("Zmm+ZPt-2Jets-JPt-bTags-HPt-DPhi-JetVeto");

  bJetOP.push_back("Loose");
  bJetOP.push_back("Medium");
  bJetOP.push_back("Tight");

  geom.push_back("StdGeom");
  geom.push_back("Phase-1");

  string hName;
  TH1F* h;
  for(unsigned int g = 0; g < geom.size(); g++) {
    vector< vector<TH1F*> > vOvMuonPt;
    hMuonPt.push_back(vOvMuonPt);
    vector< vector<TH1F*> > vOvMuonEta;
    hMuonEta.push_back(vOvMuonEta);
    vector< vector<TH1F*> > vOvMuonPhi;
    hMuonPhi.push_back(vOvMuonPhi);
    vector< vector<TH1F*> > vOvZPt;
    hZPt.push_back(vOvZPt);
    vector< vector<TH1F*> > vOvZEta;
    hZEta.push_back(vOvZEta);
    vector< vector<TH1F*> > vOvZPhi;
    hZPhi.push_back(vOvZPhi);
    vector< vector<TH1F*> > vOvZMass;
    hZMass.push_back(vOvZMass);
    vector< vector<TH1F*> > vOvJet1Pt;
    hJet1Pt.push_back(vOvJet1Pt);
    vector< vector<TH1F*> > vOvJet2Pt;
    hJet2Pt.push_back(vOvJet2Pt);
    vector< vector<TH1F*> > vOvJet3Pt;
    hJet3Pt.push_back(vOvJet3Pt);
    vector< vector<TH1F*> > vOvJetEta;
    hJetEta.push_back(vOvJetEta);
    vector< vector<TH1F*> > vOvJetPhi;
    hJetPhi.push_back(vOvJetPhi);
    vector< vector<TH1F*> > vOvHPt;
    hHPt.push_back(vOvHPt);
    vector< vector<TH1F*> > vOvHEta;
    hHEta.push_back(vOvHEta);
    vector< vector<TH1F*> > vOvHPhi;
    hHPhi.push_back(vOvHPhi);
    vector< vector<TH1F*> > vOvHMass;
    hHMass.push_back(vOvHMass);
    for(unsigned int b = 0; b < bJetOP.size(); b++) {
      vector<TH1F*> vMuonPt;
      hMuonPt[g].push_back(vMuonPt);
      vector<TH1F*> vMuonEta;
      hMuonEta[g].push_back(vMuonEta);
      vector<TH1F*> vMuonPhi;
      hMuonPhi[g].push_back(vMuonPhi);
      vector<TH1F*> vZPt;
      hZPt[g].push_back(vZPt);
      vector<TH1F*> vZEta;
      hZEta[g].push_back(vZEta);
      vector<TH1F*> vZPhi;
      hZPhi[g].push_back(vZPhi);
      vector<TH1F*> vZMass;
      hZMass[g].push_back(vZMass);
      vector<TH1F*> vJet1Pt;
      hJet1Pt[g].push_back(vJet1Pt);
      vector<TH1F*> vJet2Pt;
      hJet2Pt[g].push_back(vJet2Pt);
      vector<TH1F*> vJet3Pt;
      hJet3Pt[g].push_back(vJet3Pt);
      vector<TH1F*> vJetEta;
      hJetEta[g].push_back(vJetEta);
      vector<TH1F*> vJetPhi;
      hJetPhi[g].push_back(vJetPhi);
      vector<TH1F*> vHPt;
      hHPt[g].push_back(vHPt);
      vector<TH1F*> vHEta;
      hHEta[g].push_back(vHEta);
      vector<TH1F*> vHPhi;
      hHPhi[g].push_back(vHPhi);
      vector<TH1F*> vHMass;
      hHMass[g].push_back(vHMass);
      for(unsigned int c = 0; c < cutLevel.size(); c++) {
	hName = "MuonPt-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
	hMuonPt[g][b].push_back(h);
	hName = "MuonEta-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
	hMuonEta[g][b].push_back(h);
	hName = "MuonPhi-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 35, 0, 7.0);
	hMuonPhi[g][b].push_back(h);	
	hName = "ZPt-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
	hZPt[g][b].push_back(h);
	hName = "ZEta-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
	hZEta[g][b].push_back(h);
	hName = "ZPhi-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 35, 0, 7.0);
	hZPhi[g][b].push_back(h);	
	hName = "ZMass-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 40, 0., 200.);
	hZMass[g][b].push_back(h);	
	hName = "Jet1Pt-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
	hJet1Pt[g][b].push_back(h);
	hName = "Jet2Pt-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 30, 0, 300.);
	hJet2Pt[g][b].push_back(h);
	hName = "Jet3Pt-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 30, 0, 300.);
	hJet3Pt[g][b].push_back(h);
	hName = "JetEta-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
	hJetEta[g][b].push_back(h);
	hName = "JetPhi-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 35, 0., 7.0);
	hJetPhi[g][b].push_back(h);	
	hName = "HPt-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 30, 0., 300.);
	hHPt[g][b].push_back(h);
	hName = "HEta-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
        h = new TH1F(hName.c_str(), hName.c_str(), 35, -3.5, 3.5);
	hHEta[g][b].push_back(h);
	hName = "HPhi-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 35, 0, 7.0);
	hHPhi[g][b].push_back(h);	
	hName = "HMass-" + cutLevel[c] + "-" + bJetOP[b] + "-" + geom[g];
	h = new TH1F(hName.c_str(), hName.c_str(), 30, 0, 300.);
	hHMass[g][b].push_back(h);	
      }
    }
  }

}

void fillHists(int c, int b, int g) {
  hMuonPt[g][b][c]->Fill(muons[0].Pt());
  hMuonPt[g][b][c]->Fill(muons[1].Pt());
  hMuonEta[g][b][c]->Fill(muons[0].Eta());
  hMuonEta[g][b][c]->Fill(muons[1].Eta());
  hMuonPhi[g][b][c]->Fill(muons[0].Phi());
  hMuonPhi[g][b][c]->Fill(muons[1].Phi());
  hZPt[g][b][c]->Fill(Z.Pt());
  hZEta[g][b][c]->Fill(Z.Eta());
  hZPhi[g][b][c]->Fill(Z.Phi());
  hZMass[g][b][c]->Fill(Z.M());
  if(jets.size() > 0) {
    hJet1Pt[g][b][c]->Fill(jets[0].Pt());
    hJetEta[g][b][c]->Fill(jets[0].Eta());
    hJetPhi[g][b][c]->Fill(jets[0].Phi());
  }
  if(jets.size() > 1) {
    hJet2Pt[g][b][c]->Fill(jets[1].Pt());
    hJetEta[g][b][c]->Fill(jets[1].Eta());
    hJetPhi[g][b][c]->Fill(jets[1].Phi());
  }
  if(jets.size()>2) {
    hJet3Pt[g][b][c]->Fill(jets[2].Pt());
  }
  if(c > 1) {
    hHPt[g][b][c]->Fill(H.Pt());
    hHEta[g][b][c]->Fill(H.Eta());
    hHPhi[g][b][c]->Fill(H.Phi());
    hHMass[g][b][c]->Fill(H.M());
  }
}

void printDebugInfo(int g, int b) {
  for(unsigned int i = 0; i < bQs.size(); i++) {
    double dbq0 = deltaR(bQs[i].Eta(), bjets[g][b][0].Eta(), bQs[i].Phi(), bjets[g][b][0].Phi());
    if(dbq0 < 0.3) {
      cerr << "bQuark[" << i << "] = (" 
	   << bQs[i].Pt() << "," << bQs[i].Eta() << "," << bQs[i].Phi() << ")" 
	   << " matches bjet[" << g << "," << b << "][0] = ("
	   << bjets[g][b][0].Pt() << "," << bjets[g][b][0].Eta() << "," << bjets[g][b][0].Phi() << ")" << endl;
    }
    double dbq1 = deltaR(bQs[i].Eta(), bjets[g][b][1].Eta(), bQs[i].Phi(), bjets[g][b][1].Phi());
    if(dbq1 < 0.3) {
      cerr << "bQuark[" << i << "] = (" 
	   << bQs[i].Pt() << "," << bQs[i].Eta() << "," << bQs[i].Phi() << ")" 
	   << " matches bjet[" << g << "," << b << "][0] = ("
	   << bjets[g][b][1].Pt() << "," << bjets[g][b][1].Eta() << "," << bjets[g][b][1].Phi() << ")" << endl;
    }
  }
  for(unsigned int i = 0; i < cQs.size(); i++) {
    double dbq0 = deltaR(cQs[i].Eta(), bjets[g][b][0].Eta(), cQs[i].Phi(), bjets[g][b][0].Phi());
    if(dbq0 < 0.3) {
      cerr << "cQuark[" << i << "] = (" 
	   << cQs[i].Pt() << "," << cQs[i].Eta() << "," << cQs[i].Phi() << ")" 
	   << " matches bjet[" << g << "," << b << "][0] = ("
	   << bjets[g][b][0].Pt() << "," << bjets[g][b][0].Eta() << "," << bjets[g][b][0].Phi() << ")" << endl;
    }
    double dbq1 = deltaR(cQs[i].Eta(), bjets[g][b][1].Eta(), cQs[i].Phi(), bjets[g][b][1].Phi());
    if(dbq1 < 0.3) {
      cerr << "cQuark[" << i << "] = (" 
	   << cQs[i].Pt() << "," << cQs[i].Eta() << "," << cQs[i].Phi() << ")" 
	   << " matches bjet[" << g << "," << b << "][0] = ("
	   << bjets[g][b][1].Pt() << "," << bjets[g][b][1].Eta() << "," << bjets[g][b][1].Phi() << ")" << endl;
    }
  }
}

int main(){

  int AnalyzedEvents=0;
  int LeptonPairPreSelection=0;
  int LeptonPairSelection=0;
  int ZMassWindow=0;
  int ZptSelection=0;
  int JetPairPreSelection=0;
  int JetPairSelection=0;
  int BJetPairPreSelection[2][3] = {{0, 0, 0}, {0, 0, 0}};
  int BJetPairSelection[2][3] = {{0, 0, 0}, {0, 0, 0}};
  int HptSelection[2][3]= {{0, 0, 0}, {0, 0, 0}};
  int DphiSelection[2][3]={{0, 0, 0}, {0, 0, 0}};
  int JetVetoSelection[2][3]={{0, 0, 0}, {0, 0, 0}};
  int HMassWindow[2][3]={{0, 0, 0}, {0, 0, 0}};

  //================================================
  //event loop

  ifstream ConfigFile;
  ConfigFile.open("Config.txt");
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
    TTree *rootTree = (TTree*)file->Get("UltraFastSim");
    int nev = int(rootTree->GetEntries());
    cerr<<"Reading : "<<filename<<endl;
    cerr<<"Number of entries is : "<<nev<<endl;
    TBranch *branch = rootTree->GetBranch("UltraFastSim");
    UltraFastSim *ufs=new UltraFastSim();
    branch->SetAddress(&ufs);
    for(int i = 0; i < nev; i++){
      AnalyzedEvents++;
      if(!(AnalyzedEvents%1000))cerr<<"event number "<<AnalyzedEvents<<endl;
      rootTree->GetEvent(i);
      bQs=ufs->bQuarkList();
      cQs=ufs->cQuarkList();
      muons=ufs->muonList();
      jets=ufs->jetList();
      bjets[0][0]=ufs->bJetListLooseStdGeom();
      bjets[1][0]=ufs->bJetListLoose();
      bjets[0][1]=ufs->bJetListMediumStdGeom();
      bjets[1][1]=ufs->bJetListMedium();
      bjets[0][2]=ufs->bJetListTightStdGeom();
      bjets[1][2]=ufs->bJetListTight();
      if(muons.size()>1){
	LeptonPairPreSelection++;
	if(muons[0].Pt()>20. && muons[1].Pt()>20.){	  
	  LeptonPairSelection++;	  
	  Z=TLorentzVector(muons[0].Px(), muons[0].Py(), muons[0].Pz(), muons[0].Energy()) +
	    TLorentzVector(muons[1].Px(), muons[1].Py(), muons[1].Pz(), muons[1].Energy());
	  double Zinvmass=Z.M();
	  double Zpt=Z.Pt();
	  double Zphi=Z.Phi();
	  if(Zinvmass>70. && Zinvmass<110.){
	    for(int b=0;b<3;b++){
	      for(int g=0;g<2;g++) {
		fillHists(0, b, g);
	      }
	    }
	    ZMassWindow++;
	    // Make sure there are two jets above 30 GeV
	    // We take the top two jets only, which results in requiring that other event activity is lower
	    if(Zpt>100.){
	      for(int b=0;b<3;b++){
		for(int g=0;g<2;g++) {
		  fillHists(1, b, g);
		}
	      }
	      ZptSelection++;
	      if(jets.size()>1) {
		JetPairPreSelection++;
		if(jets[0].Pt()>50. && jets[1].Pt()>50.) {
		  H=jets[0]+jets[1];
		  for(int b=0;b<3;b++){
		    for(int g=0;g<2;g++) {
		      fillHists(2, b, g);
		    }
		  }
		  JetPairSelection++;		  
		  double Hinvmass=H.M();
		  double Hpt=H.Pt();
		  double Hphi=H.Phi();		  
		  double DphiZH=fabs(Hphi-Zphi);
		  if(DphiZH>M_PI)DphiZH = (2.*M_PI-DphiZH); 		  
		  //=======>>>>> loop over b-jet algo, and tk. geometry
		  for(int b=0;b<3;b++){
		    for(int g=0;g<2;g++) {
		      if(bjets[g][b].size()>1){
			BJetPairPreSelection[g][b]++;
			double djb0 = deltaR(jets[0].Eta(), bjets[g][b][0].Eta(), jets[0].Phi(), bjets[g][b][0].Phi());
			double djb1 = deltaR(jets[1].Eta(), bjets[g][b][1].Eta(), jets[1].Phi(), bjets[g][b][1].Phi());			
			// The top two bJets better be the top two jets
			if(djb0 < 0.01 && djb1 < 0.01) {
			  fillHists(3, b, g);
			  BJetPairSelection[g][b]++;
			  if(Hpt>120.){			    
			    fillHists(4, b, g);
			    HptSelection[g][b]++;			   
			    if(DphiZH>-1.) { // Disable 2.75 Dphi cut
			      fillHists(5, b, g);
			      DphiSelection[g][b]++;
			      if(jets[2].Pt()<5000.) {  // Disable jet veto
				fillHists(6, b, g);
				JetVetoSelection[g][b]++;
				float HinvmassMin = 105;
				float HinvmassMax = 175;
				if(OutputRootFile.find("v6") != string::npos) {
				  HinvmassMin = 90;
				  HinvmassMax = 140;
				}
				if(Hinvmass < HinvmassMax && Hinvmass > HinvmassMin) {
				  HMassWindow[g][b]++;
				  printDebugInfo(g, b);
				}//H mass window
			      }//Jet veto
			    }//H,Z dphi cut
			  }//H pt cut
			}//bjet - jet matching cut
		      }//bjet pair
		    }//loop over tk geometries
		  }//loop over b-tagging algo
		}//jet pt>30 GeV
	      }//jet.size>1
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
  cout<<" Analyzed Events          : "<<AnalyzedEvents<<endl;
  cout<<" Lepton pair preselection : "<<LeptonPairPreSelection<<endl;
  cout<<" Lepton pair    selection : "<<LeptonPairSelection<<endl;
  cout<<" Z Invariant Mass Window  : "<<ZMassWindow<<endl;
  cout<<" Z pt                     : "<<ZptSelection<<endl;
  cout<<" Jet pair preselection    : "<<JetPairPreSelection<<endl;
  cout<<" Jet pair    selection    : "<<JetPairSelection<<endl;
  for(int b=0;b<3;b++){
    for(int g=0;g<2;g++) {
      cout<<" ======================================== " <<endl;
      cout<<" Algo - "<<bJetOP[b]<<" b-tagging; "<<geom[g]<<" geometry"<<endl;
      cout<<" ======================================== " <<endl;
      cout<<" B-jet pair preselection  : "<<BJetPairPreSelection[g][b]<<endl;
      cout<<" B-jet pair    selection  : "<<BJetPairSelection[g][b]<<endl;
      cout<<" H pt                     : "<<HptSelection[g][b]<<endl;
      cout<<" Dphi(H,Z)                : "<<DphiSelection[g][b]<<endl;
      cout<<" Jet veto                 : "<<JetVetoSelection[g][b]<<endl;
      cout<<" H Invariant Mass Window  : "<<HMassWindow[g][b]<<endl;
      cout<<" ======================================== " <<endl;
    }
  }
  fout.cd();
  fout.Write();
  fout.Close();

  return 0;
}
