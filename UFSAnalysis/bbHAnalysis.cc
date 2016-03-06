
#include "bbHAnalysis.h"

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "Pythia.h"
using namespace Pythia8;

#include "UltraFastSim.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

bbHAnalysis::bbHAnalysis(TFile *o, Pythia *p, UltraFastSim *u, bool v) : outFile(o), pythia(p), ufs(u), verbose_(v), iEvent(0) {
  outFile->cd();
  tree = new TTree("bbHTree","Tree containing bbH -> tau_mu, tau_h analysis objects");
  tree->Branch("nElecs",   &nElecs,   "nElecs/I");
  tree->Branch("nMuons",   &nMuons,   "nMuons/I");
  tree->Branch("nTaus",    &nTaus,    "nTaus/I");
  tree->Branch("nRTaus",   &nRTaus,   "nRTaus/I");
  tree->Branch("nbQrks",   &nbQrks,   "nbQrks/I");
  tree->Branch("ncQrks",   &ncQrks,   "ncQrks/I");
  tree->Branch("nJets",    &nJets,    "nJets/I");
  tree->Branch("nBJets",   &nBJets,   "nBJets/I");
  tree->Branch("muonPt",   &muonPt,   "muonPt/D");
  tree->Branch("muonEta",  &muonEta,  "muonEta/D");
  tree->Branch("muonPhi",  &muonPhi,  "muonPhi/D");
  tree->Branch("tauPt",    &tauPt,    "tauPt/D");
  tree->Branch("tauEta",   &tauEta,   "tauEta/D");
  tree->Branch("tauPhi",   &tauPhi,   "tauPhi/D");
  tree->Branch("taumPt",   &taumPt,   "taumPt/D");
  tree->Branch("taumEta",  &taumEta,  "taumEta/D");
  tree->Branch("taumPhi",  &taumPhi,  "taumPhi/D");
  tree->Branch("tauhPt",   &tauhPt,   "tauhPt/D");
  tree->Branch("tauhEta",  &tauhEta,  "tauhEta/D");
  tree->Branch("tauhPhi",  &tauhPhi,  "tauhPhi/D");
  tree->Branch("rTauPt",   &rTauPt,   "rTauPt/D");
  tree->Branch("rTauEta",  &rTauEta,  "rTauEta/D");
  tree->Branch("rTauPhi",  &rTauPhi,  "rTauPhi/D");
  tree->Branch("bQrk1Pt",  &bQrk1Pt,  "bQrk1Pt/D");
  tree->Branch("bQrk1Eta", &bQrk1Eta, "bQrk1Eta/D");
  tree->Branch("bQrk1Phi", &bQrk1Phi, "bQrk1Phi/D");
  tree->Branch("bQrk2Pt",  &bQrk2Pt,  "bQrk2Pt/D");
  tree->Branch("bQrk2Eta", &bQrk2Eta, "bQrk2Eta/D");
  tree->Branch("bQrk2Phi", &bQrk2Phi, "bQrk2Phi/D");
  tree->Branch("bJet1Pt",  &bJet1Pt,  "bJet1Pt/D");
  tree->Branch("bJet1Eta", &bJet1Eta, "bJet1Eta/D");
  tree->Branch("bJet1Phi", &bJet1Phi, "bJet1Phi/D");
  tree->Branch("bJet2Pt",  &bJet2Pt,  "bJet2Pt/D");
  tree->Branch("bJet2Eta", &bJet2Eta, "bJet2Eta/D");
  tree->Branch("bJet2Phi", &bJet2Phi, "bJet2Phi/D");

  tree->Branch("bJet1PtStdGeom",  &bJet1PtStdGeom,  "bJet1PtStdGeom/D");
  tree->Branch("bJet1EtaStdGeom", &bJet1EtaStdGeom, "bJet1EtaStdGeom/D");
  tree->Branch("bJet1PhiStdGeom", &bJet1PhiStdGeom, "bJet1PhiStdGeom/D");
  tree->Branch("bJet2PtStdGeom",  &bJet2PtStdGeom,  "bJet2PtStdGeom/D");
  tree->Branch("bJet2EtaStdGeom", &bJet2EtaStdGeom, "bJet2EtaStdGeom/D");
  tree->Branch("bJet2PhiStdGeom", &bJet2PhiStdGeom, "bJet2PhiStdGeom/D");
 

  tree->Branch("gVisMass", &gVisMass, "gVisMass/D");
  tree->Branch("rVisMass", &rVisMass, "rVisMass/D");
  tree->Branch("ET",       &ET,       "ET/D");
  tree->Branch("HT",       &HT,       "HT/D");
  tree->Branch("MET",      &MET,      "MET/D");
  tree->Branch("MHT",      &MHT,      "MHT/D");

  //ml added branches for invariant mass analysis
  tree->Branch("muonE", &muonE,"muonE/D");
  tree->Branch("muonPx", &muonPx,"muonPx/D");
  tree->Branch("muonPy", &muonPy,"muonPy/D");
  tree->Branch("muonPz", &muonPz,"muonPz/D");
  tree->Branch("taumE", &taumE,"taumuonE/D");
  tree->Branch("taumPx", &taumPx,"taumPx/D");
  tree->Branch("taumPy", &taumPy,"taumPy/D");
  tree->Branch("taumPz", &taumPz,"taumPz/D");
  tree->Branch("hadE", &hadE,"hadE/D");
  tree->Branch("hadPx", &hadPx,"hadPx/D");
  tree->Branch("hadPy", &hadPy,"hadPy/D");
  tree->Branch("hadPz", &hadPz,"hadPz/D");
  tree->Branch("tauhE", &tauhE,"tauhE/D");
  tree->Branch("tauhPx", &tauhPx,"tauhPx/D");
  tree->Branch("tauhPy", &tauhPy,"tauhPy/D");
  tree->Branch("tauhPz", &tauhPz,"tauPz/D");
  
  tree->Branch("nBJetsStdGeom", &nBJetsStdGeom,"nBJetStdsGeom/I");
  
  tree->Branch("nbJetCut", &nbJetCut,"nbJetCut/I");
  tree->Branch("nbJetCutStdGeom", &nbJetCutStdGeom,"nbJetCutStdGeom/I");
  tree->Branch("nbasicCut",&nbasicCut,"nbasicCut/I");
  tree->Branch("nbasicCutStdGeom",&nbasicCutStdGeom,"nbasicCutStdGeom/I");

  tree->Branch("tau1E", &tau1E,"tau1E/D");
  tree->Branch("tau1Px", &tau1Px,"tau1Px/D");
  tree->Branch("tau1Py", &tau1Py,"tau1Py/D");
  tree->Branch("tau1Pz", &tau1Pz,"tau1Pz/D");
  tree->Branch("tau2E", &tau2E,"tau2E/D");
  tree->Branch("tau2Px", &tau2Px,"tau2Px/D");
  tree->Branch("tau2Py", &tau2Py,"tau2Py/D");
  tree->Branch("tau2Pz", &tau2Pz,"tau2Pz/D");

  tree->Branch("rTau1E", &rTau1E,"rTau1E/D");
  tree->Branch("rTau1Px", &rTau1Px,"rTau1Px/D");
  tree->Branch("rTau1Py", &rTau1Py,"rTau1Py/D");
  tree->Branch("rTau1Pz", &rTau1Pz,"rTau1Pz/D");
  tree->Branch("rTau2E", &rTau2E,"rTau2E/D");
  tree->Branch("rTau2Px", &rTau2Px,"rTau2Px/D");
  tree->Branch("rTau2Py", &rTau2Py,"rTau2Py/D");
  tree->Branch("rTau2Pz", &rTau2Pz,"rTau2Pz/D");


}

bbHAnalysis::~bbHAnalysis() {
}

bool bbHAnalysis::end() {
  outFile->cd();
  tree->Write();
  outFile->Close();
  return true;
}

bool bbHAnalysis::run() {

  if(verbose_) {   
    cout << "Number of Gen Chrgd = " << ufs->chargedHadronList().size() << endl;
    cout << "Number of Gen Neutr = " << ufs->neutralHadronList().size() << endl;
    cout << "Number of Gen Phtns = " << ufs->photonList().size() << endl;
    cout << "Number of Gen Elcns = " << ufs->electronList().size() << endl;
    cout << "Number of Gen Muons = " << ufs->muonList().size() << endl;
    cout << "Number of Gen Taus =  " << ufs->genTauList().size() << endl;
    cout << "Number of GenVisTau = " << ufs->visTauList().size() << endl;
    cout << "Number of Rec Taus =  " << ufs->tauList().size() << endl;
    cout << "Number of c Quarks =  " << ufs->cQuarkList().size() << endl;
    cout << "Number of b Quarks =  " << ufs->bQuarkList().size() << endl;
    cout << "Number of Jets =      " << ufs->jetList().size() << endl;
    cout << "Number of bJets =     " << ufs->bJetList().size() << endl;
    cout << "Number of bJetsStdGeom= " << ufs->bJetListStdGeom().size()<<endl;
    cout << endl;
  }

 

  // Store information for top muon, hadronic tau (not matching muon), top bJets (not matching tau)

  nElecs = ufs->electronList().size();
  nMuons = ufs->muonList().size();
  nTaus  = ufs->genTauList().size();
  nRTaus  = ufs->tauList().size();
  nbQrks = ufs->bQuarkList().size();
  ncQrks = ufs->cQuarkList().size();
  nJets = ufs->jetList().size();
  nBJets = ufs->bJetList().size();
  nBJetsStdGeom=ufs->bJetListStdGeom().size();
  muonPt = 0;
  muonEta = -999.;
  muonPhi = -999.;
  tauPt = 0;
  tauEta = -999.;
  tauPhi = -999.;
  taumPt = 0;
  taumEta = -999.;
  taumPhi = -999.;
  tauhPt = 0;
  tauhEta = -999.;
  tauhPhi = -999.;
  rTauPt = 0;
  rTauEta = -999.;
  rTauPhi = -999.;
  bQrk1Pt = 0;
  bQrk1Eta = -999.;
  bQrk1Phi = -999.;
  bQrk2Pt = 0;
  bQrk2Eta = -999.;
  bQrk2Phi = -999.;
  bJet1Pt = 0;
  bJet1Eta = -999.;
  bJet1Phi = -999.;
  bJet2Pt = 0;
  bJet2Eta = -999.;
  bJet2Phi = -999.;
  bJet1PtStdGeom = 0;
  bJet1EtaStdGeom  = -999.;
  bJet1PhiStdGeom  = -999.;
  bJet2PtStdGeom  = 0;
  bJet2EtaStdGeom  = -999.;
  bJet2PhiStdGeom  = -999.;
  gVisMass = -999.;
  rVisMass = -999.;
  ET = ufs->getET();
  HT = ufs->getHT();
  MET = ufs->getMET().Pt();
  MHT = ufs->getMHT().Pt();


  //ml added 
  nbJetCut=0;//passed bjet cut- 1 is passed
  nbJetCutStdGeom=0;
  nbasicCut=0;
  nbasicCutStdGeom=0;
  
  bJetPtCut=30.;
  muonE=0;
  muonPx=0;
  muonPy=0;
  muonPz=0;
  taumE=0;
  taumPx=0;
  taumPy=0;
  taumPz=0;
  
  hadE=0;
  hadPx=0;
  hadPy=0;
  hadPz=0;
  tauhE=0;
  tauhPx=0;
  tauhPy=0;
  tauhPz=0;
  
  tau1E=0;
  tau1Px=0;
  tau1Py=0;
  tau1Pz=0;
  tau2E=0;
  tau2Px=0;
  tau2Py=0;
  tau2Pz=0;
 
  rTau1E=0;
  rTau1Px=0;
  rTau1Py=0;
  rTau1Pz=0;
  rTau2E=0;
  rTau2Px=0;
  rTau2Py=0;
  rTau2Pz=0;

  hadronPt = 0;


  for(unsigned int i =0; i<ufs->chargedHadronList().size();i++){
     if(ufs->chargedHadronList()[i].Pt() > hadronPt) {
       hadE=ufs->chargedHadronList()[i].Energy();
       hadPx=ufs->chargedHadronList()[i].Px();
       hadPy=ufs->chargedHadronList()[i].Py();
       hadPz=ufs->chargedHadronList()[i].Pz();
       
       hadronPt=ufs->chargedHadronList()[i].Pt();
     }
  }

  int muonMother;
  TLorentzVector muonP;
  for(unsigned int i = 0; i < ufs->muonList().size(); i++) {
    //loops over all the muons and stores highest pt one
    if(ufs->muonList()[i].Pt() > muonPt) {
      muonMother = ufs->muonList()[i].GetFirstMother();
      muonPt = ufs->muonList()[i].Pt();
      muonEta = ufs->muonList()[i].Eta();
      muonPhi = ufs->muonList()[i].Phi();
      ufs->muonList()[i].Momentum(muonP);

      //ml added
      muonE = ufs->muonList()[i].Energy();
      muonPx = ufs->muonList()[i].Px();
      muonPy = ufs->muonList()[i].Py();
      muonPz = ufs->muonList()[i].Pz();   
    }
  }

  //this is the loop over all generated taus
  for(unsigned int i = 0; i < ufs->genTauList().size(); i++) {
    if(ufs->genTauList()[i].Pt() > tauPt) {
      tauPt = ufs->genTauList()[i].Pt();
      tauEta = ufs->genTauList()[i].Eta();
      tauPhi = ufs->genTauList()[i].Phi();

     	tau2E=tau1E;
	tau2Px=tau1Px;
	tau2Py=tau1Py;
	tau2Pz=tau1Pz;
       
	tau1E=ufs->genTauList()[i].Energy();
	tau1Px=ufs->genTauList()[i].Px();
	tau1Py=ufs->genTauList()[i].Py();
	tau1Pz=ufs->genTauList()[i].Pz();

	//	cout<<tau1E<<" "<<tau2E<<endl;
    }
  }

  TLorentzVector taumP;
  TLorentzVector tauhP;
  for(unsigned int i = 0; i < ufs->visTauList().size(); i++) {
    //loop over all visible taus that decay to lepton and hadron
    if(muonMother == ufs->visTauList()[i].GetSecondMother()) {  // We stole the second mother spot to store the pythia particle index
      if(ufs->visTauList()[i].Pt() > taumPt) {
	taumPt = ufs->visTauList()[i].Pt();
	taumEta = ufs->visTauList()[i].Eta();
	taumPhi = ufs->visTauList()[i].Phi();
	ufs->visTauList()[i].Momentum(taumP);
	//ml added
	taumE=ufs->visTauList()[i].Energy();
	taumPx=ufs->visTauList()[i].Px();
	taumPy=ufs->visTauList()[i].Py();
	taumPz=ufs->visTauList()[i].Pz();	  
      }
    }
    else {
      if(ufs->visTauList()[i].Pt() > tauhPt) {
	tauhPt = ufs->visTauList()[i].Pt();
	tauhEta = ufs->visTauList()[i].Eta();
	tauhPhi = ufs->visTauList()[i].Phi();
	ufs->visTauList()[i].Momentum(tauhP);

	tauhE=ufs->visTauList()[i].Energy();
	tauhPx=ufs->visTauList()[i].Px();
	tauhPy=ufs->visTauList()[i].Py();
	tauhPz=ufs->visTauList()[i].Pz();
      }
    }
  }
  
  TLorentzVector gTauTauP = taumP + tauhP;
  gVisMass = gTauTauP.M();


  //full tau list highest taus rtau final state (i.e. no hadrons)
  TLorentzVector rTauP;
  for(unsigned int i = 0; i < ufs->tauList().size(); i++) {
    if(ufs->tauList()[i].Pt() > rTauPt) {
      rTauPt = ufs->tauList()[i].Pt();
      rTauEta = ufs->tauList()[i].Eta();
      rTauPhi = ufs->tauList()[i].Phi();
      ufs->tauList()[i].Momentum(rTauP);

      //save two highest momentum rtaus   
	rTau2E=rTau1E;
	rTau2Px=rTau1Px;
	rTau2Py=rTau1Py;
	rTau2Pz=rTau1Pz;

	rTau1E=ufs->tauList()[i].Energy();
	rTau1Px=ufs->tauList()[i].Px();
	rTau1Py=ufs->tauList()[i].Py();
	rTau1Pz=ufs->tauList()[i].Pz();
      
    }
  }
  TLorentzVector rTauTauP = muonP + rTauP;
  rVisMass = rTauTauP.M();


  //save the highest two pt of bquarks etc
  if(nbQrks > 0) {
    bQrk1Pt = ufs->bQuarkList()[0].Pt();
    bQrk1Eta = ufs->bQuarkList()[0].Eta();
    bQrk1Phi = ufs->bQuarkList()[0].Phi();
    if(nbQrks > 1) {
      bQrk2Pt = ufs->bQuarkList()[1].Pt();
      bQrk2Eta = ufs->bQuarkList()[1].Eta();
      bQrk2Phi = ufs->bQuarkList()[1].Phi();
    }
  }

  if(nBJets > 0) {
    bJet1Pt = ufs->bJetList()[0].Pt();
    bJet1Eta = ufs->bJetList()[0].Eta();
    bJet1Phi = ufs->bJetList()[0].Phi();
    if(nBJets > 1) {
      bJet2Pt = ufs->bJetList()[1].Pt();
      bJet2Eta = ufs->bJetList()[1].Eta();
      bJet2Phi = ufs->bJetList()[1].Phi();
    }
  }

  if(nBJetsStdGeom > 0) {
    bJet1PtStdGeom = ufs->bJetListStdGeom()[0].Et();
    bJet1EtaStdGeom = ufs->bJetListStdGeom()[0].rap();
    bJet1PhiStdGeom = ufs->bJetListStdGeom()[0].phi();
    if(nBJetsStdGeom > 1) {
      bJet2PtStdGeom = ufs->bJetListStdGeom()[1].Et();
      bJet2EtaStdGeom = ufs->bJetListStdGeom()[1].rap();
      bJet2PhiStdGeom = ufs->bJetListStdGeom()[1].phi();
    }
  }
  //ml added - passes secondary pt cut - saves 1 for pass, 0 otherwise 
  if(nMuons >0 && nTaus >1 && nBJets > 1){
    nbasicCut = 1; 
     
  }
  if(nBJets>1){
    if(bJet1Pt > bJetPtCut && bJet2Pt>bJetPtCut){
      nbJetCut=1;
    }
  }

  if(nMuons>0 && nTaus>1 && nBJetsStdGeom > 1){
    nbasicCutStdGeom =1;
  }
  if(nBJetsStdGeom > 1){
    if(bJet1PtStdGeom > bJetPtCut && bJet2PtStdGeom > bJetPtCut){
      nbJetCutStdGeom=1;
    }
  }



  tree->Fill();
  iEvent++;
  return true;
}

