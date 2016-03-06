#include "ZHAnalysis.h"

using namespace std;

void ZHAnalysis::BookTree() {
  myevent_ = new myevent;
  t=new TTree("t","tree");
  t->Branch("myevent","myevent",&myevent_,256000,1);
}


bool ZHAnalysis::Run(UltraFastSim * ufs){

  vector<TParticle>          GenTaus_ = ufs->genTauList();
  vector<TParticle>          BQuarks_ = ufs->bQuarkList();
  vector<TParticle>          CQuarks_ = ufs->cQuarkList();
  vector<TParticle>          Photons_ = ufs->photonList();
  vector<TParticle>          Electrons_ = ufs->electronList();
  vector<TParticle>          Muons_ = ufs->muonList();
  vector<TParticle>          MuonsStdGeom_ = ufs->muonListStdGeom();
  vector<TParticle>          Taus_ = ufs->tauList();
  vector<TParticle>          ChargedHadrons_ = ufs->chargedHadronList();
  vector<TParticle>          NeutralHadrons_ = ufs->neutralHadronList();
  vector<TLorentzVector>     Jets_ = ufs->jetList();
  vector<TLorentzVector>     BJets_ = ufs->bJetList();
  vector<TLorentzVector>     BJetsStdGeom_ = ufs->bJetListStdGeom();

  myevent_->clear();

  for(vector<TParticle>::const_iterator it=GenTaus_.begin();it!=GenTaus_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->GenTaus).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=BQuarks_.begin();it!=BQuarks_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->BQuarks).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=CQuarks_.begin();it!=CQuarks_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->CQuarks).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=Photons_.begin();it!=Photons_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->Photons).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=Electrons_.begin();it!=Electrons_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->Electrons).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=Muons_.begin();it!=Muons_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->Muons).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=MuonsStdGeom_.begin();it!=MuonsStdGeom_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->MuonsStdGeom).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=Taus_.begin();it!=Taus_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->Taus).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=ChargedHadrons_.begin();it!=ChargedHadrons_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->ChargedHadrons).push_back(obj);
  }
  for(vector<TParticle>::const_iterator it=NeutralHadrons_.begin();it!=NeutralHadrons_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->Energy();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->NeutralHadrons).push_back(obj);
  }
  for(vector<TLorentzVector>::const_iterator it=Jets_.begin();it!=Jets_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->E();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->Jets).push_back(obj);
  }
  for(vector<TLorentzVector>::const_iterator it=BJets_.begin();it!=BJets_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->E();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->BJets).push_back(obj);
  }
  for(vector<TLorentzVector>::const_iterator it=BJetsStdGeom_.begin();it!=BJetsStdGeom_.end();it++){
    myobject obj;
    obj.px=it->Px();
    obj.py=it->Py();
    obj.pz=it->Pz();
    obj.E=it->E();
    obj.pt=it->Pt();
    obj.eta=it->Eta();
    obj.phi=it->Phi();
    (myevent_->BJetsStdGeom).push_back(obj);
  }

  t->Fill();

  return true;

}//Run 
