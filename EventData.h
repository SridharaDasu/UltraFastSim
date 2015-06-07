#ifndef EventData_H
#define EventData_H

#include <vector>

#include "TObject.h"
#include "TParticle.h"
#include "TLorentzVector.h"

class EventData : public TObject {

public:

  EventData() {;}

  virtual ~EventData() {;}

  const std::vector<TParticle>& genTauList() {return genTaus;}
  const std::vector<TParticle>& visTauList() {return visTaus;}
  const std::vector<TParticle>& bQuarkList() {return bQuarks;}
  const std::vector<TParticle>& cQuarkList() {return cQuarks;}
  const std::vector<TParticle>& photonList() {return photons;}
  const std::vector<TParticle>& electronList() {return electrons;}
  const std::vector<TParticle>& muonList() {return muons;}
  const std::vector<TParticle>& tauList() {return taus;}
  const std::vector<TParticle>& chargedHadronList() {return chargedHadrons;}
  const std::vector<TParticle>& neutralHadronList() {return neutralHadrons;}
  const std::vector<TLorentzVector>& jetList() {return jets;}

  const TLorentzVector& getMET() {return MET;}
  const TLorentzVector& getMHT() {return MHT;}

  double getET() {return ET;}
  double getHT() {return HT;}

protected:

  std::vector<TParticle> genTaus;
  std::vector<TParticle> bQuarks;
  std::vector<TParticle> cQuarks;

  std::vector<TParticle> photons;
  std::vector<TParticle> electrons;
  std::vector<TParticle> muons;
  std::vector<TParticle> visTaus;
  std::vector<TParticle> taus;
  std::vector<TParticle> chargedHadrons;
  std::vector<TParticle> neutralHadrons;

  std::vector<TLorentzVector> jets;

  TLorentzVector MET;
  TLorentzVector MHT;

  double ET;
  double HT;

  TLorentzVector primaryVertex;

  ClassDef (EventData, 1)
};

#endif
