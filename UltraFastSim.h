#ifndef UltraFastSim_H
#define UltraFastSim_H

#include <vector>

namespace Pythia8 {
  class Rndm;
  class Event;
  class Particle;
}

namespace fastjet {
  class JetDefinition;
  class ClusterSequence;
}

#include "TObject.h"
#include "TParticle.h"
#include "TLorentzVector.h"

#include "EventData.h"

class UltraFastSim : public EventData {

public:

  UltraFastSim() {;}

  virtual ~UltraFastSim() {;}

  bool run(Pythia8::Event &event, Pythia8::Rndm *r);

private:

  void clear();

  void makeJets();
  void makeTaus();
  void makeETSums();

  void setCommon(Pythia8::Particle& particle, TParticle& smearedParticle);
  void tkSmear(Pythia8::Particle& particle, TParticle& smearedParticle);
  void emSmear(Pythia8::Particle& particle, TParticle& smearedParticle);
  void hdSmear(Pythia8::Particle& particle, TParticle& smearedParticle);

};

#endif
