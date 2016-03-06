#ifndef ZHAnalysis_H
#define ZHAnalysis_H

#include <vector>
#include "TParticle.h"
#include "UltraFastSim.h"
#include "fastjet/PseudoJet.hh"
#include "TFile.h"
#include "TTree.h"
#include "myevent.h"
#include "LinkDef.h"


class ZHAnalysis {

public:

  void BookTree();
  bool Run(UltraFastSim *);

 private:

  myevent *myevent_;
  TTree *t;
  
};

#endif
