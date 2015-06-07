#include "UFSDataStore.h"

#include <string>
using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "EventData.h"

#include "LinkDef.h"

UFSDataStore::UFSDataStore(const char* n, EventData *u) : ufs(u), outFile(0), tree(0) {
  string name(n);
  name += ".root";
  outFile = new TFile(name.c_str(), "recreate");
  tree = new TTree("EventData", "UFSTree");
  tree->Branch("EventData", "EventData", &ufs, 256000, 1);
}

UFSDataStore::~UFSDataStore() {
  outFile->cd();
  outFile->Write();
  outFile->Close();
}

bool UFSDataStore::run() {
  tree->Fill();
  return true;
}
