#ifndef UFSDataStore_H
#define UFSDataStore_H

class EventData;
class TFile;
class TTree;

class UFSDataStore {
public:
  UFSDataStore(const char *name, EventData *ufs);
  ~UFSDataStore();
  bool run();
private:
  UFSDataStore();
  EventData *ufs;
  TFile *outFile;
  TTree *tree;
};

#endif
