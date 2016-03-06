#ifndef EventAnalysis_H
#define EventAnalysis_H

#include <string>
#include <map>

using namespace std;

class TH1;
class TFile;
class UltraFastSim;

class EventAnalysis {
 public:
  EventAnalysis(string histFileName);
  //  EventAnalysis(string histFileName = "histograms.root");
  virtual ~EventAnalysis() {;}
  void processEvent(UltraFastSim *ufs);
  void finalize();
 private:
  map<string, TH1*> histograms;
  TFile* file;
};

#endif
