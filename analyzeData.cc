#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <math.h>
#include <iomanip>
using namespace std;

#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "UltraFastSim.h"
#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath

#include "EventAnalysis.h"

void dump(string tag, const vector<TParticle>& particleList) {
  for(unsigned int i = 0; i < particleList.size(); i++) {
    cout << setw(12) << tag
	 << setw(4) << particleList[i].GetPdgCode()	 
	 << setw(16) << particleList[i].Pt()
	 << setw(16) << particleList[i].Eta()
	 << setw(16) << particleList[i].Phi()
	 << endl;
  }
}

void dump(string tag, const vector<TLorentzVector>& jetList) {
  for(unsigned int i = 0; i < jetList.size(); i++) {
    cout << setw(16) << tag
	 << setw(16) << jetList[i].Pt()
	 << setw(16) << jetList[i].Eta()
	 << setw(16) << jetList[i].Phi()
	 << endl;
  }
}

void dump(string tag, const TLorentzVector& object) {
  cout << setw(16) << tag
       << setw(16) << object.Pt()
       << setw(16) << object.Eta()
       << setw(16) << object.Phi()
       << endl;
}

int main(int argc, char** argv){
  ifstream ConfigFile;
  string OutputRootFile;
  cout<<"Enter histogram output file name: ";
  cin>>OutputRootFile;
  string filename;
  EventAnalysis eventAnalysis(OutputRootFile);
  cout<<"Enter input UltraFastSim file names (one per line) -- EOF when done : " << endl;
  int AnalyzedEvents=0;
  while(cin>>filename && filename!="EOF"){
    TFile *file=TFile::Open(filename.c_str());
    cout << "file = " << file << endl;
    TTree *rootTree = (TTree*)file->Get("EventData");
    cout << "rootTree = " << rootTree << endl;
    int nev = int(rootTree->GetEntries());
    cerr<<"Reading : "<<filename<<endl;
    cerr<<"Number of entries is : "<<nev<<endl;
    TBranch *branch = rootTree->GetBranch("EventData");
    UltraFastSim *ufs=new UltraFastSim();
    branch->SetAddress(&ufs);
    for(int i = 0; i < nev; i++){
      AnalyzedEvents++;
      if(!(AnalyzedEvents%1000))cerr<<"event number "<<AnalyzedEvents<<endl;
      rootTree->GetEvent(i);
      eventAnalysis.processEvent(ufs);
      // Dump event
      if(argc == 2 && strncmp(argv[1], "--dump", 3) == 0) {
	cout << "Event" << setw(5) << i << endl;
	dump("Muon", ufs->muonList());
	dump("Electron", ufs->electronList());
	dump("GenTau", ufs->genTauList());
	dump("VisTau", ufs->visTauList());
	dump("RecTau", ufs->tauList());
	dump("Photon", ufs->photonList());
	dump("bQuark", ufs->bQuarkList());
	dump("cQuark", ufs->cQuarkList());
	dump("Jet", ufs->jetList());
	dump("MET", ufs->getMET());
	dump("MHT", ufs->getMHT());
	cout << "_______________________________________________________________________________" << endl;
      }
    }//for loop on events
    delete rootTree;
    delete ufs;
    delete file;
  }//for loop on files
  eventAnalysis.finalize(); 
  cout<<" Analyzed Events          : "<<AnalyzedEvents<<endl;
  return 0;
}
