#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <math.h>
using namespace std;
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "UltraFastSim.h"
#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath

#include "EventAnalysis.h"

int main(){
  vector<TParticle> bQs;
  vector<TParticle> cQs;
  vector<TParticle> muons;
  vector<TLorentzVector> jets;
  TLorentzVector Z;
  TLorentzVector H;
  int AnalyzedEvents=0;
  ifstream ConfigFile;
  string OutputRootFile;
  cout<<"Enter histogram output file name: ";
  cin>>OutputRootFile;
  string filename;
  EventAnalysis eventAnalysis(OutputRootFile);
  cout<<"Enter input UltraFastSim file names (one per line) -- EOF when done : " << endl;
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
    }//for loop on events
    delete rootTree;
    delete ufs;
    delete file;
  }//for loop on files
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events          : "<<AnalyzedEvents<<endl;
  eventAnalysis.finalize(); 

  return 0;
}
