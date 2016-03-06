#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>

using namespace std;

int main(int argc, char **argv) {

  vector<string> bJetOP;
  bJetOP.push_back("Loose");
  bJetOP.push_back("Medium");
  bJetOP.push_back("Tight");
  vector<string> geom;
  geom.push_back("StdGeom");
  geom.push_back("Phase-1");
  int AnalyzedEvents=0;
  int LeptonPairPreSelection=0;
  int LeptonPairSelection=0;
  int ZMassWindow=0;
  int ZptSelection=0;
  int JetPairPreSelection=0;
  int JetPairSelection=0;
  int BJetPairPreSelection[2][3] = {{0, 0, 0}, {0, 0, 0}};
  int BJetPairSelection[2][3] = {{0, 0, 0}, {0, 0, 0}};
  int HptSelection[2][3]= {{0, 0, 0}, {0, 0, 0}};
  int DphiSelection[2][3]={{0, 0, 0}, {0, 0, 0}};
  int JetVetoSelection[2][3]={{0, 0, 0}, {0, 0, 0}};
  int HMassWindow[2][3]={{0, 0, 0}, {0, 0, 0}};

  fstream inFile;
  for(int i = 1; i < argc; i++) {
    cout << "Processing " << argv[i] << endl;
    inFile.open(argv[i]);
    if(inFile.fail()) exit(1);
    char line[1024];
    strcpy(line, "junk");
    while(strncmp(line, " ====", 5) != 0 && !inFile.fail()) {
      inFile.getline(line, 1024);
    }
    if(!inFile.eof()) {
      int n;
      inFile.getline(line, 1024, ':'); inFile >> n; AnalyzedEvents += n; inFile.getline(line, 1024, '\n');
      inFile.getline(line, 1024, ':'); inFile >> n; LeptonPairPreSelection += n; inFile.getline(line, 1024, '\n');
      inFile.getline(line, 1024, ':'); inFile >> n; LeptonPairSelection += n; inFile.getline(line, 1024, '\n');
      inFile.getline(line, 1024, ':'); inFile >> n; ZMassWindow += n; inFile.getline(line, 1024, '\n');
      inFile.getline(line, 1024, ':'); inFile >> n; ZptSelection += n; inFile.getline(line, 1024, '\n');
      inFile.getline(line, 1024, ':'); inFile >> n; JetPairPreSelection += n; inFile.getline(line, 1024, '\n');
      inFile.getline(line, 1024, ':'); inFile >> n; JetPairSelection += n; inFile.getline(line, 1024, '\n');
      for(int b=0;b<3;b++){
	for(int g=0;g<2;g++) {
	  inFile.getline(line, 1024);
	  inFile.getline(line, 1024);
	  inFile.getline(line, 1024);
	  inFile.getline(line, 1024, ':'); inFile >> n; BJetPairPreSelection[g][b] += n; inFile.getline(line, 1024, '\n');
	  inFile.getline(line, 1024, ':'); inFile >> n; BJetPairSelection[g][b] += n; inFile.getline(line, 1024, '\n');
	  inFile.getline(line, 1024, ':'); inFile >> n; HptSelection[g][b] += n; inFile.getline(line, 1024, '\n');
	  inFile.getline(line, 1024, ':'); inFile >> n; DphiSelection[g][b] += n; inFile.getline(line, 1024, '\n');
	  inFile.getline(line, 1024, ':'); inFile >> n; JetVetoSelection[g][b] += n; inFile.getline(line, 1024, '\n');
	  inFile.getline(line, 1024, ':'); inFile >> n; HMassWindow[g][b] += n; inFile.getline(line, 1024, '\n');
	  inFile.getline(line, 1024);
	}
      }
    }
    else
      {
	cout << "Failed processing " << argv[i] << endl;
      }
    inFile.close();
  }
  cout<<" ======================================== " <<endl;
  cout<<" Analyzed Events          : "<<AnalyzedEvents<<endl;
  cout<<" Lepton pair preselection : "<<LeptonPairPreSelection<<endl;
  cout<<" Lepton pair    selection : "<<LeptonPairSelection<<endl;
  cout<<" Z Invariant Mass Window  : "<<ZMassWindow<<endl;
  cout<<" Z pt                     : "<<ZptSelection<<endl;
  cout<<" Jet pair preselection    : "<<JetPairPreSelection<<endl;
  cout<<" Jet pair    selection    : "<<JetPairSelection<<endl;
  for(int b=0;b<3;b++){
    for(int g=0;g<2;g++) {
      cout<<" ======================================== " <<endl;
      cout<<" Algo - "<<bJetOP[b]<<" b-tagging; "<<geom[g]<<" geometry"<<endl;
      cout<<" ======================================== " <<endl;
      cout<<" B-jet pair preselection  : "<<BJetPairPreSelection[g][b]<<endl;
      cout<<" B-jet pair    selection  : "<<BJetPairSelection[g][b]<<endl;
      cout<<" H pt                     : "<<HptSelection[g][b]<<endl;
      cout<<" Dphi(H,Z)                : "<<DphiSelection[g][b]<<endl;
      cout<<" Jet Veto selection       : "<<JetVetoSelection[g][b]<<endl;
      cout<<" H Invariant Mass Window  : "<<HMassWindow[g][b]<<endl;
      cout<<" ======================================== " <<endl;
    }
  }
  return 0;
}
