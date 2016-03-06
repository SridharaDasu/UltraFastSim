#include "UNL_ZH.h"
#include <string>
#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TCollection.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TNetFile.h"

int main(int argc, char **argv) {
  if(argc != 3) {
      std::cerr << "Command syntax: " << argv[0] << " inFileList outFile" << std::endl;
      exit(2);
  }

  // Read in list of files for TChain
  std::ifstream inFile;
  inFile.open(argv[1]);
  if (!inFile) {
    std::cerr << "Can't open " << inFile << std::endl;
    exit(1);
  }
  TChain *chain = new TChain("UltraFastSim");
  while (! inFile.eof()) {
    char dummy[256];
    inFile.getline(dummy,256); // read whole line
    if (dummy[0] != '#' && dummy[0] != '\0' && dummy[0] != ' ') { // Skip lines with # as first char
      chain->AddFile(dummy);
    }
  }
  chain->ls();
  std::cout << "Reading " << chain->GetEntries() << " events..." << std::endl;

  UNL_ZH * ana = new UNL_ZH(chain);  

  // Open the output root file
  ana->SetOutputFile(argv[2]);
  ana->Begin(chain);

  TObjArray *fileElements=chain->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while ( (chEl=(TChainElement*)next()) ) {
    if (chEl->GetTitle() == "\n") continue;
    TFile *f = TFile::Open(chEl->GetTitle());
    if (!f) continue;
    std::cout << "Reading " << chEl->GetTitle() << std::endl;
    TTree *t = (TTree*)f->Get("UltraFastSim");
    t->SetNotify(ana);
    ana->Init(t);
    for (int entry=0; entry < t->GetEntries(); ++entry) {
      ana->Process(entry);
    }
    f->Close();
  }
  ana->Terminate();

    //  chain->Process(ana); // loop over everything

  // Cleanup
  delete chain;
  delete ana;
}
