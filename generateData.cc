// This simple program generates various signal processes for 
// various physics processes.
// It includes the mixing of required level of pileup.
// It then takes pythia output and runs it through very simplistic
// detector simulation using parameterizations, including making
// isolated light leptons, photons, composite tau, jet and b-tagged 
// jet objects.
// It saves the output in root format for later analysis.

#include <string>
#include <iostream>

using namespace std;
#include <math.h>

#include "Pythia8/Pythia.h"
using namespace Pythia8; 

#include "UltraFastSim.h"
#include "UFSDataStore.h"
#include "UFSFilter.h"

int poisson(double mean, Rndm *rndmPtr) {
  static double oldMean = -1;
  static double g;
  if(mean != oldMean) {
    oldMean = mean;
    if(mean == 0) {
      g = 0;
    }
    else {
      g = exp(-mean);
    }
  }    
  double em = -1;
  double t = 1;
  do {
    em++;
    t *= rndmPtr->flat();
  } while(t > g);
  return em;
}

int main(int argc, char **argv) {
  // Generator. Process selection. LHC initialization. Histogram.

  Pythia pythia;

  if(argc >= 2) pythia.readFile(argv[1]);
  else {
    cerr << "Command syntax: " 
	 << argv[0] 
	 << " <command file name> [<numberOfEvents>] [<runNumber>==<randomNumberSeed>] [<meanPileUpLevel>] [<filterFileName>]" 
	 << endl;
    exit(1);
  }

  // Get default parameters from the pythia cards file

  int   nEvents    = pythia.mode("Main:numberOfEvents");
  int     nList    = pythia.mode("Main:numberToList");
  int     nShow    = pythia.mode("Main:timesToShow");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");
  bool   showCS    = pythia.flag("Main:showChangedSettings");
  bool   showAS    = pythia.flag("Main:showAllSettings");
  bool   showCPD   = pythia.flag("Main:showChangedParticleData");
  bool   showAPD   = pythia.flag("Main:showAllParticleData");

  // Overwrite with command line options, when appropriate

  if(argc >= 3) nEvents = atoi(argv[2]);
  int runNumber = 0;
  if(argc >= 4) runNumber = atoi(argv[3]);
  float meanPileupEventCount = 0;
  if(argc >= 5) meanPileupEventCount = atof(argv[4]);
  UFSFilter* filter = 0;
  if(argc == 6) filter = new UFSFilter(argv[5]);

  // Initialize pythia

  pythia.init();

  const int beamA = pythia.info.idA();
  const int beamB = pythia.info.idB();
  const double cmEnergy = pythia.info.eCM();

  // List settings.
  if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  Rndm * rndmPtr = &pythia.rndm;
  rndmPtr->init(runNumber);

  // Pileup pythia for proton beams

  Pythia *pileupPythia;
  if(beamA == 2212 && beamB == 2212) {
    pileupPythia = new Pythia();
    pileupPythia->init();
    // Use beamA = beamB = 2212 and cmEnergy = Hard Interaction cmEnergy
    pileupPythia->readString("SoftQCD:minBias = on");
    cout << "beam A, B, s" << beamA << ", " << beamB << ", " << cmEnergy << endl;
    pileupPythia->readString("Beams:idA = 2212");
    pileupPythia->readString("Beams:idB = 2212");
    std::stringstream ss;
    ss << "Beams:eCM = " << cmEnergy;
    pileupPythia->readString(ss.str());
    pileupPythia->rndm.init(runNumber);
  }

  // Setup analysis class
  string hFileName(argv[1]);
  hFileName.erase(hFileName.rfind("."));  

  // Ultra Fast Simulator

  UltraFastSim ufs;
  UFSDataStore dataStore(hFileName.c_str(), &ufs);

  // Begin event loop
  for (int iEvent = 0, iAbort = 0; iEvent < nEvents; ) {

    // Generate event. Skip if error. Unless, errors are too many.
    if (!pythia.next()) {
      iAbort++;
      if(iAbort > nAbort) {
	cerr << "Too many aborted events = " << nAbort << endl;
	exit(1);
      }
      continue;
    }

    // Filter if requested
    if(filter != 0) {
      if(!filter->doFilter(&pythia)) continue;
    }

    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
    }

    // Add pileup for proton beams

    if(beamA == 2212) {
      int pileupEventCount = 0;
      if(meanPileupEventCount > 0) pileupEventCount = poisson(meanPileupEventCount, rndmPtr);
      for (int puEvent = 0; puEvent < pileupEventCount; ) {
	if(!pileupPythia->next()) continue;
	double vx = rndmPtr->gauss()*2.; // Mean beam spread is 2 mm in x-y and 75 mm in z 
	double vy = rndmPtr->gauss()*2.;
	double vz = rndmPtr->gauss()*75.;
	for(int i = 0; i < pileupPythia->event.size(); i++) {
	  Particle& particle = pileupPythia->event[i];
	  particle.xProd(vx+particle.xProd());
	  particle.yProd(vy+particle.yProd());
	  particle.zProd(vz+particle.zProd());
	  if(particle.status() > 0) {
	    if(particle.isVisible()) {
	      if(particle.pT() > 0.5) {
		pythia.event.append(particle);
	      }
	    }
	  }
	}
	puEvent++;
      }
    }

    // Ultra fast simulation of detector effects
    if(!ufs.run(pythia.event, rndmPtr))
      {
	cerr << "Ultra fast simulation failed - aborting" << endl;
	exit(2);
      }

    // Store data

    if(!dataStore.run())
      {
	cerr << "Failed to store data" << endl;
	exit(3);
      }

    // End of event loop. Statistics. Histogram. Done.
    if(!(iEvent % 100)) cout << "Processed event " << iEvent << endl;
    iEvent++;
    
  }

  // Save output

  pythia.stat();
  if(beamA == 2212) {
    pileupPythia->stat();
  }
  if(filter != 0) delete filter;

  return 0;

}
