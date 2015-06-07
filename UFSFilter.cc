#include "UFSFilter.h"

#include <stdio.h>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

UFSFilter::UFSFilter(const char *name) : trueCount(0), falseCount(0) {
  FILE* file = fopen(name, "r");
  UFSFilterData *data = new UFSFilterData();
  while(fscanf(file, "%d %d %d %f %f %f %f", 
	       &data->id,
	       &data->minCount, &data->maxCount, 
	       &data->ptMin, &data->ptMax, 
	       &data->etaMin, &data->etaMax) == 7) {
    filterData.push_back(data);
    data = new UFSFilterData();
  }
  fclose(file);
}

UFSFilter::~UFSFilter() {
  std::cout << "UFSFilter: trueCount = " << trueCount << std::endl;
  std::cout << "UFSFilter: falseCount = " << falseCount << std::endl;
}

bool UFSFilter::doFilter(Pythia* pythia) {

  for(int j = 0; j < filterData.size(); j++)filterData[j]->count = 0;

  for(int i = 0; i < pythia->event.size(); i++) {
    Particle& particle = pythia->event[i];
    for(int j = 0; j < filterData.size(); j++) {
      if(particle.id() == filterData[j]->id &&
	 particle.eta() > filterData[j]->etaMin &&
	 particle.eta() < filterData[j]->etaMax &&
	 particle.pT() > filterData[j]->ptMin &&
	 particle.pT() < filterData[j]->ptMax) {
	if(abs(filterData[j]->id) == 15){
	  if (particle.status() > 0)
	    return false;
	  else 
	    filterData[j]->count++;
	}
	else filterData[j]->count++;
      }
    }
  }
  int nFilt=0;
  for(int j = 0; j < filterData.size(); j++) {
    if(filterData[j]->count >= filterData[j]->minCount &&
       filterData[j]->count <= filterData[j]->maxCount) {
      nFilt++;
    }
  }  
  if (nFilt == filterData.size()) {
    trueCount++;
    return true;
  }
  else {
    falseCount++;
    return false;
  }
}
