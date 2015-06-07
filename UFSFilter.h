#ifndef UFSFilter_H
#define UFSFilter_H

#include <vector>

namespace Pythia8 {
  class Pythia;
}

class UFSFilterData {
 public:
  UFSFilterData() {;}
  ~UFSFilterData() {;}
  int id;
  int minCount;
  int maxCount;
  float ptMin;
  float ptMax;
  float etaMin;
  float etaMax;
  int count;
};

class UFSFilter {
public:
  UFSFilter(const char *name);
  ~UFSFilter();
  bool doFilter(Pythia8::Pythia* pythia);
private:
  UFSFilter();
  std::vector<UFSFilterData*> filterData;
  int trueCount;
  int falseCount;
};

#endif
