CXX           = g++
CXXFLAGS      = $(OPT2) -pipe -Wall -W -Woverloaded-virtual -g 
LD            = g++
LDFLAGS       = $(OPT2)
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)

PYINCLUDES = -I$(PYTHIA8DATA)/../include
PYLDLIBS = -L$(PYTHIA8DATA)/../lib -lpythia8 -llhapdfdummy

FASTJET ?= $(CMSSW_RELEASE_BASE)/../../../external/fastjet/3.0.1/
FASTJETINCLUDES = -I$(FASTJET)/include
FASTJETLDLIBS = -L$(FASTJET)/lib/ -lfastjet

ROOTINCLUDES = -I$(ROOTSYS)/include
ROOTLDLIBS = `$(ROOTSYS)/bin/root-config --cflags --glibs` -lEG

INCLUDES = -I. -I../UltraFastSim $(PYINCLUDES) $(FASTJETINCLUDES) $(ROOTINCLUDES)
LDLIBS = -L. -lUFSTauTauAnalysis -L../UltraFastSim -lUltraFastSim $(PYLDLIBS) $(FASTJETLDLIBS) $(ROOTLDLIBS)

all : analyzeData
	@echo Built $@

clean :
	rm -rf analyzeData libUFSTauTauAnalysis.a

analyzeData : analyzeData.cc ../UltraFastSim/eventdict.cc libUFSTauTauAnalysis.a
	$(LD) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o $@
	@echo Built $@

libUFSTauTauAnalysis.a : libUFSTauTauAnalysis.a(EventAnalysis.o)
	@echo Built $@

analyzeData.cc : EventAnalysis.h

EventAnalysis.cc : EventAnalysis.h

.cc.a :
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^
	ar rv $@ $*.o
	rm -f $*.o
	@echo Compiled $@
