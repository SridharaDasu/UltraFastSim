CXX           = g++
CXXFLAGS      = $(OPT2) -pipe -Wall -W -Woverloaded-virtual -g -std=c++11
LD            = g++
LDFLAGS       = $(OPT2)
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)

PYINCLUDES ?= -I`pythia8-config --includedir`
PYLDLIBS ?= -L`pythia8-config --libdir` -lpythia8

FASTJET ?= $(ROOTSYS)/../../../external/fastjet/3.1.0-odfocd/
FASTJETINCLUDES ?= -I$(FASTJET)/include
FASTJETLDLIBS ?= -L$(FASTJET)/lib/ -lfastjet

ROOTINCLUDES = -I$(ROOTSYS)/include
ROOTLDLIBS = `$(ROOTSYS)/bin/root-config --cflags --glibs` -lEG

INCLUDES = -I. $(PYINCLUDES) $(FASTJETINCLUDES) $(ROOTINCLUDES)
LDLIBS = -L. -lUltraFastSim $(PYLDLIBS) $(FASTJETLDLIBS) $(ROOTLDLIBS)

all : generateData analyzeData libUltraFastSim.a
	@echo "Built $@"

clean :
	rm -rf analyzeData generateData eventdict.cc eventdict.h libUltraFastSim.a

generateData : generateData.cc  eventdict.cc libUltraFastSim.a
	$(LD) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o $@
	@echo Built $@

analyzeData : analyzeData.cc  eventdict.cc libUltraFastSim.a
	$(LD) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o $@
	@echo Built $@

libUltraFastSim.a : libUltraFastSim.a(UltraFastSim.o) libUltraFastSim.a(UFSDataStore.o) libUltraFastSim.a(EventAnalysis.o) libUltraFastSim.a(UFSFilter.o)
	@echo Built $@

UltraFastSim.cc : UltraFastSim.h

UFSDataStore.cc : UFSDataStore.h

analyzeData.cc : EventAnalysis.h

EventAnalysis.cc : EventAnalysis.h

UFSFilter.cc : UFSFilter.h

.cc.a :
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^
	ar rv $@ $*.o
	rm -f $*.o
	@echo Compiled $@

eventdict.cc : UltraFastSim.h
	rootcling -f eventdict.cc -c EventData.h LinkDef.h
