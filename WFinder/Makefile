CXX           = g++
CXXFLAGS      = $(OPT2) -m32 -pipe -Wall -W -Woverloaded-virtual -g
LD            = g++
LDFLAGS       = $(OPT2) -m32 -bind_at_load
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)

ROOTINCLUDES = -I$(ROOTSYS)/include
ROOTLDLIBS = `$(ROOTSYS)/bin/root-config --cflags --glibs` -lEG

UFSDIR = $(HOME)/CMSSW_4_2_6/UserCode/dasu/UltraFastSim
UFSINCLUDES = -I$(UFSDIR)

INCLUDES = $(ROOTINCLUDES) $(UFSINCLUDES)
LDLIBS = -L. $(ROOTLDLIBS)

all : ZFinder WFinder
	@echo Built $@

clean :
	rm -rf eventdict.cc eventdict.h ZFinder WFinder

ZFinder : ZFinder.C eventdict.cc 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o ZFinder

WFinder : WFinder.C eventdict.cc 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o WFinder

eventdict.cc :
	rootcint -f eventdict.cc -c $(UFSDIR)/EventData.h $(UFSDIR)/LinkDef.h
