CXX           = g++
CXXFLAGS      = $(OPT2) -pipe -Wall -W -Woverloaded-virtual -g
LD            = g++
LDFLAGS       = $(OPT2) 
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)

ROOTINCLUDES = -I$(ROOTSYS)/include
ROOTLDLIBS = `$(ROOTSYS)/bin/root-config --cflags --glibs` -lEG

UFSDIR = ../UltraFastSim
UFSINCLUDES = -I$(UFSDIR)

INCLUDES = $(ROOTINCLUDES) $(UFSINCLUDES)
LDLIBS = -L. $(ROOTLDLIBS)

all : ZFinder
	@echo Built $@

clean :
	rm -rf eventdict.cc eventdict.h ZFinder

ZFinder : ZFinder.C eventdict.cc 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o ZFinder

eventdict.cc :
	rootcint -f eventdict.cc -c $(UFSDIR)/EventData.h $(UFSDIR)/LinkDef.h
