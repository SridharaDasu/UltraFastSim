CXX           = g++
CXXFLAGS      = $(OPT2) -m32 -pipe -Wall -W -Woverloaded-virtual -g
LD            = g++
LDFLAGS       = $(OPT2) -m32 -bind_at_load
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)

PYINCLUDES = -I$(CMS_PATH)/slc5_ia32_gcc434/external/pythia8/135/include/
PYLDLIBS = -L$(CMS_PATH)/slc5_ia32_gcc434/external/pythia8/135/lib -lpythia8 -llhapdfdummy

FASTJETINCLUDES = -I$(CMS_PATH)/slc5_ia32_gcc434/external/fastjet/2.4.0/include
FASTJETLDLIBS = -L$(CMS_PATH)/slc5_ia32_gcc434/external/fastjet/2.4.0/lib -lfastjet

ROOTINCLUDES = -I$(CMS_PATH)/slc5_ia32_gcc434/lcg/root/5.22.00d-cms19/include
ROOTLDLIBS = -L$(CMS_PATH)/slc5_ia32_gcc434/lcg/root/5.22.00d-cms19/lib -lPhysics -lEG -lMathCore -lHist -lTree -lCore -lNet -lRIO

INCLUDES = $(PYINCLUDES) $(FASTJETINCLUDES) $(ROOTINCLUDES)
LDLIBS = -L. -lUltraFastSim $(PYLDLIBS) $(FASTJETLDLIBS) $(ROOTLDLIBS)

all : generateData ZHOfflineAnalysis  BBAOfflineAnalysis
	@echo Built $@

clean :
	rm -rf generateData generateData.dSYM ZHOfflineAnalysis ZHOfflineAnalysis.dSYM   BBAOfflineAnalysis BBAOfflineAnalysis.dSYM eventdict.cc eventdict.h libUltraFastSim.a

generateData : generateData.cc  eventdict.cc libUltraFastSim.a
	$(LD) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -lUltraFastSim -o $@
	@echo Built $@

libUltraFastSim.a : libUltraFastSim.a(UltraFastSim.o) libUltraFastSim.a(bTagData.o) libUltraFastSim.a(muonTagData.o) libUltraFastSim.a(UFSDataStore.o)
	@echo Built $@

UltraFastSim.cc : UltraFastSim.h

bTagData.cc : bTagData.h

muonTagData.cc : muonTagData.h

UFSDataStore.cc : UFSDataStore.h

.cc.a :
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^
	ar rv $@ $*.o
	rm -f $*.o
	@echo Compiled $@

ZHOfflineAnalysis : ZHOfflineAnalysis.C eventdict.cc 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o ZHOfflineAnalysis

BBAOfflineAnalysis : BBAOfflineAnalysis.C eventdict.cc 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o BBAOfflineAnalysis

eventdict.cc :
	rootcint -f eventdict.cc -c UltraFastSim.h LinkDef.h
