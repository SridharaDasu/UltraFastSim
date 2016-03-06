// Load Sridhara's class into an active ROOT session
{
  gSystem->Load("/opt/osg/app/cmssoft/cms/slc5_ia32_gcc434/external/fastjet/2.4.0/lib/libfastjet.so");
  gSystem->AddIncludePath("-I/opt/osg/app/cmssoft/cms/slc5_ia32_gcc434/external/pythia8/142/include -I/opt/osg/app/cmssoft/cms/slc5_ia32_gcc434/external/fastjet/2.4.0/include");
  gSystem->CompileMacro("UltraFastSim.cc","k");
}
