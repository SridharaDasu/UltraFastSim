#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <math.h>
using namespace std;
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TDCacheFile.h"
#include "TSystem.h"
#include "UltraFastSim.h"
#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath

double deltaR(double eta1, double eta2, double phi1, double phi2) 
{
  double dphi = fabs(phi1-phi2);
  if(dphi > M_PI) 
    {
      dphi = (2*M_PI - dphi);
    }
  double deta = fabs(eta1-eta2);
  return sqrt(dphi*dphi + deta*deta);
}
double deltaphi(double phi1, double phi2) 
{
  double dphi = fabs(phi1-phi2);
  if(dphi > M_PI) 
    {
      dphi = (2*M_PI - dphi);
    }
  return dphi;
}
double CorrectPhi(double phi, double x, double y)
{
  if(phi>0 && x<0 && y<0)phi+=M_PI;
  if(phi<0 && x>0 && y<0)phi=(2*M_PI-fabs(phi));
  if(phi<0 && x<0 && y>0)phi=(M_PI-fabs(phi));
  if(y<0)phi-=(2*M_PI);//without this you have 0<phi<2pi, but with this
  return phi;//you get -pi<phi<pi
}

int main(int argc, char **argv){
  if(argc!=2) {cout<<"error - need Config*.txt argument "<<argc<<endl; exit(1);}

  vector<TParticle> muons;
  vector<TParticle> vistaus;
  vector<TParticle> gentaus;
  vector<TParticle> chadrons;
  vector<TParticle>  rtaus;

const int ntypes = 6;
  vector<TLorentzVector> bjets[ntypes];
  int mhPreSelection[ntypes]={0};
  int mhSelection[ntypes]={0};
  int mhAddSelection[ntypes]={0};
  
  int mtPreSelection[ntypes]={0};
  int mtSelection[ntypes]={0};
   int mtAddSelection[ntypes]={0};

  int ttPreSelection[ntypes]={0};
  int ttSelection[ntypes]={0};
  int ttAddSelection[ntypes]={0};

  int mh_lmPreSelection=0;
  int mh_lmSelection[4]={0};
  int mt_lmSelection[4]={0};
  
  int nlomepre[2]={0};

  int mh_lmSelectionP[4]={0};
  int mt_lmSelectionP[4]={0};
  
  int AnalyzedEvents=0;
 
  //================================================
  //event loop

  
  ifstream ConfigFile;
  ConfigFile.open(argv[1]);
  if(!ConfigFile){cout<<"Unable to read the config file "<<argv[1];exit(1);}
  
  string OutputRootFile;
  double CrossSection;
  double Lumi;
  // bool Phase1;

  ConfigFile>>OutputRootFile;
  ConfigFile>>CrossSection;
  ConfigFile>>Lumi;
 // ConfigFile>>Phase1;

  //file should be created before tree to avoid memory resident trees
  TFile fout(OutputRootFile.c_str(),"recreate");
  TH1F* hptmuons = new TH1F("Muons Pt","Muons Pt",2000,0.,2000.);
  TH1F* hetamuons = new TH1F("Muons Eta","Muons Eta",100,-5.,5.);
  TH1F* hphimuons = new TH1F("Muons Phi","Muons Phi",100,-7.,7.);

  
  ///////////////////ml added histograms
	TH1F* hbjet_mult[ntypes];
	 hbjet_mult[0] = new TH1F("B-jets Mult. Loose","B-jets Mult. Loose",10,0.,10.);
	 hbjet_mult[1] = new TH1F("B-jets Mult. StdGeom Loose","B-jets Mult. StdGeom Loose",10,0.,10.);
	 hbjet_mult[2] = new TH1F("B-jets Mult. Medium","B-jets Mult. Medium",10,0.,10.);
	 hbjet_mult[3] = new TH1F("B-jets Mult. StdGeom Medium","B-jets Mult. StdGeom Medium",10,0.,10.);
	 hbjet_mult[4] = new TH1F("B-jets Mult. Tight","B-jets Mult. Tight",10,0.,10.);
	 hbjet_mult[5] = new TH1F("B-jets Mult. StdGeom Tight","B-jets Mult. StdGeom Tight",10,0.,10.);
    

	TH1F* hetbjets[ntypes];
	hetbjets[0] = new TH1F("B-jets Et Loose","B-jets Et Loose",2000,0.,2000.);
	hetbjets[1] = new TH1F("B-jets Et StdGeom Loose","B-jets Et StdGeom Loose",2000,0.,2000.);
	hetbjets[2] = new TH1F("B-jets Et Medium","B-jets Et Medium",2000,0.,2000.);
	hetbjets[3] = new TH1F("B-jets Et StdGeom Medium","B-jets Et StdGeom Medium",2000,0.,2000.);
	hetbjets[4] = new TH1F("B-jets Et Tight","B-jets Et Tight",2000,0.,2000.);
	hetbjets[5] = new TH1F("B-jets Et StdGeom Tight","B-jets Et StdGeom Tight",2000,0.,2000.);
  
	TH1F* hetabjets[ntypes];
	hetabjets[0] = new TH1F("B-jets Eta Loose","B-jets Eta Loose",2000,-5.,5.);
	hetabjets[1] = new TH1F("B-jets Eta StdGeom Loose","B-jets Eta StdGeom Loose",2000,-5.,5.);
	hetabjets[2] = new TH1F("B-jets Eta Medium","B-jets Eta Medium",2000,-5.,5.);
	hetabjets[3] = new TH1F("B-jets Eta StdGeom Medium","B-jets Eta StdGeom Medium",2000,-5.,5.);
	hetabjets[4] = new TH1F("B-jets Eta Tight","B-jets Eta Tight",2000,-5.,5.);
	hetabjets[5] = new TH1F("B-jets Eta StdGeom Tight","B-jets Eta StdGeom Tight",20,-5.,5.);

	
	
	TH1F* hphibjets[ntypes];
	hphibjets[0] = new TH1F("B-jets phi Loose","B-jets phi Loose",2000,-7.,7.);
	hphibjets[1] = new TH1F("B-jets phi StdGeom Loose","B-jets phi StdGeom Loose",2000,-7.,7.);
	hphibjets[2] = new TH1F("B-jets phi Medium","B-jets phi Medium",2000,-7.,7.);
	hphibjets[3] = new TH1F("B-jets phi StdGeom Medium","B-jets phi StdGeom Medium",2000,-7.,7.);
	hphibjets[4] = new TH1F("B-jets phi Tight","B-jets phi Tight",2000,-7.,7.);
	hphibjets[5] = new TH1F("B-jets phi StdGeom Tight","B-jets phi StdGeom Tight",2000,-7.,7.);

	
const int nsel = 3;
//***//
	TH1F* mh_hinvmass[ntypes][nsel];
	mh_hinvmass[0][0]= new TH1F("H invmass preselection from Muon and Hadron Loose","H invmass preselection from Muon and Hadron Loose",2000,0.,2000.);
   mh_hinvmass[1][0]= new TH1F("H invmass preselection from Muon and Hadron stdgeom Loose","H invmass preselection from Muon and Hadron stdgeom Loose",2000,0.,2000.);
	mh_hinvmass[2][0]= new TH1F("H invmass preselection from Muon and Hadron Medium","H invmass preselection from Muon and Hadron Medium",2000,0.,2000.);
   mh_hinvmass[3][0]= new TH1F("H invmass preselection from Muon and Hadron stdgeom Medium","H invmass preselection from Muon and Hadron stdgeom Medium",2000,0.,2000.);
	mh_hinvmass[4][0]= new TH1F("H invmass preselection from Muon and Hadron Tight","H invmass preselection from Muon and Hadron Tight",2000,0.,2000.);
   mh_hinvmass[5][0]= new TH1F("H invmass preselection from Muon and Hadron stdgeom Tight","H invmass preselection from Muon and Hadron stdgeom Tight",2000,0.,2000.);

	mh_hinvmass[0][1]= new TH1F("H invmass from Muon and Hadron Loose","H invmass from Muon and Hadron Loose",2000,0.,2000.);
   mh_hinvmass[1][1]= new TH1F("H invmass from Muon and Hadron stdgeom Loose","H invmass from Muon and Hadron stdgeom Loose",2000,0.,2000.);
	mh_hinvmass[2][1]= new TH1F("H invmass from Muon and Hadron Medium","H invmass from Muon and Hadron Medium",2000,0.,2000.);
   mh_hinvmass[3][1]= new TH1F("H invmass from Muon and Hadron stdgeom Medium","H invmass from Muon and Hadron stdgeom Medium",2000,0.,2000.);
	mh_hinvmass[4][1]= new TH1F("H invmass from Muon and Hadron Tight","H invmass from Muon and Hadron Tight",2000,0.,2000.);
   mh_hinvmass[5][1]= new TH1F("H invmass from Muon and Hadron stdgeom Tight","H invmass from Muon and Hadron stdgeom Tight",2000,0.,2000.);

	mh_hinvmass[0][2]= new TH1F("H invmass selection from Muon and Hadron Loose","H invmass from Muon and Hadron Loose",2000,0.,2000.);
   mh_hinvmass[1][2]= new TH1F("H invmass selection from Muon and Hadron stdgeom Loose","H invmass from Muon and Hadron stdgeom Loose",2000,0.,2000.);
	mh_hinvmass[2][2]= new TH1F("H invmass selection from Muon and Hadron Medium","H invmass from Muon and Hadron Medium",2000,0.,2000.);
   mh_hinvmass[3][2]= new TH1F("H invmass selection from Muon and Hadron stdgeom Medium","H invmass from Muon and Hadron stdgeom Medium",2000,0.,2000.);
	mh_hinvmass[4][2]= new TH1F("H invmass selection from Muon and Hadron Tight","H invmass from Muon and Hadron Tight",2000,0.,2000.);
   mh_hinvmass[5][2]= new TH1F("H invmass selection from Muon and Hadron stdgeom Tight","H invmass from Muon and Hadron stdgeom Tight",2000,0.,2000.);

//higgs mass from mu+tauh
	TH1F* mutau_hinvmass[ntypes][nsel];
	 mutau_hinvmass[0][0]= new TH1F("H invmass preselection from Tauh+Muon Loose","H invmass preselection TM Loose",2000,0.,2000.);
	 mutau_hinvmass[1][0]= new TH1F("H invmass preselection from Tauh+Muon stdgeom Loose","H invmass preselection TM stdgeom Loose",2000,0.,2000.);
	 mutau_hinvmass[2][0]= new TH1F("H invmass preselection from Tauh+Muon Medium","H invmass preselection TM Medium",2000,0.,2000.);
	 mutau_hinvmass[3][0]= new TH1F("H invmass preselection from Tauh+Muon stdgeom Medium","H invmass preselection TM stdgeom Medium",2000,0.,2000.);
	 mutau_hinvmass[4][0]= new TH1F("H invmass preselection from Tauh+Muon Tight","H invmass preselection TM Tight",2000,0.,2000.);
	 mutau_hinvmass[5][0]= new TH1F("H invmass preselection from Tauh+Muon stdgeom Tight","H invmass preselection TM stdgeom Tight",2000,0.,2000.);

	 mutau_hinvmass[0][1]= new TH1F("H invmass from Tauh+Muon Loose","H invmass TM Loose",2000,0.,2000.);
	 mutau_hinvmass[1][1]= new TH1F("H invmass from Tauh+Muon stdgeom Loose","H invmass TM stdgeom Loose",2000,0.,2000.);
	 mutau_hinvmass[2][1]= new TH1F("H invmass from Tauh+Muon Medium","H invmass TM Medium",2000,0.,2000.);
	 mutau_hinvmass[3][1]= new TH1F("H invmass from Tauh+Muon stdgeom Medium","H invmass TM stdgeom Medium",2000,0.,2000.);
	 mutau_hinvmass[4][1]= new TH1F("H invmass from Tauh+Muon Tight","H invmass TM Tight",2000,0.,2000.);
	 mutau_hinvmass[5][1]= new TH1F("H invmass from Tauh+Muon stdgeom Tight","H invmass TM stdgeom Tight",2000,0.,2000.);

	 mutau_hinvmass[0][2]= new TH1F("H invmass selection from Tauh+Muon Loose","H invmass selection TM Loose",2000,0.,2000.);
	 mutau_hinvmass[1][2]= new TH1F("H invmass selection from Tauh+Muon stdgeom Loose","H invmass selection TM stdgeom Loose",2000,0.,2000.);
	 mutau_hinvmass[2][2]= new TH1F("H invmass selection from Tauh+Muon Medium","H invmass selection TM Medium",2000,0.,2000.);
	 mutau_hinvmass[3][2]= new TH1F("H invmass selection from Tauh+Muon stdgeom Medium","H invmass selection TM stdgeom Medium",2000,0.,2000.);
	 mutau_hinvmass[4][2]= new TH1F("H invmass selection from Tauh+Muon Tight","H invmass selection TM Tight",2000,0.,2000.);
	 mutau_hinvmass[5][2]= new TH1F("H invmass selection from Tauh+Muon stdgeom Tight","H invmass selection TM stdgeom Tight",2000,0.,2000.);


//higgs mass from taus
	TH1F* tau_hinvmass[ntypes][nsel];
	tau_hinvmass[0][0]= new TH1F("H invmass preselection from Taus Loose","H invmass preselection Taus Loose",2000,0.,2000.);
	tau_hinvmass[1][0]= new TH1F("H invmass preselection from Taus stdgeom Loose","H invmass preselection Taus stdgeom Loose",2000,0.,2000.);
	tau_hinvmass[2][0]= new TH1F("H invmass preselection from Taus Medium","H invmass preselection Taus Medium",2000,0.,2000.);
	tau_hinvmass[3][0]= new TH1F("H invmass preselection from Taus stdgeom Medium","H invmass preselection Taus stdgeom Medium",2000,0.,2000.);
	tau_hinvmass[4][0]= new TH1F("H invmass preselection from Taus Tight","H invmass preselection Taus Tight",2000,0.,2000.);
	tau_hinvmass[5][0]= new TH1F("H invmass preselection from Taus stdgeom Tight","H invmass preselection Taus stdgeom Tight",2000,0.,2000.);

 	tau_hinvmass[0][1]= new TH1F("H invmass from Taus Loose","H invmass Taus Loose",2000,0.,2000.);
	tau_hinvmass[1][1]= new TH1F("H invmass from Taus stdgeom Loose","H invmass Taus stdgeom Loose",2000,0.,2000.);
	tau_hinvmass[2][1]= new TH1F("H invmass from Taus Medium","H invmass Taus Medium",2000,0.,2000.);
	tau_hinvmass[3][1]= new TH1F("H invmass from Taus stdgeom Medium","H invmass Taus stdgeom Medium",2000,0.,2000.);
	tau_hinvmass[4][1]= new TH1F("H invmass from Taus Tight","H invmass Taus Tight",2000,0.,2000.);
	tau_hinvmass[5][1]= new TH1F("H invmass from Taus stdgeom Tight","H invmass Taus stdgeom Tight",2000,0.,2000.);

	tau_hinvmass[0][2]= new TH1F("H invmass selection from Taus Loose","H invmass selection Taus Loose",2000,0.,2000.);
	tau_hinvmass[1][2]= new TH1F("H invmass selection from Taus stdgeom Loose","H invmass selection Taus stdgeom Loose",2000,0.,2000.);
	tau_hinvmass[2][2]= new TH1F("H invmass selection from Taus Medium","H invmass selection Taus Medium",2000,0.,2000.);
	tau_hinvmass[3][2]= new TH1F("H invmass selection from Taus stdgeom Medium","H invmass selection Taus stdgeom Medium",2000,0.,2000.);
	tau_hinvmass[4][2]= new TH1F("H invmass selection from Taus Tight","H invmass selection Taus Tight",2000,0.,2000.);
	tau_hinvmass[5][2]= new TH1F("H invmass selection from Taus stdgeom Tight","H invmass selection Taus stdgeom Tight",2000,0.,2000.);


	//muon plots
	TH1F* mhptmuons[3];
	mhptmuons[0] = new TH1F("Muons Pt -highest Pt Loose","Muons Pt -highest Pt Loose",2000,0.,2000.);
	mhptmuons[1] = new TH1F("Muons Pt -highest Pt Medium","Muons Pt -highest Pt Medium",2000,0.,2000.);
	mhptmuons[2] = new TH1F("Muons Pt -highest Pt Tight","Muons Pt -highest Pt Tight",2000,0.,2000.);
	
	TH1F* mhetamuons[3];
	mhetamuons[0] = new TH1F("Muons Eta -highest Pt Loose","Muons Eta -highest Pt Loose",100,-5.,5.);
	mhetamuons[1] = new TH1F("Muons Eta -highest Pt Medium","Muons Eta -highest Pt Medium",100,-5.,5.);
	mhetamuons[2] = new TH1F("Muons Eta -highest Pt Tight","Muons Eta -highest Pt Tight",100,-5.,5.);
	
	TH1F* mhphimuons[3];
	mhphimuons[0] = new TH1F("Muons Phi -highest Pt Loose","Muons Phi -highest Pt Loose",100,-7.,7.);
	mhphimuons[1] = new TH1F("Muons Phi -highest Pt Medium","Muons Phi -highest Pt Medium",100,-7.,7.);
	mhphimuons[2] = new TH1F("Muons Phi -highest Pt Tight","Muons Phi -highest Pt Tight",100,-7.,7.);

	
	//taus
	TH1F* mhpttaus[3];
	mhpttaus[0] = new TH1F("Taus Pt -highest Pt Loose","Tau Pt -highest Pt Loose",2000,0.,2000.);
	mhpttaus[1] = new TH1F("Taus Pt -highest Pt Medium","Tau Pt -highest Pt Medium",2000,0.,2000.);
	mhpttaus[2] = new TH1F("Taus Pt -highest Pt Tight","Tau Pt -highest Pt Tight",2000,0.,2000.);
   
	TH1F* mhetataus[3];
	mhetataus[0] = new TH1F("Taus Eta -highest Pt Loose","Tau Eta -highest Pt Loose",100,-5.,5.);
	mhetataus[1] = new TH1F("Taus Eta -highest Pt Medium","Tau Eta -highest Pt Medium",100,-5.,5.);
	mhetataus[2] = new TH1F("Taus Eta -highest Pt Tight","Tau Eta -highest Pt Tight",100,-5.,5.);
   
	TH1F* mhphitaus[3];
	mhphitaus[0] = new TH1F("Taus Phi -highest Pt Loose","Tau Phi -highest Pt Loose",100,-7.,7.);
	mhphitaus[1] = new TH1F("Taus Phi -highest Pt Medium","Tau Phi -highest Pt Medium",100,-7.,7.);
	mhphitaus[2] = new TH1F("Taus Phi -highest Pt Tight","Tau Phi -highest Pt Tight",100,-7.,7.);
   
   //hadrons
   TH1F* mhpthadrons[3];
   mhpthadrons[0] = new TH1F("Hadrons Pt -highest Pt Loose","Hadrons Pt -highest Pt Loose",2000,0.,2000.);
   mhpthadrons[1] = new TH1F("Hadrons Pt -highest Pt Medium","Hadrons Pt -highest Pt Medium",2000,0.,2000.);
   mhpthadrons[2] = new TH1F("Hadrons Pt -highest Pt Tight","Hadrons Pt -highest Pt Tight",2000,0.,2000.);
   
   TH1F* mhetahadrons[3];
   mhetahadrons[0] = new TH1F("Hadrons Eta -highest Pt Loose","Hadrons Eta -highest Pt Loose",100,-5.,5.);
   mhetahadrons[1] = new TH1F("Hadrons Eta -highest Pt Medium","Hadrons Eta -highest Pt Medium",100,-5.,5.);
   mhetahadrons[2] = new TH1F("Hadrons Eta -highest Pt Tight","Hadrons Eta -highest Pt Tight",100,-5.,5.);
   
   TH1F* mhphihadrons[3];
   mhphihadrons[0] = new TH1F("Hadrons Phi -highest Pt Loose","Hadrons Phi -highest Pt Loose",100,-7.,7.);
   mhphihadrons[1] = new TH1F("Hadrons Phi -highest Pt Medium","Hadrons Phi -highest Pt Medium",100,-7.,7.);
   mhphihadrons[2] = new TH1F("Hadrons Phi -highest Pt Tight","Hadrons Phi -highest Pt Tight",100,-7.,7.);
   
   //r tauus
   TH1F* mhptrtaus[3];
   mhptrtaus[0] = new TH1F("Tauh Pt -highest Pt Loose","Tauh Pt -highest Pt Loose",2000,0.,2000.);
   mhptrtaus[1] = new TH1F("Tauh Pt -highest Pt Medium","Tauh Pt -highest Pt Medium",2000,0.,2000.);
   mhptrtaus[2] = new TH1F("Tauh Pt -highest Pt Tight","Tauh Pt -highest Pt Tight",2000,0.,2000.);
   
   TH1F* mhetartaus[3];
   mhetartaus[0] = new TH1F("Tauh Eta -highest Pt Loose","Tauh Eta -highest Pt Loose",100,-5.,5.);
   mhetartaus[1] = new TH1F("Tauh Eta -highest Pt Medium","Tauh Eta -highest Pt Medium",100,-5.,5.);
   mhetartaus[2] = new TH1F("Tauh Eta -highest Pt Tight","Tauh Eta -highest Pt Tight",100,-5.,5.);
   
   TH1F* mhphirtaus[3];
   mhphirtaus[0] = new TH1F("Tauh Phi -highest Pt Loose","Tauh Phi -highest Pt Loose",100,-7.,7.);
   mhphirtaus[1] = new TH1F("Tauh Phi -highest Pt Medium","Tauh Phi -highest Pt Medium",100,-7.,7.);
   mhphirtaus[2] = new TH1F("Tauh Phi -highest Pt Tight","Tauh Phi -highest Pt Tight",100,-7.,7.);
   
   //additional  loose/tight medium/tight
TH1F* mh_hinvmass_lm[4];

	mh_hinvmass_lm[0]= new TH1F("H invmass from Muon and Hadron Loose-Tight","H invmass from Muon and Hadron Loose-Tight",2000,0.,2000.);
	mh_hinvmass_lm[1]= new TH1F("H invmass from Muon and Hadron stdgeom Loose-Tight","H invmass from Muon and Hadron stdgeom Loose-Tight",2000,0.,2000.);
	mh_hinvmass_lm[2]= new TH1F("H invmass from Muon and Hadron Medium-Tight","H invmass from Muon and Hadron Medium-Tight",2000,0.,2000.);
	mh_hinvmass_lm[3]= new TH1F("H invmass from Muon and Hadron stdgeom Medium-Tight","H invmass from Muon and Hadron stdgeom Medium-Tight",2000,0.,2000.);

  
	TH1F* tau_hinvmass_lm[4];
 	tau_hinvmass_lm[0]= new TH1F("H invmass from Taus Loose-Tight","H invmass Taus Loose-Tight",2000,0.,2000.);
	tau_hinvmass_lm[1]= new TH1F("H invmass from Taus stdgeom Loose-Tight","H invmass Taus stdgeom Loose-Tight",2000,0.,2000.);
	tau_hinvmass_lm[2]= new TH1F("H invmass from Taus Medium-Tight","H invmass Taus Medium-Tight",2000,0.,2000.);
	tau_hinvmass_lm[3]= new TH1F("H invmass from Taus stdgeom Medium-Tight","H invmass Taus stdgeom Medium-Tight",2000,0.,2000.);
	
	TH1F* mutau_hinvmass_lm[4];
	mutau_hinvmass_lm[0]= new TH1F("H invmass from Tauh+Muon Loose-Tight","H invmass TM Loose-Tight",2000,0.,2000.);
	 mutau_hinvmass_lm[1]= new TH1F("H invmass from Tauh+Muon stdgeom Loose-Tight","H invmass TM stdgeom Loose-Tight",2000,0.,2000.);
	 mutau_hinvmass_lm[2]= new TH1F("H invmass from Tauh+Muon Medium-Tight","H invmass TM Medium-Tight",2000,0.,2000.);
	 mutau_hinvmass_lm[3]= new TH1F("H invmass from Tauh+Muon stdgeom Medium-Tight","H invmass TM stdgeom Medium-Tight",2000,0.,2000.);


	 //post selection
	 TH1F* mh_hinvmass_lmP[4];

	 mh_hinvmass_lmP[0]= new TH1F("H invmass from Muon and Hadron Loose-Tight Post","H invmass from Muon and Hadron Loose-Tight Post",2000,0.,2000.);
	mh_hinvmass_lmP[1]= new TH1F("H invmass from Muon and Hadron stdgeom Loose-Tight Post","H invmass from Muon and Hadron stdgeom Loose-Tight Post",2000,0.,2000.);
	mh_hinvmass_lmP[2]= new TH1F("H invmass from Muon and Hadron Medium-Tight Post","H invmass from Muon and Hadron Medium-Tight Post",2000,0.,2000.);
	mh_hinvmass_lmP[3]= new TH1F("H invmass from Muon and Hadron stdgeom Medium-Tight Post","H invmass from Muon and Hadron stdgeom Medium-Tight Post",2000,0.,2000.);

  
	TH1F* tau_hinvmass_lmP[4];
 	tau_hinvmass_lmP[0]= new TH1F("H invmass from Taus Loose-Tight Post","H invmass Taus Loose-Tight Post",2000,0.,2000.);
	tau_hinvmass_lmP[1]= new TH1F("H invmass from Taus stdgeom Loose-Tight Post","H invmass Taus stdgeom Loose-Tight Post",2000,0.,2000.);
	tau_hinvmass_lmP[2]= new TH1F("H invmass from Taus Medium-Tight Post","H invmass Taus Medium-Tight Post",2000,0.,2000.);
	tau_hinvmass_lmP[3]= new TH1F("H invmass from Taus stdgeom Medium-Tight Post","H invmass Taus stdgeom Medium-Tight Post",2000,0.,2000.);
	
	TH1F* mutau_hinvmass_lmP[4];
	mutau_hinvmass_lmP[0]= new TH1F("H invmass from Tauh+Muon Loose-Tight Post","H invmass TM Loose-Tight Post",2000,0.,2000.);
	 mutau_hinvmass_lmP[1]= new TH1F("H invmass from Tauh+Muon stdgeom Loose-Tight Post","H invmass TM stdgeom Loose-Tight Post",2000,0.,2000.);
	 mutau_hinvmass_lmP[2]= new TH1F("H invmass from Tauh+Muon Medium-Tight Post","H invmass TM Medium-Tight Post",2000,0.,2000.);
	 mutau_hinvmass_lmP[3]= new TH1F("H invmass from Tauh+Muon stdgeom Medium-Tight Post","H invmass TM stdgeom Medium-Tight Post",2000,0.,2000.);


		
	string filename;
	vector<string> fns;//file name container
	int SumEvents=0;
	while(ConfigFile>>filename && filename!="EOF"){
	  fns.push_back(filename);
	  cout<<"checking "<<filename<<endl;
	  TFile *file= TFile::Open(filename.c_str());
	  TTree *rootTree = (TTree*)file->Get("UltraFastSim");
	  int nev = int(rootTree->GetEntries());
	  cout<<"It has "<<nev<<" events."<<endl;
	  SumEvents+=nev;
	  delete rootTree;
	  delete file;
	}
	
  // ========== >>>> Loop over files starts here 
  for(vector<string>::iterator iterstr=fns.begin();iterstr!=fns.end();iterstr++){
    
    cout<<"Reading "<<*iterstr<<endl;
    TFile *file= TFile::Open((*iterstr).c_str());
    UltraFastSim *ufs=new UltraFastSim();
    TTree *rootTree = (TTree*)file->Get("UltraFastSim");
    int nev = int(rootTree->GetEntries());
    cout<<"number of entries is : "<<nev<<endl;
    TBranch *branch = rootTree->GetBranch("UltraFastSim");
    branch->SetAddress(&ufs);
  

    for(int i = 0; i < nev; i++){ 
      AnalyzedEvents++;
      if(!(AnalyzedEvents%1000))cout<<"event number "<<AnalyzedEvents<<endl;
      
      rootTree->GetEvent(i);
      muons=ufs->muonList();
      gentaus=ufs->genTauList(); //use gen taus to reconstruct mass
      //vistaus=ufs->visTauList(); //use vis taus to check mothers if needed
      rtaus=ufs->tauList();
      chadrons=ufs->chargedHadronList();
      
		bjets[0]=ufs->bJetListLoose();
		bjets[1]=ufs->bJetListLooseStdGeom();
		bjets[2]=ufs->bJetListMedium();
 		bjets[3]=ufs->bJetListMediumStdGeom();
		bjets[4]=ufs->bJetListTight();
		bjets[5]=ufs->bJetListTightStdGeom();

      
      for(vector<TParticle>::iterator itm=muons.begin();itm!=muons.end();itm++){
		hptmuons->Fill(itm->Pt()); 
		hetamuons->Fill(itm->Eta()); 
		hphimuons->Fill(itm->Phi()); 
      }

		for(int i=0; i<6 ; i++){
		  for(vector<TLorentzVector>::iterator itj=bjets[i].begin();itj!=bjets[i].end();itj++){
			hetbjets[i]->Fill(itj->Pt()); 
			hetabjets[i]->Fill(itj->Eta()); 
			hphibjets[i]->Fill(itj->Phi()); 
		  }
		}
       


      
      //need to sort hadrons and muons into highest pt
      double tau1E=0.;
      double tau1Px=0.;
      double tau1Py=0.;
      double tau1Pz=0.;
      double tau1Pt=0.;
      double tau1Eta=-999.;
      double tau1Phi=-999.;
      double tau2E=0.;
      double tau2Px=0.;
      double tau2Py=0.;
      double tau2Pz=0.;
      double tau2Pt=0.;
      double tau2Eta=-999.;
      double tau2Phi=-999.;
      double tauPt = 0.;
      //first taus
      for(vector<TParticle>::iterator itt=gentaus.begin();itt!=gentaus.end();itt++) {
	if(itt->Pt() > tauPt) {
	  tauPt = itt->Pt();
		  
	  tau2E=tau1E;
	  tau2Px=tau1Px;
	  tau2Py=tau1Py;
	  tau2Pz=tau1Pz;
	  tau2Pt=tau1Pt;
	  tau2Eta=tau1Eta;
	  tau2Phi=tau1Phi;
	  
	tau1E= itt->Energy();
	tau1Px=itt->Px();
	tau1Py=itt->Py();
	tau1Pz=itt->Pz();
	tau1Pt=itt->Pt();
	tau1Eta=itt->Eta();
	tau1Phi=itt->Phi();
	}
      }
      
      double muonPt=0.;
      double muonE=0.;
      double muonPx=0.;
      double muonPy=0.;
      double muonPz=0.;
      double muonEta=-999.;
      double muonPhi=-999.;
      
      //muons 
      for(vector<TParticle>::iterator imm=muons.begin();imm!=muons.end();imm++) { 
	//loops over all the muons and stores highest pt one
		if(imm->Pt() > muonPt) {
		  muonPt = imm->Pt();
		  muonEta= imm->Eta();
		  muonPhi=imm->Phi();
		  muonE = imm->Energy();
		  muonPx = imm->Px();
		  muonPy = imm->Py();
		  muonPz = imm->Pz();   
		}
      }

      double dphi = 0.;
      double maxdphi=0.;
      int gensel = 0; 
      //find gentau associated to muon
      for(int i = 0; i<gentaus.size(); i++){
	dphi = deltaphi(muonPhi,gentaus[i].Phi());
	if(dphi>maxdphi)
	  { 
	    maxdphi = dphi;
	    gensel = i;
	  }
	     
      }


      double hadronE=0.;
      double hadronPx=0.;
      double hadronPy=0.;
      double hadronPz=0.;
      double hadronPt=0.;
      double hadronEta=-999.;
      double hadronPhi=-999.;
      //chadrons
      for(vector<TParticle>::iterator ic=chadrons.begin();ic!=chadrons.end();ic++){
	if(ic->Pt() > hadronPt) {
	  hadronE = ic->Energy();
	  hadronPx= ic->Px();
	  hadronPy= ic->Py();
	  hadronPz= ic->Pz();
	  hadronEta=ic->Eta();
	  hadronPhi=ic->Phi();
	 
	  hadronPt= ic->Pt();
	}
      }

 double rtauE=0.;
 double rtauPx=0.;
 double rtauPy=0.;
 double rtauPz=0.;
 double rtauPt=0.;
 double rtauPhi=-999.;
 double rtauEta=-999.;
 
 //vistau - ie hadron tau
 for(vector<TParticle>::iterator irt=rtaus.begin();irt!=rtaus.end();irt++) {
   if(irt->Pt() > rtauPt) {
     rtauPt = irt->Pt();
     rtauEta=irt->Eta();
     rtauPhi=irt->Phi();
     
     rtauE= irt->Energy();
     rtauPx=irt->Px();
     rtauPy=irt->Py();
     rtauPz=irt->Pz();  
   }
 }

for(int i =0; i< 3; i++){
  if(muons.size()>0 &&  gentaus.size() >1 && bjets[i].size()>1){
    if(muonPt>0.){
      mhptmuons[i]->Fill(muonPt);
      mhetamuons[i]->Fill(muonEta);
      mhphimuons[i]->Fill(muonPhi);
    }
      
    if(tau1Pt >0. && tau2Pt>0.){
      mhpttaus[i]->Fill(tau1Pt);
      mhpttaus[i]->Fill(tau2Pt);
      mhetataus[i]->Fill(tau1Eta);
      mhetataus[i]->Fill(tau2Eta);
      mhphitaus[i]->Fill(tau1Phi);
      mhphitaus[i]->Fill(tau2Phi);
    }
    if(hadronPt>0){
      mhpthadrons[i]->Fill(hadronPt);
      mhetahadrons[i]->Fill(hadronEta);
      mhphihadrons[i]->Fill(hadronPhi);
    }
    if(rtauPt>0){
      mhptrtaus[i]->Fill(rtauPt);
      mhetartaus[i]->Fill(rtauEta);
      mhphirtaus[i]->Fill(rtauPhi);
    }
  }
}
     
      double mhHinvmass=0.;
      double mtHinvmass=0;
      double tauHinvmass=0.;
      
	for(int i=0; i<6;i++){
	  hbjet_mult[i]->Fill(bjets[i].size());
	}
 

if(muons.size()>0 &&  gentaus.size() >1 ){ //basic cuts
  mhHinvmass=(pow((muonE+hadronE),2)-
	      pow((muonPx+hadronPx),2)-
	      pow((muonPy+hadronPy),2)-
	      pow((muonPz+hadronPz),2));
  if(mhHinvmass>0.) mhHinvmass=sqrt(mhHinvmass);
  else mhHinvmass = 0.;
    //two taus//
  tauHinvmass=(pow((tau1E+tau2E),2)-
	       pow((tau1Px+tau2Px),2)-
	       pow((tau1Py+tau2Py),2)-
	       pow((tau1Pz+tau2Pz),2));
   if(tauHinvmass>0.)tauHinvmass=sqrt(tauHinvmass);
   else tauHinvmass = 0.;

  mtHinvmass=(pow((muonE+gentaus[gensel].Energy()),2)-
	      pow((muonPx+gentaus[gensel].Px()),2)-
	      pow((muonPy+gentaus[gensel].Py()),2)-
	      pow((muonPz+gentaus[gensel].Pz()),2));
  if(mtHinvmass>0.)mtHinvmass=sqrt(mtHinvmass);
  else mtHinvmass = 0.;

for(int i=0; i<6; i++){ //loop over types of loose/medium/tight squared selections of btagging
  if(bjets[i].size() > 1){
    if(mhHinvmass>0.){
      mhPreSelection[i]++;
      mh_hinvmass[i][0]->Fill(mhHinvmass);
    }
    
    if(tauHinvmass>0.){
      ttPreSelection[i]++;
      tau_hinvmass[i][0]->Fill(tauHinvmass);
    }
    //tau + muon   
    if(mtHinvmass>0.){
      mtPreSelection[i]++;
      mutau_hinvmass[i][0]->Fill(mtHinvmass);
    }
    
    
	//additional constraints
    if(bjets[i][0].Pt()>30. && bjets[i][1].Pt()>30.){
      if(mhHinvmass>0.){
	mhSelection[i]++;
	mh_hinvmass[i][1]->Fill(mhHinvmass);
      }
    
    if(tauHinvmass>0.){
      ttSelection[i]++;
      tau_hinvmass[i][1]->Fill(tauHinvmass);
    }
    //tau + muon   
    if(mtHinvmass>0.){
       mtSelection[i]++;
       mutau_hinvmass[i][1]->Fill(mtHinvmass);
    }
    
      
      if(muonPt>50. && rtauPt > 50.){
	 if(mhHinvmass>0.){
	   mhAddSelection[i]++;
	   mh_hinvmass[i][2]->Fill(mhHinvmass);
	 }
	 
	 if(tauHinvmass>0.){
	   ttAddSelection[i]++;
	   tau_hinvmass[i][2]->Fill(tauHinvmass);
	 }
	 //tau + muon   
	 if(mtHinvmass>0.){
	   mtAddSelection[i]++;
	   mutau_hinvmass[i][2]->Fill(mtHinvmass);
	 }
	
      }
    }//close additional constraints      
	
    
  }//bjets>1
 }//6 types

//loose/medium -tight
  if(bjets[4].size()>0){
       nlomepre[0]++;
       if(bjets[0].size()>(bjets[4].size())){ //tight ph1
	   mh_lmSelection[0]++;
	   mh_hinvmass_lm[0]->Fill(mhHinvmass); //loose/tight ph1
	   tau_hinvmass_lm[0]->Fill(tauHinvmass);
	   if(gentaus.size()>0){//tau+muon
	     mt_lmSelection[0]++;
	     mutau_hinvmass_lm[0]->Fill(mtHinvmass);
	   }
	   
	   if( (bjets[4][0].Pt()>bjets[0][1].Pt() && bjets[0][1].Pt()>30.) ||(bjets[0][0].Pt()>bjets[4][1].Pt() && bjets[4][1].Pt()>30.) || (bjets[4][0].Pt()==bjets[0][1].Pt() && bjets[0].size()>1 && bjets[0][2].Pt()> 30.) ){
	       mh_lmSelectionP[0]++;
	       mh_hinvmass_lmP[0]->Fill(mhHinvmass); //loose/tight ph1
	       tau_hinvmass_lmP[0]->Fill(tauHinvmass);
	       if(gentaus.size()>0){//tau+muon
		 mt_lmSelectionP[0]++;
		 mutau_hinvmass_lmP[0]->Fill(mtHinvmass);
	       }
	     
	   }  
       }
       
       if(bjets[2].size()>bjets[4].size()){    
	 mh_lmSelection[2]++; //med tight ph1
	 mh_hinvmass_lm[2]->Fill(mhHinvmass);
	 tau_hinvmass_lm[2]->Fill(tauHinvmass);
	 mt_lmSelection[2]++;
	 mutau_hinvmass_lm[2]->Fill(mtHinvmass);
	   
	 
	
	 if( (bjets[4][0].Pt()>bjets[2][1].Pt() && bjets[2][1].Pt()>30.) ||(bjets[2][0].Pt()>bjets[4][1].Pt() && bjets[4][1].Pt()>30.) || (bjets[4][0].Pt()==bjets[2][1].Pt() && bjets[2].size()>1 && bjets[2][2].Pt()> 30.) ){
	   mh_lmSelectionP[2]++;
	   mh_hinvmass_lmP[2]->Fill(mhHinvmass); //loose/tight ph1
	   tau_hinvmass_lmP[2]->Fill(tauHinvmass);
	   mt_lmSelectionP[2]++;
	   mutau_hinvmass_lmP[2]->Fill(mtHinvmass);
	     
	   
	 } 
       }
  }
 
     
  // below is to do Loose/Medium with Tight for stdg and ph1. 
     //std

 mh_lmPreSelection++;    
     if(bjets[5].size()>0){
       nlomepre[1]++;
       if(bjets[1].size()>(bjets[5].size())){ //tight ph1
	   mh_lmSelection[1]++;
	   mh_hinvmass_lm[1]->Fill(mhHinvmass); //loose/tight ph1
	   tau_hinvmass_lm[1]->Fill(tauHinvmass);
	   if(gentaus.size()>0){//tau+muon
	     mt_lmSelection[1]++;
	     mutau_hinvmass_lm[1]->Fill(mtHinvmass);
	   }
	   
	   if( (bjets[5][0].Pt()>bjets[1][1].Pt() && bjets[1][1].Pt()>30.) ||(bjets[1][0].Pt()>bjets[5][1].Pt() && bjets[5][1].Pt()>30.) || (bjets[5][0].Pt()==bjets[1][1].Pt() && bjets[1].size()>1 && bjets[1][2].Pt()> 30.) ){
	       mh_lmSelectionP[1]++;
	       mh_hinvmass_lmP[1]->Fill(mhHinvmass); //loose/tight ph1
	       tau_hinvmass_lmP[1]->Fill(tauHinvmass);
	       if(gentaus.size()>0){//tau+muon
		 mt_lmSelectionP[1]++;
		 mutau_hinvmass_lmP[1]->Fill(mtHinvmass);
	       }
	     
	   }  
       }
   
       if(bjets[3].size()>bjets[5].size()){    
	 mh_lmSelection[3]++; //med tight ph1
	 mh_hinvmass_lm[3]->Fill(mhHinvmass);
	 tau_hinvmass_lm[3]->Fill(tauHinvmass);
	 if(gentaus.size()>0){//tau+muon
	   mt_lmSelection[3]++;
	   mutau_hinvmass_lm[3]->Fill(mtHinvmass);
	 }
	 
	
	 if( (bjets[5][0].Pt()>bjets[3][1].Pt() && bjets[3][1].Pt()>30.) ||(bjets[3][0].Pt()>bjets[5][1].Pt() && bjets[5][1].Pt()>30.) || (bjets[5][0].Pt()==bjets[3][1].Pt() && bjets[3].size()>1 && bjets[3][2].Pt()> 30.) ){
	   mh_lmSelectionP[3]++;
	   mh_hinvmass_lmP[3]->Fill(mhHinvmass); //loose/tight ph1
	   tau_hinvmass_lmP[3]->Fill(tauHinvmass);
	   if(gentaus.size()>0){//tau+muon
	     mt_lmSelectionP[3]++;
	     mutau_hinvmass_lmP[3]->Fill(mtHinvmass);
	     } 
	   
	 } 
       }
      }
  //Lmt end loop


    }//basic
    
    
  }
  delete rootTree;
  delete ufs;
  delete file;
}//for loop on files

    TString typetitle[ntypes]={"Phase 1 Loose","StdGeom Loose","Phase 1 Medium","StdGeom Medium","Phase 1 Tight","StdGeom Tight"};
  cout<<"====================================== " <<endl;
  cout<<" Analyzed Events               : "<<AnalyzedEvents<<endl;
  cout<<"Lumi                    : "<<Lumi<<endl;
  cout<< "Cross-Section          : "<< CrossSection <<endl;

	for(int i=0; i<ntypes; i++){
	  cout<<"==="<< typetitle[i]<<"====" <<endl;
	 //<<<<<<< BBAOfflineAnalysis.C
	
	 cout<<" ==================================== "<<endl;
	 cout<<" Muon+Hadron         preselection      : "<<mhPreSelection[i]<<endl;
	 cout<<" Muon+Hadron    	    selection  : "<<mhSelection[i]<<endl;
	 cout<<" Muon+Hadron additional selection      : "<<mhAddSelection[i]<<endl;
	 
	 cout<<" Muon+Tauh         preselection         : "<<mtPreSelection[i]<<endl;
	 cout<<" Muon+Tauh     	    selection          : "<<mtSelection[i]<<endl;
	 cout<<" Muon+Tauh additional selection        : "<<mtAddSelection[i]<<endl;
	 
	 cout<<" Taus         preselection             : "<<mtPreSelection[i]<<endl;
	 cout<<" Taus     	    selection          : "<<mtSelection[i]<<endl;
	 cout<<" Taus additional selection             : "<<mtAddSelection[i]<<endl;
	 
	 
	}


	cout<<"=================="<<endl;
	cout<<"mh_lmPreSelection :"<<mh_lmPreSelection<<endl;

	cout<<" Tight Bjet.size Ph1 >0            : "<<nlomepre[0]<<endl;
	cout<<" Tight Bjet.size Std >0            : "<<nlomepre[1]<<endl;	
	
	cout<<" M+H Loose-Tight               :"<<mh_lmSelection[0]<<endl; 
	cout<<" M+H Loose-Tight StdGeom       :"<< mh_lmSelection[1]<<endl;
	cout<<" M+H Medium-Tight              :"<<mh_lmSelection[2]<<endl; 
	cout<<" M+H Medium-Tight StdGeom      :"<< mh_lmSelection[3]<<endl;

	cout<<" M+Tau Loose-Tight             :"<<mt_lmSelection[0]<<endl; 
	cout<<" M+Tau Loose-Tight StdGeom     :"<< mt_lmSelection[1]<<endl;
	cout<<" M+Tau Medium-Tight            :"<<mt_lmSelection[2]<<endl; 
	cout<<" M+Tau Medium-Tight StdGeom    :"<< mt_lmSelection[3]<<endl;
	
	cout<<" genTaus Loose-Tight =MH       :"<<mh_lmSelection[0]<<endl; 
	cout<<" genTaus Loose-Tight StdGeom   :"<< mh_lmSelection[1]<<endl;
	cout<<" genTaus Medium-Tight          :"<<mh_lmSelection[2]<<endl; 
	cout<<" genTaus Medium-Tight StdGeom  :"<< mh_lmSelection[3]<<endl;

	cout<<" ====================================== " <<endl;

	cout<<" M+H Loose-Tight P             :"<<mh_lmSelectionP[0]<<endl; 
	cout<<" M+H Loose-Tight StdGeom P     :"<< mh_lmSelectionP[1]<<endl;
	cout<<" M+H Medium-Tight     P        :"<<mh_lmSelectionP[2]<<endl; 
	cout<<" M+H Medium-Tight StdGeom P    :"<< mh_lmSelectionP[3]<<endl;

	cout<<" M+Tau Loose-Tight       P     :"<<mt_lmSelectionP[0]<<endl; 
	cout<<" M+Tau Loose-Tight StdGeom P   :"<< mt_lmSelectionP[1]<<endl;
	cout<<" M+Tau Medium-Tight        P   :"<<mt_lmSelectionP[2]<<endl; 
	cout<<" M+Tau Medium-Tight StdGeom P  :"<< mt_lmSelectionP[3]<<endl;
	
	cout<<" genTaus Loose-Tight         P  :"<<mh_lmSelectionP[0]<<endl; 
	cout<<" genTaus Loose-Tight StdGeom P  :"<< mh_lmSelectionP[1]<<endl;
	cout<<" genTaus Medium-Tight         P :"<<mh_lmSelectionP[2]<<endl; 
	cout<<" genTaus Medium-Tight StdGeom P :"<< mh_lmSelectionP[3]<<endl;
	
	fout.cd();
	fout.Write();
	fout.Close();

	return 0;
}
						       
