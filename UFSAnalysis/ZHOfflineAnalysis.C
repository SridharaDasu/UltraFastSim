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

int main(){

  vector<TParticle> muons;
  vector<TLorentzVector> bjets[6];


  int AnalyzedEvents=0;
  int LeptonPairPreSelection=0;
  int LeptonPairSelection=0;
  int ZMassWindow=0;
  int BJetPairPreSelection[6]={0};
  int BJetPairSelection[6]={0};
  int ZptSelection[6]={0};
  int HptSelection[6]={0};
  int DphiSelection[6]={0};
  int HMassWindow[6]={0};

  double Weight=1.;

  //================================================
  //event loop

  ifstream ConfigFile;
  ConfigFile.open("Config.txt");
  if(!ConfigFile){cout<<"Unable to read the config file";exit(1);}
  
  string OutputRootFile;
  double CrossSection;
  double Lumi;

  ConfigFile>>OutputRootFile;
  ConfigFile>>CrossSection;
  ConfigFile>>Lumi;

  //file should be created before tree to avoid memory resident trees
  TFile fout(OutputRootFile.c_str(),"recreate");

  TH1F* hptmuons  = new TH1F("Muons Pt","Muons Pt",300,0.,300.);
  TH1F* hetamuons = new TH1F("Muons Eta","Muons Eta",100,-5.,5.);
  TH1F* hphimuons = new TH1F("Muons Phi","Muons Phi",100,-7.,7.);
  TH1F* hzinvmass= new TH1F("Z invmass","Z invmass",300,0.,300.);
  TH1F* hphiz= new TH1F("Z Phi","Z Phi",100,-7.,7.);

  TH1F* hbjet_mult[6];
  hbjet_mult[0] = new TH1F("B-jets Mult. Loose StdGeom","B-jets Mult. Loose StdGeom",10,0.,10.);
  hbjet_mult[1] = new TH1F("B-jets Mult. Loose Phase 1","B-jets Mult. Loose Phase 1",10,0.,10.);
  hbjet_mult[2] = new TH1F("B-jets Mult. Medium StdGeom","B-jets Mult. Medium StdGeom",10,0.,10.);
  hbjet_mult[3] = new TH1F("B-jets Mult. Medium Phase 1","B-jets Mult. Medium Phase 1",10,0.,10.);
  hbjet_mult[4] = new TH1F("B-jets Mult. Tight StdGeom","B-jets Mult. Tight StdGeom",10,0.,10.);
  hbjet_mult[5] = new TH1F("B-jets Mult. Tight Phase 1","B-jets Mult. Tight Phase 1",10,0.,10.);

  TH1F* hetbjets[6];
  hetbjets[0] = new TH1F("B-jets Et Loose StdGeom","B-jets Et Loose StdGeom",300,0.,300.);
  hetbjets[1] = new TH1F("B-jets Et Loose Phase 1","B-jets Et Loose Phase 1",300,0.,300.);
  hetbjets[2] = new TH1F("B-jets Et Medium StdGeom","B-jets Et Medium StdGeom",300,0.,300.);
  hetbjets[3] = new TH1F("B-jets Et Medium Phase 1","B-jets Et Medium Phase 1",300,0.,300.);
  hetbjets[4] = new TH1F("B-jets Et Tight StdGeom","B-jets Et Tight StdGeom",300,0.,300.);
  hetbjets[5] = new TH1F("B-jets Et Tight Phase 1","B-jets Et Tight Phase 1",300,0.,300.);

  TH1F* hetabjets[6];
  hetabjets[0] = new TH1F("B-jets Eta Loose StdGeom","B-jets Eta Loose StdGeom",100,-5.,5.);
  hetabjets[1] = new TH1F("B-jets Eta Loose Phase 1","B-jets Eta Loose Phase 1",100,-5.,5.);
  hetabjets[2] = new TH1F("B-jets Eta Medium StdGeom","B-jets Eta Medium StdGeom",100,-5.,5.);
  hetabjets[3] = new TH1F("B-jets Eta Medium Phase 1","B-jets Eta Medium Phase 1",100,-5.,5.);
  hetabjets[4] = new TH1F("B-jets Eta Tight StdGeom","B-jets Eta Tight StdGeom",100,-5.,5.);
  hetabjets[5] = new TH1F("B-jets Eta Tight Phase 1","B-jets Eta Tight Phase 1",100,-5.,5.);

  TH1F* hphibjets[6];
  hphibjets[0] = new TH1F("B-jets Phi Loose StdGeom","B-jets Phi Loose StdGeom",100,-7.,7.);
  hphibjets[1] = new TH1F("B-jets Phi Loose Phase 1","B-jets Phi Loose Phase 1",100,-7.,7.);
  hphibjets[2] = new TH1F("B-jets Phi Medium StdGeom","B-jets Phi Medium StdGeom",100,-7.,7.);
  hphibjets[3] = new TH1F("B-jets Phi Medium Phase 1","B-jets Phi Medium Phase 1",100,-7.,7.);
  hphibjets[4] = new TH1F("B-jets Phi Tight StdGeom","B-jets Phi Tight StdGeom",100,-7.,7.);
  hphibjets[5] = new TH1F("B-jets Phi Tight Phase 1","B-jets Phi Tight Phase 1",100,-7.,7.);

  TH1F* hptz[6];
  hptz[0] = new TH1F("Z Pt Loose StdGeom","Z Pt Loose StdGeom",300,0.,300.);
  hptz[1] = new TH1F("Z Pt Loose Phase 1","Z Pt Loose Phase 1",300,0.,300.);
  hptz[2] = new TH1F("Z Pt Medium StdGeom","Z Pt Medium StdGeom",300,0.,300.);
  hptz[3] = new TH1F("Z Pt Medium Phase 1","Z Pt Medium Phase 1",300,0.,300.);
  hptz[4] = new TH1F("Z Pt Tight StdGeom","Z Pt Tight StdGeom",300,0.,300.);
  hptz[5] = new TH1F("Z Pt Tight Phase 1","Z Pt Tight Phase 1",300,0.,300.);

  TH1F* hpth[6];
  hpth[0] = new TH1F("H Pt Loose StdGeom","H Pt Loose StdGeom",300,0.,300.);
  hpth[1] = new TH1F("H Pt Loose Phase 1","H Pt Loose Phase 1",300,0.,300.);
  hpth[2] = new TH1F("H Pt Medium StdGeom","H Pt Medium StdGeom",300,0.,300.);
  hpth[3] = new TH1F("H Pt Medium Phase 1","H Pt Medium Phase 1",300,0.,300.);
  hpth[4] = new TH1F("H Pt Tight StdGeom","H Pt Tight StdGeom",300,0.,300.);
  hpth[5] = new TH1F("H Pt Tight Phase 1","H Pt Tight Phase 1",300,0.,300.);

  TH1F* hhinvmass[6];
  hhinvmass[0] = new TH1F("H invmass preselection Loose StdGeom","H invmass preselection Loose StdGeom",300,0.,300.);
  hhinvmass[1] = new TH1F("H invmass preselection Loose Phase 1","H invmass preselection Loose Phase 1",300,0.,300.);
  hhinvmass[2] = new TH1F("H invmass preselection Medium StdGeom","H invmass preselection Medium StdGeom",300,0.,300.);
  hhinvmass[3] = new TH1F("H invmass preselection Medium Phase 1","H invmass preselection Medium Phase 1",300,0.,300.);
  hhinvmass[4] = new TH1F("H invmass preselection Tight StdGeom","H invmass preselection Tight StdGeom",300,0.,300.);
  hhinvmass[5] = new TH1F("H invmass preselection Tight Phase 1","H invmass preselection Tight Phase 1",300,0.,300.);

  TH1F* hhinvmass2[6];
  hhinvmass2[0] = new TH1F("H invmass Loose StdGeom","H invmass Loose StdGeom",300,0.,300.);
  hhinvmass2[1] = new TH1F("H invmass Loose Phase 1","H invmass Loose Phase 1",300,0.,300.);
  hhinvmass2[2] = new TH1F("H invmass Medium StdGeom","H invmass Medium StdGeom",300,0.,300.);
  hhinvmass2[3] = new TH1F("H invmass Medium Phase 1","H invmass Medium Phase 1",300,0.,300.);
  hhinvmass2[4] = new TH1F("H invmass Tight StdGeom","H invmass Tight StdGeom",300,0.,300.);
  hhinvmass2[5] = new TH1F("H invmass Tight Phase 1","H invmass Tight Phase 1",300,0.,300.);

  TH1F* hphih[6];
  hphih[0] = new TH1F("H Phi Loose StdGeom","H Phi Loose StdGeom",100,-7.,7.);
  hphih[1] = new TH1F("H Phi Loose Phase 1","H Phi Loose Phase 1",100,-7.,7.);
  hphih[2] = new TH1F("H Phi Medium StdGeom","H Phi Medium StdGeom",100,-7.,7.);
  hphih[3] = new TH1F("H Phi Medium Phase 1","H Phi Medium Phase 1",100,-7.,7.);
  hphih[4] = new TH1F("H Phi Tight StdGeom","H Phi Tight StdGeom",100,-7.,7.);
  hphih[5] = new TH1F("H Phi Tight Phase 1","H Phi Tight Phase 1",100,-7.,7.);

  TH1F* hdphiZH[6];
  hdphiZH[0] = new TH1F("Dphi(Z,H) Loose StdGeom","Dphi(Z,H) Loose StdGeom",100,-7.,7.);
  hdphiZH[1] = new TH1F("Dphi(Z,H) Loose Phase 1","Dphi(Z,H) Loose Phase 1",100,-7.,7.);
  hdphiZH[2] = new TH1F("Dphi(Z,H) Medium StdGeom","Dphi(Z,H) Medium StdGeom",100,-7.,7.);
  hdphiZH[3] = new TH1F("Dphi(Z,H) Medium Phase 1","Dphi(Z,H) Medium Phase 1",100,-7.,7.);
  hdphiZH[4] = new TH1F("Dphi(Z,H) Tight StdGeom","Dphi(Z,H) Tight StdGeom",100,-7.,7.);
  hdphiZH[5] = new TH1F("Dphi(Z,H) Tight Phase 1","Dphi(Z,H) Tight Phase 1",100,-7.,7.);


  //Just to know the total statistics of the root files
  //This is needed for the event normalization
  string filename;
  vector<string> fns;//file name container
  int SumEvents=0;
  while(ConfigFile>>filename && filename!="EOF"){
    fns.push_back(filename);
    cout<<"checking "<<filename<<endl;
    TFile *file=TFile::Open(filename.c_str());
    TTree *rootTree = (TTree*)file->Get("UltraFastSim");
    int nev = int(rootTree->GetEntries());
    cout<<"It has "<<nev<<" events."<<endl;
    SumEvents+=nev;
    delete rootTree;
    delete file;
  }

  cout<<"Will run on total "<<SumEvents<<" events ... "<<endl;
  Weight=CrossSection*Lumi/SumEvents;
  cout<<"Event weight is "<<Weight<<endl;

  // ========== >>>> Loop over files starts here 
  for(vector<string>::iterator iterstr=fns.begin();iterstr!=fns.end();iterstr++){
    
    cout<<"Reading "<<*iterstr<<endl;
    TFile *file=TFile::Open((*iterstr).c_str());
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
      bjets[0]=ufs->bJetListLooseStdGeom();
      bjets[1]=ufs->bJetListLoose();
      bjets[2]=ufs->bJetListMediumStdGeom();
      bjets[3]=ufs->bJetListMedium();
      bjets[4]=ufs->bJetListTightStdGeom();
      bjets[5]=ufs->bJetListTight();

      for(vector<TParticle>::iterator itm=muons.begin();itm!=muons.end();itm++){
	hptmuons->Fill(itm->Pt(),Weight); 
	hetamuons->Fill(itm->Eta(),Weight); 
	hphimuons->Fill(itm->Phi(),Weight); 
      }
      for(int i=0;i<6;i++){
	for(vector<TLorentzVector>::iterator itj=bjets[i].begin();itj!=bjets[i].end();itj++){
	  hetbjets[i]->Fill(itj->Pt(),Weight); 
	  hetabjets[i]->Fill(itj->Eta(),Weight); 
	  hphibjets[i]->Fill(itj->Phi(),Weight); 
	}
      }
      
      double Zinvmass=0;
      double Zpt=0;
      double Zphi=0;

      double Hinvmass[6]={0};
      double Hpt[6]={0};
      double Hpx[6]={0};
      double Hpy[6]={0};
      double Hphi[6]={0};
      double DphiZH[6]={0};
      
      if(muons.size()>1){
	
	LeptonPairPreSelection++;
	
	if(muons[0].Pt()>20. && muons[1].Pt()>20.){
	  
	  LeptonPairSelection++;
	  
	  for(int i=0;i<6;i++)hbjet_mult[i]->Fill(bjets[i].size());

	  Zinvmass=(pow((muons[0].Energy()+muons[1].Energy()),2)-
		    pow((muons[0].Px()+muons[1].Px()),2)-
		    pow((muons[0].Py()+muons[1].Py()),2)-
		    pow((muons[0].Pz()+muons[1].Pz()),2));
	  if(Zinvmass>0)Zinvmass=sqrt(Zinvmass);
	  else Zinvmass=0;
	  Zpt=sqrt(pow((muons[0].Px()+muons[1].Px()),2)+pow((muons[0].Py()+muons[1].Py()),2));
	  
	  double Zpx=0;
	  double Zpy=0;
	  Zpx=muons[0].Px()+muons[1].Px();
	  Zpy=muons[0].Py()+muons[1].Py();
	  Zphi=atan(Zpy/Zpx);
	  Zphi=CorrectPhi(Zphi,Zpx,Zpy);
	  
	  hzinvmass->Fill(Zinvmass,Weight);
	  hphiz->Fill(Zphi,Weight);      
	  
	  if(Zinvmass>70. && Zinvmass<110.){
	    
	    ZMassWindow++;
	    
	    //=======>>>>> loop over b-jet algo, and tk. geometry
	    for(int i=0;i<6;i++){
	      
	      if(bjets[i].size()>1){
		
		BJetPairPreSelection[i]++;
		
		if(bjets[i][0].Pt()>30. && bjets[i][1].Pt()>30.){
		  
		  BJetPairSelection[i]++;
		  
		  Hinvmass[i]=(pow((bjets[i][0].Energy()+bjets[i][1].Energy()),2)-
			       pow((bjets[i][0].Px()+bjets[i][1].Px()),2)-
			       pow((bjets[i][0].Py()+bjets[i][1].Py()),2)-
			       pow((bjets[i][0].Pz()+bjets[i][1].Pz()),2));
		  if(Hinvmass[i]>0)Hinvmass[i]=sqrt(Hinvmass[i]);
		  else Hinvmass[i]=0;
		  
		  Hpt[i]=sqrt(pow((bjets[i][0].Px()+bjets[i][1].Px()),2)
			      +pow((bjets[i][0].Py()+bjets[i][1].Py()),2));
		  
		  Hpx[i]=bjets[i][0].Px()+bjets[i][1].Px();
		  Hpy[i]=bjets[i][0].Py()+bjets[i][1].Py();
		  Hphi[i]=atan(Hpy[i]/Hpx[i]);
		  Hphi[i]=CorrectPhi(Hphi[i],Hpx[i],Hpy[i]);
		  
		  hhinvmass[i]->Fill(Hinvmass[i],Weight); 
		  hphih[i]->Fill(Hphi[i],Weight);
		  
		  DphiZH[i]=fabs(Hphi[i]-Zphi);
		  if(DphiZH[i]>M_PI)DphiZH[i] = (2.*M_PI-DphiZH[i]); 
		  
		  hptz[i]->Fill(Zpt,Weight);
		  if(Zpt>150.){
		    
		    ZptSelection[i]++;
		    
		    hpth[i]->Fill(Hpt[i],Weight);
		    if(Hpt[i]>150.){
		      
		      HptSelection[i]++;
		      
		      hdphiZH[i]->Fill(DphiZH[i],Weight);
		      if(DphiZH[i]>2.5){
			
			DphiSelection[i]++;
			
			hhinvmass2[i]->Fill(Hinvmass[i],Weight); 
			if(Hinvmass[i]<140. && Hinvmass[i]>100.){
			  
			  HMassWindow[i]++;
			  
			}// H mass window
		      }//dphi cut
		    }//h pt cut
		  }//z pt cut
		}//bjet et > 30 GeV
	      }//bjet pair
	    }//loop over b-tagging algo and tk geometries
	  }//Z mass window
	}//muon pt>20 GeV
      }//muon.size==2
    }//for loop on events
    delete rootTree;
    delete ufs;
    delete file;
  }//for loop on files
  
  string algoname[6];
  algoname[0]="Loose b-tagging, Tk. Std. Geom.";
  algoname[1]="Loose b-tagging, Tk. Phase 1";
  algoname[2]="Medium b-tagging, Tk. Std. Geom.";
  algoname[3]="Medium b-tagging, Tk. Phase 1";
  algoname[4]="Tight b-tagging, Tk. Std. Geom.";
  algoname[5]="Tight b-tagging, Tk. Phase 1";

  for(int i=0;i<6;i++){
    cout<<" ====================================== " <<endl;
    cout<<" Algo : "<<algoname[i]<<endl;
    cout<<" Analyzed Events                  : "<<AnalyzedEvents<<endl;
    cout<<" ====================================== " <<endl;
    cout<<" Lepton pair preselection : "<<LeptonPairPreSelection<<endl;
    cout<<" Lepton pair    selection : "<<LeptonPairSelection<<endl;
    cout<<" Z Invariant Mass Window  : "<<ZMassWindow<<endl;
    cout<<" B-jet  pair preselection : "<<BJetPairPreSelection[i]<<endl;
    cout<<" B-jet  pair    selection : "<<BJetPairSelection[i]<<endl;
    cout<<" Z pt                     : "<<ZptSelection[i]<<endl;
    cout<<" H pt                     : "<<HptSelection[i]<<endl;
    cout<<" Dphi(H,Z)                : "<<DphiSelection[i]<<endl;
    cout<<" H Invariant Mass Window  : "<<HMassWindow[i]<<endl;
    cout<<" Total efficiency     (%) : "<<100*double(HMassWindow[i])/AnalyzedEvents<<endl;
    cout<<" Number of Events         : "<<double(HMassWindow[i])/AnalyzedEvents*CrossSection*Lumi<<endl;
    cout<<" ====================================== " <<endl;
  }
  fout.cd();
  fout.Write();
  fout.Close();

  return 0;
}
