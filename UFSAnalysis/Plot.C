#include <string>

void Plot()
{

  char fName[1024];

  std::string zhName("ZHmmbb-v5-plots.root");
  std::string zzName("ZZmmbb-v5-plots.root");
  std::string zjName("Zmm-v5-plots.root");

  // Make plots

  TFile zh(zhName.c_str());
  TFile zj(zjName.c_str());
  TFile zz(zzName.c_str());

  TH1F *zhmuonpt=zh.Get("Muons Pt");
  TH1F *zjmuonpt=zj.Get("Muons Pt");
  TH1F *zzmuonpt=zz.Get("Muons Pt");

  zjmuonpt->SetLineColor(kGreen);
  zzmuonpt->SetLineColor(kOrange+3);
  zhmuonpt->SetLineColor(kRed);

  zjmuonpt->Scale(1./zjmuonpt->Integral());
  zzmuonpt->Scale(1./zzmuonpt->Integral());
  zhmuonpt->Scale(1./zhmuonpt->Integral());

  TCanvas c1;

  zhmuonpt->Draw();
  zzmuonpt->Draw("same");
  zjmuonpt->Draw("same");

  TLegend *lg=new TLegend(0.5,0.5,0.7,0.7);
  lg->SetFillColor(kWhite);
  lg->AddEntry(zhmuonpt,"ZH, m(H)=120 GeV","l");
  lg->AddEntry(zjmuonpt,"Z+jets","l");
  lg->AddEntry(zzmuonpt,"ZZ","l");
  lg->Draw();

  c1.SaveAs("muonPt.png");;

  /////////////////////////////

  TH1F *zhzinvmass=zh.Get("Z invmass");
  TH1F *zjzinvmass=zj.Get("Z invmass");
  TH1F *zzzinvmass=zz.Get("Z invmass");

  zjzinvmass->SetLineColor(kGreen);
  zzzinvmass->SetLineColor(kOrange+3);
  zhzinvmass->SetLineColor(kRed);

  zjzinvmass->Scale(1./zjzinvmass->Integral());
  zzzinvmass->Scale(1./zzzinvmass->Integral());
  zhzinvmass->Scale(1./zhzinvmass->Integral());

  TCanvas c5;

  zhzinvmass->Draw();
  zzzinvmass->Draw("same");
  zjzinvmass->Draw("same");

  TLegend *lg5=new TLegend(0.5,0.5,0.7,0.7);
  lg5->SetFillColor(kWhite);
  lg5->AddEntry(zhzinvmass,"ZH, m(H)=120 GeV","l");
  lg5->AddEntry(zjzinvmass,"Z+jets","l");
  lg5->AddEntry(zzzinvmass,"ZZ","l");
  lg5->Draw();

  c5.SaveAs("ZMass.png");;

  /////////////////////////////

  vector<string> geom;
  geom.push_back("Phase 1");
  geom.push_back("StdGeom");

  vector<string> bJetOP;
  bJetOP.push_back("Loose");
  bJetOP.push_back("Medium");
  bJetOP.push_back("Tight");

  vector<string> hNames;
  hNames.push_back("B-jets Et ");
  hNames.push_back("B-jets Eta ");
  hNames.push_back("B-jets Phi ");
  hNames.push_back("B-jets Mult. ");
  hNames.push_back("H invmass preselection ");
  hNames.push_back("H Phi ");
  hNames.push_back("Z Pt ");
  hNames.push_back("H Pt ");
  hNames.push_back("Dphi(Z,H) ");
  hNames.push_back("H invmass ");
  
  for(int i = 0; i < geom.size(); i++) {
    for(int j = 0; j < bJetOP.size(); j++) {
      for(int k = 0; k < hNames.size(); k++) {
	
	cout << hNames[k] << endl;

	string hName(hNames[k].c_str());
	hName += bJetOP[j];
	hName += " ";
	hName += geom[i];

	cout << hName << endl;
	
	TH1F *zhHist=zh.Get(hName.c_str());
	TH1F *zjHist=zj.Get(hName.c_str());
	TH1F *zzHist=zz.Get(hName.c_str());
	
	if(zhHist == 0 || zjHist == 0 || zzHist == 0) {
	  cout << "Failed to find " << hName << endl;
	  break;
	}

	if(i > 7) {
	  zhHist.Rebin(10);
	  zjHist.Rebin(10);
	  zzHist.Rebin(10);
	}
	
	zjHist->SetLineColor(kGreen);
	zzHist->SetLineColor(kOrange);
	zhHist->SetLineColor(kRed);
	//zjHist->SetFillColor(kYellow);
	zzHist->SetFillColor(kBlue);
	zhHist->SetFillColor(kBlack);
	
	// zjHist->Scale(1./zjHist->Integral());
	// zzHist->Scale(1./zzHist->Integral());
	// zhHist->Scale(1./zhHist->Integral());
	
	TCanvas *c = new TCanvas();
	
	zhHist->Draw();
	zzHist->Draw("same");
	zjHist->Draw("same");
	
	TLegend *lg2=new TLegend(0.5,0.5,0.7,0.7);
	lg2->SetFillColor(kWhite);
	lg2->AddEntry(zhHist,"ZH, m(H)=120 GeV","l");
	lg2->AddEntry(zjHist,"Z+jets","l");
	lg2->AddEntry(zzHist,"ZZ","l");
	lg2->Draw();
	
	strcpy(fName, hName.c_str());
	for(int l = 0; l < strlen(fName); l++) {
	  if(fName[l] == ' ') fName[l] = '-';
	}
	strcat(fName, ".png");
	c->SaveAs(fName);

	TCanvas *cStack = new TCanvas();
	THStack *stack = new THStack("Stack", hName.c_str());
	stack->Add(zjHist);
	stack->Add(zzHist);
	stack->Add(zhHist);
	gPad->SetLogy(1);
	stack->Draw();
	lg2->Draw();
	fName[strlen(fName)-4] = 0;
	strcat(fName, "-stack.png");
	cStack->SaveAs(fName);
	gPad->SetLogy(0);

      }
    }
  }
  
}
