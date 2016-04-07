#include "TTigressAnalysis.h"

#include <TFile.h>
#include <TLegend.h>
#include<TLegendEntry.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TArrow.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>

#include<fstream>
#include <cmath>
#include <algorithm>    // std::count
#include <vector>       

ClassImp(TTigressAnalysis)

std::string TTigressAnalysis::histfile = "";
std::string TTigressAnalysis::nndcfile = "";

TH3S *TTigressAnalysis::hexcgamgam = NULL;
TH3S *TTigressAnalysis::hexcthcmgam = NULL;

TH2F *TTigressAnalysis::hgamgam = NULL;
TH2F *TTigressAnalysis::hexcgam = NULL;
TH1D *TTigressAnalysis::hgam = NULL;
TH1D *TTigressAnalysis::hexc = NULL;

Int_t TTigressAnalysis::gambinsz = 4;
Int_t TTigressAnalysis::excbinsz = 40;
Bool_t TTigressAnalysis::addback = true;		

////////////////////////////////////////////////////////////////////////////////

Bool_t TTigressAnalysis::verbose = true;
TList *TTigressAnalysis::list = 0;
std::vector<double> TTigressAnalysis::energy;
std::vector<int> TTigressAnalysis::states;
Double_t TTigressAnalysis::ExcSig = 200.0;
Double_t TTigressAnalysis::NSig = 1.5;


TH1D *TTigressAnalysis::hint = NULL;
TH2F *TTigressAnalysis::htrans = NULL;
TH2F *TTigressAnalysis::hgams = NULL;
TH2F *TTigressAnalysis::hseq = NULL;
Int_t TTigressAnalysis::nlines = 0;
Int_t TTigressAnalysis::nsequences = 0;
Int_t TTigressAnalysis::nstates = 0;
TF1 *TTigressAnalysis::TigSigma = NULL;
TF1 *TTigressAnalysis::TigEfficiency = NULL;





TTigressAnalysis::TTigressAnalysis(const char *fname)	{	
	LoadHistos(fname); // read histograms
	
	InitGammasLevels(); // load nndc stuff
}

TTigressAnalysis::~TTigressAnalysis()	{	}

void TTigressAnalysis::Print(Option_t *opt) {

	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
	printf("\n\n\t____TTigressAnalysis____");

	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n\n");
		
}

void TTigressAnalysis::Clear(Option_t *opt) {

	histfile.clear();
	hexcgamgam = NULL;
	hexcgam = NULL;
	hgamgam = NULL;
	hgam = NULL;		
	hexc = NULL;			
	
	gambinsz = 1;
	excbinsz = 1;	
	addback = true;
}

void TTigressAnalysis::LoadHistos(const char *fname){
	
	TFile *f1 = new TFile(fname,"READ");
	printf("\n\t Reading Histos File ' %s ' :\n",fname);
	
	if(f1->IsOpen()){
		histfile.assign(fname);
    hexcgamgam  = (TH3S*)f1->Get("ExcGamGam_dp"); 
    hexcthcmgam = (TH3S*)f1->Get("ExcGamThetaCmSmooth_dp"); // fixed!
		hexcgam 	  = (TH2F*)f1->Get("ExcGam_dp");       
    hgamgam 	  = (TH2F*)f1->Get("GamGam_dp");
    hgam			  = (TH1D*)f1->Get("Gam_dp");
    hexc 			  = (TH1D*)f1->Get("Exc_dp");  
    
    printf("\n\t Loaded Histograms :-\n");
    if(hexcgamgam) hexcgamgam->Print();  		
    if(hexcthcmgam) hexcthcmgam->Print();  		    
    if(hexcgam) hexcgam->Print();  		
    if(hgamgam) hgamgam->Print();  		
    if(hgam) hgam->Print();  		
    if(hexc) hexc->Print();  		
	} else {
		histfile = "";
		return;
	}

}

TH1D *TTigressAnalysis::Gammas(Double_t emin, Double_t emax){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
		
	TH1D *h = (TH1D*)hgam->Clone("Gammas");
	h->Rebin(gambinsz/h->GetBinWidth(0));	
  //h->SetTitle(Form("%s; Gamma Energy [keV]; Counts / %i keV",hgam->GetTitle(),gambinsz));	
	h->GetXaxis()->SetRangeUser(emin,emax);
	return h;	
}

TH1D *TTigressAnalysis::GamGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);
	
  TH2F *h2;
  if(exc_lo>=0.0 && exc_hi>exc_lo){ // set exc energy gate if included
  	Int_t zp[2] = {hexcgamgam->GetZaxis()->FindBin(exc_lo), hexcgamgam->GetZaxis()->FindBin(exc_hi)};
  	printf("\n Excitation energy range : %.2f-%.2f keV [ bins %i - %i ]\n",exc_lo,exc_hi,zp[0],zp[1]);
  	hexcgamgam->GetZaxis()->SetRange(zp[0],zp[1]);
  // produce an excitation energy (Z axis) gated gamma gamma matrix
  	h2 = (TH2F*)hexcgamgam->Project3D("yx");   
  } else 
  	h2 = hgamgam;
  
  TH1D *hpx[4], *hpy[4], *hpt[4];    
	Double_t peakszx, peakszy, bgloszx, bgloszy, bghiszx, bghiszy;
  // peak x & y & tot
  hpx[0] = TH2Proj(h2,'x',emin,emax,peakszy); 
  hpy[0] = TH2Proj(h2,'y',emin,emax,peakszx);  
  hpt[0] = TH1Sum(hpx[0],hpy[0]);
  // bglo x & y & tot
  hpx[1] = TH2Proj(h2,'x',bg0,bg1,bgloszy);
  hpy[1] = TH2Proj(h2,'y',bg0,bg1,bgloszx);  
  hpt[1] = TH1Sum(hpx[1],hpy[1],peakszy/bgloszy,peakszx/bgloszx);
  //hpt[1]->Scale(0.5);
  // bghi x & y & tot  
  hpx[2] = TH2Proj(h2,'x',bg2,bg3,bghiszy);
  hpy[2] = TH2Proj(h2,'y',bg2,bg3,bghiszx);  
  hpt[2] = TH1Sum(hpx[2],hpy[2],peakszy/bghiszy,peakszx/bghiszx);
 // hpt[1]->Scale(0.5);
    
	// bgtot  
  hpt[3] = TH1Sum(hpt[1],hpt[2]);
  // peak - bgtot
  TH1D *hg = TH1Sum(hpt[0],hpt[3],1,-0.5);
  hg->Scale(0.5);
	hg->Rebin(gambinsz/hg->GetBinWidth(0));
  hg->SetNameTitle("GatedGammas",Form("Gammma Energy Coincident With Gated Gam %sEnergy; Gamma energy [keV]; Counts / %i keV",exc_hi>0?"And Exc ":"",gambinsz));
      
	return hg;
}	


TH1D *TTigressAnalysis::ExcGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}

	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);

	TH1D *hp[4];
	TH2F *h2 = hexcgam;
	Double_t peaksz, bglosz, bghisz;
  // peak
  hp[0] = TH2Proj(h2,'y',emin,emax,peaksz);	
	// bglo  
  hp[1] = TH2Proj(h2,'y',bg0,bg1,bglosz);		
	// bghi  
  hp[2] = TH2Proj(h2,'y',bg2,bg3,bghisz);	
  // bgtot	
  hp[3] = TH1Sum(hp[1],hp[2],peaksz/bglosz,peaksz/bghisz);
  
  // peak - bgtot
  TH1D *he = TH1Sum(hp[0],hp[3],1,-0.5);
       	
	he->SetLineColor(1);
	he->GetYaxis()->SetTitleOffset(1.3);
	he->Rebin(excbinsz/he->GetBinWidth(0));
  he->SetNameTitle("GatedProtons",Form("Excitation Energy Coincident With Gated Gamma Energy; Excitation Energy [keV]; Counts / %i keV",excbinsz));

	return he;
}

TH2F *TTigressAnalysis::ExcThetaGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	// so far it always uses dp in cm frame
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);

	TH2F *hp[4];
	TH3S *h3 = hexcthcmgam;
	// x = thcm
	// y = gam
	// z = exc

	Double_t peaksz, bglosz, bghisz;
  // peak
  hp[0] = TH3Proj(h3,"zx",emin,emax,peaksz);	 
	// bglo  
  hp[1] = TH3Proj(h3,"zx",bg0,bg1,bglosz);		
	// bghi  
  hp[2] = TH3Proj(h3,"zx",bg2,bg3,bghisz);	
  // bgtot	
  hp[3] = TH2Sum(hp[1],hp[2],peaksz/bglosz,peaksz/bghisz);
  
  // peak - bgtot
  TH2F *het = TH2Sum(hp[0],hp[3],1,-0.5);
 	
	het->RebinY(excbinsz/het->GetYaxis()->GetBinWidth(0));
  het->SetNameTitle("ExcThetaCm",Form("Excitation Energy Versus Theta Cm; Theta Cm [Deg]; Excitation Energy [keV]"));

	return het;
}


void TTigressAnalysis::SetBackgroundLims(Double_t emin, Double_t emax, Double_t &bg0, Double_t &bg1, Double_t &bg2, Double_t &bg3){
	// default behaviour for setting background region
	
	if(bg0 && bg1 && bg2 && bg3){
		return; // all background regions have been set
	}else if(bg0 && bg1 && !bg2 && !bg3){ // only use lower background
		bg2 = emax;
		bg3 = emax;
		printf("\n\t Error :  A background must be specified above and below the peak !\n\n");
	} else if(!bg0 && !bg1 && bg2 && bg3){ // only use upper background
		bg0 = emin;
		bg1 = emin;
		printf("\n\t Error :  A background must be specified above and below the peak !\n\n");
	} else {
		Double_t ewid = emax-emin;
		bg0 = emin - 0.5*ewid;
		bg1 = emin;
		bg2 = emax;
		bg3 = emax + 0.5*ewid;		
	}
	
	return;	
}

TCanvas *TTigressAnalysis::AnalyzeGammas(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	TCanvas *c = 0;
	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return c;
	}
	
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
      
  c = new TCanvas("GammaAnalysis","Gamma Coincidence Analysis",0,0,1200,800);
  c->Divide(2,2);
  
  ////////////////////////////////////////////////////////////////////////////////////

  c->cd(1);
  // DRAW GAMMA SINGLES (FITTED) ON PAD1  
  hgam->SetTitle("Gamma Singles");  
  hgam->GetXaxis()->SetRangeUser(bg0-20,bg3+80);
  
  hgam->GetYaxis()->SetTitle(Form("Counts / %.0f keV",hgam->GetBinWidth(0)));
  hgam->GetYaxis()->SetTitleOffset(1.3);
  hgam->DrawCopy();
  
  // selected peak and background region boxes
  TBox *bp = new TBox(emax,hgam->GetMaximum()*.01,emin,hgam->GetMaximum()*1.05);
  bp->SetFillColor(3);
  bp->SetFillStyle(3002); 
  TBox *b1 = new TBox(bg0,hgam->GetMaximum()*.01,bg1,hgam->GetMaximum()*1.05);
  b1->SetFillColor(2);
  b1->SetFillStyle(3002);
  TBox *b2 = new TBox(bg2,hgam->GetMaximum()*.01,bg3,hgam->GetMaximum()*1.05);
  b2->SetFillColor(2);
  b2->SetFillStyle(3002);   
    
  bp->Draw();
  b1->Draw();
  b2->Draw();

  FitPeakExcludeRange(hgam,emin,emax,bg0,bg1,bg2,bg3);

  TPaveText *bstat = new TPaveText(bg3+10,hgam->GetMinimum()*1.4,bg3+100,hgam->GetMaximum()*1.05);
  bstat->SetLineColor(1);
  bstat->SetShadowColor(0);
  bstat->AddText("Fit Results :");
  bstat->AddLine(.0,.85,1.,.85);
  TF1 *func = hgam->GetFunction("gauss_linbg_exc");
  for(int fn=0; fn<5; fn++)
  	bstat->AddText(Form("%s : %6.2f +/- %6.2f",func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn)));
 
 // bstat->AddText(Form("Chi2/DOF : %8.2f",func->GetChisquare()/func->GetNDF()));
  bstat->Draw();
 	  
////////////////////////////////////////////////////////////////////////////////////
  c->cd(3);

  TH1D *hexcgated = ExcGated(emin,emax,bg0,bg1,bg2,bg3); 
  hexcgated->GetXaxis()->SetRangeUser(0,6000);
  hexcgated->DrawCopy();    

  TBox *bpex = new TBox(exc_lo,hexcgated->GetMinimum()*1.05,exc_hi,hexcgated->GetMaximum()*1.05);
  bpex->SetFillColor(3);
  bpex->SetFillStyle(3002);  

  if(exc_lo>=0.0 && exc_hi>exc_lo)
  	  bpex->Draw();

  ////////////////////////////////////////////////////////////////////////////////////
  c->cd(2);
  TH1D *hgamgated = GamGated(emin,emax,bg0,bg1,bg2,bg3); 
  hgamgated->GetXaxis()->SetRangeUser(0,1500);
  hgamgated->DrawCopy();  
  
  ////////////////////////////////////////////////////////////////////////////////////  
  c->cd(4);
	if(exc_lo>=0.0 && exc_hi>exc_lo){
		TH1D *hexcgamgated = GamGated(emin,emax,bg0,bg1,bg2,bg3,exc_lo,exc_hi); 
		hexcgamgated->GetXaxis()->SetRangeUser(0,1500);  
		hexcgamgated->DrawCopy();
	} else {
		TH2F *hexcthcm = ExcThetaGated(emin,emax,bg0,bg1,bg2,bg3); 
		hexcthcm->DrawCopy("colz");
	}
  ////////////////////////////////////////////////////////////////////////////////////  
  
	return c;
}

void TTigressAnalysis::FitPeakExcludeRange(TH1 *hist, Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){ 

	Double_t gBgConstant, gBgSlope, gContent, gMean, gSigma, gBinW, gChi2pNDF;

	TF1 *fitfunc = new TF1("gauss_linbg_exc",gaus_lbg_exc, bg0, bg3, 9);
	fitfunc->SetLineColor(3);
	fitfunc->SetLineWidth(1);
	// *** Obtaining and specifying the start values for the fit ***
	gContent = hist->Integral(hist->FindBin(emin),hist->FindBin(emax)); 
	gMean    = 0.5 * ( emax + emin);  
	gSigma   = 0.3 * ( emax - emin); 
	gBinW 		= hist->GetBinWidth(1);
	//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);

	fitfunc->SetParameters(0, 0, gContent, gMean, gSigma); 

	fitfunc->SetParName(0,"BgConstant");
	fitfunc->SetParName(1,"BgSlope   ");
	fitfunc->SetParName(2,"Content   "); fitfunc->SetParLimits(2,0,gBinW*gContent);
	fitfunc->SetParName(3,"Mean      "); fitfunc->SetParLimits(3,emin,emax);
	fitfunc->SetParName(4,"Sigma     "); fitfunc->SetParLimits(4,0,30.);	

	fitfunc->SetParName(5,"ExcLo Min "); fitfunc->FixParameter(5,bg1);
	fitfunc->SetParName(6,"ExcLo Max "); fitfunc->FixParameter(6,emin);
	fitfunc->SetParName(7,"ExcHi Min "); fitfunc->FixParameter(7,emax);
	fitfunc->SetParName(8,"ExcHi Man "); fitfunc->FixParameter(8,bg2);         

	hist->Fit(fitfunc,"RQ","SAME");

	gBgConstant = fitfunc->GetParameter(0);
	gBgSlope    = fitfunc->GetParameter(1);
	gContent    = fitfunc->GetParameter(2)/gBinW;
	gMean       = fitfunc->GetParameter(3);
	gSigma      = fitfunc->GetParameter(4);	
	gChi2pNDF   = fitfunc->GetChisquare() / fitfunc->GetNDF();

}

Double_t TTigressAnalysis::gaus_lbg_exc(Double_t *x, Double_t *par){
/*
	// excludes the region between backgrounds and peak
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
  par[5]	 bglo max point
  par[6]	 peak min point
	par[7]	 peak max point
  par[8]	 bghi min point      
*/	 
  static Double_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);  
  Double_t arg = (x[0] - par[3])/(sqrt2*par[4]);

  Double_t fitval = par[0] + x[0]*par[1] + par[2]/(sqrt2pi*par[4]) * exp(-arg*arg);
  
	 if((x[0]>par[5] &&  x[0]<par[6]) || (x[0]>par[7] &&  x[0]<par[8])){
		 	TF1::RejectPoint();
	 		return fitval;
	 }
	   
   return fitval;
}

TH1D *TTigressAnalysis::TH1Sum(TH1D *ha, TH1D *hb, Double_t sca, Double_t scb){
	TH1D *h = (TH1D*) ha->Clone(Form("%s_%s",ha->GetName(),hb->GetName()));
	h->Scale(sca);
	h->Add(hb,scb);
	return h;
}

TH2F *TTigressAnalysis::TH2Sum(TH2F *ha, TH2F *hb, Double_t sca, Double_t scb){
	TH2F *h = (TH2F*) ha->Clone(Form("%s_%s",ha->GetName(),hb->GetName()));
	h->Scale(sca);
	h->Add(hb,scb);
	return h;
}

TH1D *TTigressAnalysis::TH2Proj(TH2F *h, char ax, Double_t minval, Double_t maxval, Double_t &sz){

	TAxis *axis;
	if(ax=='x')
		axis =	h->GetYaxis();
	else
		axis =	h->GetXaxis();
	
  int b[2] = {axis->FindBin(minval), axis->FindBin(maxval)};
	double c[2] = {axis->GetBinCenter(b[0]), axis->GetBinCenter(b[1])};
	printf("\n\t vals[] = {%.1f, %.1f}  b[] = {%i, %i}, c[] = {%.1f, %.1f}",minval,maxval,b[0],b[1],c[0],c[1]);
	// actual gate width
	sz = c[1] - c[0];
	
	if(ax=='x')
		return (TH1D*) h->ProjectionX(Form("%s_px_%iTo%i",h->GetName(),b[0],b[1]),b[0],b[1]);
	else 
		return (TH1D*) h->ProjectionY(Form("%s_py_%iTo%i",h->GetName(),b[0],b[1]),b[0],b[1]);
}

TH2F *TTigressAnalysis::TH3Proj(TH3S *h, std::string str, Double_t minval, Double_t maxval, Double_t &sz){

	TAxis *axis;	
	if(str.find('x')==std::string::npos)
		axis = h->GetXaxis();
	else if(str.find('y')==std::string::npos)
		axis = h->GetYaxis();
	else if(str.find('z')==std::string::npos)
		axis = h->GetZaxis();
					
  int b[2] = {axis->FindBin(minval), axis->FindBin(maxval)};
	double c[2] = {axis->GetBinCenter(b[0]), axis->GetBinCenter(b[1])};
//	printf("\n\tstr = %s .. vals[] = {%.1f, %.1f}  b[] = {%i, %i}, c[] = {%.1f, %.1f}",str.c_str(),minval,maxval,b[0],b[1],c[0],c[1]);

	axis->SetRange(b[0],b[1]);
	// actual gate width
	sz = c[1] - c[0];
	
	TH2F *h2 = (TH2F*) h->Project3D(str.c_str());
	
	TH2F *h2safe = (TH2F*)h2->Clone("SafeCopy");
	h2safe->SetName(Form("%s_py_%iTo%i",h->GetName(),b[0],b[1]));
	return h2safe;
}

////////////////////////////////////////////////////////////////////////////////
// DEMO - FOR PETER
// 		RunMe()
// 		PrintStates()
// 		PrintStates(1200,2200)
//		DrawDecayMats(6,false)
//		DrawDecayMats(6,true)
// 		DrawGammas(6,true)
//		DrawGammasExcRange(1200,1800)
//		DrawGammasExcRange(1200,1800,300)
//		DrawLevelsFromGamma(414,1200,1800,300)
//		DrawGammasExcRangeGamGate(414,1200,1800,300)
// 		ShowData(1200,1800);
//		ShowDataGamGate --! NOT YET IMPLEMENTED

// NOTES
// DONE = determine addback efficiency curve from Eu152 and apply it to THEORY
// DONE = excitation energy & gamma peak width increase with energy. get these calibrations.

// Gate on exc energy range produce gamma spectrum. Compare to data and eliminate what we dont see experimentally
// normalize theory to specific peak (direct to ground if possible)

// read out and read in a table of intensities instead of assuming a uniform distribution
// Make a function : SetStateIntensity(int state, double intensity)
// RemoveState(eng)!!

// first point is to remove levels which are clearly not populated
// for ShowData second plot show chisquare between data and theory

// extract population intensities from gamma gated spectra too.
// coincidence efficiency needs to be known so that theory is decreased accordingly
// use 1229 keV [414 & 815 peaks] for gam-gam coincidence efficiency?

// ACCOUNTS FOR THE FACT THAT THE REAL DATA HAS GAMMAS FROM STATES OUTSIDE OF THE SELECTED RANGE
// BY INCLUDING ALL STATES THAT ARE UP TO ±1.5 SIGMA BEYOND THE SPECIFIED EXC RANGE AND
// WEIGHTING THEIR GAMMA SPECTRA BY THE COUNTS EXPECTED WITHIN THE EXC RANGE [MIN. ~7%]

// IDENTIFY IF STATES ARE *ABSENT*, 0 OR 1.
// INCLUDE GAMMA GAMMA COINCIDENCES
// 

////////////////////////////////////////////////////////////////////////////////

void TTigressAnalysis::InitGammasLevels(Int_t nmax, Bool_t verb){
	gStyle->SetOptStat(0);
	gStyle->SetTitleXOffset(1.3);
	gStyle->SetTitleYOffset(1.5);	
	
	SetVerbose(verb);
	LoadLevelsGammas("Sr96_LevelsGammas.txt",nmax);
	
	// with appropriate peak widths we can identify multiple adjacent peaks	
	TigSigma = new TF1("func_sig","[0]+[1]*TMath::Power(x,[2])",0,4000); // power law?
	TigSigma->SetParameters(-4.74,0.679,0.409); // from fitting singles data	
	TigSigma->SetNpx(4000);
	TigSigma->SetTitle("TIGRESS Resolution Curve; Energy [keV]; Peak Sigma [keV]");
		
	TigEfficiency = new TF1("func_eff","[0]*pow(10.,[1]*log10(x)+[2]*pow(log10(x),2.)+[3]*pow(1/x,2.))",0,4000);
	TigEfficiency->SetParameters(1.77634e+06,0.394454,-0.145119,-848.751);	
	TigEfficiency->SetParameter(0,TigEfficiency->GetParameter(0)/TigEfficiency->GetMaximum());
	TigEfficiency->SetNpx(4000);
	TigEfficiency->SetTitle("TIGRESS Efficiency Curve; Energy [keV]; Relative Efficiency"); 
//	ShowData();

}

void TTigressAnalysis::LoadLevelsGammas(std::string fname, int nmax){

	ClearVars();
	
	std::ifstream infile(fname.c_str());	
	if(!infile.is_open()){
		printf("\nFAILED TO OPEN FILE :  %s\n\n",fname.c_str());
		return;
	}
	
	nndcfile.assign(fname);
	
	nlines = std::count(std::istreambuf_iterator<char>(infile),std::istreambuf_iterator<char>(), '\n')+1;
	infile.close();
	infile.open(fname.c_str());	
	
	TH2F *htmp = new TH2F("TmpMat","TmpMat",nlines,0,nlines,nlines,0,nlines);
	TH2F *htmp2 = new TH2F("TmpMat2","TmpMat2",nlines,0,nlines,10000,0,10000);

	Double_t initial, final, intensity, egamma, maxgam=0;
	int statenum;
	nstates = 0;
	printf("\n Reading Input File %s:-\n",fname.c_str());
	if(verbose) printf("\n\t   State Energy\tGamma Energy\tIntensity[%%]\tFinal State [num]");
	while(infile.good()){
		
		infile >> initial >> egamma >> intensity >> final;

//		if(intensity<intensity_min)
//			continue;	
		
		GetIndex(initial,true); // add new state, but don't duplicate		

		statenum = GetIndex(final)+1;
		
		if(verbose) printf("\t\t%5.1f\t%8.1f\t%9.1f\t%7.1f   [%i]\n",initial,egamma,intensity,final,statenum);
		
		if(statenum==0){ // final state doesn't match any state
			printf("\n\t! Error :  Final state %.1f not valid !\n\n",final);
			continue;
		}
				
		htmp->SetBinContent(nstates,statenum,intensity);	// fill intensity-transition matrix	
		if(maxgam<egamma)
			maxgam = egamma;
			
		htmp2->SetBinContent(nstates,egamma,intensity);  // fill gammas histogram
		
		if(nstates==nmax){
			printf("\n..............................................................\n\t\tSTOPPED READING [reached %i states] ",nstates);			
			break;
		}		
	}
	printf("\n\n\t --- Complete! [eof = %i  fail = %i  bad = %i]\n\n",infile.eof(),infile.fail(),infile.bad());
	
	infile.close();
	
	TObject *obj = gROOT->FindObjectAny("TransitionMat");
	if(obj)
		obj->Delete();
	htrans = new TH2F("TransitionMat","",nstates,0,nstates,nstates,0,nstates);
	htrans->SetTitle("Transition Intensity Matrix; Initial State; Final State");
	
	for(int i=1; i<=nstates; i++) // go to the maximum row
		for(int j=1; j<=i; j++) // go up to the diagonal
			htrans->SetBinContent(i,j,htmp->GetBinContent(i,j));
	
	htrans->Scale(1./htrans->GetMaximum()); // normalize to max strength 1

	
	obj = gROOT->FindObjectAny("GammasMat");
	if(obj)
		obj->Delete();
	hgams = new TH2F("GammasMat","",nstates,0,nstates,(Int_t)maxgam+100,0,maxgam+100);
	hgams->SetTitle("Gammas Emitted Matrix; State; Gamma Energy [keV]");

	for(int i=1; i<=nstates; i++){ // go to the maximum row
		TH1D *hprojy = htmp2->ProjectionY("",i,i);
		for(int j=hprojy->FindFirstBinAbove(); j<=hprojy->FindLastBinAbove(); j++){ // go up to the diagonal
			if(htmp2->GetBinContent(i,j))
				hgams->SetBinContent(i,j,htmp2->GetBinContent(i,j));
		}
	}	
	hgams->Scale(1./hgams->GetMaximum()); // normalize to max strength 1
		

	PrintStates();

	return;
}

std::vector<int> TTigressAnalysis::PrintStates(Double_t emin, Double_t emax){

	std::vector<int> vals;
	if(!energy.size()){
		printf("\n\n No States Have Been Loaded !!\n\n");
		return vals;
	}
	printf("\n\n List of States:-");

	for(int i=0; i<(Int_t)energy.size(); i++){
		if(energy.at(i)<emin)
			continue;
					
		if(energy.at(i)>emax)
			break;
						
		vals.push_back(i+1);	
		printf("\n\tState %2i =   %6.1f keV",i+1,energy.at(i));
	}
	printf("\n\n");
	return vals;
}

Int_t TTigressAnalysis::PrintCascades(Int_t from_state, Bool_t printeng){

	BuildDecayScheme(from_state);

	printf("\n\n List of Decays from state %2i %s:-",from_state,printeng?"{With Gamma Energy}":"");
	TH1D *htmp;	
	for(int j=1; j<=nsequences; j++){
	
		htmp = hseq->ProjectionX("",j,j);
		printf("\n\tCascade %2i = ",j);
		
		if(printeng) 
			printf(" %6.1f ",energy.at((Int_t)htmp->GetBinContent(1)-1));
		else 
			printf("  %2i",(Int_t)htmp->GetBinContent(1));
			
		for(int i=2; i<=htmp->FindLastBinAbove(); i++){
			if(printeng) {
				printf(" -{ %6.1f }->  %6.1f ",energy.at((Int_t)htmp->GetBinContent(i-1)-1) - 
				energy.at((Int_t)htmp->GetBinContent(i)-1),	energy.at((Int_t)htmp->GetBinContent(i)-1));
			} else
				printf("  -> %2i",(Int_t)htmp->GetBinContent(i));
		}
		if(printeng)printf(" keV * * *");
		
	}

	printf("\n\n");
	
	return nsequences;

}

void TTigressAnalysis::ClearVars(){
	
	nndcfile.clear();
	list = new TList();
	hseq = 0;
	hint = 0;
	nsequences = 0;

	states.clear();	
	energy.clear();
	htrans = 0;	
	hgams  = 0;
	nlines = 0;
	nstates = 0;	

}


void TTigressAnalysis::FixIntensities(Double_t emin, Double_t emax, Int_t level, Double_t strength){

	static TH1D *hdata=0, *htmp, *hbg;
	static THStack *hstack;
	static TList *hlist;
	static Double_t eminval=0.0, emaxval=-10000.0;
	if(eminval!=emin || emaxval!=emax){
		hdata=0;
		eminval = emin;
		emaxval = emax;
	}
	static std::vector<double> intensities, total, str;
	static std::vector<int> substates;
	static TCanvas *c;	
	c = (TCanvas*) gROOT->FindObject("StackCanvas");

	if(!hdata || !c){
	
		if(!hint)
			SetIntensity(0,100.0);

		total.clear();
		intensities.clear();	
		substates.clear();
		str.clear();	
		
	// get theory spectrum with uniform population of all states	
	// this	spectrum is made up of delta function peaks
		TH1D *h = DrawGammasExcRange(emin,emax);
		if(!h){
			printf("\n\t! Error : Check energy range [emin > emax]\n\n");
			return;
		}		
						
		SetVerbose(false);
		if(!htrans){
			printf("\n\t! Error : No intensity matrix\n\n");
			return;
		}
		const char *fname = "Results_GamGam_GoodStrips.root";
		TFile *f = new TFile(fname,"READ");
		if(!f){
			printf("\n\t! Error : File '%s' could not be read\n\n",fname);
			return;
		}

		const char *hname = "ExcGam";	
		TH2F *hexcgam = (TH2F*)f->Get(hname);
		if(!hexcgam){
			printf("\n\t! Error : Histogram '%s' could not be read\n\n",hname);
			return;
		}
		
		Int_t binlo = hexcgam->GetYaxis()->FindBin(emin);
		Int_t binhi = hexcgam->GetYaxis()->FindBin(emax);	
		// use excitation energy gate on gammas
		TH1D *hpeak = hexcgam->ProjectionX("peak",binlo,binhi);
		/*
		Int_t gatesize = fabs(0.5*(binhi-binlo));
		TH1D *hlo = hexcgam->ProjectionX("bglo",binlo-1-gatesize,binlo-1);
		TH1D *hhi = hexcgam->ProjectionX("bghi",binhi+1,binhi+1+gatesize);
		// subtract background		
		hpeak->Print();
		hlo->Print();	
		hhi->Print();			
		
		hpeak->Add(hlo,-1.0);
		hpeak->Add(hhi,-1.0);
		*/
		
		hdata = new TH1D("GammaSinglesData","",4000,0,4000);
		hdata->GetYaxis()->SetTitleOffset(1.5);
		hdata->SetLineColorAlpha(kBlack,0.3);
		
		// I WANT THIS CORRECTED ELSEWHERE -- NOT IN THIS CODE
		Int_t binshift = -3; // correct the energies so they line up with theory	
		for(int i=1; i<hpeak->GetNbinsX()-binshift; i++)
			hdata->SetBinContent(i,hpeak->GetBinContent(i+binshift));
		
 		hdata->Rebin();			// REBINNING DATA 
		//hdata->Add(hdata->ShowBackground(20,"0N"),-1.0);
		hbg = (TH1D*)hdata->ShowBackground(22,"0N");
		hbg->SetLineColorAlpha(kRed,0.5);
		hdata->SetTitle(Form("Gamma Singles From Data With Exc Range %.1f - %.1f keV; Gamma Energy [keV]; Counts / %.0f keV",emin,emax,hdata->GetBinWidth(0)));
		
		TH1D *htmp2;
		hlist = new TList;	
		std::string s;		
		
		// make realistic copies of all the spectra
		for(int i=0; i<(Int_t)states.size(); i++){
			htmp = (TH1D*)list->FindObject(Form("GammaSpectrum_State%i_Strength",states.at(i)));
			if(!htmp)
				continue;
			
			s.assign(htmp->GetTitle());
			str.push_back(atof(s.substr(s.length()-6,5).c_str()));
									
			htmp2 = MakeRealistic(htmp);
			htmp2->SetLineColor(htmp->GetLineColor());
			htmp2->SetLineColorAlpha(htmp2->GetLineColor(),str.back()/100.0);
		//	htmp2->SetFillStyle(0);
			htmp2->SetLineWidth(2);
			
 			htmp2->Rebin(htmp2->GetNbinsX()/hdata->GetNbinsX());
						
			hlist->Add(htmp2);
			substates.push_back(states.at(i)); // state numbers
			total.push_back(htmp2->Integral()); // total counts
		//	SetIntensity(states.at(i),200.0); // fill histogram
		}			
		intensities.assign(substates.size(),200.0); // intensities
		
				
		c = new TCanvas("StackCanvas","",1100,500);		
		c->Divide(2,1);
		c->cd(1);		
	}
		
	if(strength<0)
		printf("\n\tError :  Negative State Population [%.1f] Is Not Allowed \n",strength);
	
	
	Int_t indx = level-substates.at(0);	
	if(indx>=(Int_t)substates.size()){
		printf("\n\tError :  State '%i' out of range [%i - %i]\n",level,substates.at(0),substates.back());
	} else {
	  printf("\n\tSET :  State %i [indx %i] Strength = %.1f",level,indx,strength);
		//intensities.at(indx) = strength;		// modify the scale factor for one spectrum
		SetIntensity(level,strength);
	}

	TLegend *leg = new TLegend(0.5,0.5,1.0,0.92);
	leg->SetFillStyle(1001);

	hstack = new THStack(Form("RealisticGammaStacked_%.1fTo%.1f",emin,emax),"");
	hstack->SetTitle("Gamma Spectrum From Each State; Gamma Energy [keV]; Intensity [arb]");
	hstack->Add(hbg);
	
	Double_t strng;
	for(int i=0; i<(Int_t)substates.size(); i++){
	
		htmp = (TH1D*)hlist->At(i);
	//	printf("\n state %i :  total [%.1f]\t hist [%.1f]\t integral [%.1f]\t scale [%.1f]",substates.at(i),total.at(i),hint->GetBinContent(substates.at(i)),htmp->Integral(),hint->GetBinContent(substates.at(i))*total.at(i)/htmp->Integral());
		strng=hint->GetBinContent(substates.at(i));
		if(strng)
			htmp->Scale(strng*total.at(i)/htmp->Integral());
		else
			htmp->Scale(0);
			
		hstack->Add(htmp);

		leg->AddEntry(htmp,Form("State %2i [%5.1f keV @ %2.f%%]: Strength %.1f",substates.at(i),energy.at(substates.at(i)-1),str.at(i),hint->GetBinContent(substates.at(i))),"l");	
	}		
	printf("\n\n");
	
	c->cd(1);
	hstack->Draw("hist "); // noclear?
	hdata->Draw("same");		
	
	leg->Draw();
	c->cd(2);
	hint->GetXaxis()->SetRangeUser(substates.at(0)-2,substates.back()+1);
	hint->Draw();

}

void TTigressAnalysis::SetIntensity(Int_t state, Double_t strength){

	if(!hint){	
		hint = new TH1D("StateIntensites","",nstates,0,nstates);
		hint->SetTitle("State Intensities Set Using Data; State Number; Populated Strength");
		for(int i=1; i<=nstates; i++)
			hint->SetBinContent(i,strength);
	} else if(!state){ // fill all bins with given strength
		for(int i=1; i<=nstates; i++)
			hint->SetBinContent(i,strength);
	}
	
	hint->SetBinContent(state,strength);
	return;
}	

void TTigressAnalysis::ReadIntensities(const char *fname){

	std::ifstream infile(fname);
	if(!infile.is_open()){
		printf("\n\t! Error : Could not locate file ' %s '\n\n",fname);
		return;	
	}
	
	SetIntensity(0,0);// initiate the histogram
	Int_t bin;
	Double_t eng, val;
	while(infile.good()){
		infile >> bin >> eng >> val;
		hint->SetBinContent(bin,val);
	}
	printf("\n\n\tRead in file ' %s ' and set TH1D*'hint' intensities.\n\n",fname);
}

void TTigressAnalysis::WriteIntensities(const char *fname){

	if(!hint){
		printf("\n\t Error :  At least one intensity must be set!");
		return;
	}

	std::ofstream ofile(fname);
	
	for(int i=1; i<hint->GetNbinsX(); i++)
		ofile << i << "\t" << energy.at(i-1) << "\t" << hint->GetBinContent(i) << "\n";

	printf("\n\n\tWrote Contents of TH1D*'hint' to file!\n\n");
	ofile.close();
}

TH1D *TTigressAnalysis::DrawGammas(Int_t from_state, Bool_t draw){
	
	TH1D *hgam = (TH1D*)list->FindObject(Form("GammaSpectrum_State%i",from_state));
	if(hgam){
		printf("\n %s already exists. Using this histogram.\n",hgam->GetName());
		return hgam;		
	}
	
	printf("\n* * * * * * * * * * * * * * * * * * * * * * ");	
	printf("\nDrawing gammas from state %i :-\n",from_state);	
				
	if(!BuildDecayScheme(from_state))
		return 0;
	
	hgam = new TH1D(Form("GammaSpectrum_State%i",from_state),"",4000.0,0,4000.0);
	hgam->SetTitle(Form("Gammas Emitted From %.1f keV State [%i]; Gamma Energy [keV]; Intensity [# of gammas per state];",energy.at(from_state-1),from_state));
	hgam->SetLineColor(kRed);
	hgam->SetLineWidth(2);
	hgam->GetYaxis()->SetTitleOffset(1.5);
	if(list->FindObject(hgam->GetName()))
		list->Remove(list->FindObject(hgam->GetName()));
	list->Add(hgam);
	
	Double_t egam, strength, intensity;
	TH1D *htmp, *htmp2;
	Int_t initial, final; 
	
	for(int j=1; j<=nsequences; j++){
	
		htmp = hseq->ProjectionX("",j,j);
		intensity = 1.0;
		if(verbose) printf("\n Drawing Cascade %2i :-",j);
		for(int i=1; i<htmp->FindLastBinAbove(); i++){

			initial = (Int_t)htmp->GetBinContent(i);
			final = (Int_t)htmp->GetBinContent(i+1);
			if(verbose) printf("\n\tDecay %2i -> %2i : ",initial,final);
			
			// remember that vector index = bin number -1
			egam = energy.at(initial-1)-energy.at(final-1);	
			
			// strength of transition [expressed as a ratio to 100%]
			strength = htrans->GetBinContent(initial,final);
			// project out all transition strengths from this state so that we can sum over all strengths
			htmp2 = htrans->ProjectionY("transtmp",initial,initial);		
			// sum all transition strengths to normalize and convert this transition to a probability
			strength /= htmp2->Integral(1,initial);
			
			intensity*=strength;
			if(verbose) printf(" Gam. Eng. = %4.1f \t Int. = %4.3f \t Tot. Int. = %4.3f",egam,strength,intensity);
			
			hgam->Fill(egam,intensity);
		}
	}
	if(verbose) printf("\n\t --- Complete!\n\n");
	//gPad->SetLogy();
	
	if(draw){
		TCanvas *c = new TCanvas("GammaSpec");
		TH1D *h = (TH1D*)hgam->Clone(Form("%s_Clone",htmp->GetName()));				
		h->Draw();
	}
//	hgam = 0;
		
//	h->Scale(1.0/h->Integral()); // nomalize!
	
	return hgam;
}

TH1D *TTigressAnalysis::DrawGammasExcRange(Double_t emin, Double_t emax){
	
	
	TH1D *h = (TH1D*)list->FindObject(Form("GammaSpectrumFromMultiStates_%.1fTo%.1f",emin,emax));
	if(h){
		printf("\n %s already exists. Using this histogram.\n",h->GetName());
		return h;		
	}
	
	
// excitation energy resolution means that states outside of emin-emax range will still be 
// in this window so the theory should do the same in order to reproduce what we see.	
// Include states which are outside of emin-emax range, but only counts their contribution in this range	
// a peak at ±1.5 sigma from the emin-emax range limits would contribute a maximum of ~6.7%
	Double_t emin_ext=emin, emax_ext=emax;	
	Bool_t extend_excrng = true;	
	if(extend_excrng){ 
		emin_ext -= NSig*ExcSig;
		emax_ext += NSig*ExcSig;
	}
	

	printf("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");	
	printf("\n\nDrawing gammas from states in%s range %.1f - %.1f keV :-",extend_excrng?Form(" *EXTENDED [+/- %.1f keV]*",NSig*ExcSig):"",emin_ext,emax_ext);	
		
	states = PrintStates(emin_ext,emax_ext);
	if(!states.size()){
		printf("\n\t! Error : No states in this range\n\n");
		return 0;
	}

	TH1D *htmp;
	THStack *hstack = new THStack(Form("GammaStacked_%.1fTo%.1f%s",emin,emax,extend_excrng?"_Extended":""),"");
	hstack->SetTitle("Gamma Spectrum From Each State; Gamma Energy [keV]; Intensity [# of gammas per state]");
	
	h = new TH1D(Form("GammaSpectrumFromMultiStates_%.1fTo%.1f%s",emin,emax,extend_excrng?"_Extended":""),
							Form("Gammas Emitted From %5.1f-%5.1f keV State [%i-%i]",energy.at(states.at(0)-1),energy.at(states.back()-1),states.at(0),states.back())
							,4000.0,0,4000.0);
	h->SetTitle(Form("%s; Gamma Energy [keV]; Intensity [# of gammas per state];",h->GetTitle()));
	h->SetLineColor(kRed);
	h->SetLineWidth(2);
	
	
	TF1 *func = new TF1("gaus","gaus",0,8000.0);
	func->SetParameters(1.0/(sqrt(2*3.1415926)*ExcSig),0.0,ExcSig);		
	
	Double_t strength;
	TLegend *leg = new TLegend(0.6,0.4,0.98,0.92);
	leg->SetFillStyle(1001);
	
	TH1D *hgam;
	TObject *obj;	
	Int_t n=1;
	for(int i=0; i<(Int_t)states.size(); i++){
		hgam = DrawGammas(states.at(i));
		if(!hgam)
			continue;

		// calculate the contribution from each state is its integral intensity over excitation gate range 
		// take each state to be a normalized gaussian 
		func->SetParameter(1,energy.at(states.at(i)-1));					
		// rescale this spectrum to account for expected contribution in emin-emax range
		strength = func->Integral(emin,emax);
		
		obj = list->FindObject(Form("%s_Strength",hgam->GetName()));
		if(obj)
			list->Remove(obj);
			
		htmp = (TH1D*) hgam->Clone(Form("%s_Strength",hgam->GetName()));
		htmp->Scale(strength);
		list->Add(htmp);

		// keep track of the strength 
		htmp->SetTitle(Form("%s @ Strength %4.1f%%",htmp->GetTitle(),strength*100));
		
		// all low [out of range] state gammas are black and all high [o.o.r.] gammas are blue
		if(energy.at(states.at(i)-1)<emin)
			htmp->SetLineColorAlpha(kBlack,strength);		
		else if(energy.at(states.at(i)-1)>emax)
			htmp->SetLineColorAlpha(kBlue,strength);		
		else{
		 	if(n==10)
		 		n=20;
			htmp->SetLineColor(n++);
		}	
		
		h->Add(htmp);
		hstack->Add(htmp);
	//	leg->AddEntry(htmp,Form("State %2i, %6.1f keV [%2.f%%]",states.at(i),energy.at(states.at(i)-1),strength*100),"lp");		
		leg->AddEntry(htmp,Form("State %2i [%5.1f keV @ %2.f%%]",states.at(i),energy.at(states.at(i)-1),strength*100.0),"l");
	}
	
	obj = list->FindObject(h->GetName());
	if(obj)
		list->Remove(obj);
	list->Add(h);	

	obj = list->FindObject(hstack->GetName());
	if(obj)
		list->Remove(obj);
	list->Add(hstack);	

	TCanvas *obj2 = (TCanvas*)gROOT->FindObject("CanvasEnergyRange");
	if(obj2)
		obj2->Close();
		
	TCanvas *c = new TCanvas("CanvasEnergyRange","Canvas",1000,500);
	c->Divide(2,1);
	c->cd(1);
	THStack *hstack2 = (THStack*)hstack->Clone(Form("%s_Clone",hstack->GetName()));	
	hstack2->Draw();
	hstack2->GetYaxis()->SetTitleOffset(1.5);	
	leg->Draw();
	
	c->cd(2);
	TH1D *hh = (TH1D*)h->Clone(Form("%s_Clone",h->GetName()));	
	hh->DrawCopy();
	hh->GetYaxis()->SetTitleOffset(1.5);	
	
	return h;
}

TCanvas *TTigressAnalysis::DrawDecayMats(Int_t from_state, Bool_t engaxis){
	
	BuildDecayScheme(from_state);

	TCanvas *c = new TCanvas("Canvas","Canvas",1100,500);
	c->Divide(2,1);
	c->cd(1);

	TH2F *hint2 = (TH2F*)htrans->Clone(Form("%s_Copy",htrans->GetName()));
	hint2->GetXaxis()->SetNdivisions(from_state);
	hint2->GetXaxis()->SetRangeUser(0,from_state);
	hint2->GetYaxis()->SetNdivisions(from_state);
	hint2->GetYaxis()->SetRangeUser(0,from_state);	
	gPad->SetGrid();
	hint2->Draw("colz");
	
	TF1 *ff = new TF1("pol1","pol1",0,from_state);
	ff->SetParameters(0,1);
	ff->SetLineColor(1);
	ff->SetLineWidth(1);	
	ff->Draw("same");
	
	c->cd(2);
	//DrawGammas(from_state)->Draw();
	THStack *hdecay = GetDecayScheme(from_state,engaxis);
	hdecay->Draw("nostack"); // first draw is necessary
	hdecay->SetTitle(Form("Cascades From State %i; Cascade Number; State %s",from_state,engaxis?"Energy [keV]":"Number"));
	hdecay->GetYaxis()->SetTitleOffset(1.5);
	hdecay->GetXaxis()->SetNdivisions(nsequences);
	gPad->SetGridx();	
	if(!engaxis)hdecay->GetYaxis()->SetNdivisions(from_state+1);	
	hdecay->GetYaxis()->SetRangeUser(1,from_state);
	gPad->Modified();
	gPad->Update();
	
	DrawTransitions(from_state,engaxis);
	
	return c;
}

int TTigressAnalysis::GetIndex(Double_t val, bool add_element){ // compares a value to contents of vector to see if it exists as an element

// std::count(energy.begin(),energy.end(),energy.size())

	for(int i=0; i<(Int_t)energy.size(); i++){
		if(fabs(val-energy.at(i))<1.0) // within 1 keV means same state
			return i;
	}
	
	// if element hasn't been found add it in
 	if(add_element){
 		energy.push_back(val);
		nstates++; 		
 		if(verbose) printf("\n State %2i :  \r",nstates);
 		return energy.size()-1;
 	}
	return -1;
}

void TTigressAnalysis::AddSeq(std::vector<int> states){

	if(!hseq){
		TObject *obj = gROOT->FindObjectAny("CascadeMat");
		if(obj)
			obj->Delete();	
		hseq = new TH2F("CascadeMat","",nstates,0,nstates,50,0,50);
		hseq->SetTitle("Cascade To Ground State; Series of States; Cascade Number");	
		hseq->SetLineWidth(2);
	}
	nsequences++;
	
	if(verbose)	printf("\n\t* Added Cascade %3i :  [ ",nsequences);
	for(int i=0; i<(Int_t)states.size(); i++){
		if(verbose) printf( "%i ",states.at(i));
		hseq->SetBinContent(i+1,nsequences,states.at(i));
		hseq->SetBinError(i+1,nsequences,0.0001);
	}
	if(verbose) printf("]");
		
	return;
}

int TTigressAnalysis::BuildDecayScheme(int from_state){

	hseq = 0;
	nsequences = 0;	
	
	if(!htrans){
		printf("\n\t! Error : Read Input File First !\n\n");
		return 0;
	} else if(from_state>nstates){
		printf("\n\t! Error: Maximum allowed state = % i!\n\n",nstates);
		return 0;
	}
		
	std::vector<int>s; // index to show the current column on each row
	
	s.push_back(from_state);
	int m = from_state, i=0;
	
	if(verbose) printf("\n Building Cascade from state %i..",from_state);
	for(i=1; i<=m; i++){
			
		if(i==from_state)
			break;
		
		if(i==m){ // end of loop for row m, so go back to previous row
			
//			printf("\n \t*[ %i , %i ]*\t end of sequence :  ",m,i);		
//			for(int j=0;j<s.size(); j++) printf(" %i ",s.at(j)); 
			
			i = s.at(s.size()-1); // s[end] is previous column [the -1 is immediately incremented to zero..]
			m	= s.at(s.size()-2); // s[end-1] is previous row
			s.pop_back(); // remove last element
			
//			printf(" -> ");
//			for(int j=0;j<s.size(); j++) printf(" %i ",s.at(j)); 		
//			printf("\n\tProceeding to [ %i , %i] ...\n",m,i+1);
			
			continue;
		}
					
		//printf("\n [ %i , %i ] =  %.0f\t sequence length = %lu",m,i,htrans->GetBinContent(m,i),s.size());		
					
		if(htrans->GetBinContent(m,i)>0){
			s.push_back(i);
			if(i==1) // ground state has been reached
				AddSeq(s);
			m = i; // set current row
			i = 0; // restart loop on next row	(0 will be incremented immediately to 1)		
			continue;	
		}
				
	}
	if(verbose) printf("\n\t --- Complete!\n");
	
	return nsequences;
}

THStack *TTigressAnalysis::GetDecayScheme(Int_t from_state, Bool_t engaxis){
	
	if(!hseq)
		if(!BuildDecayScheme(from_state))
			return 0;			

	TH2F *hseq2 = ShrinkToFit(from_state,hseq);
	
	TH1D *hproj;
	THStack *hdecay = new THStack(Form("DecayScheme_State%i",from_state),Form("Cascades From State %i",from_state));
	
	Int_t binmax = hseq2->ProjectionX()->FindLastBinAbove();
	for(int i=1; i<=binmax; i++){
		hproj = hseq2->ProjectionY(Form("Step%i",i),i,i);
		
		if(engaxis) ConvertToEnergyAxis(hproj);
		
		if(i==1){
			hproj->SetLineWidth(4);		
			hproj->SetLineColor(kBlue);
		}
		else{
			hproj->SetLineWidth(2);
			hproj->SetLineColor(kBlack);
		}

		hdecay->Add(hproj);
	}

	return hdecay;
}

void TTigressAnalysis::ConvertToEnergyAxis(TH1D *h){

	for(int j=1; j<=h->FindLastBinAbove();j++){
		if(h->GetBinContent(j)==0)
			continue;
		h->SetBinContent(j,energy.at((Int_t)h->GetBinContent(j)-1));
	}
	return;
}

TH2F *TTigressAnalysis::ShrinkToFit(int from_state, TH2F *h){

	Int_t maxbinx = h->ProjectionX()->FindLastBinAbove();
	Int_t maxbiny = h->ProjectionY()->FindLastBinAbove();	

	TH2F *h2 = new TH2F(Form("%s_State%i",h->GetName(),from_state),"",maxbinx,0,maxbinx,maxbiny,0,maxbiny);
	h2->SetTitle(Form("%s; %s; %s",h->GetTitle(),h->GetXaxis()->GetTitle(),h->GetYaxis()->GetTitle()));
	h2->SetLineColor(h->GetLineColor());
	h2->SetLineWidth(h->GetLineWidth());
	
	for(int ix=1; ix<=maxbinx; ix++){
		for(int iy=1; iy<=maxbiny; iy++){
			h2->SetBinContent(ix,iy,h->GetBinContent(ix,iy));	
			h2->SetBinError(ix,iy,h->GetBinError(ix,iy));	
		}
	}
	return h2;
}

void TTigressAnalysis::DrawTransitions(int from_state, Bool_t engaxis){
		
	TArrow *arrow;
	TPaveText *pt;
	TH2F *hseq2 = ShrinkToFit(from_state,hseq);
	TH1D *h;
	Double_t dh=0.1, ypos;
	if(engaxis)
		dh=25;	
	
	for(int i=1; i<=nsequences; i++){
		h = hseq2->ProjectionX("",i,i);
		if(engaxis) ConvertToEnergyAxis(h);		
			
		for(int j=1; j<=h->FindLastBinAbove(); j++){
			if(!engaxis && h->GetBinContent(j+1)==0)
				continue;
			arrow = new TArrow(i-0.5,h->GetBinContent(j)-dh,i-0.5,h->GetBinContent(j+1)+dh,0.01,"|>");
			arrow->SetAngle(40);
			arrow->SetFillColor(kRed);
			arrow->SetFillStyle(3001);
			arrow->SetLineColorAlpha(kRed,0.25);
			arrow->Draw();
			if(engaxis){
				ypos = 0.5*(h->GetBinContent(j)+h->GetBinContent(j+1));
				pt = new TPaveText(i-0.75,ypos-dh,i-0.25,ypos+dh,"nb");
				pt->AddText(Form("%6.1f keV",h->GetBinContent(j)-h->GetBinContent(j+1)));
				pt->SetFillColor(0);
				pt->SetTextColor(1);
//				pt->SetTextSize(1);
				pt->Draw();
			}
		}
	}
	
	return;
}

TH1D *TTigressAnalysis::MakeRealistic(TH1D *hgam){

	// converts delta function spectrum into sigma(E) and efficiency corrected spectrum
	TH1D *htmp = new TH1D(Form("%s_Realistic",hgam->GetName()),hgam->GetTitle(),hgam->GetNbinsX(),0,4000);		
	htmp->GetXaxis()->SetTitle(hgam->GetXaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(hgam->GetYaxis()->GetTitle());	
	
	Double_t x, val, sig, eff;
	TF1 *func = new TF1("gauspeak","gaus",0,4000);
	func->SetNpx(4000);
		
	// integrate the counts under the data peaks 
	for(int i=1; i<hgam->GetNbinsX(); i++){
		
		val = hgam->GetBinContent(i);
		if(!val)
			continue;
		
		// define a gaussian centered at bin center with sigma = sig and integral = hcalc->GetBinContent(i)		
		x = hgam->GetBinCenter(i);
		sig = TigSigma->Eval(x);
		eff = TigEfficiency->Eval(x);		
		
		func->SetParameters(eff*val/(sqrt(2*3.14159)*sig),x,sig);
		// now get this histogram and use it as a filter over the data
		htmp->Add((TH1D*)func->GetHistogram());
	}
	
	return htmp;
}

