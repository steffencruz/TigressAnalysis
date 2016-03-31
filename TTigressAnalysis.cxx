#include "TTigressAnalysis.h"

#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

#include <cmath>

ClassImp(TTigressAnalysis)

std::string TTigressAnalysis::file = "";

TH3F *TTigressAnalysis::hexcgamgam = NULL;
TH2F *TTigressAnalysis::hgamgam = NULL;
TH2F *TTigressAnalysis::hexcgam = NULL;

TH1D *TTigressAnalysis::hgam = NULL;
TH1D *TTigressAnalysis::hexc = NULL;

Int_t TTigressAnalysis::gambinsz = 1;
Int_t TTigressAnalysis::excbinsz = 1;
Bool_t TTigressAnalysis::addback = true;		


TTigressAnalysis::TTigressAnalysis(std::string fname)	{	
	LoadHistos(fname.c_str());
}

TTigressAnalysis::~TTigressAnalysis()	{	}

void TTigressAnalysis::Reset(Option_t *opt){
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
	
	if(f1->IsOpen()){
		file.assign(fname);
    hexcgamgam = (TH3F*)f1->Get("ExcGamGam");  
		hexcgam 	 = (TH2F*)f1->Get("ExcGam");        
    hgamgam 	 = (TH2F*)f1->Get("GamGam");
    hgam			 = (TH1D*)f1->Get("Gam");
    hexc 			 = (TH1D*)hexcgam->ProjectionY("Exc");    		
	} else {
		file = "";
		return;
	}

}

TH1D *TTigressAnalysis::GamGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	if(!file.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);
	
  TH1D *hpx[5], *hpy[5], *hpt[5];
  TH2F *h2;
  if(exc_lo>=0.0 && exc_hi>exc_lo){ // set exc energy gate if included
  	Int_t zp[2] = {hexcgamgam->GetZaxis()->FindBin(exc_lo), hexcgamgam->GetZaxis()->FindBin(exc_hi)};
  	printf("\n Excitation energy range : %.2f-%.2f keV [ bins %i - %i ]\n",exc_lo,exc_hi,zp[0],zp[1]);
  	hexcgamgam->GetZaxis()->SetRange(zp[0],zp[1]);
  // produce an excitation energy (Z axis) gated gamma gamma matrix
  	h2 = (TH2F*)hexcgamgam->Project3D("yx");   
  } else {
  	h2 = hgamgam;
  }
  
  TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
  // project peak onto x and y axis  
  int xp[2] = {ax->FindBin(emin), ax->FindBin(emax)};
  int yp[2] = {ay->FindBin(emin), ay->FindBin(emax)};
  hpx[0] = (TH1D*)h2->ProjectionX("hpx0",yp[0],yp[1]);  
  hpy[0] = (TH1D*)h2->ProjectionY("hpy0",xp[0],xp[1]);  
  
  // get axis ranges emax project background region onto x and y
  Int_t xl[2] = {ax->FindBin(bg0), ax->FindBin(bg1)};
  Int_t xh[2] = {ax->FindBin(bg2), ax->FindBin(bg3)};
  Int_t yl[2] = {ay->FindBin(bg0), ay->FindBin(bg1)};  
  Int_t yh[2] = {ay->FindBin(bg2), ay->FindBin(bg3)};
  ////////////////////////////////////////////////////////////////////////////////////
  
  // actual bin centers
  Double_t enpeak[2]= {ax->GetBinCenter(xp[0]),ax->GetBinCenter(xp[1])};
  Double_t enleft[2]= {ax->GetBinCenter(xl[0]),ax->GetBinCenter(xl[1])};
  Double_t enrigh[2]= {ax->GetBinCenter(xh[0]),ax->GetBinCenter(xh[1])};

// determine widths of each region so that the background can be weighted appropriately
  Double_t width[3];
  width[0] = ax->GetBinCenter(xp[1])-ax->GetBinCenter(xp[0]);
  width[1] = ax->GetBinCenter(xl[1])-ax->GetBinCenter(xl[0]);
  width[2] = ax->GetBinCenter(xh[1])-ax->GetBinCenter(xh[0]);
  
  ////////////////////////////////////////////////////////////////////////////////////

	// add x and y peak 1D projections together
  hpt[0] = (TH1D*)hpx[0]->Clone("hpt0");
  hpt[0]->Add(hpy[0],1); // total peak = x projection + y projection

	// add x and y BG_LOW 1D projections together
  hpx[1] = (TH1D*)h2->ProjectionX("hpx1",yl[0],yl[1]);
  hpx[1]->Scale(0.5*width[0]/width[1]);
  
  hpy[1] = (TH1D*)h2->ProjectionY("hpy1",xl[0],xl[1]);
  hpy[1]->Scale(0.5*width[0]/width[1]);
	// total low background
  hpt[1] = (TH1D*)hpx[1]->Clone("hpt1");
  hpt[1]->Add(hpy[1],1); 
  
	// add x and y BG_HIGH 1D projections together  
  hpx[2] = (TH1D*)h2->ProjectionX("hpx2",yh[0],yh[1]);
  hpx[2]->Scale(0.5*width[0]/width[2]);
  
  hpy[2] = (TH1D*)h2->ProjectionY("hpy2",xh[0],xh[1]);
  hpy[2]->Scale(0.5*width[0]/width[2]);
	// total high background
  hpt[2] = (TH1D*)hpx[2]->Clone("hpt1");
  hpt[2]->Add(hpy[2],1);

  
	// subtract weighted backgrounds emin gated gamma gamma peak 
  hpx[3] = (TH1D*)hpx[0]->Clone("hpx3");
  hpx[3]->Add(hpx[1],-1);
  hpx[3]->Add(hpx[2],-1);
  
  hpy[3] = (TH1D*)hpy[0]->Clone("hpy3");
  hpy[3]->Add(hpy[1],-1);
  hpy[3]->Add(hpy[2],-1);

  hpt[3] = (TH1D*)hpt[0]->Clone("GammaGatedGammas");
  hpt[3]->SetTitle(Form("%s; Gamma Energy [keV]; Counts / %.0f keV",hpt[3]->GetTitle(),ax->GetBinWidth(0)));
  hpt[3]->Add(hpt[1],-1); 
  hpt[3]->Add(hpt[2],-1);
  
    
  // also add together x and y projections for total spectra
  hpx[4] = (TH1D*)hpx[1]->Clone("hpx4");
  hpx[4]->Add(hpx[2],1);
  hpy[4] = (TH1D*)hpy[1]->Clone("hpy4");
  hpy[4]->Add(hpy[2],1);
  hpt[4] = (TH1D*)hpt[1]->Clone("hpt4");
  hpt[4]->Add(hpt[2],1);
   ////////////////////////////////////////////////////////////////////////////////////
 
  for(Int_t i=0;i<4;i++){
    for(Int_t b=xp[0];b<xp[1];b++){
      hpt[i]->SetBinContent(b,hpt[i]->GetBinContent(b)*0.5);
    }
  }
  
 return hpt[3];
}

TH1D *TTigressAnalysis::ExcGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

	if(!file.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);

	TH1D *hpye[4];
  TAxis *ax = hexcgam->GetXaxis();
 
  // excitation energy doesn't need all these x and y projections as we only have one gamma axis
  hpye[0] = (TH1D*)hexcgam->ProjectionY("hpye0",ax->FindBin(emin),ax->FindBin(emax));   
  
  hpye[1] = (TH1D*)hexcgam->ProjectionY("hpye1",ax->FindBin(bg0),ax->FindBin(bg1));
  hpye[1]->Scale(0.5*(emax-emin)/(bg1-bg0)); // 0.5 is because we take half background emin below and half emin above
  
  hpye[2] = (TH1D*)hexcgam->ProjectionY("hpye2",ax->FindBin(bg2),ax->FindBin(bg3));
  hpye[2]->Scale(0.5*(emax-emin)/(bg3-bg2));
  
  hpye[3] = (TH1D*)hpye[0]->Clone("Gamma Gated Protons");
  hpye[3]->Add(hpye[1],-1);
  hpye[3]->Add(hpye[2],-1);  
  
  hpye[3]->SetLineColor(1);
  hpye[3]->GetYaxis()->SetTitleOffset(1.3);
  
  return hpye[3];

}

void TTigressAnalysis::SetBackgroundLims(Double_t emin, Double_t emax, Double_t &bg0, Double_t &bg1, Double_t &bg2, Double_t &bg3){
	// default behaviour for background
	
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
	} else if(!bg0 && !bg1 && !bg2 && !bg3){
		Double_t ewid = emax-emin;
		bg0 = emin - 0.5*ewid;
		bg1 = emin;
		bg2 = emax;
		bg3 = emax + 0.5*ewid;		
	}

	return;	
}

void TTigressAnalysis::AnalyzeGammas(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	if(!file.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return;
	}
	LoadHistos(file.c_str());
	
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
      
  TCanvas *c = new TCanvas("c","c",0,0,1200,800);
  c->Divide(2,2);
  
  ////////////////////////////////////////////////////////////////////////////////////

  c->cd(1);
  // DRAW GAMMA SINGLES (FITTED) ON PAD1  
  hgam->SetTitle("Gamma Singles");  
  hgam->Rebin(2.0/hgam->GetBinWidth(0));
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
  hexcgated->Rebin(excbinsz/hexcgated->GetBinWidth(0));
  hexcgated->GetXaxis()->SetRangeUser(0,6000);
  hexcgated->SetNameTitle("GammaGatedProtons",Form("Excitation Energy Coincident With Gated Gamma Energy; Excitation Energy [keV]; Counts / %.0f keV",hexcgated->GetBinWidth(0)));
  hexcgated->DrawCopy();    

  TBox *bpex = new TBox(exc_lo,hexcgated->GetMinimum()*1.05,exc_hi,hexcgated->GetMaximum()*1.05);
  bpex->SetFillColor(3);
  bpex->SetFillStyle(3002);  

  if(exc_lo>=0.0 && exc_hi>exc_lo)
  	  bpex->Draw();


  ////////////////////////////////////////////////////////////////////////////////////
  c->cd(2);
  TH1D *hgamgated = GamGated(emin,emax,bg0,bg1,bg2,bg3); 
  hgamgated->Rebin(gambinsz/hgamgated->GetBinWidth(0));
  hgamgated->GetXaxis()->SetRangeUser(0,1500);
  hgamgated->SetNameTitle("GammaGatedGammas",Form("Gammma Energy Coincident With Gated Gamma Energy; Gamma energy [keV]; Counts / %.0f keV",hgamgated->GetBinWidth(0)));
  hgamgated->DrawCopy();  
  
  ////////////////////////////////////////////////////////////////////////////////////  
  c->cd(4);
  TH1D *hexcgamgated = GamGated(emin,emax,bg0,bg1,bg2,bg3,exc_lo,exc_hi); 
  hexcgamgated->Rebin(gambinsz/hexcgamgated->GetBinWidth(0));
  hexcgamgated->GetXaxis()->SetRangeUser(0,1500);  
  hexcgamgated->SetLineColor(1);
  hexcgamgated->SetNameTitle("GammaExcGatedGammas",Form("Gammma Energy Coincident With Gated Gamma And Excitation Energy; Gamma energy [keV]; Counts / %.0f keV",hexcgamgated->GetBinWidth(0))); 
  hexcgamgated->DrawCopy();
    
  ////////////////////////////////////////////////////////////////////////////////////  
  
	return;
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

	fitfunc->SetParameters(0, 0, gSigma, gContent, gMean); 

	fitfunc->SetParName(0,"BgConstant");
	fitfunc->SetParName(1,"BgSlope   ");
	fitfunc->SetParName(2,"Content   "); fitfunc->SetParLimits(2,0,gContent);
	fitfunc->SetParName(3,"Mean      "); fitfunc->SetParLimits(3,emin,emax);
	fitfunc->SetParName(4,"Sigma     "); fitfunc->SetParLimits(4,0,30.);	

	fitfunc->SetParName(5,"BgLoMin   "); fitfunc->FixParameter(5,bg1);
	fitfunc->SetParName(6,"BgLoMax   "); fitfunc->FixParameter(6,emin);
	fitfunc->SetParName(7,"BgHiMin   "); fitfunc->FixParameter(7,emax);
	fitfunc->SetParName(8,"BgHiMax   "); fitfunc->FixParameter(8,bg2);         

	hist->Fit(fitfunc, "R", "SAME");

	gBgConstant = fitfunc->GetParameter(0);
	gBgSlope    = fitfunc->GetParameter(1);
	gContent    = fitfunc->GetParameter(2)/gBinW;
	gMean       = fitfunc->GetParameter(3);
	gSigma      = fitfunc->GetParameter(4);	
	gChi2pNDF   = fitfunc->GetChisquare() / fitfunc->GetNDF();

	printf("\n      Chi Square: %f\n",fitfunc->GetChisquare());
	printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc->GetParError(2));
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
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);  
  Double_t arg = (x[0] - par[3])/(sqrt2*par[4]);

  Double_t fitval = par[0] + x[0]*par[1] + par[2]/(sqrt2pi*par[4]) * exp(-arg*arg);
  
	 if((x[0]>par[5] &&  x[0]<par[6]) || (x[0]>par[7] &&  x[0]<par[8])){
		 	TF1::RejectPoint();
	 		return fitval;
	 }
	   
   return fitval;
}

