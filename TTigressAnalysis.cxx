#include "TTigressAnalysis.h"

#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
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

Int_t TTigressAnalysis::binwidth = 1;
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
	
	binwidth = 1;
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

TH1D *TTigressAnalysis::GammasGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	if(!file.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
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
  
  // get axis ranges to project background region onto x and y
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
  hpt[3]->SetTitle(Form("Gam_ExcGated; Gamma Energy [keV]; Counts / %.0f keV",hpt[3]->GetTitle(),hpt[3]->GetXaxis()->GetBinWidth(0)));
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

