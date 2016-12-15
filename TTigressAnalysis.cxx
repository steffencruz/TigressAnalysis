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
#include <TVirtualFitter.h>

#include<fstream>
#include <cmath>
#include <algorithm>    // std::count
#include <vector>     

#include <TTigress.h>  


ClassImp(TTigressAnalysis)

std::string TTigressAnalysis::histfile = "";
std::string TTigressAnalysis::reaction = "";
std::string TTigressAnalysis::nndcfile = "";
std::string TTigressAnalysis::detsel   = "";

TH3S *TTigressAnalysis::hegg                = NULL;
TH3S *TTigressAnalysis::hexcgamgam          = NULL;
TH3S *TTigressAnalysis::hexcgamgam_sel[3]   = {NULL,NULL,NULL};

TH3S *TTigressAnalysis::hegt                = NULL;
TH3S *TTigressAnalysis::hexcgamthtig        = NULL;
TH3S *TTigressAnalysis::hexcgamthtig_sel[3] = {NULL,NULL,NULL};

TH3S *TTigressAnalysis::hexcthcmgam         = NULL;

TH2F *TTigressAnalysis::heg                 = NULL;
TH2F *TTigressAnalysis::hexcgam             = NULL;
TH2F *TTigressAnalysis::hexcgam_sel[3]      = {NULL,NULL,NULL};

TH2F *TTigressAnalysis::hgg                 = NULL;
TH2F *TTigressAnalysis::hgamgam             = NULL;
TH2F *TTigressAnalysis::hgamgam_sel[3]      = {NULL,NULL,NULL};

TH1D *TTigressAnalysis::he                  = NULL;
TH1D *TTigressAnalysis::hexc                = NULL;
TH1D *TTigressAnalysis::hexc_sel[3]         = {NULL,NULL,NULL};

TH1D *TTigressAnalysis::hg                  = NULL;
TH1D *TTigressAnalysis::hgam                = NULL;
TH1D *TTigressAnalysis::hgam_sel[3]         = {NULL,NULL,NULL};


Double_t TTigressAnalysis::gambinsz         = 1.0;
Double_t TTigressAnalysis::excbinsz         = 40.0;
Bool_t TTigressAnalysis::addback            = true;		

////////////////////////////////////////////////////////////////////////////////

Bool_t TTigressAnalysis::verbose = true;
TList *TTigressAnalysis::list = 0;
std::vector<double> TTigressAnalysis::energies;
std::vector<int> TTigressAnalysis::states;
std::map<int,std::vector<int> > TTigressAnalysis::cascades;

Double_t TTigressAnalysis::ExcSig = 180.0;
Double_t TTigressAnalysis::NSig = 1.5;


TH1D *TTigressAnalysis::hint = NULL;
TH2F *TTigressAnalysis::htrans = NULL;
TH2F *TTigressAnalysis::hgams = NULL;
TH2F *TTigressAnalysis::hseq = NULL;
Int_t TTigressAnalysis::nlines = 0;
Int_t TTigressAnalysis::nsequences = 0;
Int_t TTigressAnalysis::ncascades = 0;
Int_t TTigressAnalysis::nstates = 0;
TF1 *TTigressAnalysis::fTigSigma = NULL;
TF1 *TTigressAnalysis::fTigEff = NULL;
TF1 *TTigressAnalysis::fTigEffLow = NULL;
TF1 *TTigressAnalysis::fTigEffUp = NULL;
TGraphErrors *TTigressAnalysis::gTigEff = NULL;



TTigressAnalysis::TTigressAnalysis()	{	
}

TTigressAnalysis::~TTigressAnalysis()	{	}

Bool_t TTigressAnalysis::Init(){
	
// gStyle->SetPalette(1);
//  gStyle->SetOptStat(0);
//  gStyle->SetTitleOffset(1.5,"Y");
  SetEfficiencyCurve();
      
	LoadHistos(Form("%s/Results_RedwoodMats.root",TIGDIR),"dp"); // read histograms	
	return InitLevelsGammas(100,false); // load nndc stuff
}

void TTigressAnalysis::Print(Option_t *opt) {

	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
	printf("\n\n\t____TTigressAnalysis____");
  printf("\n\n\t Reaction   :- %s",histfile.c_str());
  printf("\n\t Histo File   :- %s",histfile.c_str());
  printf("\n\t NNDC File    :- %s",nndcfile.c_str());

  printf("\n\n\t Gam bin size :- %.1f",gambinsz);
  printf("\n\t Exc bin size :- %.1f",excbinsz);
  printf("\n\t Addback      :- %s",addback?"TRUE":"FALSE");
  printf("\n\t Verbose      :- %s",verbose?"TRUE":"FALSE");
	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n\n");
		
}

void TTigressAnalysis::Clear(Option_t *opt) {

	histfile.clear();
	reaction.clear();
	nndcfile.clear();

	hexcgamgam = NULL;
	hexcthcmgam = NULL;
	hexcgamthtig = NULL;
	hexcgam = NULL;
	hgamgam = NULL;
	hgam = NULL;		
	hexc = NULL;			

  for(int i=0; i<3; i++){
    hexcgamgam_sel[i] = NULL;
    hexcgamthtig_sel[i] = NULL;
    hexcgam_sel[i] = NULL;
    hgamgam_sel[i] = NULL; 
    hexc_sel[i] = NULL;     
    hgam_sel[i] = NULL;
  }
	
	gambinsz = 1.;
	excbinsz = 1.;	
	addback = true;
	
}

Bool_t TTigressAnalysis::LoadHistos(const char *fname, const char *reac){
	
	if(strcmp(fname,histfile.c_str())!=0 || strcmp(reac,reaction.c_str())!=0){
		
		TFile *f1 = new TFile(fname,"READ");	
		printf("\n Reading Histos File ' %s ' and searching for ' %s ' :\n",fname,reac);				
			
		if(!f1->IsOpen()){
			printf("\n\t Error :  Couldn't open file!\n");
			histfile = "";
			return false;
		}
		
		// this assigns errors to all new histograms
	//	TH1::SetDefaultSumw2();

    // read file and set histograms
    if(f1->IsOpen()){			
	  	histfile.assign(fname);
	  	
      hexcgamgam          = (TH3S*)f1->Get(Form("%s/ExcGamGam_%s",reac,reac)); 
      hexcgamgam_sel[0]   = (TH3S*)f1->Get(Form("%s/UQ/ExcGamGamUQ_%s",reac,reac)); 
      hexcgamgam_sel[1]   = (TH3S*)f1->Get(Form("%s/UB/ExcGamGamUB_%s",reac,reac)); 
      hexcgamgam_sel[2]   = (TH3S*)f1->Get(Form("%s/DB/ExcGamGamDB_%s",reac,reac)); 

      hexcthcmgam         = (TH3S*)f1->Get(Form("%s/ExcGamThetaCm_%s",reac,reac));
      
      hexcgamthtig        = (TH3S*)f1->Get(Form("%s/ExcGamTigressTheta_%s",reac,reac));     
      hexcgamthtig_sel[0] = (TH3S*)f1->Get(Form("%s/UQ/ExcGamTigressThetaUQ_%s",reac,reac));     
      hexcgamthtig_sel[1] = (TH3S*)f1->Get(Form("%s/UB/ExcGamTigressThetaUB_%s",reac,reac));     
      hexcgamthtig_sel[2] = (TH3S*)f1->Get(Form("%s/DB/ExcGamTigressThetaDB_%s",reac,reac));     

      hexcgam 	          = (TH2F*)f1->Get(Form("%s/ExcGam_%s",reac,reac));   
      hexcgam_sel[0] 	    = (TH2F*)f1->Get(Form("%s/UQ/ExcGamUQ_%s",reac,reac));       
      hexcgam_sel[1] 	    = (TH2F*)f1->Get(Form("%s/UB/ExcGamUB_%s",reac,reac));       
      hexcgam_sel[2] 	    = (TH2F*)f1->Get(Form("%s/DB/ExcGamDB_%s",reac,reac));       

      hgamgam 	          = (TH2F*)f1->Get(Form("%s/GamGam_%s",reac,reac));  
      hgamgam_sel[0]      = (TH2F*)f1->Get(Form("%s/UQ/GamGamUQ_%s",reac,reac));       
      hgamgam_sel[1]      = (TH2F*)f1->Get(Form("%s/UB/GamGamUB_%s",reac,reac));       
      hgamgam_sel[2]      = (TH2F*)f1->Get(Form("%s/DB/GamGamDB_%s",reac,reac));   

      hexc 			         = (TH1D*)f1->Get(Form("%s/Exc_%s",reac,reac));  
      hexc_sel[0]        = (TH1D*)f1->Get(Form("%s/UQ/ExcUQ_%s",reac,reac));       
      hexc_sel[1]        = (TH1D*)f1->Get(Form("%s/UB/ExcUB_%s",reac,reac));       
      hexc_sel[2]        = (TH1D*)f1->Get(Form("%s/DB/ExcDB_%s",reac,reac));           

      hgam			         = (TH1D*)f1->Get(Form("%s/Gam_%s",reac,reac));
      hgam_sel[0]        = (TH1D*)f1->Get(Form("%s/UQ/GamUQ_%s",reac,reac));       
      hgam_sel[1]        = (TH1D*)f1->Get(Form("%s/UB/GamUB_%s",reac,reac));       
      hgam_sel[2]        = (TH1D*)f1->Get(Form("%s/DB/GamDB_%s",reac,reac));               
    }
    
    // print status and set flags
    if(f1->IsOpen()){
      printf("\n Loaded Basic Histograms :-\n");
      if(hexcgamgam)   printf("\t-> (TH3S*) %s\n",hexcgamgam->GetName()); 
      if(hexcthcmgam)  printf("\t-> (TH3S*) %s\n",hexcthcmgam->GetName());  
      if(hexcgamthtig) printf("\t-> (TH3S*) %s\n",hexcgamthtig->GetName()); 	
      if(hexcgam)      printf("\t-> (TH2F*) %s\n",hexcgam->GetName());   	
      if(hgamgam)      printf("\t-> (TH2F*) %s\n",hgamgam->GetName()); 
      if(hexc)         printf("\t-> (TH1D*) %s\n",hexc->GetName());   
      if(hgam)         printf("\t-> (TH1D*) %s\n",hgam->GetName()); 

      printf("\n Loaded Addditional Histograms :-\n");
      for(int i=0; i<3; i++){
        printf(" + Detector region %i :\n",i);
        if(hexcgamgam_sel[i])   printf("\t+  (TH3S*) %s\n",hexcgamgam_sel[i]->GetName()); 
        if(hexcgamthtig_sel[i]) printf("\t+  (TH3S*) %s\n",hexcgamthtig_sel[i]->GetName()); 
        if(hexcgam_sel[i])      printf("\t+  (TH2F*) %s\n",hexcgam_sel[i]->GetName()); 
        if(hgamgam_sel[i])      printf("\t+  (TH2F*) %s\n",hgamgam_sel[i]->GetName()); 
        if(hexc_sel[i])         printf("\t+  (TH1D*) %s\n",hexc_sel[i]->GetName()); 
        if(hgam_sel[i])         printf("\t+  (TH1D*) %s\n",hgam_sel[i]->GetName()); 
      }
    }
    
    // by default the histograms cover all of SHARC  [not a specific range]
    if(!DetectorSelector("all")||!hexcthcmgam)
      return false;
     
    reaction.assign(reac);  
      			  		      
	} else

	printf("\n\t File ' %s ' Has Been Loaded !\n",fname);
	return true;
}

Bool_t TTigressAnalysis::DetectorSelector(std::string opt){

  if(!IncludeRegion(-1)){
    printf("\n Looks like there are no available histograms. Please load histograms.\n\n");  
    return false;
  }  
  
  printf("\n Detector Selector : '%s'\n",opt.c_str());
  const unsigned long npos = std::string::npos;  
  if(opt.find("all")!=npos || opt.size()==0){ // use all angles hist
    printf("\n Histograms have been set to default ranges [all angles]\n\n");  
    return true;
  }
  detsel.clear();
  detsel.assign(opt);
// keep titles and formatting but get rid of data  
  ResetDefaults();
  
  if(opt.find("uqqq")!=npos)
    if(!IncludeRegion(0)) return false;

  if(opt.find("ubox")!=npos)
    if(!IncludeRegion(1)) return false;

  if(opt.find("dbox")!=npos)
    if(!IncludeRegion(2)) return false;  
  
  return true;    
}

Bool_t TTigressAnalysis::IncludeRegion(Int_t indx){

  if( indx<0 ){
    if(!hexcgamgam||!hexcgamthtig||!hexcgam||!hgamgam||!hgam||!hexc)
      return false;
    else {
      CloneDefaults();
      return true;
    }
  }

  if(!hexcgamgam_sel[indx] || !hexcgam_sel[indx] || !hgamgam_sel[indx] || !hexc_sel[indx] || !hgam_sel[indx]){
    printf("\n\t Error :  Could not locate all of the histograms for this region. Returning to default\n"); 
    IncludeRegion(-1);
    return false;
  }
  
  printf("\t-> hexcgamgam   += (TH3S*) %s\n",hexcgamgam_sel[indx]->GetName()); 
  printf("\t-> hexcgamthtig += (TH3S*) %s\n",hexcgamthtig_sel[indx]->GetName()); 
  printf("\t-> hexcgam      += (TH2F*) %s\n",hexcgam_sel[indx]->GetName());  
  printf("\t-> hgamgam      += (TH2F*) %s\n",hgamgam_sel[indx]->GetName());  
  printf("\t-> hexc         += (TH1D*) %s\n",hexc_sel[indx]->GetName());  
  printf("\t-> hgam         += (TH1D*) %s\n",hgam_sel[indx]->GetName()); 

  hegg->Add(hexcgamgam_sel[indx]);      
  hegt->Add(hexcgamthtig_sel[indx]);      
  heg->Add(hexcgam_sel[indx]);   
  hgg->Add(hgamgam_sel[indx]);      
  he->Add(hexc_sel[indx]);      
  hg->Add(hgam_sel[indx]);    
  
  return true;
}

TH1D *TTigressAnalysis::Gam(Double_t exc_lo, Double_t exc_hi){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
		
  TH1D *h;
  if(exc_lo>=0.0 && exc_hi>exc_lo){ // set exc energy gate if included
		Double_t gatesz;
  // produce an excitation energy (Y axis) gated gamma spectrum  			
		h = (TH1D*)TH2Proj(hexcgam,'x',exc_lo,exc_hi,gatesz); 
		TAxis *yax = heg->GetYaxis();
  	Int_t yp[2] = {yax->FindBin(exc_lo), yax->FindBin(exc_hi)};
		exc_lo = yax->GetBinCenter(yp[0]);
		exc_hi = yax->GetBinCenter(yp[1]);
  	printf("\n Gam:  Excitation energy range : %.2f-%.2f keV [ bins %i - %i ]\n",exc_lo,exc_hi,yp[0],yp[1]);
		
  	h->SetTitle(Form("Gamma Singles Gated on Excitation Energy Range %.1f - %.1f keV",exc_lo,exc_hi));
  } else 
  	h = (TH1D*)hg->Clone("GammaEnergy");

  Double_t reb = gambinsz/h->GetBinWidth(1);
  if(reb>1)
	  h->Rebin((Int_t)reb);	

  h->SetNameTitle("GammaEnergy",Form("%s; E_{#gamma} [keV]; Counts / %.1f keV",h->GetTitle(),h->GetBinWidth(0)));	

	return h;	
}

TH1D *TTigressAnalysis::GamGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);
	
  TH2F *h2 = GamGam(exc_lo,exc_hi);
  
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
  TH1D *h = TH1Sum(hpt[0],hpt[3],1,-0.5);
  h->Scale(0.5); // accounts for double filling of gam-gam matrix
  
  Double_t reb = gambinsz/h->GetBinWidth(1);
  if(reb>1)
	  h->Rebin((Int_t)reb);
	  
  h->SetNameTitle("GatedGammaEnergy",Form("Gammma Energy Coincident With Gated Gam %sEnergy; E_{#gamma} [keV]; Counts / %.1f keV",exc_hi>0?"& Exc ":"",h->GetBinWidth(0)));
      
	return h;
}	

TH2F *TTigressAnalysis::GamGam(Double_t exc_lo, Double_t exc_hi){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
  TH2F *h2;
  if(exc_lo>=0.0 && exc_hi>exc_lo){ // set exc energy gate if included
		TAxis *zax = hegg->GetZaxis();  	
  	Int_t zp[2] = {zax->FindBin(exc_lo), zax->FindBin(exc_hi)};
		exc_lo = zax->GetBinCenter(zp[0]);
		exc_hi = zax->GetBinCenter(zp[1]);  	
  	printf("\n GamGam :  Excitation energy range : %.2f-%.2f keV [ bins %i - %i ]\n",exc_lo,exc_hi,zp[0],zp[1]);
  	hegg->GetZaxis()->SetRange(zp[0],zp[1]);
  // produce an excitation energy (Z axis) gated gamma gamma matrix
  	h2 = (TH2F*)hegg->Project3D("yx");   
    // manually set error
  	for(int i1=1; i1<=h2->GetNbinsX(); i1++)
      for(int i2=1; i2<=h2->GetNbinsY(); i2++)
        h2->SetBinError(i1,i2,sqrt((double)h2->GetBinContent(i1,i2)));
    
  } else 
  	h2 = (TH2F*)hgg->Clone("GammaGamma");

/*
  Double_t reb = gambinsz/h->GetBinWidth(1);
  if(reb>1)
	  h2->Rebin((Int_t)reb);	
*/
  h2->SetNameTitle("GammaGamma",Form("%s; E_{#gamma} [keV]; E_{#gamma} [keV]; Counts / %.1f keV",h2->GetTitle(),h2->GetYaxis()->GetBinWidth(1)));	

	return h2;	
}


TH1D *TTigressAnalysis::ExcGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}

	TH1D *h;
	if(!emin && !emax)
		h = (TH1D*) he->Clone("ExcEnergy");
	else{
		SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);

		TH1D *hp[4];
		TH2F *h2 = heg;
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
		h = TH1Sum(hp[0],hp[3],1,-0.5);
  }    
   	
	h->SetLineColor(1);
	h->GetYaxis()->SetTitleOffset(1.3);
  Double_t reb = excbinsz/h->GetBinWidth(1);
  if(reb>1)	
	  h->Rebin(reb);
	  
  if(!emin && !emax)	  
    h->SetNameTitle("ExcEnergy",Form("Excitation Energy; E_{exc} [keV]; Counts / %.1f keV",excbinsz));
  else
    h->SetNameTitle("GatedExcEnergy",Form("Excitation Energy Coincident With Gated Gamma Energy; E_{exc} [keV]; Counts / %.1f keV",excbinsz));

	return h;
}

TH2F *TTigressAnalysis::ExcThetaGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
	TH2F *hp[4], *het;
	if(emin>=0 && emax>emin){
		
		SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);	
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
		het = TH2Sum(hp[0],hp[3],1,-0.5);
  } else{
  	Double_t tmpsz;
  	het = TH3Proj(hexcthcmgam,"zx",0,4000.0,tmpsz);	 
 	}
 	
  Double_t reb = excbinsz/het->GetYaxis()->GetBinWidth(1);
  if(reb>1)	
	  het->RebinY(reb); 	

	het->SetMinimum(0);  
  het->SetNameTitle("GatedExcVsThetaCm",Form("Excitation Energy Versus Theta Cm; #theta_{CM} [#circ]; E_{exc} [keV]"));

	return het;
}

TH2F *TTigressAnalysis::ExcGamGated(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
		
	TH2F *hp[4], *h2;
	if(emin>=0 && emax>emin){

		SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);	
		TH3S *h3 = hegg;
		// x = gam
		// y = gam
		// z = exc		
		
		Double_t peaksz, bglosz, bghisz;
		// peak
		hp[0] = TH3Proj(h3,"zy",emin,emax,peaksz);	 
		// bglo  
		hp[1] = TH3Proj(h3,"zy",bg0,bg1,bglosz);		
		// bghi  
		hp[2] = TH3Proj(h3,"zy",bg2,bg3,bghisz);	
		// bgtot	
		hp[3] = TH2Sum(hp[1],hp[2],peaksz/bglosz,peaksz/bghisz);
	
		// peak - bgtot
		h2 = TH2Sum(hp[0],hp[3],1,-0.5);
	} else 
		h2 = (TH2F*) heg->Clone("ExcVsGamma");

  Double_t rebx = gambinsz/h2->GetXaxis()->GetBinWidth(1); 	
  if(rebx>1)	
	  h2->RebinX(rebx);
	  
  Double_t reby = excbinsz/h2->GetYaxis()->GetBinWidth(1);
  if(reby>1)	
	  h2->RebinY(reby);  	

	h2->SetMinimum(0);
  h2->SetNameTitle("ExcVsGamma",Form("Excitation Energy Versus Gamma Energy; E_{#gamma} [keV]; E_{exc} [keV]"));

  if(emin>=0 && emax>emin)
    h2->SetTitle(Form("Gated %s",h2->GetTitle()));
    
	return h2;	
}

void TTigressAnalysis::SetBackgroundLims(Double_t emin, Double_t emax, Double_t &bg0, Double_t &bg1, Double_t &bg2, Double_t &bg3){
	// default behaviour for setting background region
	
	if(!emin && !emax)
		return;	
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
	printf("\n Background Limits Automatically Set :-");
	printf("\n\t-> PEAK :  %.1f - %.1f keV \n\t-> BGLO :  %.1f - %.1f keV \n\t-> BGHI :  %.1f - %.1f keV\n",emin,emax,bg0,bg1,bg2,bg3);
	return;	
}

void TTigressAnalysis::SetPeakLims(Double_t egam, Double_t &emin, Double_t &emax){

  if(!fTigSigma || fTigSigma->Eval(egam)<1.0){
    printf("\n\t Error :  Peak limits could not be set for egam = %.1f keV \n\n",egam);
    return;
  }

  emin = egam - 2.0*fTigSigma->Eval(egam);
  emax = egam + 2.0*fTigSigma->Eval(egam);
	
	printf("\n Peak Limits Automatically Set For %.0f keV Gamma:-",egam);
	printf("\n\t-> PEAK :  %.1f - %.1f keV \n",emin,emax);  
	
	return;
}

TCanvas *TTigressAnalysis::QuickAnalysis(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	TCanvas *canv = 0;
	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return canv;
	}
		
  canv = new TCanvas("QuickAnalysis","Gamma-Exc Analysis",0,0,800,600);
  canv->cd();
  
  TPad *p1 = new TPad("p1","p1",0,0.5,1,1);
  p1->Draw(); 
  p1->SetBottomMargin(0.0);

  TPad *p2 = new TPad("p2","p2",0,0,1,0.5);
  p2->Draw();  
  p2->SetTopMargin(0.01);
  p2->SetLogz(1);  
  
  Double_t xmax = 2000.;
  if(exc_hi>xmax)
    xmax = exc_hi+400.0;
  ////////////////////////////////////////////////////////////////////////////////////

  p1->cd();   
  // DRAW GAMMA SINGLES (FITTED, AND IN EXC RANGE) ON PAD1  
  TH1D *hgamx;
  if(emin && emax){  
    hgamx = GamGated(emin,emax,bg0,bg1,bg2,bg3,exc_lo,exc_hi);
    hgamx->GetXaxis()->SetRangeUser(0,4000.0);  
  } else
    hgamx = Gam(exc_lo,exc_hi);

  hgamx->SetTitleSize(0.05,"Y");
    
  TAxis *xax = hgamx->GetXaxis();  
  xax->SetRangeUser(0,xmax);  
  
  TAxis *yax = hgamx->GetYaxis();
  yax->CenterTitle();
  yax->SetTitleOffset(1.0);  
  
  hgamx->Draw("hist");  
    
  ////////////////////////////////////////////////////////////////////////////////////
  p2->cd();   
  TH2F *hexcgamx = ExcGamGated(emin,emax,bg0,bg1,bg2,bg3);
  hexcgamx->SetMinimum(1);
  hexcgamx->SetTitleSize(0.05,"X");
  hexcgamx->SetTitleSize(0.05,"Y");
  hexcgamx->SetTitle("");
    
  xax = hexcgamx->GetXaxis();
  xax->SetRangeUser(0,xmax);  
  yax->SetTitleOffset(1.0);
  xax->CenterTitle();  
  
  yax = hexcgamx->GetYaxis();
  yax->SetRangeUser(exc_lo,exc_hi+50.0);
  yax->SetTitleOffset(1.0);  
  yax->CenterTitle();
  
  hexcgamx->Draw("colz");
  
  TF1 *fun = new TF1("Direct_Population","pol1",0,4000);
  fun->SetParameters(0,1.0);
  fun->SetLineWidth(2.0);
  fun->SetLineStyle(2);
  fun->Draw("same");

  ////////////////////////////////////////////////////////////////////////////////////  
  printf("\n\tMade Quick Analysis Plots!\n\n");
	return canv;

}

TCanvas *TTigressAnalysis::AnalyzeGammas(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	TCanvas *c = 0;
	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return c;
	}
	
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);	
	
  c = new TCanvas("GammaAnalysis","Gamma Coincidence Analysis",0,0,1200,800);
  c->Divide(2,2);
  
  ////////////////////////////////////////////////////////////////////////////////////

  c->cd(1);
  // DRAW GAMMA SINGLES (FITTED, AND IN EXC RANGE) ON PAD1  
  TH1D *hgamx = Gam(exc_lo,exc_hi);
  FitPeakStats(hgamx,emin,emax,bg0,bg1,bg2,bg3);
 	 
////////////////////////////////////////////////////////////////////////////////////
  c->cd(3);

  TH1D *hexcgated = ExcGated(emin,emax,bg0,bg1,bg2,bg3); 
  hexcgated->GetXaxis()->SetRangeUser(0,6000);
  hexcgated->SetLineColor(kBlue);   
  hexcgated->Draw();    

  TBox *bpex = new TBox(exc_lo,hexcgated->GetMinimum()*1.05,exc_hi,hexcgated->GetMaximum()*1.05);
  bpex->SetFillColor(3);
  bpex->SetFillStyle(3002);  

  if(exc_lo>=0.0 && exc_hi>exc_lo)
  	  bpex->Draw();

  ////////////////////////////////////////////////////////////////////////////////////
  c->cd(2);
  TH1D *hexcgamgated = GamGated(emin,emax,bg0,bg1,bg2,bg3,exc_lo,exc_hi); 
   hexcgamgated->Sumw2();
  hexcgamgated->GetXaxis()->SetRangeUser(0,1500); 
  hexcgamgated->SetLineColor(kBlue); 
  hexcgamgated->Draw("");

  
  ////////////////////////////////////////////////////////////////////////////////////  
  c->cd(4);
  TH2F *hexcthcm = ExcThetaGated(emin,emax,bg0,bg1,bg2,bg3); 
  hexcthcm->Draw("colz");

  ////////////////////////////////////////////////////////////////////////////////////  
  printf("\n\tMade Analysis Plots!\n\n");
	return c;
}

TH1D *TTigressAnalysis::GamAngCorr(Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Double_t exc_lo, Double_t exc_hi){

	TH2F *h2 = GamAngCorrMat(exc_lo,exc_hi);
	if(!h2)
		return 0;
	
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);
	
	TH1D *hp[4];	
	Double_t peaksz, bglosz, bghisz;
  // peak
  hp[0] = TH2Proj(h2,'x',emin,emax,peaksz);	
	// bglo  
  hp[1] = TH2Proj(h2,'x',bg0,bg1,bglosz);		
	// bghi  
  hp[2] = TH2Proj(h2,'x',bg2,bg3,bghisz);	
  // bgtot	
  hp[3] = TH1Sum(hp[1],hp[2],peaksz/bglosz,peaksz/bghisz);
  
  // peak - bgtot
  TH1D *htot = TH1Sum(hp[0],hp[3],1,-0.5);
  htot->SetName(Form("CountsVsTigTheta_Gam%.0fTo%.0f_Exc%.0fTo%.0f",emin,emax,exc_lo,exc_hi));
  
  // rebin data to remove segments and return to crystal positions
  // divide through by basic hit pattern (normalized so max=1)
  TH1D *ht = NormalizeThetaHits(htot,0);
	
	/*
	TH1D *htmp = h2->ProjectionX();	
	Double_t tot, val;
	for(int i=htmp->FindFirstBinAbove(); i<=htmp->FindLastBinAbove(); i++){
		tot = h2->ProjectionY("",i,i)->Integral();
		if(!tot)
			continue;
			
		val = ht->GetBinContent(i);
		ht->SetBinContent(i,val/tot); // normalize: ratio to total counts in slice
		ht->SetBinError(i,sqrt(val)/tot); // set errors
	}

	ht->Scale(1/ht->GetMaximum());	
	*/
	const char *msg = "";
	if(exc_lo>0 && exc_hi>exc_lo)
		msg = Form("And %.1f-%.1f keV Exc. Energy Range",exc_lo,exc_hi);
	ht->SetTitle(Form("Angular Correlations For %.1f-%.1f keV Gamma Gate %s",emin,emax,msg));

	return ht;	
}

TH2F *TTigressAnalysis::GamAngCorrMat(Double_t exc_lo, Double_t exc_hi){

	if(!histfile.size()){
		printf("\n\t Error :  No file has been loaded! \n\n");
		return 0;
	}
	
	// x = gam
	// y = thetatig
	// z = exc
			
  TH2F *h2;
	TAxis *zax = hegt->GetZaxis();
	Double_t gatesz;
  
  if(exc_lo>=0.0 && exc_hi>exc_lo){ // set exc energy gate if included
		h2 = (TH2F*)TH3Proj(hegt,"yx",exc_lo,exc_hi,gatesz); 
  	Int_t zp[2] = {zax->FindBin(exc_lo), zax->FindBin(exc_hi)};
		exc_lo = zax->GetBinCenter(zp[0]);
		exc_hi = zax->GetBinCenter(zp[1]);
  	printf("\n GamAngCorrMat:  Excitation energy range : %.2f-%.2f keV [ bins %i - %i ]\n",exc_lo,exc_hi,zp[0],zp[1]);
		
  	h2->SetTitle(Form("Gammas Versus Tigress Theta Gated on Excitation Energy Range %.1f - %.1f keV",exc_lo,exc_hi));
  } else {
		exc_lo = zax->GetBinCenter(1);
		exc_hi = zax->GetBinCenter(hegt->GetNbinsZ());
		h2 = (TH2F*)TH3Proj(hegt,"yx",exc_lo,exc_hi,gatesz); 
	}
	
  Double_t reby = gambinsz/h2->GetYaxis()->GetBinWidth(1);
  if(reby>1)	
	  h2->RebinY(reby);  	
	
  h2->SetNameTitle("GamVsThetaTig",Form("%s; #theta_{TIGRESS} [#circ]; E_{#gamma} [keV]",h2->GetTitle()));	

	return h2;		
}

TH1D *TTigressAnalysis::GetHitPattern(Int_t detmin, Int_t detmax, Int_t crymin, Int_t crymax, Int_t segmin, Int_t segmax){
  
  TH1D *hhit = new TH1D("TigHitPattern","TIGRESS Hit Pattern",180,0,180);
  
  for(int d=detmin; d<=detmax; d++)
    for(int c=crymin; c<=crymax; c++)
      for(int s=segmin; s<=segmax; s++)
        hhit->Fill(TTigress::GetPosition(d,c,s).Theta()*TMath::RadToDeg());
  
  return hhit;
}

TH1D *TTigressAnalysis::NormalizeThetaHits(TH1D *htheta, Int_t segval){
  
  TH1D *hhit;
  if(segval<0)
    hhit = GetHitPattern(5,16,0,3,0,7);
  else 
    hhit = GetHitPattern(5,16,0,3,segval,segval);
  
  std::vector<int> binvals;
  // build hit pattern without segments
  for(int i=hhit->FindFirstBinAbove(); i<=hhit->FindLastBinAbove(); i++){
    if(hhit->GetBinContent(i)==0)
      continue;
    binvals.push_back(i);
  }
  
  TH1D *h = (TH1D*) htheta->Clone(Form("%s_NoSeg",htheta->GetName()));
  h->Reset(); // remove all data
  Double_t val;
  Int_t bin_near, bin_diff;
  for(int i=htheta->FindFirstBinAbove(); i<=htheta->FindLastBinAbove(); i++){
    
    val = htheta->GetBinContent(i);
    if(!val)
      continue; 
    
    bin_near = 0.;
    bin_diff = 100.0;
    // determine where to round the theta values to by using hit pattern
    for(int j=0; j<binvals.size(); j++){
      
      if(abs(i-binvals.at(j))<bin_diff){
    //    printf("\n i=%3i, j=%i\t htheta->GetBinCenter(i) = %.1f,  thetavals.at(j) = %.1f, diff = %.1f deg",i,j,htheta->GetBinCenter(i),thetavals.at(j),fabs(htheta->GetBinCenter(i)-thetavals.at(j)));      
        bin_near = binvals.at(j);
        bin_diff = abs(i-binvals.at(j));
      }
    }
   
   // put the content of htheta->Bin(i) into the nearest nonzero theta
    h->Fill(bin_near-1,val); 
  }
 
  // scale the hit pattern
  hhit->Scale(1/hhit->GetMaximum()); 
  h->Divide(hhit);
  
  return h;
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
	 //  printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);

	fitfunc->SetParameters(0, 0, gContent, gMean, gSigma); 
    
	fitfunc->SetParName(0,"BgConstant");
	fitfunc->SetParName(1,"BgSlope");
	fitfunc->SetParName(2,"Area");          fitfunc->SetParLimits(2,0,gBinW*gContent);
	fitfunc->SetParName(3,"Mean");          fitfunc->SetParLimits(3,emin,emax);
	fitfunc->SetParName(4,"Sigma");         fitfunc->SetParLimits(4,0,500.);	
	fitfunc->SetParName(5,"ExcludeLo Min"); fitfunc->FixParameter(5,bg1);
	fitfunc->SetParName(6,"ExcludeLo Max"); fitfunc->FixParameter(6,emin);
	fitfunc->SetParName(7,"ExcludeHi Min"); fitfunc->FixParameter(7,emax);
	fitfunc->SetParName(8,"ExcludeHi Max"); fitfunc->FixParameter(8,bg2);         

	hist->Fit(fitfunc,"QR","SAME");

	gBgConstant = fitfunc->GetParameter(0);
	gBgSlope    = fitfunc->GetParameter(1);
	gContent    = fitfunc->GetParameter(2)/gBinW;
//	printf("\n\t gContent = %.2f\n",gContent);
	//fitfunc->SetParameter(2,gContent);
	//fitfunc->SetParError(2,fitfunc->GetParError(2)/gBinW);
	gMean       = fitfunc->GetParameter(3);
	gSigma      = fitfunc->GetParameter(4);	
	gChi2pNDF   = fitfunc->GetChisquare() / fitfunc->GetNDF();

}

void TTigressAnalysis::FitBgExcludeRange(TH1 *hist, Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Bool_t quad_fit){ 

	TF1 *fitfunc = new TF1("pol_bg_exc",pol_bg_exc, bg0, bg3, 7);
	fitfunc->SetLineColor(2);
	fitfunc->SetLineWidth(1);

	fitfunc->SetParameters(0, 0, 0); 

	fitfunc->SetParName(0,"BgConstant"); 
	fitfunc->SetParName(1,"BgSlope");
	fitfunc->SetParName(2,"BgQuad");
	if(!quad_fit) fitfunc->FixParameter(2,0.0);
	
	fitfunc->SetParName(3,"ExcludeLo Min "); fitfunc->FixParameter(3,bg1);
	fitfunc->SetParName(4,"ExcludeLo Max "); fitfunc->FixParameter(4,emin);
	fitfunc->SetParName(5,"ExcludeHi Min "); fitfunc->FixParameter(5,emax);
	fitfunc->SetParName(6,"ExcludeHi Max "); fitfunc->FixParameter(6,bg2);         

	hist->Fit(fitfunc,"QR","SAME");
}


Double_t TTigressAnalysis::FitPeakStats(TH1 *hist, Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3, Bool_t bg_only, Bool_t quad_fit){

  if(!hist)
    return 0;
	SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);	
  
  TAxis *xax = hist->GetXaxis();
  xax->SetRangeUser(bg0-20,bg3+80);
  
  Double_t bin = hist->GetBinWidth(1);  
  hist->GetYaxis()->SetTitle(Form("Counts / %.0f keV",bin));
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->DrawCopy();
  
  // selected peak and background region boxes
  TBox *bp = new TBox(emax,hist->GetMaximum()*.01,emin,hist->GetMaximum()*1.05);
  bp->SetFillColor(3);
  bp->SetFillStyle(3002); 
  TBox *b1 = new TBox(bg0,hist->GetMaximum()*.01,bg1,hist->GetMaximum()*1.05);
  b1->SetFillColor(2);
  b1->SetFillStyle(3002);
  TBox *b2 = new TBox(bg2,hist->GetMaximum()*.01,bg3,hist->GetMaximum()*1.05);
  b2->SetFillColor(2);
  b2->SetFillStyle(3002);   
    
  bp->Draw();
  b1->Draw();
  b2->Draw();
  //printf("\n\n\t\t [%.0f-%.0f] <%.0f-%.0f> [%.0f-%.0f]\n\n",bg0,bg1,emin,emax,bg2,bg3);
	
  TPaveText *bstat = new TPaveText(bg3+10,hist->GetMinimum()*1.4,bg3+100,hist->GetMaximum()*1.05);
  bstat->SetLineColor(1);
  bstat->SetShadowColor(0);
  bstat->AddText("Fit Results :");
  bstat->AddLine(.0,.85,1.,.85);
  
  TF1 *func;
  Double_t counts, cnterr;
  if(!bg_only) {
    
    FitPeakExcludeRange(hist,emin,emax,bg0,bg1,bg2,bg3);  
    func = hist->GetFunction("gauss_linbg_exc");
    
    counts = func->GetParameter(2)/bin;
    cnterr = func->GetParError(2)/bin;  

    printf("\n\n ____ peak+lbg fit :- ___");  
    printf("\n\n -> cnt = %.1f +/- %.1f\n\n",counts,cnterr);  
    for(int fn=0; fn<5; fn++){
      bstat->AddText(Form("%10s : %6.4f +/- %6.4f",func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn)));
      printf("\n Par %i.  %10s : %6.4f +/- %6.4f",fn,func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn));  
    }                
  } else {
    Double_t tot, tot_err, bg, bg_err;  
    
    FitBgExcludeRange(hist,emin,emax,bg0,bg1,bg2,bg3,quad_fit);    
    
    func = hist->GetFunction("pol_bg_exc");  
    tot = hist->IntegralAndError(xax->FindBin(emin),xax->FindBin(emax),tot_err);
    bg = func->Integral(emin,emax)/hist->GetBinWidth(1);
    bg_err = func->IntegralError(emin,emax);  
    
    counts = tot - bg;
    cnterr = sqrt(pow(tot_err/tot,2.0)+pow(bg_err/bg,2.0))*counts;   
    
    printf("\n\n ____ bg pol%i fit :- ___",quad_fit?2:1);
    printf("\n\n -> tot = %.1f +/- %.1f",tot,tot_err);
    printf("\n -> bg  = %.1f +/- %.1f",bg,bg_err);    
    printf("\n -> cnt = %.1f +/- %.1f\n",counts,cnterr);   
    
    Int_t fnmax = 2;
    if(quad_fit)
      fnmax=3;    
    for(int fn=0; fn<fnmax; fn++){
      bstat->AddText(Form("%10s : %6.4f +/- %6.4f",func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn)));
      printf("\n Par %i.  %10s : %6.4f +/- %6.4f",fn,func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn));  
    }    
  }
  bstat->AddText(Form("\n* COUNTS  : %6.2f +/- %6.2f *",counts,cnterr));  
  printf("\n * COUNTS  : %6.2f +/- %6.2f *\n\n",counts,cnterr);
   
  printf(" _______________________ \n\n");

 // bstat->AddText(Form("Chi2/DOF : %8.2f",func->GetChisquare()/func->GetNDF()));
  bstat->Draw();

  return counts;
}

// Normalizes theory curve to data in region emin-emax so that observed peaks can be compared to theory
// ** emin and emax etc are for the normalization, not the gamma gate **
TH1D *TTigressAnalysis::FitTheory(TH1D *hdata, Double_t exc, Double_t egam, Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

  TCanvas *c;
  if(!hdata){
    c = new TCanvas;
    if(!egam)
      hdata = Gam(exc-400.0,exc+400.0);
    else{
      Double_t sig = Resolution(egam);
      hdata = GamGated(egam-3.5*sig,egam+3.5*sig,0,0,0,0,exc-400.0,exc+400.0);
    }
  }
  // make sure peak and background range is set
  if(!emin || !emax)
    SetPeakLims(egam,emin,emax);
  SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);  
  
  // get counts in data peak
  FitPeakExcludeRange(hdata,emin,emax,bg0,bg1,bg2,bg3);
  TF1 *func = hdata->GetFunction("gauss_linbg_exc");
  for(int fn=0; fn<5; fn++)
  	printf("\n\t%s : %6.2f +/- %6.2f",func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn));
 
  Double_t bin = hdata->GetBinWidth(1);
  Double_t counts_data = func->GetParameter(2)/bin, cnterr_data = func->GetParError(2)/bin;
  printf("\n\t** COUNTS  : %6.2f +/- %6.2f **  \n- - - - - - - - - - - - - - - - - - - - -\n\n",counts_data,cnterr_data);      

  // make theory photo-peak
	Int_t from_state = GetStateIndex(exc,false)+1;  
	printf("\n\t EXC = %.2f keV ---> FROM_STATE = %i\n",exc,from_state);
	if(from_state<=0)
	  return hdata;
	  
  TH1D *hthry_tmp = DrawGammas(from_state,egam);
  if(hthry_tmp->GetEntries()==0){
    printf("\n\t No gamma-rays could be found in coincidence with %.1f keV.\n",egam);
    return hdata;
  }
  printf("\n\t GAM \t INTENSITY");
  Double_t val=0.0;
  for(int i=hthry_tmp->FindFirstBinAbove();i<=hthry_tmp->FindLastBinAbove(); i++){
    val = hthry_tmp->GetBinContent(i);
    if(val)
      printf("\n\t-> %.1f\t%.1e",hthry_tmp->GetBinCenter(i),val);
  }
  printf("\n\n");
  // add realistic width and apply efficiency
  TH1D *hthry = MakeRealistic(hthry_tmp,true);
  hthry->SetLineColor(kRed);
  
  Double_t reb = bin/hthry->GetBinWidth(1);
  if(reb>1)
    hthry->Rebin(reb);  
  TAxis *xax = hthry->GetXaxis();
  
  // now get theory counts and find scaling factor
  Double_t counts_thry = hthry->Integral(xax->FindBin(emin),xax->FindBin(emax));
  Double_t scale_thry = counts_data/counts_thry, relerr_thry = cnterr_data/counts_data;
  Double_t scalerr_thry = relerr_thry*scale_thry;
  
  hthry->Scale(scale_thry);    
  printf("\n\t** SCALE  : %6.2f +/- %6.2f **  \n- - - - - - - - - - - - - - - - - - - - -\n\n",scale_thry,scalerr_thry);     
  
  if(egam){
    Double_t eff=Efficiency(egam), eff_err=EfficiencyError(egam);
    scale_thry/=eff;
    scalerr_thry = sqrt(pow(relerr_thry,2)+pow(eff_err/eff,2))*scale_thry;
    printf("\n\t Efficiency correction for %.1f keV gamma gate = %.2f\n",egam,1/eff);
  } 
  printf("\n\t** SCALE  : %6.2f +/- %6.2f **  \n- - - - - - - - - - - - - - - - - - - - -\n\n",scale_thry,scalerr_thry);     
  
  if(c){
    hdata->Draw("hist");
    hthry->Draw("same");
  }
  
  return hthry;
}

Double_t TTigressAnalysis::GetPopulationStrength(TH1D *hist, Double_t exc, Double_t egam, Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){
  
  if(!hist){
    new TCanvas;
    printf("\n\n\t Making histogram for excitation energy range :  %.1f - %.1f keV.\n\n",exc-400.0,exc+400.0);
    hist = Gam(exc-400.0,exc+400.0);
  }
  
  if(!emin || !emax)
    SetPeakLims(egam,emin,emax);  
  AnalyzeGammas(emin,emax,bg0,bg1,bg2,bg3,exc-400.0,exc+400.0);
  
  TF1 *func = hist->GetFunction("gauss_linbg_exc");
  for(int fn=0; fn<5; fn++)
  	printf("\n\t%s : %6.2f +/- %6.2f",func->GetParName(fn),func->GetParameter(fn),func->GetParError(fn));
 
  Double_t bin = hist->GetBinWidth(1);
  Double_t counts = func->GetParameter(2)/bin, cnterr = func->GetParError(2)/bin;
  printf("\n\t** COUNTS  : %6.2f +/- %6.2f **  \n- - - - - - - - - - - - - - - - - - - - -\n\n",counts,cnterr);  
  
  Double_t eff = fTigEff->Eval(egam), efferr = eff*fTigEff->GetParError(0)/fTigEff->GetParameter(0);
  Double_t br = BranchingRatio(exc,egam);

  if(!fTigEff || eff<0 || !br)
    return 0;
  
  Double_t strength = counts/(br*eff*0.01);  
  Double_t strerr = strength*(cnterr/counts+efferr/eff);
  
  printf("\n\n\t Eff @ ~%.1f keV = %.2f +/- %.2f %%, Branching ratio = %.2f",egam,eff,efferr,br);
  printf("\n\t -> Strength = %.2f +/- %.2f\n\n",strength,strerr);
   
  return strength;
}

TList *TTigressAnalysis::CheckDopplerCorrection(Double_t egam, TH2F *h2, Double_t emin, Double_t emax, Double_t bg0, Double_t bg1, Double_t bg2, Double_t bg3){

  TList *list = new TList;
  
  if(!emin && !emax)
    SetPeakLims(egam,emin,emax);
  
  if(!h2)
    h2 = GamAngCorrMat(0,5000);
  if(!h2)
    return list;
        
  list->Add(h2);
//  h2->RebinX();
  h2->RebinY();
  printf("\n Setting Gam axis to %.0f keV per Bin\n\n",h2->GetYaxis()->GetBinWidth(1));
  
  Double_t therr = h2->GetXaxis()->GetBinWidth(1)/2;  
    
  TGraphErrors *gdop = new TGraphErrors();
  gdop->SetNameTitle(Form("OldDop_Egam%.0f",egam),"Doppler Corrected Peak; #theta_{LAB} [#circ]; Peak Mean [keV];");
  
  TGraphErrors *graw = new TGraphErrors();
  graw->SetNameTitle(Form("RawGam_Egam%.0f",egam),"Un-Corrected Peak; #theta_{LAB} [#circ]; Peak Mean [keV];");
  graw->SetMarkerColor(kRed);
  graw->SetLineColor(kRed); 
  
  TGraphErrors *gdop2 = new TGraphErrors();
  gdop2->SetNameTitle(Form("NewDop_Egam%.0f",egam),"Fitted Doppler Peak; #theta_{LAB} [#circ]; Peak Mean [keV];"); 
  gdop2->SetMarkerColor(kBlue);
  gdop2->SetLineColor(kBlue);
  
  SetBackgroundLims(emin,emax,bg0,bg1,bg2,bg3);	
  TAxis *yax = h2->GetYaxis();
  yax->SetRangeUser(bg0,bg3);

  TH1D *htmp = h2->ProjectionX(), *h1;
  list->Add(htmp);
  
  Double_t theta, peak_mean, peak_err, cnts, cnts_rng;
  Double_t beta = 0.0953, raw_eng, raw_eng_err; // for improving doppler
  Double_t d2r = TMath::DegToRad();
  TF1 *func;// = new TF1("gauss_linbg_exc","pol1(0)+gaus(2)",bg0,bg3);
  Double_t mincnt = 30; // minimum total counts in a spectrum range bg0-bg3
    
  // loop over theta slices and fit the gamma peak for each slice
  for(int i=htmp->FindFirstBinAbove(mincnt); i<=htmp->FindLastBinAbove(mincnt); i++){
    
    theta = h2->GetXaxis()->GetBinCenter(i);
    h1 = h2->ProjectionY(Form("Gam_Theta%.0f",theta),i,i);
    // root can't fit without proper bin errors. Sumw2 doesn't seem to be setting this right
    for(int j=1; j<=h1->GetNbinsX(); j++)
      h1->SetBinError(j,sqrt(h1->GetBinContent(j)));
      
   // h1->SetTitle(Form("Gamma Spectrum For #theta_{TIG} = %.1f#circ; E_{#gamma} [keV]; Counts",theta));
    list->Add(h1);

    cnts = h1->Integral();
    if(cnts<mincnt){
    //  printf("\n%i\t Histogram did not have enough counts.. %.0f",i,cnts);
      continue;
    }
    cnts_rng = h1->Integral(h1->GetXaxis()->FindBin(emin),h1->GetXaxis()->FindBin(emax));    
    if(verbose) printf("\n%2i\t Histogram Counts = %.0f total and %.0f in range [%.0f - %.0f]..",i,cnts,cnts_rng,emin,emax);
    
    FitPeakExcludeRange(h1,emin,emax,bg0,bg1,bg2,bg3);  
//    h1->Fit(func,"R","",bg0,bg3);
    
    func = h1->GetFunction("gauss_linbg_exc");
    if(!func){
 //     printf("\n%i\t Function was not found.\n",i);
      continue;
    }
    
    peak_mean = func->GetParameter(3);
    peak_err  = func->GetParError(3);
    if(verbose){ printf("\n  \t Theta = %.0f deg.\n",theta);
      printf("\t\t\tAmpl = %.2f +/- %.2f keV\n",func->GetParameter(2),func->GetParError(2));
      printf("\t\t\tMean = %.2f +/- %.2f keV\n",peak_mean,peak_err);
      printf("\t\t\tSigm = %.2f +/- %.2f keV\n",func->GetParameter(4),func->GetParError(4));
      printf("\t\t\tChi2 / NDF = %.2f / %.2d = %.2f\n",func->GetChisquare(),func->GetNDF(),func->GetChisquare()/func->GetNDF());
    }
    gdop->SetPoint(gdop->GetN(),theta,peak_mean);
    gdop->SetPointError(gdop->GetN()-1,therr,peak_err);
    
    // remove AVERAGE doppler correction using 'beta'
    raw_eng = peak_mean*sqrt(1-beta*beta)/(1-beta*cos(theta*d2r));
    raw_eng_err = peak_err/peak_mean*raw_eng;
    
    graw->SetPoint(graw->GetN(),theta,raw_eng);
    graw->SetPointError(graw->GetN()-1,therr,raw_eng_err);  
  }
  
  std::string fexpr = "[1]*sqrt(1-[0]*[0])/(1-[0]*cos(x*3.1415/180))";
 // std::string fexpr = Form("%.2f*sqrt(1-[0]*[0])/(1-[0]*cos(x*3.1415/180))",egam);
  TF1 *fCorr = new TF1("DopplerCorrection",fexpr.c_str(),0,180);
  fCorr->SetParameter(0,beta);
  fCorr->SetParameter(1,egam);
  fCorr->SetParLimits(1,emin,emax);
  
  TCanvas *c3 = new TCanvas("c3","c3");
  graw->Draw("AP");  
  graw->Fit(fCorr);
  
  Double_t betafit = fCorr->GetParameter(0), corr, edop, edop_err;
  
  printf("\n\n Reconstructing doppler correction..\n");
  for(int i=0; i<graw->GetN(); i++){
    theta = graw->GetX()[i];
    therr = graw->GetEX()[i];
    raw_eng = graw->GetY()[i];
    raw_eng_err = graw->GetEY()[i];    
    corr  = 1/sqrt(1-betafit*betafit)*(1-betafit*cos(theta*d2r));
    edop = raw_eng*corr;
    edop_err = raw_eng_err*corr;
    if(verbose) printf("\n%.2i\t Theta = %.0f deg.\tCorr = %4f\tEdop = %.2f +/- %.2f keV\n",i,theta,corr,edop,edop_err);
    gdop2->SetPoint(i,theta,edop);
    gdop2->SetPointError(i,therr,edop_err);
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  gdop->Draw("AP");
//  TCanvas *c2 = new TCanvas("c2","c2");
  gdop2->Draw("same P");
  
  list->Add(gdop);
  list->Add(graw);  
  list->Add(gdop2);   
  
  printf("\n\n COMPLETE \n\n");
  return list;
}

TF1 *TTigressAnalysis::CorrelationFunction(Double_t par0, Double_t par1, Double_t par2){
	
	TF1 *func = new TF1("CorrFunc","[0] + [1]*pow(cos(x*TMath::DegToRad()),2.0) + [2]*pow(cos(x*TMath::DegToRad()),4.0)",0,180);
	func->SetParameters(par0,par1,par2);
			
	return func;
}


TCanvas *TTigressAnalysis::SetEfficiencyCurve(const char *efname, Double_t engabs, Double_t abseff, Double_t abserr){
  // MAKE THE EFFICIENCY CURVE FIT ERROR SIMPLY THE VARIANCE OF THE DATA ABOUT THE FIT

  const char *absmsg, *addmsg;
  if(engabs>0 && abseff>0 && abserr>=0){
    absmsg = "Absolute ";
    addmsg = "Addback ";
  } else {
    absmsg = " ";
    addmsg = " ";    
  }
  printf("\n\n Setting %s%sEfficiency Curve :\n",absmsg,addmsg);
    
  gTigEff = new TGraphErrors();    
  gTigEff->SetNameTitle("EffData",Form("TIGRESS %s%sEfficiency Curve; E_{#gamma} [keV]; %s%sEfficiency [%%]",absmsg,addmsg,absmsg,addmsg));

  std::ifstream infile(efname);
  Double_t eng, eff, err;
  while(infile.good()){
    infile >> eng >> eff >> err;

    gTigEff->SetPoint(gTigEff->GetN(),eng,eff);
    gTigEff->SetPointError(gTigEff->GetN()-1,0,err);
  }
  infile.close();
  if(gTigEff->GetN()<3){
    printf("\t -> Insufficient data in %s [%i points]!\n\n",efname,gTigEff->GetN());
    gTigEff->Clear();
    TCanvas *cfail = new TCanvas();
    return cfail;
  } else printf("\t -> Complete [added %i points]\n",gTigEff->GetN());

	// Relative Efficiency		
	fTigEff = new TF1("func_eff","[0]*pow(10.,[1]*log10(x)+[2]*pow(log10(x),2.)+[3]*pow(1/x,2.))",0,4000);
	fTigEff->SetNpx(2000);
	fTigEff->SetLineWidth(2);
	fTigEff->SetParameters(35,-0.14,0.052,-2308);	
	
//	fTigEff->SetParameter(0,fTigEff->GetParameter(0)/fTigEff->GetMaximum());
  fTigEff->SetNameTitle("EffFit",Form("TIGRESS %s%sEfficiency Curve; E_{#gamma} [keV]; %s%sEfficiency [%%]",absmsg,addmsg,absmsg,addmsg));
  
   // fit data
  gTigEff->Fit(fTigEff,"QM");  
  Double_t *xx = gTigEff->GetX(), *yy = gTigEff->GetY(), *ye = gTigEff->GetEY();
  
  Double_t chi2=0;
  for(int i=0; i<gTigEff->GetN();i++)
    chi2 += pow(yy[i]-fTigEff->Eval(xx[i]),2.0)/ye[i];
  // if(verbose) printf("\n\n\n \tChi2 = %f, Chi2/NDF = %f\n\n",chi2,chi2/(gTigEff->GetN()-1));   
  
  TCanvas *canvas = new TCanvas("EfficiencyCurve","EfficiencyCurve",800,500); 
  canvas->SetGrid();
  gTigEff->Draw("AP");  
  
  // calculate fit error using only absolute scaling... approx, but still helpful
  TF1 *flo = (TF1*)fTigEff->Clone("flower_abs"); 
  flo->SetLineWidth(1);flo->SetLineStyle(2); flo->SetLineColor(3);
  TF1 *fup = (TF1*)fTigEff->Clone("fupper_abs");
  fup->SetLineWidth(1);fup->SetLineStyle(2); fup->SetLineColor(3); 

  // use absolute scaling factor uncertainty to make confidence band
  Double_t conf = fTigEff->GetParError(0);
  flo->SetParameter(0,fTigEff->GetParameter(0)-conf);
  fup->SetParameter(0,fTigEff->GetParameter(0)+conf);  
//  printf("\n\tval = %f +/- %f\n\n",fTigEff->GetParameter(0),conf);
  
  if(engabs<=0 || abseff<=0 || abserr<0)
    return canvas;
    
  // get the scaling factor for the fit to reproduce the required point
  Double_t abscale = abseff/fTigEff->Eval(engabs);
  
  // now apply absolute scaling
  for(int i=0; i<gTigEff->GetN(); i++){  
    gTigEff->SetPoint(i,xx[i],yy[i]*abscale);      
    gTigEff->SetPointError(i,0.0,ye[i]*abscale); 
  }
  
  // refit absolute data
  fTigEff->SetParameter(0,fTigEff->GetParameter(0)*abscale);
  gTigEff->Fit(fTigEff,"RQEM","",0,4000);
  
  int ngr = gTigEff->GetN();
  TGraphErrors *grint = new TGraphErrors(ngr);
  for(int i=0; i<ngr; i++)
    grint->SetPoint(i, gTigEff->GetX()[i], 0);

  //Compute the confidence intervals at the x points of the created graph
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint);
  TGraph *gup = new TGraph(ngr);
  TGraph *glo = new TGraph(ngr);  
  double val;
  for(int i=0; i<ngr; i++){
  // val is the sum of fit error with scaling error
    val = grint->GetY()[i]*sqrt(pow(grint->GetEY()[i]/grint->GetY()[i],2.0)+pow(abserr/abseff,2.0));
    gup->SetPoint(i,gTigEff->GetX()[i],grint->GetY()[i]+val);
    glo->SetPoint(i,gTigEff->GetX()[i],grint->GetY()[i]-val);
  }


  // make confidence band larger due to uncertainty on absolute scaling
 // conf = fTigEff->GetParameter(0)*abserr/abseff;
  conf = fTigEff->GetParError(0);
  
  // MAKE THE EFFICIENCY CURVE FIT ERROR SIMPLY THE VARIANCE OF THE DATA ABOUT THE FIT

  fTigEffLow = (TF1*)fTigEff->Clone("flower_abs");  
  fTigEffLow->SetLineWidth(1);fTigEffLow->SetLineStyle(3);
 // flo->SetParameter(0,fTigEff->GetParameter(0)-conf);
 // flo->SetParameter(1,fTigEff->GetParameter(1)-fTigEff->GetParError(1));
  
  fTigEffUp = (TF1*)fTigEff->Clone("fupper_abs");
  fTigEffUp->SetLineWidth(1);fTigEffUp->SetLineStyle(3); 
 // fup->SetParameter(0,fTigEff->GetParameter(0)+conf);   
 // fup->SetParameter(1,fTigEff->GetParameter(1)+fTigEff->GetParError(1));
  gup->Fit(fTigEffUp,"RQEM","",0,4000);
  glo->Fit(fTigEffLow,"RQEM","",0,4000);

  //fTigEff->Draw("e3");
  gTigEff->Draw("AP");
//  gup->Draw("same P");
//  glo->Draw("same P");  
  fTigEffLow->Draw("same");
  fTigEffUp->Draw("same");
 
   //Now the "grint" graph contains function values as its y-coordinates
   //and confidence intervals as the errors on these coordinates
   //Draw the graph, the function and the confidence intervals 
 //  grint->SetLineColor(kGreen);
 //  grint->SetFillColor(kGreen);
 //  grint->SetFillStyle(1);
 //  grint->Draw("e3 same");
 //  grint->Print();
   
   
  return canvas;
}

Double_t TTigressAnalysis::Efficiency(Double_t eng){

	if(!fTigEff){
		printf("\n\t Warning :  Tigress Efficiency has not been set!\n\n");
		return 0.0;
	}
	return fTigEff->Eval(eng)*0.01; 
}

Double_t TTigressAnalysis::EfficiencyError(Double_t eng){
  
  Double_t eff = Efficiency(eng);
  if(eff==0)
    return 0.0;
     
//    return eff*fTigEff->GetParError(0)/fTigEff->GetParameter(0); 
  return 0.5*(fTigEffUp->Eval(eng)-fTigEffLow->Eval(eng))*0.01;
//  Double_t relerr = 0.0906;
//  return Efficiency(eng)*relerr;
}


void TTigressAnalysis::SetResolutionCurve(){
  /*
  TGraphErrors *g = new TGraphErrors();
  g->SetNameTitle("TigressResolutionGraph","TIGRESS Resolution; E_#gamma [keV]; #sigma [keV]");
  g->SetPoint(g->GetN(),414 ,	3.369	); g->SetPointError(g->GetN()-1,10.0, 0.106);
  g->SetPoint(g->GetN(),815	, 5.235	); g->SetPointError(g->GetN()-1,10.0, 0.043);
  g->SetPoint(g->GetN(),978	, 5.911	); g->SetPointError(g->GetN()-1,10.0, 0.177);
  g->SetPoint(g->GetN(),1035, 5.849	); g->SetPointError(g->GetN()-1,10.0, 0.238);
  g->SetPoint(g->GetN(),1180,	6.618	); g->SetPointError(g->GetN()-1,10.0, 0.440);
  g->SetPoint(g->GetN(),1305,	6.706	); g->SetPointError(g->GetN()-1,10.0, 0.387);
  g->SetPoint(g->GetN(),1400,	6.784	); g->SetPointError(g->GetN()-1,10.0, 0.902);
  g->SetPoint(g->GetN(),1500,	7.992	); g->SetPointError(g->GetN()-1,10.0, 1.194);
  g->SetPoint(g->GetN(),1762,	10.741); g->SetPointError(g->GetN()-1,10.0,	1.303);
  g->SetPoint(g->GetN(),2084,	11.473); g->SetPointError(g->GetN()-1,10.0,	0.611);
  g->SetPoint(g->GetN(),3500,	19.641); g->SetPointError(g->GetN()-1,10.0,	2.348);
  */
	// with appropriate peak widths we can identify multiple adjacent peaks	
	fTigSigma = new TF1("func_sig","[0]+[1]*TMath::Power(x,[2])",0,4000); // power law?
	fTigSigma->SetParameters(-4.74,0.679,0.409); // from fitting singles data	
	fTigSigma->SetNpx(4000);
	fTigSigma->SetTitle("TIGRESS Resolution Curve; E_{#gamma} [keV]; #sigma [keV]");

  return;
}

TF1 *TTigressAnalysis::GetResolutionCurve(){

  if(!fTigSigma)
    SetResolutionCurve();
  
 // TF1 *fun =  (TF1*) fTigSigma->Clone("TigResolution");
  return fTigSigma;
}

Double_t TTigressAnalysis::Resolution(Double_t eng, Bool_t FWHM){

	if(!fTigSigma){
		printf("\n\t Warning :  Tigress Resolution has not been set!\n\n");
		return 0.0;
	}
	Double_t sig = fTigSigma->Eval(eng);
	if(FWHM)
	  sig*=2.355;
  return sig;
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

Double_t TTigressAnalysis::pol_bg_exc(Double_t *x, Double_t *par){
/*
	// excludes the region between backgrounds and peak
  par[0]   background constant
  par[1]   background slope
  par[1]   background quad
  par[2]	 bglo max point
  par[3]	 peak min point
	par[4]	 peak max point
  par[5]	 bghi min point      
*/	   
  Double_t fitval = par[0] + x[0]*par[1] + x[0]*x[0]*par[2];

  // no peak, just background in bg regions
  if(x[0]>par[3] && x[0]<par[6]){
    TF1::RejectPoint();
    return fitval;
  }
	   
  return fitval;
}

TH1D *TTigressAnalysis::TH1Sum(TH1D *ha, TH1D *hb, Double_t sca, Double_t scb){
	
	if(!ha->GetSumw2N()) ha->Sumw2();
	if(!hb->GetSumw2N()) hb->Sumw2();
	  	
	TH1D *h = (TH1D*) ha->Clone(Form("%s_%s",ha->GetName(),hb->GetName()));
	h->Scale(sca);
	h->Add(hb,scb);
	return h;
}

TH2F *TTigressAnalysis::TH2Sum(TH2F *ha, TH2F *hb, Double_t sca, Double_t scb){

	if(!ha->GetSumw2N()) ha->Sumw2();
	if(!hb->GetSumw2N()) hb->Sumw2();
	
	TH2F *h = (TH2F*) ha->Clone(Form("%s_%s",ha->GetName(),hb->GetName()));
	h->Scale(sca);
	h->Add(hb,scb);
	return h;
}

TH1D *TTigressAnalysis::TH2Proj(TH2F *h, char ax, Double_t minval, Double_t maxval, Double_t &sz){
  
	if(!h->GetSumw2N()) h->Sumw2();
  
	TAxis *axis;
	if(ax=='x')
		axis =	h->GetYaxis();
	else
		axis =	h->GetXaxis();
	
  int b[2] = {axis->FindBin(minval), axis->FindBin(maxval)};
	double c[2] = {axis->GetBinCenter(b[0]), axis->GetBinCenter(b[1])};
//	printf("\n\t vals[] = {%.1f, %.1f}  b[] = {%i, %i}, c[] = {%.1f, %.1f}",minval,maxval,b[0],b[1],c[0],c[1]);
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
	if(!h2->GetSumw2N()) h2->Sumw2();
	
	TH2F *h2safe = (TH2F*)h2->Clone("SafeCopy");
	h2safe->SetName(Form("%s_py_%iTo%i",h->GetName(),b[0],b[1]));
	return h2safe;
}


////////////////////////////////////////////////////////////////////////////////
// NOTES

// Gate on exc energy range produce gamma spectrum. Compare to data and eliminate what we dont see experimentally
// normalize theory to specific peak (direct to ground if possible)

// ACCOUNTS FOR THE FACT THAT THE REAL DATA HAS GAMMAS FROM STATES OUTSIDE OF THE SELECTED RANGE
// BY INCLUDING ALL STATES THAT ARE UP TO ±1.5 SIGMA BEYOND THE SPECIFIED EXC RANGE AND
// WEIGHTING THEIR GAMMA SPECTRA BY THE COUNTS EXPECTED WITHIN THE EXC RANGE [MIN. ~7%]

////////////////////////////////////////////////////////////////////////////////

Bool_t TTigressAnalysis::InitLevelsGammas(Int_t nmax, Bool_t verb){
//	gStyle->SetOptStat(0);
//	gStyle->SetTitleXOffset(1.3);
//	gStyle->SetTitleYOffset(1.5);	
	
	SetVerbose(verb);
	Bool_t success = LoadLevelsGammas(Form("%s/Sr96_LevelsGammas.txt",TIGDIR),nmax);
	
	if(!fTigSigma) SetResolutionCurve();
  if(!gTigEff)   SetEfficiencyCurve();
  
	return success;
}

Bool_t TTigressAnalysis::LoadLevelsGammas(std::string fname, int nmax){

	ClearVars();
	
	std::ifstream infile(fname.c_str());	
	if(!infile.is_open()){
		printf("\nFAILED TO OPEN FILE :  %s\n\n",fname.c_str());
		return false;
	}
	
	nndcfile.assign(fname);
	
	nlines = std::count(std::istreambuf_iterator<char>(infile),std::istreambuf_iterator<char>(), '\n')+1;
	infile.close();
	infile.open(fname.c_str());	
	
	GetRid("TmpMat");
	GetRid("TmpMat2");
	TH2F *htmp = new TH2F("TmpMat","TmpMat",nlines,0,nlines,nlines,0,nlines);
	TH2F *htmp2 = new TH2F("TmpMat2","TmpMat2",nlines,0,nlines,10000,0,10000);

	Double_t initial, final, intensity, egamma, maxgam=0;
	int statenum;
	nstates = 0;
	ncascades = 0;
	
	printf("\n Reading Input File ' %s ' :\n",fname.c_str());
	if(verbose) printf("\n\t   State Energy\tGamma Energy\tIntensity[%%]\tFinal State [num]");
	while(infile.good()){
		
		infile >> initial >> egamma >> intensity >> final;

		GetStateIndex(initial,true); // add new state, but don't duplicate		

		statenum = GetStateIndex(final)+1;
		
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
	printf("\n\t --- Complete! [eof = %i  fail = %i  bad = %i]\n\n",infile.eof(),infile.fail(),infile.bad());
	infile.close();
	
	GetRid("TransitionMat");
	htrans = new TH2F("TransitionMat","",nstates,0,nstates,nstates,0,nstates);
	htrans->SetTitle("Transition Intensity Matrix; Initial State; Final State");
	
	for(int i=1; i<=nstates; i++) // go to the maximum row
		for(int j=1; j<=i; j++) // go up to the diagonal
			htrans->SetBinContent(i,j,htmp->GetBinContent(i,j));
	
	htrans->Scale(1./htrans->GetMaximum()); // normalize to max strength 1

	GetRid("GammasMat");
	hgams = new TH2F("GammasMat","",nstates,0,nstates,6000.0,0,6000.0);
	hgams->SetTitle("Gammas Emitted Matrix; State; E_{#gamma} [keV]");

	for(int i=1; i<=nstates; i++){ // go to the maximum row
		TH1D *hprojy = htmp2->ProjectionY("",i,i);
		for(int j=hprojy->FindFirstBinAbove(); j<=hprojy->FindLastBinAbove(); j++){ // go up to the diagonal
			if(htmp2->GetBinContent(i,j))
				hgams->SetBinContent(i,j,htmp2->GetBinContent(i,j));
		}
	}	
	hgams->Scale(1./hgams->GetMaximum()); // normalize all transitions to max strength 1
		
	for(int i=1; i<=nstates; i++)
		BuildDecayScheme(i);
				
	return true;
}

std::vector<int> TTigressAnalysis::PrintStates(Double_t emin, Double_t emax){

	std::vector<int> vals;
	if(!energies.size()){
		printf("\n\n No States Have Been Loaded !!\n\n");
		return vals;
	}
	if(!hint) SetIntensity(0,1);
	
	printf("\n\n List of States:-");

	for(int i=0; i<(Int_t)energies.size(); i++){
		if(energies.at(i)<emin)
			continue;
					
		if(energies.at(i)>emax)
			break;
						
		vals.push_back(i+1);	
		printf("\n\tState %2i =   %6.1f keV",i+1,energies.at(i));
		if(verbose) printf("\t\t[ Pop. Strength = %.3E ]",hint->GetBinContent(i));
	}
	printf("\n\n");
	return vals;
}

std::vector<int>  TTigressAnalysis::PrintCascades(Int_t from_state, Double_t egam, Bool_t printeng){

	std::string msgam = "";
	if(egam) msgam = Form("with %.1f keV gamma gate ",egam);
	printf("\n\n List of Decays from state %2i %s%s:-",from_state,msgam.c_str(),printeng?"{With Gamma Energy}":"");
	std::vector<int> indx = GetCascadeIndex(from_state,egam);
	Int_t ncas = indx.size();
	int j=0;
	
	for(int k=0; k<ncas; k++){

		j = indx.at(k);
		printf("\n\tCascade %2i = ",k+1);
		
		if(printeng) 
			printf(" %6.1f ",energies.at(cascades[j].at(0)-1));
		else 
			printf("  %2i",cascades[j].at(0));
			
		for(int i=1; i<cascades[j].size(); i++){
			if(printeng) {
				printf(" -{ %6.1f }->  %6.1f ",	energies.at(cascades[j].at(i-1)-1)-energies.at(cascades[j].at(i)-1)
																			,	energies.at(cascades[j].at(i)-1));
			} else
				printf("  -> %2i",cascades[j].at(i));
		}
		if(printeng)printf(" keV * * *");
	}
	
	printf("\n\n");
	return indx;
}

Double_t TTigressAnalysis::BranchingRatio(Double_t state_eng, Double_t egam){
// This allows gates on transitions such as B->C in the cascade A->B->C

	Int_t from_state = GetStateIndex(state_eng,false);
	if(!BuildDecayScheme(from_state+1))
		return 0;
	
	Double_t eng1 = energies.at(from_state);	
  TH1D *hg = DrawGammas(from_state+1); // uses gamma spectrum
  //hg->Draw();
  Double_t val = hg->Integral((int)egam-1,(int)egam+1);
	if(!val){
	  printf("\n\t Error :  The transition %.1f keV -> %.1f keV was not found in file.\n\n",state_eng,state_eng-egam);
	  return 0;
  }  
  
	if(verbose)printf("\n\t %.1f keV gamma is produced from %.1f keV state cascades with total intensity = %.2e\n",egam,eng1,val);
  
  return val;  
}		

void TTigressAnalysis::ClearVars(){
	
	nndcfile.clear();
	list = new TList();
	hseq = 0;
	hint = 0;
	nsequences = 0;

	states.clear();	
	energies.clear();
	htrans = 0;	
	hgams  = 0;
	nlines = 0;
	nstates = 0;	

}

void TTigressAnalysis::FixIntensities(Double_t emin, Double_t emax, Double_t egam, Int_t level, Double_t strength){

	static TH1D *hdata=0, *htmp, *hbg;
	static THStack *hstack;
	static TList *hlist;
	static Double_t eminval=0.0, emaxval=-10000.0, egamset;
	if(eminval!=emin || emaxval!=emax){
		hdata=0;
		eminval = emin;
		emaxval = emax;
	}
	static std::vector<double> intensities, total, str;
	static std::vector<int> substates;
	static TCanvas *c;	
	c = (TCanvas*) gROOT->FindObject("StackCanvas");

	if(!hdata || !c || egamset!=egam){
	
		if(!hint)
			SetIntensity(0,100.0);

		total.clear();
		intensities.clear();	
		substates.clear();
		str.clear();	
		
	// get theory spectrum with uniform population of all states	
	// this	spectrum is made up of delta function peaks
		TH1D *h = DrawGammasGated(emin,emax,egam);
		if(!h)
			return;
		
		egamset = egam;		
		Double_t effcorr = 1.0;
		if(egam){
			// perform basic gamma gate with subtraction regions
			hdata = GamGated(egam-3*fTigSigma->Eval(egam),egam+3*fTigSigma->Eval(egam),0,0,0,0,emin,emax);
			hdata->SetTitle(Form("Gamma Coincidences [w. %.1f keV] From Data With Exc Range %.1f - %.1f keV; E_{#gamma} [keV]; Counts / %.0f keV",egam,emin,emax,hdata->GetBinWidth(0)));
			effcorr = Efficiency(egam);
		}	else {
			hdata = Gam(emin,emax);
			hdata->Add(hdata->ShowBackground(20,"0N"),-1.0);
			hdata->SetTitle(Form("Gamma Singles From Data With Exc Range %.1f - %.1f keV; E_{#gamma} [keV]; Counts / %.0f keV",emin,emax,hdata->GetBinWidth(0)));
		}
		hdata->SetLineColorAlpha(kBlack,0.3);		

		TH1D *htmp2;
		hlist = new TList;	
		std::string s;		
				
		// make realistic copies of all the spectra
		for(int i=0; i<(Int_t)states.size(); i++){

			htmp = (TH1D*)list->FindObject(Form("GammaSpectrum_State%i%s_Strength",states.at(i),egam>0?Form("_GamGate%.1f",egam):""));
			if(!htmp)
				continue;
			
			s.assign(htmp->GetTitle());
			str.push_back(atof(s.substr(s.length()-6,5).c_str()));
									
			htmp2 = MakeRealistic(htmp,true);
			htmp2->Scale(effcorr); // coincidences lower the efficiency further
						
			htmp2->SetLineColor(htmp->GetLineColor());
			htmp2->SetLineColorAlpha(htmp2->GetLineColor(),str.back()/100.0);
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
	hstack->SetTitle("Gamma Spectrum From Each State; E_{#gamma}^{thry} [keV]; Intensity [arb]");
	//hstack->Add(hbg);
	
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

		leg->AddEntry(htmp,Form("State %2i [%5.1f keV @ %2.f%%]: Strength %.1f",
		substates.at(i),energies.at(substates.at(i)-1),str.at(i),hint->GetBinContent(substates.at(i))),"l");	
	}		
	printf("\n\n");
	
	c->cd(1);
	hdata->Draw();			
	hstack->Draw("hist same"); // noclear?
	
	leg->Draw();
	c->cd(2);
	hint->GetXaxis()->SetRangeUser(substates.at(0)-2,substates.back()+1);
	hint->Draw();

}

Bool_t TTigressAnalysis::SetIntensity(Int_t state, Double_t strength){

  if(!energies.size()){
    printf("\n\t! Error : No states have been loaded yet. Intensities cannot be set now.\n\n"); 
    return false;
  }

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
	return true;
}	

Bool_t TTigressAnalysis::ReadIntensities(const char *fname){

	std::ifstream infile(fname);
	if(!infile.is_open()){
		printf("\n\t! Error : Could not locate file ' %s '\n\n",fname);
		return false;	
	}
	
	Bool_t success = SetIntensity(0,0);// initiate the histogram
	if(!success)
	  return false;
	  
	Int_t bin;
	Double_t eng, val;
	while(infile.good()){
		infile >> bin >> eng >> val;
		hint->SetBinContent(bin,val);
	}
	printf("\n\n\tRead in file ' %s ' and set TH1D*'hint' intensities.\n\n",fname);
	return true;
}

void TTigressAnalysis::WriteIntensities(const char *fname){

	if(!hint){
		printf("\n\t Error :  At least one intensity must be set!");
		return;
	}

	std::ofstream ofile(fname);
	
	for(int i=1; i<hint->GetNbinsX(); i++)
		ofile << i << "\t" << energies.at(i-1) << "\t" << hint->GetBinContent(i) << "\n";

	printf("\n\n\tWrote Contents of TH1D*'hint' to file!\n\n");
	ofile.close();
}

TH1D *TTigressAnalysis::DrawGammas(Int_t from_state, Double_t egam){
	
	const char *name = Form("GammaSpectrum_State%i%s",from_state,egam>0?Form("_GamGate%.1f",egam):"");
	TH1D *hgcalc = (TH1D*)list->FindObject(name);
	if(hgcalc){
		if(verbose)printf("\n %s already exists. Using this histogram.\n",hgcalc->GetName());
		return hgcalc;		
	}
	
	if(verbose)printf("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * *");	
	if(verbose)printf("\nDrawing gammas from state %i [%.1f keV] :-\n",from_state,energies.at(from_state-1));	
				
	if(!BuildDecayScheme(from_state))
		return 0;
	
	GetRid(name);
	hgcalc = new TH1D(name,"",4000.0,0,4000.0);
	const char *gnamefull = Form("With %.1f keV Gamma Gate",egam);	
	hgcalc->SetTitle(Form("Gammas Emitted From %.1f keV State [%i] %s; E_{#gamma}^{thry}[keV]; Intensity [# of #gamma's emitted per state];",energies.at(from_state-1),from_state,egam>0?gnamefull:""));
	hgcalc->SetLineColor(kRed);
	hgcalc->SetLineWidth(2);
	list->Add(hgcalc);
	
	std::vector<int> indx = GetCascadeIndex(from_state,egam);
	Int_t ncas = indx.size(), j=0, initial, final;
	Double_t gam, strength, intensity, scale;
		
	GetRid("tmp");
	TH1D *htmpi, *htmpf, *htmpg = new TH1D("tmp","",4000.0,0,4000.0);
	
	for(int k=0; k<ncas; k++){

		j	= indx.at(k);
		intensity = 1.0;
		if(verbose) printf("\n Drawing Cascade %2i :-",k+1);
		for(int i=0; i<(Int_t)cascades[j].size()-1; i++){

			initial = cascades[j].at(i);
			final = cascades[j].at(i+1);
			if(verbose) printf("\n\tDecay %2i -> %2i : ",initial,final);
			
			// remember that vector index = bin number -1
			gam = energies.at(initial-1)-energies.at(final-1);			
			// strength of transition [expressed as a ratio to 100%]
			strength = htrans->GetBinContent(initial,final);
			// project out all transition strengths from this state so that we can sum over all strengths
			htmpi = htrans->ProjectionY("transtmpi",initial,initial);				
			// sum all transition strengths to normalize and convert this transition to a probability
			strength /= htmpi->Integral(1,initial);
			
			// the transition strength is shared equally between each cascade
			htmpf = htrans->ProjectionY("transtmpf",final,final);				
			if(htmpf->GetEntries())
			  strength /= htmpf->GetEntries(); // reduce the strength of a specific cascade
			
			intensity*=strength;
			if(verbose) printf(" Gam. Eng. = %4.1f\tInt. = %5.4f \t Tot. Int. = %5.4f",gam,strength,intensity);
	
			htmpg->SetBinContent(gam,intensity);
		}
		scale = 1.0;
		if(egam>0){
			Int_t bin = floor(egam);
			scale = htmpg->GetBinContent(bin); // rescale by gamma intensity
			if(verbose){
				printf("\n\t---- Gamma Spectrum Contents ----");
				for(int m=htmpg->FindFirstBinAbove(); m<=htmpg->FindLastBinAbove(); m++){
					if(htmpg->GetBinContent(m))
						printf("\n\t  Bin = %i, Content = %.4f",m,htmpg->GetBinContent(m));
				}		
				printf("\n\tGAMMA GATE:  Bin = %i, Content = %.2f\n",bin,scale);		
			}
			htmpg->SetBinContent(bin,0.0); // we don't see gated transition
		}			
		hgcalc->Add(htmpg,scale);
		htmpg->Reset();
	}
	
	if(verbose){
		printf("\n\t---- TOTAL Spectrum Contents ----");
		Double_t tot = 0.0;
		for(int m=hgcalc->FindFirstBinAbove(); m<=hgcalc->FindLastBinAbove(); m++){
				if(hgcalc->GetBinContent(m)){
					printf("\n\t  Bin = %4i, Content = %.4f",m,hgcalc->GetBinContent(m));
			    tot+=hgcalc->GetBinContent(m);
        }
			}	
    printf("\n\t Total Content = %.4f",tot);			
		printf("\n\t --- Complete!\n\n");
	}		
	
	return hgcalc;
}

TH1D *TTigressAnalysis::DrawGammasGated(Double_t emin, Double_t emax, Double_t egam){

	const char *name = Form("GammaSpectra%s%s",emin>=0&&emax>emin?Form("_Exc%.1fTo%.1f",emin,emax):"",egam>0?Form("_Gam%.1f",egam):"");
	TH1D *h = (TH1D*)list->FindObject(name);
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
	
	if(emin<0 && emax<0){
		emin_ext = 0.0;
		emax_ext = energies.back()+500.0;
		extend_excrng = false;
	}
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
	hstack->SetTitle(Form("Gamma Spectrum From Selected States %s; E_{#gamma}^{thry} [keV]; Intensity [# of #gamma's per state]",egam>0?Form("With %.1f keV Gam Gate",egam):""));
	
	h = new TH1D(name,"",4000.0,0,4000.0);
	h->SetTitle(Form("Gammas Emitted From %5.1f-%5.1f keV State [%i-%i], %s; E_{#gamma}^{thry} [keV]; Intensity [# of #gamma's per state];",
			energies.at(states.at(0)-1),energies.at(states.back()-1),states.at(0),states.back(),egam>0?Form("With %.1f keV Gam Gate",egam):""));
	h->SetLineColor(kRed);
	h->SetLineWidth(2);
	
	
	TF1 *func = new TF1("gaus","gaus",0,8000.0);
	func->SetParameters(1.0/(sqrt(2*3.1415926)*ExcSig),0.0,ExcSig);		
	
	Double_t strength;
	TLegend *leg = new TLegend(0.6,0.4,0.98,0.92);
	leg->SetFillStyle(1001);
	
	TH1D *hgcalc;
	Int_t n=1;
	for(int i=0; i<(Int_t)states.size(); i++){
		hgcalc = DrawGammas(states.at(i),egam);
		if(!hgcalc)
			continue;

		// calculate the contribution from each state is its integral intensity over excitation gate range 
		// take each state to be a normalized gaussian 
		func->SetParameter(1,energies.at(states.at(i)-1));					
		// rescale this spectrum to account for expected contribution in emin-emax range
		strength = func->Integral(emin,emax);
		
		GetRid(Form("%s_Strength",hgcalc->GetName()));
		htmp = (TH1D*) hgcalc->Clone(Form("%s_Strength",hgcalc->GetName()));
		htmp->Scale(strength);
		list->Add(htmp);

		// keep track of the strength 
		htmp->SetTitle(Form("%s @ Strength %4.1f%%",htmp->GetTitle(),strength*100));
		
		// all low [out of range] state gammas are black and all high [o.o.r.] gammas are blue
		if(energies.at(states.at(i)-1)<emin || energies.at(states.at(i)-1)>emax)
			htmp->SetLineColorAlpha(kBlack,strength);		
		else{
		 	if(n==9)
		 		n=19;
			htmp->SetLineColor(++n);
		}	
		
		h->Add(htmp);
		hstack->Add(htmp);
	//	leg->AddEntry(htmp,Form("State %2i, %6.1f keV [%2.f%%]",states.at(i),energies.at(states.at(i)-1),strength*100),"lp");		
		leg->AddEntry(htmp,Form("State %2i [%5.1f keV @ %2.f%%]",states.at(i),energies.at(states.at(i)-1),strength*100.0),"l");
	}

	GetRid(h->GetName(),false);
	list->Add(h);	

	GetRid(hstack->GetName(),false);
	list->Add(hstack);	

	GetRid("CanvasEnergyRange");
	TCanvas *c = new TCanvas("CanvasEnergyRange","Canvas",1100,500);
	c->Divide(2,1);
	c->cd(1);
	GetRid(Form("%s_Clone",hstack->GetName()));	
	THStack *hstack2 = (THStack*)hstack->Clone(Form("%s_Clone",hstack->GetName()));	
	hstack2->Draw();
	leg->Draw();
	
	c->cd(2);
	GetRid(Form("%s_Clone",h->GetName()));		
	TH1D *hh = (TH1D*)h->Clone(Form("%s_Clone",h->GetName()));	
	hh->DrawCopy();
	
	return h;
}

TH2F *TTigressAnalysis::DrawExcGam(Double_t egam, Bool_t use_int){

	const char *name = Form("ExcGam%s%s",egam>0?Form("_Gam%.1f",egam):"",use_int?"_UsingSetIntensities":"");
	TH2F *h = (TH2F*)list->FindObject(name);
	if(h){
		printf("\n %s already exists. Using this histogram.\n",h->GetName());
		return h;		
	}

	printf("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");	
	printf("\n\nDrawing excitation energy verus gamma energy :-");	
	
	Double_t emax = energies.back()+500.0;	
	states = PrintStates(0.0,emax);
	if(!states.size()){
		printf("\n\t! Error : No states in this range\n\n");
		return 0;
	}

	h = new TH2F(name,"",2000,0,4000.0,emax/40.0,0.0,emax);
	h->SetTitle(Form("Excitation Versus Gamma %s; E_{#gamma}^{thry} [keV]; E_{exc} [keV];",
	                egam>0?Form("With %.1f keV Gam Gate",egam):""));		

	TAxis *yax = h->GetYaxis();
	Double_t binexwid = yax->GetBinWidth(0), exeng, energy, scale=1.0, val=1.0;	
	Int_t binexlo, binexhi;
	TF1 *func = new TF1("gaus","gaus",0,8000.0);
	func->SetParameters(1.0/(sqrt(2*3.1415926)*ExcSig),0.0,ExcSig);		
		
	TH1D *hgcalc;
	printf("\n\t Progress...\r");
	for(int i=0; i<(Int_t)states.size(); i++){
    printf("\t Progress... %4.2f %%\r",(Double_t)i/(Double_t)states.size()*100.0);     
	  fflush(stdout);
		hgcalc = DrawGammas(states.at(i),egam);
		if(!hgcalc)
			continue;
    energy = energies.at(states.at(i)-1);
    hgcalc = MakeRealistic(hgcalc); // include realistic gamma peak width

    if(use_int && hint)
      scale = hint->GetBinContent(states.at(i));

		// take each state to be a normalized gaussian 
		func->SetParameter(1,energy);	     
    binexlo = yax->FindBin(energy-2.5*ExcSig);  
    binexhi = yax->FindBin(energy+2.5*ExcSig);
    for(int k=binexlo; k<=binexhi; k++){
      exeng = yax->GetBinCenter(k);
      for(int j=hgcalc->FindFirstBinAbove(); j<=hgcalc->FindLastBinAbove(); j++){
        val = hgcalc->GetBinContent(j)*func->Eval(exeng)*binexwid;
        if(val)
          h->Fill(hgcalc->GetBinCenter(j),exeng,val*scale);
      }
    }
	}
  printf("\t Progress... 100.00%%     COMPLETE\n\n");
	list->Add(h);	

	return h;
}

TCanvas *TTigressAnalysis::DrawDecayMats(Int_t from_state, Bool_t engaxis){
	
	TCanvas *c = new TCanvas("Canvas","Canvas",1100,500);
	
	if(!htrans || !cascades.size()){
		printf("\n\tError :  Read Input File First!\n");
		return c;
	}
	
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
	hdecay->SetTitle(Form("Cascades From State %i : %.1f keV ; Cascade Number; State %s",from_state,energies.at(from_state-1),engaxis?"E_{#gamma}^{thry} [keV]":"Number"));
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

Int_t TTigressAnalysis::GetStateIndex(Double_t val, bool add_element){ // compares a value to contents of vector to see if it exists as an element
	for(int i=0; i<(Int_t)energies.size(); i++){
		if(fabs(val-energies.at(i))<1.0) // within 1 keV means same state
			return i;
	}
	
	// if element hasn't been found we can add it in
 	if(add_element){
 		energies.push_back(val);
		nstates++; 		
 		if(verbose) printf("\n State %2i :  \r",nstates);
 		return energies.size()-1;
 	}
	return -1;
}

std::vector<int> TTigressAnalysis::GetCascadeIndex(Int_t from_state, Double_t egam){

	std::vector<int> indx;
	int s1 = 0, s2 = 0;
	if(egam>0){
		TH1D *htmp;
		Double_t binval;
		Int_t bin = round(egam);
		for(int i=1; i<=from_state; i++){
			htmp = hgams->ProjectionY("",i,i);
			/*
			for(int k=htmp->FindFirstBinAbove(); k<=htmp->FindLastBinAbove(); k++){
				binval = htmp->GetBinContent(k);
				if(binval) printf("\n\t state %i @ %.1f keV [bin %i] = %.1f",i,htmp->GetBinCenter(k),k,binval);
			}
			*/
			binval = htmp->GetBinContent((int)egam);
		
			if(binval){
				s2 = i;
				s1 = GetStateIndex(energies.at(s2-1)-egam)+1; // indx + 1 = state number 
				//printf("\n\t State %i @ %.1f keV [bin %i] = %.1f",i,egam,htmp->GetXaxis()->FindBin(egam),binval);				
				if(verbose)printf("\n\t >--{ Found transition : State %i -> %i }--> ",s2,s1);
			}
		}
		if(!s2 || s1<0){
			if(verbose)printf("\n\t No Gamma @ %.1f keV found in cascades from state %i!\n",egam,from_state);
			return indx;
		} 
	}
	
	for(int i=0; i<cascades.size(); i++){
		if(cascades[i].at(0) == from_state){
			if(s1 && s2){
				// search through cascade for s2->s1 transition
				for(int j=0; j<(Int_t)cascades[i].size()-1; j++){
					if(cascades[i].at(j) == s2){ // if s2 is found then s1 MUST be next
					 if(cascades[i].at(j+1) ==s1)
					 		indx.push_back(i);
					 else 
					 	break;
					}
				}					 	
			} else{
				indx.push_back(i);
			}
		}
	}
	
	return indx;		
}

Int_t TTigressAnalysis::BuildDecayScheme(int from_state){

	hseq = 0;
	nsequences = 0;	
	
	if(!htrans){
		printf("\n\t! Error : Read Input File First !\n\n");
		return 0;
	} else if(from_state>nstates){
		printf("\n\t! Error: Maximum allowed state = % i!\n\n",nstates);
		return 0;
	}
	
	// make a sequence hist for this state and see if there are already cascades built
	if(MakeSequenceHist(from_state)) // if cascade is already built then just return 
		return nsequences;
	
	std::vector<int> cas;  // index to show the current column on each row
	cas.push_back(from_state);
	int m = from_state, i=0;
	
	if(verbose) printf("\n Building Cascades from state %i..",from_state);
	for(i=1; i<=m; i++){
			
		if(i==from_state)
			break;
		
		if(i==m){ // end of loop for row m, so go back to previous row
			i = cas.at(cas.size()-1); // s[end] is previous column [the -1 is immediately incremented to zero..]
			m	= cas.at(cas.size()-2); // s[end-1] is previous row
			cas.pop_back(); // remove last element
			continue;
		}
					
		if(htrans->GetBinContent(m,i)>0){
			cas.push_back(i);
			if(i==1) // ground state has been reached
				AddCascade(cas);
			m = i; // set current row
			i = 0; // restart loop on next row	(0 will be incremented immediately to 1)		
			continue;	
		}	
	}
	
	MakeSequenceHist(from_state);
	if(verbose) printf("\n\t --- Complete!\n");
	
	return nsequences;
}

void TTigressAnalysis::AddCascade(std::vector<int> states){

	cascades[ncascades] = states;
	ncascades++;
	
	if(verbose){
		printf("\n\t* Added Cascade %3i :  [ ",ncascades);
		for(int i=0; i<(Int_t)states.size(); i++)
			printf( "%i ",states.at(i));
		printf("]");
	}
	return;
}

Int_t TTigressAnalysis::MakeSequenceHist(Int_t from_state){
	
	std::vector<int> cas = GetCascadeIndex(from_state); 
	nsequences = cas.size();
	if(!nsequences)
		return 0;
	if(verbose) printf("\n\t State %i has %i cascades.\n",from_state,nsequences);	
	
	Int_t maxnstates=0;
	for(int i=0; i<nsequences; i++){
		///printf("\n\t Cascade %i :  ",cas.at(i));
		for(int j=0; j<cascades[cas.at(i)].size(); j++)
//printf(" -> %i",cascades[cas.at(i)].at(j));
		if(cascades[cas.at(i)].size()>maxnstates)
			maxnstates = cascades[cas.at(i)].size();
	}//printf("\n maxnstates = %i\n\n",maxnstates);
	
	const char *name = Form("CascadeMat_State%i",from_state);	
	GetRid(name);
	hseq = new TH2F(name,"",maxnstates,0,maxnstates,nsequences,0,nsequences);
	hseq->SetTitle("Cascade To Ground State; Series of States; Cascade Number");	
	hseq->SetLineWidth(2);
	
	Int_t j=0;
	for(int k=0; k<nsequences; k++){
		j = cas.at(k);
		for(int i=0; i<(Int_t)cascades[j].size(); i++){
		//	printf("\n\t Bin [ %i , %i ] = %i",i+1,k+1,cascades[j].at(i));
			hseq->SetBinContent(i+1,k+1,cascades[j].at(i));
			hseq->SetBinError(i+1,k+1,0.0001);
		}
	}
	
	return nsequences;
}

THStack *TTigressAnalysis::GetDecayScheme(Int_t from_state, Bool_t engaxis){
	
	MakeSequenceHist(from_state);	
	
	TH1D *hproj;
	THStack *hdecay = new THStack(Form("DecayScheme_State%i",from_state),Form("Cascades From State %i",from_state));
	
	Int_t binmax = hseq->ProjectionX()->FindLastBinAbove();
	for(int i=1; i<=binmax; i++){
		hproj = hseq->ProjectionY(Form("Step%i",i),i,i);
		
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
		h->SetBinContent(j,energies.at((Int_t)h->GetBinContent(j)-1));
	}
	return;
}

void TTigressAnalysis::DrawTransitions(int from_state, Bool_t engaxis){
		
	TArrow *arrow;
	TPaveText *pt;
	TH1D *h;
	Double_t dh=0.1, ypos;
	if(engaxis)
		dh=50;	
	
	for(int i=1; i<=nsequences; i++){
		h = hseq->ProjectionX("",i,i);
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
				pt = new TPaveText(i-0.95,ypos-dh,i-0.05,ypos+dh,"nb");
				TText *text = pt->AddText(Form("%6.1f keV",h->GetBinContent(j)-h->GetBinContent(j+1)));
				pt->SetFillColor(0);
				pt->SetTextColorAlpha(kBlack,0.5);
				text->SetTextAngle(25.);
		//		pt->SetTextSize(2);
				pt->Draw();
			}
		}
	}
	
	return;
}

TH1D *TTigressAnalysis::MakeRealistic(TH1D *hgcalc, Bool_t abseff){

	// converts delta function spectrum into sigma(E) and efficiency corrected spectrum
	GetRid(Form("%s_Realistic",hgcalc->GetName()));
	TH1D *htmp = new TH1D(Form("%s_Realistic",hgcalc->GetName()),hgcalc->GetTitle(),hgcalc->GetNbinsX(),0,4000);		
	htmp->GetXaxis()->SetTitle(hgcalc->GetXaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(hgcalc->GetYaxis()->GetTitle());	
	
	Double_t x, val, sig, eff;
	TF1 *func = new TF1("gauspeak","gaus",0,4000);
	func->SetNpx(4000);
		
	// integrate the counts under the data peaks 
	for(int i=1; i<hgcalc->GetNbinsX(); i++){
		
		val = hgcalc->GetBinContent(i);
		if(!val)
			continue;
		
		// define a gaussian centered at bin center with sigma = sig and integral = hcalc->GetBinContent(i)		
		x = hgcalc->GetBinCenter(i);
		sig = fTigSigma->Eval(x);
		eff = Efficiency(x);		
		
		func->SetParameters(eff*val/(sqrt(2*3.14159)*sig),x,sig);
		// now get this histogram and use it as a filter over the data
		htmp->Add((TH1D*)func->GetHistogram());
	}
	
	return htmp;
}

void TTigressAnalysis::GetRid(const char *name, Bool_t delete_all){
	TObject *obj;

	obj = list->FindObject(name);
	if(obj) list->Remove(obj);

	if(delete_all){
		obj = gROOT->FindObjectAny(name);
		if(!obj) return;
		/*if(strcmp(obj->ClassName(),"TCanvas")==0){
			TCanvas *c = (TCanvas*)obj;
			c->Close();
		}
		else */if(obj) obj->Delete();
	}
	
}

void TTigressAnalysis::CloneDefaults(){
  hegg = (TH3S*) hexcgamgam->Clone("ExcGamGam");
  hegt = (TH3S*) hexcgamthtig->Clone("ExcGamThetaTig");
  heg  = (TH2F*) hexcgam->Clone("ExcGam");
  hgg  = (TH2F*) hgamgam->Clone("GamGam");
  he   = (TH1D*) hexc->Clone("Exc");
  hg   = (TH1D*) hgam->Clone("Gam");
}

void TTigressAnalysis::ResetDefaults(){
  hegg->Reset(); 
  hegt->Reset(); 
  heg->Reset(); 
  hgg->Reset(); 
  he->Reset(); 
  hg->Reset(); 
}
