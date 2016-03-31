#ifndef TTIGRESSANALYSIS_H
#define TTIGRESSANALYSIS_H

#include <TObject.h>
#include <TList.h>

#include <TH3F.h>
#include <TH2F.h>
#include <TH1D.h>

#include <string>

class TTigressAnalysis 	{
	
	public:

		TTigressAnalysis(std::string fname = "");
		~TTigressAnalysis();
	
		static void Print(Option_t * = "");			
		static void Clear(Option_t * = "");		
		static void Reset(Option_t * = "");
		static void LoadHistos(const char *fname);
				
		TH1D *Gammas(Double_t emin=0.0, Double_t emax=4000.0);
		TH1D *GammasGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);

	private: 
	
		static std::string file;

		static TH3F *hexcgamgam;
		static TH2F *hgamgam;
		static TH2F *hexcgam;
		static TH1D *hgam;
		static TH1D *hexc;
		
		static Int_t binwidth;
		static Bool_t addback;		
		        
	ClassDef(TTigressAnalysis,0)
};


#endif
