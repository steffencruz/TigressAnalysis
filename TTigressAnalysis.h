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
		
		static void SetGamBinSz(Int_t wid) { gambinsz = wid; }
		static void SetExcBinSz(Int_t wid) { excbinsz = wid; }
				
		TH1D *Gammas(Double_t emin=0.0, Double_t emax=4000.0);
		TH1D *GamGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);

		TH1D *ExcGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);
											
		void AnalyzeGammas(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);
													
		void FitPeakExcludeRange(TH1 *hist, 
											Double_t emin, Double_t emax, 
											Double_t bg0, Double_t bg1, 
											Double_t bg2, Double_t bg3);

		static void SetBackgroundLims(Double_t emin, Double_t emax, 
										Double_t &bg0, Double_t &bg1, 
										Double_t &bg2, Double_t &bg3);
											
	private: 
	
		// 
		static Double_t gaus_lbg_exc(Double_t *x, Double_t *par);

		
		static std::string file;

		static TH3F *hexcgamgam;
		static TH2F *hgamgam;
		static TH2F *hexcgam;
		static TH1D *hgam;
		static TH1D *hexc;
		
		static Int_t gambinsz;
		static Int_t excbinsz;	
		static Bool_t addback;		
		        
	ClassDef(TTigressAnalysis,0)
};


#endif
