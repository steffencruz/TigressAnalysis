#ifndef TTIGRESSANALYSIS_H
#define TTIGRESSANALYSIS_H

#include<TObject.h>
#include<RTypes.h>
#include<TList.h>

#include<TH3F.h>
#include<TH2F.h>
#include<TH1D.h>
#include<THStack.h>
#include<TCanvas.h>

#include <string>

// Reads a root file which contains gamma-gamma-excitation energy matrices and carries out coincidence analyses
// * exc gamma gamma
// * gamma gamma 
// * exc gamma 

class TTigressAnalysis 	{
	
	public:

		TTigressAnalysis(const char *fname = "");
		~TTigressAnalysis();
	
		static void Print(Option_t * = "");			
		static void Clear(Option_t * = "");		
		static void LoadHistos(const char *fname);
		
		static void SetGamBinSz(Int_t wid) { gambinsz = wid; }
		static void SetExcBinSz(Int_t wid) { excbinsz = wid; }
		static Int_t GetGamBinSz() { return gambinsz; }
		static Int_t GetExcBinSz() { return excbinsz; }
						
		static TH1D *Gammas(Double_t emin=0.0, Double_t emax=4000.0);
		static TH1D *GamGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);

		static TH1D *ExcGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);
											
		static void AnalyzeGammas(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);
													
		static void FitPeakExcludeRange(TH1 *hist, 
											Double_t emin, Double_t emax, 
											Double_t bg0, Double_t bg1, 
											Double_t bg2, Double_t bg3);

		static void SetBackgroundLims(Double_t emin, Double_t emax, 
										Double_t &bg0, Double_t &bg1, 
										Double_t &bg2, Double_t &bg3);


	private: 
	
		static Double_t gaus_lbg_exc(Double_t *x, Double_t *par);

		static std::string histfile;
		static TH3F *hexcgamgam;
		static TH2F *hgamgam,*hexcgam;
		static TH1D *hgam,*hexc;
		
		static Int_t gambinsz;
		static Int_t excbinsz;	
		static Bool_t addback;	

		
		
		//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//
	public:
		
		static void InitGammasLevels(Int_t nmax=100, Bool_t verb=true);

		static void LoadLevelsGammas(std::string fname="Sr96_LevelsGammas.txt", Int_t nmax = 100);
	
	  static std::vector<int> PrintStates(Double_t emin=0.0, Double_t emax = 10000.0);
	  static Int_t PrintCascades(Int_t from_state, Bool_t printeng=false);
		static void ClearVars();		

	  static void SetVerbose(bool verb) { verbose = verb;}

	  static void FixIntensities(Double_t emin=0.0, Double_t emax=10000.0, Int_t level=1, Double_t strength=50.0);
	  static void SetIntensity(Int_t state, Double_t strength);
	  static void ReadIntensities(const char *fname = "SavedIntensities.txt");
	  static void WriteIntensities(const char *fname = "SavedIntensities.txt");

	  static TH1D *DrawGammas(Int_t from_state, Bool_t draw=false);
	  static TH1D *DrawGammasExcRange(Double_t emin=0.0, Double_t emax = 10000.0);
	  static TCanvas *DrawDecayMats(Int_t from_state, Bool_t engaxis = true);
	  
									
	private:
	
		static Int_t GetIndex(Double_t val, bool add_element=false);
	
	  static void AddSeq(std::vector<int> states);
	  static Int_t BuildDecayScheme(Int_t from_state);	  
	  static THStack *GetDecayScheme(Int_t from_state, Bool_t engaxis = true);

	  static void ConvertToEnergyAxis(TH1D *h);
	  static TH2F *ShrinkToFit(int from_state, TH2F *h); // resizes histograms so axis ranges are relevant
	  static void DrawTransitions(Int_t from_state, Bool_t engaxis);	  
	  static TH1D *MakeRealistic(TH1D *hgam);

		static std::string nndcfile;
		
		static Bool_t verbose;
		static TList *list;

		static Int_t nlines; 
		static Int_t nsequences;
		static Int_t nstates;		
		static Double_t ExcSig;
		static Double_t NSig;			
		static std::vector<double> energy;
		static std::vector<int> states;
		
		static TH2F *htrans, *hgams, *hseq;
		static TH1D *hint;
		static TF1 *TigSigma, *TigEfficiency;
		        
	ClassDef(TTigressAnalysis,0)
};


#endif
