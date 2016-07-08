#ifndef TTIGRESSANALYSIS_H
#define TTIGRESSANALYSIS_H

#include<TObject.h>
#include<RTypes.h>
#include<TList.h>

#include<TH3S.h>
#include<TH2F.h>
#include<TH1D.h>
#include<THStack.h>
#include<TCanvas.h>

#include <string>
#include <map>

// Reads a root file which contains gamma-gamma-excitation energy matrices and carries out coincidence analyses
// WARNING : THE TH3F
// * exc gamma gamma
// * gamma gamma 
// * exc gamma 

// gamma gamma coincidences are now done on theory spectra too.
// there is no efficiency applied so far.
// it should be as simple as scaling the final spectrum by TigEff->Eval(egam)

class TTigressAnalysis 	{
	
	public:

		TTigressAnalysis(void);
		~TTigressAnalysis();
	
		static Bool_t Init();
		static void Print(Option_t * = "");			
		static void Clear(Option_t * = "");		
		static Bool_t LoadHistos(const char *fname, const char *reac="dp");
		
		static void SetGamBinSz(Int_t wid) { gambinsz = wid; }
		static void SetExcBinSz(Int_t wid) { excbinsz = wid; }
		static Int_t GetGamBinSz() { return gambinsz; }
		static Int_t GetExcBinSz() { return excbinsz; }
						
		static TH1D *Gam(Double_t exmin=-1.0, Double_t exmax=-1.0);
											
		static TH1D *GamGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);

		static TH1D *ExcGated(Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);

		static TH2F *ExcGamGated(Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);
																																	
		static TH2F *ExcThetaGated(Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);
																						
		static TCanvas *AnalyzeGammas(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);

		static TH1D *GamAngCorr(Double_t emin, Double_t emax,
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,		
											Double_t exmin=-1.0, Double_t exmax=-1.0);
											
		static TH2F *GamAngCorrMat(Double_t exmin=-1.0, Double_t exmax=-1.0);													
													
		static void FitPeakExcludeRange(TH1 *hist, 
											Double_t emin, Double_t emax, 
											Double_t bg0, Double_t bg1, 
											Double_t bg2, Double_t bg3);

		static void SetBackgroundLims(Double_t emin, Double_t emax, 
										Double_t &bg0, Double_t &bg1, 
										Double_t &bg2, Double_t &bg3);

		static Double_t Efficiency(Double_t eng);
		
		static TF1 *CorrelationFunction(Double_t par0=1.0,Double_t par1=1.0,Double_t par2=1.0);
		
	private: 

		static TH1D *TH1Sum(TH1D *ha, TH1D *hb, Double_t sca=1.0, Double_t scb=1.0);
		static TH2F *TH2Sum(TH2F *ha, TH2F *hb, Double_t sca=1.0, Double_t scb=1.0);

		static TH1D *TH2Proj(TH2F *h, char ax, Double_t minval, Double_t maxval, Double_t &sz);
		static TH2F *TH3Proj(TH3S *h, std::string ax, Double_t minval, Double_t maxval, Double_t &sz);

		static Double_t gaus_lbg_exc(Double_t *x, Double_t *par);

		static std::string histfile;
		static std::string reaction;
		static TH3S *hexcgamgam,*hexcthcmgam,*hexcgamthtig;
		static TH2F *hgamgam,*hexcgam;
		static TH1D *hgam,*hexc;
		
		static Int_t gambinsz;
		static Int_t excbinsz;	
		static Bool_t addback;	

		
		
		//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//	//
	public:
		
		static Bool_t InitLevelsGammas(Int_t nmax=100, Bool_t verb=true);

		static Bool_t LoadLevelsGammas(std::string fname="Sr96_LevelsGammas.txt", Int_t nmax = 100);
	
	  static std::vector<int> PrintStates(Double_t emin=0.0, Double_t emax=10000.0);
	  static std::vector<int> PrintCascades(Int_t from_state, Double_t egam=0.0, Bool_t printeng=false);
		static void ClearVars();		

	  static TCanvas *DrawDecayMats(Int_t from_state, Bool_t engaxis = true);
	  static void SetVerbose(bool verb) { verbose = verb;}

	  static void FixIntensities(Double_t exmin=0.0, Double_t exmax=10000.0, Double_t egam=0.0, Int_t level=1, Double_t strength=50.0);
	  static void SetIntensity(Int_t state, Double_t strength);
	  static void ReadIntensities(const char *fname = "SavedIntensities.txt");
	  static void WriteIntensities(const char *fname = "SavedIntensities.txt");

	  static TH1D *DrawGammas(Int_t from_state, Double_t egam=0.0);
		static TH1D *DrawGammasGated(Double_t exmin=-1.0, Double_t exmax=-1.0, Double_t egam=-1.0);
	  static TH1D *MakeRealistic(TH1D *hgam, Bool_t abseff = false);
									
	private:
	
		static Int_t GetStateIndex(Double_t val, bool add_element=false);
		static std::vector<int> GetCascadeIndex(Int_t from_state, Double_t egam=0.0);
	
	  static Int_t BuildDecayScheme(Int_t from_state);	  	
	  static void AddCascade(std::vector<int> states);
		static Int_t MakeSequenceHist(Int_t from_state);
	  static THStack *GetDecayScheme(Int_t from_state, Bool_t engaxis = true);

	  static void ConvertToEnergyAxis(TH1D *h);
	  static void DrawTransitions(Int_t from_state, Bool_t engaxis);	  
	  
	  static void GetRid(const char *name, Bool_t delete_all = true); // removes object

		static std::string nndcfile;
		
		static Bool_t verbose;
		static TList *list;

		static Int_t nlines; 
		static Int_t nsequences;
		static Int_t ncascades;
		static Int_t nstates;		
		static Double_t ExcSig;
		static Double_t NSig;			
		static std::vector<double> energies;
		static std::vector<int> states;
		static std::map<int,std::vector<int> > cascades;
		
		static TH2F *htrans, *hgams, *hseq;
		static TH1D *hint;
		static TF1 *TigSigma, *TigEfficiency;
		        
	ClassDef(TTigressAnalysis,0)
};


#endif
