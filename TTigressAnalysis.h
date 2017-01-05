#ifndef TTIGRESSANALYSIS_H
#define TTIGRESSANALYSIS_H

#include<TObject.h>
#include<RTypes.h>
#include<TList.h>

#include<TH3S.h>
#include<TH2F.h>
#include<TH1D.h>
#include<TF1.h>
#include<TGraphErrors.h>
#include<THStack.h>
#include<TCanvas.h>

#include <string>
#include <map>

#ifndef TIGDIR
#define TIGDIR "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/TigressAnalysis"
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// _____________ TTigressAnalysis : a gamma-particle analysis class ____________
// 
//  TTigressAnalysis is used to analyze data by applying [in almost any combination] :-
//   1. Gamma-gates
//   2. Particle-gates
//   3. Angular-gates
// 
//  _____ Input ____
//  Requires a specific which is produced by MakeExcGamThetaMats.C
//  This file contains [but is not limited to] the following main data structures :-
//  * TH3S  exc-gamma-gamma
//  * TH2F  exc-gamma-thetacm 
//  * TH2F  exc-gamma 
//  * TH2F  gamma-gamma 
//  * TH1D  exc
//  * TH1D  gamma
//
// _____ Optional additions _____
//  Same as above, but for specific regions of SHARC [UQQQ, UBOX and DBOX].
//  Names and conventions will be taken care of by MakeExcGamThetaMats.C
//
//
// _____ Known Bugs _____
// * Background subtraction must include a region both above and below the peak 
// Intensities in calculated spectra may be incorrect.  
// States which decay through multiple cascades are adding the same transitions multiple 
// times, once for each cascade, which results in an overestimation of final intensity.
//      -  FIXED [16th Nov 2016]
//
//    Made by Steffen Cruz [2016]
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


class TTigressAnalysis 	{
	
	public:

		TTigressAnalysis(void);
		~TTigressAnalysis();
	
		static Bool_t Init();
		static void Print(Option_t * = "");			
		static void Clear(Option_t * = "");		
		static Bool_t LoadHistos(const char *fname, const char *reac="dp", const char *detsel="all");
		
		static void SetGamBinSz(Int_t wid) { gambinsz = wid; }
		static void SetExcBinSz(Int_t wid) { excbinsz = wid; }
		static Int_t GetGamBinSz() { return gambinsz; }
		static Int_t GetExcBinSz() { return excbinsz; }
		
		static Bool_t DetectorSelector(std::string opt="all");
    static std::string SelectedDetector(void){ return detsel; }
						
		static TH1D *Gam(Double_t exmin=-1.0, Double_t exmax=-1.0);
											
		static TH1D *GamGated(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);

		static TH2F *GamGam(Double_t exmin=-1.0, Double_t exmax=-1.0);

		static TH1D *ExcGated(Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);

		static TH2F *ExcGamGated(Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);
																																	
		static TH2F *ExcThetaGated(Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);

    static TCanvas *QuickAnalysis(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);	
																	
		static TCanvas *AnalyzeGammas(Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Double_t exmin=-1.0, Double_t exmax=-1.0);												

    static TH1D *FitTheory(TH1D *hdata, Double_t exc, Double_t egam,
                      Double_t emin=0.0, Double_t emax=0.0, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);

		static Double_t FitPeakStats(TH1 *hist, 
											Double_t emin, Double_t emax, 
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,
											Bool_t bg_only=false,Bool_t quad_fit=false);
													
		static void FitPeakExcludeRange(TH1 *hist, 
											Double_t emin, Double_t emax, 
											Double_t bg0, Double_t bg1, 
											Double_t bg2, Double_t bg3);

		static void FitBgExcludeRange(TH1 *hist, 
											Double_t emin, Double_t emax, 
											Double_t bg0, Double_t bg1, 
											Double_t bg2, Double_t bg3,
											Bool_t quad_fit=false);												

		static void SetBackgroundLims(Double_t emin, Double_t emax, 
								  		Double_t &bg0, Double_t &bg1, 
									  	Double_t &bg2, Double_t &bg3);

    static void SetPeakLims(Double_t egam, Double_t &emin, Double_t &emax);
  
    static Double_t GetPopulationStrength(TH1D *hist, Double_t exc, Double_t egam,
                                              Double_t emin=0.0, Double_t emax=0.0, 
                                              Double_t bg0=0.0, Double_t bg1=0.0, 
                                              Double_t bg2=0.0, Double_t bg3=0.0);
    
    static TList *CheckDopplerCorrection(Double_t egam, TH2F *h2=NULL,
                      Double_t emin=0.0, Double_t emax=0.0,
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0);
    
    static void SetResolutionCurve();
    
    static TF1 *GetResolutionCurve();
    
    static TCanvas *SetEfficiencyCurve(const char *efname=Form("%s/AddbackEfficiencyData.txt",TIGDIR),
                    Double_t engabs=815.0, Double_t effabs=7.33, Double_t abserr=0.04);
		
		static TGraphErrors *GetEfficiencyCurve(void){ return gTigEff; }
		
		static Double_t Efficiency(Double_t eng);
		
		static Double_t EfficiencyError(Double_t eng);	
		
    static Double_t Resolution(Double_t eng, Bool_t FWHM=false);
		
		static TH1D *GamAngCorr(Double_t emin, Double_t emax,
											Double_t bg0=0.0, Double_t bg1=0.0, 
											Double_t bg2=0.0, Double_t bg3=0.0,		
											Double_t exmin=-1.0, Double_t exmax=-1.0);
											
		static TH2F *GamAngCorrMat(Double_t exmin=-1.0, Double_t exmax=-1.0);				
		
    static TH1D *GetHitPattern(Int_t detmin=5, Int_t detmax=16, Int_t crymin=0, Int_t crymax=3, Int_t segmin=0, Int_t segmax=7);

		static TH1D *NormalizeThetaHits(TH1D *htheta, Int_t segval=0);

		static TF1 *CorrelationFunction(Double_t par0=1.0,Double_t par1=1.0,Double_t par2=1.0);
		
	private: 

		static TH1D *TH1Sum(TH1D *ha, TH1D *hb, Double_t sca=1.0, Double_t scb=1.0);
		static TH2F *TH2Sum(TH2F *ha, TH2F *hb, Double_t sca=1.0, Double_t scb=1.0);

		static TH1D *TH2Proj(TH2F *h, char ax, Double_t minval, Double_t maxval, Double_t &sz);
		static TH2F *TH3Proj(TH3S *h, std::string ax, Double_t minval, Double_t maxval, Double_t &sz);

		static Double_t gaus_lbg_exc(Double_t *x, Double_t *par);
		static Double_t pol_bg_exc(Double_t *x, Double_t *par);

    static Bool_t IncludeRegion(Int_t indx);
    static void CloneDefaults();
    static void ResetDefaults();    
    
		static std::string histfile;
		static std::string reaction;
		static std::string detsel;
		
		// general histograms [all SHARC]
		static TH3S *hexcgamgam, *hexcthcmgam, *hexcgamthtig;
		static TH2F *hexcgam, *hgamgam;
		static TH1D *hgam, *hexc;		
    // detector specific histos [UQQQ, UBOX, DBOX]
		static TH3S *hexcgamgam_sel[3];
		static TH3S *hexcgamthtig_sel[3];
		static TH2F *hexcgam_sel[3];
		static TH2F *hgamgam_sel[3];		
		static TH1D *hexc_sel[3];
		static TH1D *hgam_sel[3];
		// detector selected histograms
		static TH3S *hegg, *hegt;
		static TH2F *hgg, *heg;
		static TH1D *hg, *he;

		static Double_t gambinsz;
		static Double_t excbinsz;	
		static Bool_t addback;
		
		static TF1 *fTigSigma, *fTigEff, *fTigEffLow, *fTigEffUp;
		static TGraphErrors *gTigEff;			

		
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
	  static Bool_t SetIntensity(Int_t state, Double_t strength);
	  static Bool_t ReadIntensities(const char *fname = "SavedIntensities.txt");
	  static void WriteIntensities(const char *fname = "SavedIntensities.txt");

	  static TH1D *DrawGammas(Int_t from_state, Double_t egam=0.0);
		static TH1D *DrawGammasGated(Double_t exmin=-1.0, Double_t exmax=-1.0, Double_t egam=-1.0);
    static TH2F *DrawExcGam(Double_t egam=-1.0, Bool_t use_int=false);	
	  static TH1D *MakeRealistic(TH1D *hgam, Bool_t abseff = false);
				
		static TH1D *DrawIntensitites(){ return hint; }				
		static Double_t BranchingRatio(Double_t state_eng, Double_t egam);			
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
		        
	ClassDef(TTigressAnalysis,0)
};


#endif
