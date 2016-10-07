//g++ -O2 -std=c++0x MakeExcGamThetaMats.C -L/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis -lSharcAnalysis -I/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis  `grsi-config --cflags --all-libs --root` -o egtmats.out
#include<stdio.h>

#include<TMath.h>
#include<Rtypes.h>

#include<TFile.h>
#include<TList.h>
#include<TChain.h>
#include<TH2F.h>
#include<TH3S.h>
#include<TGraph.h>
#include<TCutG.h>

#include"TSharcAnalysis.h"
#include"TReaction.h"

Int_t MakeExcGamThetaMats(Bool_t all_types=false, Bool_t veasy=false){

	TChain *chain = new TChain("EasyTree");
	
	//chain->Add("$DATADIR/sapling*");
	
	if(veasy)
		chain->Add("$worktime/Jul22_2015/Veasy/veasy*");
	else 
		chain->Add("$DATADIR/redwood*");
		
	chain->ls();	
			
	TReaction *r[3];	
	r[0] = new TReaction("sr95","d","p","sr96",510.9,0,true);  	
	printf("\n**  Made Reaction: %s",r[0]->GetNameFull());
	r[1] = new TReaction("sr95","p","p","sr95",510.9,0,true); 
	printf("\n**  Made Reaction: %s",r[1]->GetNameFull());
	r[2] = new TReaction("sr95","d","d","sr95",510.9,0,true); 
	printf("\n**  Made Reaction: %s\n",r[2]->GetNameFull());

	TSharcAnalysis::SetTarget(0.0,0.0,0.0,4.5,"cd2",0.5,0.5,0.5);
	const char *fname = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";
	TList *acclist = TSharcAnalysis::GetAcceptanceList(r[0],fname);

	// set tree branch addresses
	UInt_t det=0,front=0,back=0,type=0,addsize=0;
	Double_t thetalab=0,thetacm=0,phicorr=0,denergy=0,ekin=0,exc=0,phase=0;
	Double_t tadd[10], eadd[10], thetadd[10];	
	
	chain->SetBranchAddress("det",&det);
	chain->SetBranchAddress("front",&front);
	chain->SetBranchAddress("back",&back);
	
	chain->SetBranchAddress("type",&type);	
	chain->SetBranchAddress("phase",&phase);	
				
	chain->SetBranchAddress("thetalab",&thetalab);
	chain->SetBranchAddress("thetacm",&thetacm);

	chain->SetBranchAddress("denergy",&denergy);
	chain->SetBranchAddress("ekin",&ekin);
	chain->SetBranchAddress("exc",&exc);

	chain->SetBranchAddress("addsize",&addsize);		
	chain->SetBranchAddress("eadd",&eadd);
	chain->SetBranchAddress("tadd",&tadd);
	chain->SetBranchAddress("thetadd",&thetadd);
	
	// for different reactions
  TH3S *hexcgamgam[3], *hexcgamthetacm[3], *hexcgamthetalab[3], *hexcgamthetatig[3];//, *hexcgamthetacm2[3], *hexcgamthetalab2[3];
  TH2F *hgamgam[3], *hexcgam[3], *hexcthetacm[3], *hgamthetatig[3], *hthetatheta[3], *hphasetadd[3], *htaddtadd[3];
	TH1D *hgam[3], *hexc[3], *hexcg[3];
	
	// for different sharc sections
	TH3S *hexcgamgam_sec[3][3];	
	TH2F *hexcgam_sec[3][3], *hgamgam_sec[3][3];
	TH1D *hexc_sec[3][3], *hgam_sec[3][3];
	
	TList *list[3];
	TGraph *glab2cm[3];
  
  TH2F *hkin = new TH2F("KinVsThetaLab","Kinematics For ^{95}Sr; Theta Lab [deg]; Energy [keV]",720,0,180,3000,0,30000);
  TH2F *hdengcm = new TH2F("DelVsThetaCm","ElasticKinematicsCm; Theta Cm [deg]; Delta Energy [keV]",720,0,180,1000,0,10000);    
    
	printf("\n\n Making empty histograms .. %4.1f%%\r",0);	
  fflush(stdout);  
  const char *opt, *sec;
  for(int i=0; i<3; i++){
  	
  	list[i] = new TList;
  
  	if(i==0)
  		opt = "dp";
  	if(i==1)
  		opt = "pp";	
  	if(i==2)
  		opt = "dd";
  		
  	// useful for gamma gamma analysis	-- binning is large and gamma energy range is small
		hexcgamgam[i] = new TH3S(Form("ExcGamGam_%s",opt),"",1000,0,4000,1000,0,4000,350,0,7000); 
		hexcgamgam[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Gamma Energy '%s'; Gamma energy [keV]; Gamma energy [keV]; Excitation energy [keV]",opt));
		list[i]->Add(hexcgamgam[i]);				

		hexcgam[i] = new TH2F(Form("ExcGam_%s",opt),"",2000,0,4000,700,0,7000);
		hexcgam[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy '%s'; Gamma energy [keV]; Excitation energy [keV]",opt));	  		
  	list[i]->Add(hexcgam[i]);	
		
		hgamgam[i] = new TH2F(Form("GamGam_%s",opt),"",2000,0,4000,2000,0,4000);
		hgamgam[i]->SetTitle(Form("Gamma Energy vs. Gamma Energy '%s'; Gamma energy [keV]; Gamma energy [keV]",opt));	  		
  	list[i]->Add(hgamgam[i]);
  	
		hexc[i] = new TH1D(Form("Exc_%s",opt),"",800,-1000,7000);
		hexc[i]->SetTitle(Form("Excitation Energy '%s'; Excitation energy [keV];",opt));	  		
  	list[i]->Add(hexc[i]);
  	  	
		hgam[i] = new TH1D(Form("Gam_%s",opt),"",4000,0,4000);
		hgam[i]->SetTitle(Form("Gamma Energy '%s'; Gamma energy [keV];",opt));	  		
  	list[i]->Add(hgam[i]);	  	
	
  	// sectional matrices for DBOX, UBOX & UQQQ [dp only]  	
		for(int j=0; j<3; j++){
		
			if(j==0)
				sec = "UQ";
			if(j==1)
				sec = "UB";	
			if(j==2)
				sec = "DB";		
				
			printf(" Making empty histograms .. %4.1f%%\r",(double)(3*i+j)/9.0*100.0);	
			fflush(stdout);
			
			hexcgamgam_sec[i][j] = new TH3S(Form("ExcGamGam%s_%s",sec,opt),"",1000,0,4000,1000,0,4000,350,0,7000); 
			hexcgamgam_sec[i][j]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Gamma Energy [%s] '%s'; Gamma energy [keV]; Gamma energy [keV]; Excitation energy [keV]",sec,opt));
			list[i]->Add(hexcgamgam_sec[i][j]);		
			
			hexcgam_sec[i][j] = new TH2F(Form("ExcGam%s_%s",sec,opt),"",2000,0,4000,700,0,7000);
			hexcgam_sec[i][j]->SetTitle(Form("Excitation Energy vs. Gamma Energy [%s] '%s'; Gamma energy [keV]; Excitation energy [keV]",sec,opt));	  		
			list[i]->Add(hexcgam_sec[i][j]);	
		
			hgamgam_sec[i][j] = new TH2F(Form("GamGam%s_%s",sec,opt),"",2000,0,4000,2000,0,4000);
			hgamgam_sec[i][j]->SetTitle(Form("Gamma Energy vs. Gamma Energy [%s] '%s'; Gamma energy [keV]; Gamma energy [keV]",sec,opt));	  		
			list[i]->Add(hgamgam_sec[i][j]);
		
			hexc_sec[i][j] = new TH1D(Form("Exc%s_%s",sec,opt),"",800,-1000,7000);
			hexc_sec[i][j]->SetTitle(Form("Excitation Energy [%s] '%s'; Excitation energy [keV];",sec,opt));	  		
			list[i]->Add(hexc_sec[i][j]);
				
			hgam_sec[i][j] = new TH1D(Form("Gam%s_%s",sec,opt),"",4000,0,4000);
			hgam_sec[i][j]->SetTitle(Form("Gamma Energy [%s] '%s'; Gamma energy [keV];",sec,opt));	  		
			list[i]->Add(hgam_sec[i][j]);	  	
			

  	}			
  	
		hexcg[i] = new TH1D(Form("Exc_AllGam_%s",opt),"",800,-1000,7000);
		hexcg[i]->SetTitle(Form("Excitation Energy '%s' [with all gammas]; Excitation energy [keV];",opt));	  		
  	list[i]->Add(hexcg[i]);	   		 
  	  	
  	// pixel center theta plots		
  	hexcgamthetacm[i] = new TH3S(Form("ExcGamThetaCm_%s",opt),"",180,0,180,2000,0,4000,800,-1000,7000);
  	hexcgamthetacm[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Theta Cm For '%s'; Theta Cm [deg]; Gamma Energy [keV]; Excitation Energy [keV]",opt));
  	list[i]->Add(hexcgamthetacm[i]);

		hexcthetacm[i] = new TH2F(Form("ExcThetaCm_%s",opt),"",180,0,180,800,-1000,7000);
		hexcthetacm[i]->SetTitle(Form("Excitation Energy vs. Theta Cm '%s'; Theta Cm [deg]; Excitation energy [keV]",opt));	  		
  	list[i]->Add(hexcthetacm[i]);	  	
  	
  	hexcgamthetalab[i] = new TH3S(Form("ExcGamThetaLab_%s",opt),"",180,0,180,2000,0,4000,800,-1000,7000);
  	hexcgamthetalab[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Theta Lab For '%s'; Theta Lab [deg]; Gamma Energy [keV]; Excitation Energy [keV]",opt));
  	list[i]->Add(hexcgamthetalab[i]);	

  	hexcgamthetatig[i] = new TH3S(Form("ExcGamTigressTheta_%s",opt),"",180,0,180,2000,0,4000,800,-1000,7000);
  	hexcgamthetatig[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. TIGRESS Theta For '%s'; Theta Lab [deg]; Gamma Energy [keV]; Excitation Energy [keV]",opt));
  	list[i]->Add(hexcgamthetatig[i]);  
  	
  	hgamthetatig[i] = new TH2F(Form("GamTigressTheta_%s",opt),"",180,0,180,2000,0,4000);
  	hgamthetatig[i]->SetTitle(Form("Gamma Energy vs. TIGRESS Theta For '%s'; Theta Lab [deg]; Gamma Energy [keV]",opt));
  	list[i]->Add(hgamthetatig[i]);    	  	
  			
		// secondary hists to check everything is working
  	hthetatheta[i] = new TH2F(Form("ThetaLabThetaCm_%s",opt),"",180,0,180,180,0,180);
  	hthetatheta[i]->SetTitle(Form("Theta Lab vs. Theta Cm For '%s'; Theta Cm [deg]; Theta Lab [deg]",opt));
		list[i]->Add(hthetatheta[i]);

		hphasetadd[i] = new TH2F(Form("PhaseTAdd_%s",opt),"",400,0,400,700,0,7);
		hphasetadd[i]->SetTitle(Form("Phase vs. TAdd For '%s'; TIGRESS time [ns]; RF Phase [ ~ SHARC time]",opt));
		list[i]->Add(hphasetadd[i]); 
	
		htaddtadd[i] = new TH2F(Form("TAddTAdd_%s",opt),"",400,0,400,400,0,400);
		htaddtadd[i]->SetTitle(Form("TAdd vs. TAdd For '%s'; TIGRESS time [ns]; TIGRESS time [ns]",opt));
		list[i]->Add(htaddtadd[i]);
		
		glab2cm[i] = r[i]->ThetaVsTheta();
		list[i]->Add(glab2cm[i]);
	}
	
	printf(" Making empty histograms .. %4.1f%%      COMPLETE \n\n",100.0);	
			
/*		
	TFile cutfile("$worktime/Jul22_2015/Veasy/timecuts.root","READ");
	// TIMECUT1 = first band [	way better!! ]	
	// TIMECUT2 = first and second band	combined
	TCutG *timecut = (TCutG*) cutfile.Get("TIMECUT1");
	//list->Add(timecut);
*/
	Double_t tmin = 230.0, tmax = 265.0; // large time window

	// these should be suitable for all trees
	TFile cutfile("$DATADIR/TimeCuts.root","READ");
	TCutG* time_dbox = (TCutG*) cutfile.Get("TIME_DBOX");
	TCutG* time_ubox = (TCutG*) cutfile.Get("TIME_UBOX");
	TCutG* time_ubox2 = (TCutG*) cutfile.Get("TIME_UBOX2");
	TCutG* time_uqqq = (TCutG*) cutfile.Get("TIME_UQQQ");
	TCutG* timecut;
	
	TH1D *hgam1 = new TH1D("gam_notime","gam_notime",4000,0,4000);
	TH1D *hgam2 = new TH1D("gam_time","gam_time",4000,0,4000);	
	
	printf("\n\n Filling histograms:-\n\n");
	Int_t nentries=chain->GetEntries(), sec_indx;
	for(int n=0; n<nentries; n++){
	
		chain->GetEntry(n);
		
		if(n%250000==0){
			printf("\t processing event %8i/%i\t[%5.2f%%]\r",n,nentries,(Double_t)n/nentries*100.0);
			fflush(stdout);
		}
		
		// fill hist excluding dead strips and other options
		if(TSharcAnalysis::BadStrip(det,front,-1) || TSharcAnalysis::BadStrip(det,-1,back))
			continue;							
		
		if(veasy){
			if(det==5 && n<23e6)
				continue;		
		// for veasy trees
			thetalab = TSharcAnalysis::RandomizeThetaLab(det,front,back);
			thetacm = glab2cm[type]->Eval(thetalab);
		}	
		
		hkin->Fill(thetalab,ekin);
	
		if(type==0 && exc>7000)//|| thetalab<60)) // gets rid of inelastic 95Sr stuff in (d,p)
			continue;		
			
		if(det<9)
			sec_indx = 2;	
		else if(det<13)
			sec_indx = 1;			
		else if(det>=13)
			sec_indx = 0;	
						
		hexc[type]->Fill(exc);	
		// sectional excitaion spectra
		hexc_sec[type][sec_indx]->Fill(exc);	
				
		// fill the gamma energy = 0 bin with everything [including events not coincident with TIGRESS]				
		hexcthetacm[type]->Fill(thetacm,exc);	

		hexcgamthetalab[type]->Fill(thetalab,0.0,exc);
		hexcgamthetacm[type]->Fill(thetacm,0.0,exc);	

		hdengcm->Fill(thetacm,denergy);
		
		hthetatheta[type]->Fill(thetacm,thetalab);
		/*
		if(det<=8)
			timecut = time_dbox;
		else if(det>8 && det<=12)
			timecut = time_ubox2;
		else if(det>12 && det<=16)	
			timecut = time_uqqq;	
		*/
		if(addsize>4) // only use first 3 gammas
			addsize=3;
			
		// we require a gamma energy for these hists
		for(int i=0; i<(int)addsize; i++){
			for(int j=i+1; j<(int)addsize; j++){
				htaddtadd[type]->Fill(tadd[i],tadd[j]);
				
//				if(!timecut->IsInside(tadd[i],phase) || eadd[i]<50.0) 
//					continue;

				if(tadd[i]<tmin || tadd[i]>tmax || eadd[i]<50.0) 
					continue;								
								
//				if(!timecut->IsInside(tadd[j],phase) || eadd[j]<50.0) 
//					continue;			

				if(tadd[j]<tmin || tadd[j]>tmax || eadd[j]<50.0) 
					continue;									
					
				hexcgamgam[type]->Fill(eadd[i],eadd[j],exc);
				hexcgamgam[type]->Fill(eadd[j],eadd[i],exc);
													
				hgamgam[type]->Fill(eadd[i],eadd[j]);
				hgamgam[type]->Fill(eadd[j],eadd[i]);
				
				// fill in detector section plots
				hexcgamgam_sec[type][sec_indx]->Fill(eadd[i],eadd[j],exc);				
				hexcgamgam_sec[type][sec_indx]->Fill(eadd[j],eadd[i],exc);		
				
				hgamgam_sec[type][sec_indx]->Fill(eadd[i],eadd[j]);				
				hgamgam_sec[type][sec_indx]->Fill(eadd[j],eadd[i]);				

			}
			hphasetadd[type]->Fill(tadd[i],phase);
			
			if(exc>800 && exc<1600 && type==0)
				hgam1->Fill(eadd[i]);
				
//			if(!timecut->IsInside(tadd[i],phase) || eadd[i]<50.0) // gamma must be > 50 keV 
//				continue;

			if(tadd[i]<tmin || tadd[i]>tmax || eadd[i]<50.0) 
				continue;							
			
			if(exc>800 && exc<1600 && type==0)
				hgam2->Fill(eadd[i]);

			hexcgam[type]->Fill(eadd[i],exc);		
			hexcg[type]->Fill(exc);									
			hgam[type]->Fill(eadd[i]);
			
			hexcgam_sec[type][sec_indx]->Fill(eadd[i],exc);		
			hgam_sec[type][sec_indx]->Fill(eadd[i]);			
			
			hgamthetatig[type]->Fill(thetadd[i],eadd[i]);

			hexcgamthetalab[type]->Fill(thetalab,eadd[i],exc);
			hexcgamthetacm[type]->Fill(thetacm,eadd[i],exc);		
			hexcgamthetatig[type]->Fill(thetadd[i],eadd[i],exc);
		}
		
	}
	
	hexcgamthetalab[0]->GetYaxis()->SetRangeUser(0,3);
	TH2F *h2 = (TH2F*)hexcgamthetalab[0]->Project3D("zx");
	TH1D *h = (TH1D*)h2->ProjectionX("px",30,70);
	list[0]->Add(h2);
	list[0]->Add(h);
	

  printf("\n\t Complete!\n\n Now writing results to file:-\n");
	fflush(stdout);
	
	TFile *file;
	file = new TFile("Results_ExcGamThetaMats_RedwoodDets.root","RECREATE");	
	/*
	if(veasy) 
		file = new TFile("Results_ExcGamThetaMats_Veasy.root","RECREATE");
	else
		file = new TFile("Results_ExcGamThetaMats_Redwood.root","RECREATE");		
	*/
	
	file->SetCompressionLevel(4);
	
	file->mkdir("dp");
	file->cd("dp");	
	list[0]->Write();
	
	if(all_types){	
		file->cd("/");	
		file->mkdir("pp");
		file->cd("pp");	
		list[1]->Write();	
	
		file->cd("/");	
		file->mkdir("dd");
		file->cd("dd");		
		list[2]->Write();	
	}	

	file->cd("/");		
	file->mkdir("Acceptances");
	file->cd("Acceptances");
	acclist->Write();

	file->cd("/");			
	hkin->Write();
	hdengcm->Write();
	
	hgam1->Write();
	hgam2->Write();	

	printf("\n\t Made File ' %s '\n\n",file->GetName());
	
	file->Close();
	
  return 0;		
}
