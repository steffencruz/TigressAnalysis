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

Int_t MakeRedwoodMats(std::string outfilename = "Results_RedwoodMats.root", Bool_t all_types=false){

	TChain *chain = new TChain("EasyTree");
  chain->Add("$DATADIR/redwood_Sr95*");		
  printf("\n* Loaded:	%i redwood trees *\n",chain->GetNtrees());
  
	TReaction *r[3];	
	r[0] = new TReaction("sr95","d","p","sr96",510.9,0,true);  	
	printf("\n**  Made Reaction: %s\n",r[0]->GetNameFull());
	if(all_types){
    r[1] = new TReaction("sr95","p","p","sr95",510.9,0,true); 
    printf("\n**  Made Reaction: %s\n",r[1]->GetNameFull());
    r[2] = new TReaction("sr95","d","d","sr95",510.9,0,true); 
    printf("\n**  Made Reaction: %s\n\n",r[2]->GetNameFull());
	}

  
	TSharcAnalysis::SetTarget(0.0,0.0,0.0,4.5,"cd2",0.5,0.5,0.5);
	const char *fname = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";	
	TList *acclist = TSharcAnalysis::GetAcceptanceList(r[0],fname);
	
	TSharcAnalysis::Print();

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
	TList *list[3];	 // list for each reaction, list[0]="dp", list[1]="pp", list[2]="dd"		
	
  TH3S *hexcgamgam[3], *hexcgamthetacm[3], *hexcgamthetalab[3], *hexcgamthetatig[3];//, *hexcgamthetacm2[3], *hexcgamthetalab2[3];
  TH2F *hgamgam[3], *hexcgam[3], *hexcthetacm[3], *hgamthetatig[3], *hthetatheta[3], *hphasetadd[3], *htaddtadd[3];
	TH1D *hgam[3], *hexc[3], *hexcg[3];
	
	
	// for different sharc sections
	TList *list_sec[3][3];	 // list for each reaction and SHARC section combination		
	
	TH3S *hexcgamgam_sec[3][3], *hexcgamthetatig_sec[3][3];	
	TH2F *hexcgam_sec[3][3], *hgamgam_sec[3][3], *hgamthetatig_sec[3][3], *hphasetadd_sec[3][3];
	TH1D *hexc_sec[3][3], *hgam_sec[3][3];
	
  TH2F *hkin = new TH2F("KinVsThetaLab","Kinematics For ^{95}Sr; #theta_{LAB} [#circ]; Energy [keV]",720,0,180,3000,0,30000);
  TH2F *hdengcm = new TH2F("DelVsThetaCm","ElasticKinematicsCm; #theta_{CM} [#circ]; Delta Energy [keV]",720,0,180,1000,0,10000);    
    
	printf("\n\n Making empty histograms .. %4.1f%%\r",0.0);	
  fflush(stdout);  
  const char *opt, *sec;
  
  int imax = 3;
  if(!all_types) // limit this to dp 
    imax=1;
  
  for(int i=0; i<imax; i++){
  	
  	list[i] = new TList;
  
  	if(i==0)
  		opt = "dp";
  	if(i==1)
  		opt = "pp";	
  	if(i==2)
  		opt = "dd";
  	  	
  	// useful for gamma gamma analysis	-- binning is large and gamma energy range is small
		hexcgamgam[i] = new TH3S(Form("ExcGamGam_%s",opt),"",1000,0,4000,1000,0,4000,350,0,7000); 
		hexcgamgam[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Gamma Energy '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]; E_{exc} [keV]",opt));
		list[i]->Add(hexcgamgam[i]);				

  	hexcgamthetatig[i] = new TH3S(Form("ExcGamTigressTheta_%s",opt),"",180,0,180,2000,0,4000,800,-1000,7000);
  	hexcgamthetatig[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. TIGRESS Theta For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",opt));
  	list[i]->Add(hexcgamthetatig[i]);  
  	
  	hgamthetatig[i] = new TH2F(Form("GamTigressTheta_%s",opt),"",180,0,180,2000,0,4000);
  	hgamthetatig[i]->SetTitle(Form("Gamma Energy vs. TIGRESS Theta For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]",opt));
  	list[i]->Add(hgamthetatig[i]);    	 

		hexcgam[i] = new TH2F(Form("ExcGam_%s",opt),"",2000,0,4000,700,0,7000);
		hexcgam[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy '%s'; E_{#gamma} [keV]; E_{exc} [keV]",opt));	  		
  	list[i]->Add(hexcgam[i]);	
		
		hgamgam[i] = new TH2F(Form("GamGam_%s",opt),"",2000,0,4000,2000,0,4000);
		hgamgam[i]->SetTitle(Form("Gamma Energy vs. Gamma Energy '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]",opt));	  		
  	list[i]->Add(hgamgam[i]);
  	
		hexc[i] = new TH1D(Form("Exc_%s",opt),"",800,-1000,7000);
		hexc[i]->SetTitle(Form("Excitation Energy '%s'; E_{exc} [keV];",opt));	  		
  	list[i]->Add(hexc[i]);
  	  	
		hgam[i] = new TH1D(Form("Gam_%s",opt),"",4000,0,4000);
		hgam[i]->SetTitle(Form("Gamma Energy '%s'; E_{#gamma} [keV];",opt));	  		
  	list[i]->Add(hgam[i]);	  	
	
  	// sectional matrices for DBOX, UBOX & UQQQ [dp only]  	
		for(int j=0; j<3; j++){
		  if(i>0)
		    break;
		
      list_sec[i][j] = new TList;		
		
			if(j==0)
				sec = "UQ";
			if(j==1)
				sec = "UB";	
			if(j==2)
				sec = "DB";		
				
			printf(" Making empty histograms .. %4.1f%%\r",(double)(imax*i+j)/(3*imax)*100.0);	
			fflush(stdout);
			
			hexcgamgam_sec[i][j] = new TH3S(Form("ExcGamGam%s_%s",sec,opt),"",1000,0,4000,1000,0,4000,350,0,7000); 
			hexcgamgam_sec[i][j]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Gamma Energy [%s] '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]; E_{exc} [keV]",sec,opt));
			list_sec[i][j]->Add(hexcgamgam_sec[i][j]);	
			
      hexcgamthetatig_sec[i][j] = new TH3S(Form("ExcGamTigressTheta%s_%s",sec,opt),"",180,0,180,2000,0,4000,800,-1000,7000);
      hexcgamthetatig_sec[i][j]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. TIGRESS Theta [%s] For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",sec,opt));
      list_sec[i][j]->Add(hexcgamthetatig_sec[i][j]); 				

      hgamthetatig_sec[i][j] = new TH2F(Form("GamTigressTheta%s_%s",sec,opt),"",180,0,180,2000,0,4000);
      hgamthetatig_sec[i][j]->SetTitle(Form("Gamma Energy vs. TIGRESS Theta For [%s] '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]",sec,opt));
      list_sec[i][j]->Add(hgamthetatig_sec[i][j]);    
      			
			hexcgam_sec[i][j] = new TH2F(Form("ExcGam%s_%s",sec,opt),"",2000,0,4000,700,0,7000);
			hexcgam_sec[i][j]->SetTitle(Form("Excitation Energy vs. Gamma Energy [%s] '%s'; E_{#gamma} [keV]; E_{exc} [keV]",sec,opt));	  		
			list_sec[i][j]->Add(hexcgam_sec[i][j]);	
		
			hgamgam_sec[i][j] = new TH2F(Form("GamGam%s_%s",sec,opt),"",2000,0,4000,2000,0,4000);
			hgamgam_sec[i][j]->SetTitle(Form("Gamma Energy vs. Gamma Energy [%s] '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]",sec,opt));	  		
			list_sec[i][j]->Add(hgamgam_sec[i][j]);
		
			hexc_sec[i][j] = new TH1D(Form("Exc%s_%s",sec,opt),"",800,-1000,7000);
			hexc_sec[i][j]->SetTitle(Form("Excitation Energy [%s] '%s'; E_{exc} [keV];",sec,opt));	  		
			list_sec[i][j]->Add(hexc_sec[i][j]);
				
			hgam_sec[i][j] = new TH1D(Form("Gam%s_%s",sec,opt),"",4000,0,4000);
			hgam_sec[i][j]->SetTitle(Form("Gamma Energy [%s] '%s'; E_{#gamma} [keV];",sec,opt));	  		
			list_sec[i][j]->Add(hgam_sec[i][j]);	  	 		
			
      hphasetadd_sec[i][j] = new TH2F(Form("PhaseTAdd%s_%s",sec,opt),"",400,0,400,700,0,7);
      hphasetadd_sec[i][j]->SetTitle(Form("Phase vs. TAdd [%s] For '%s'; TIGRESS time [ns]; RF Phase [ ~ SHARC time]",sec,opt));
      list_sec[i][j]->Add(hphasetadd_sec[i][j]); 				
  	}			
  	
		hexcg[i] = new TH1D(Form("Exc_AllGam_%s",opt),"",800,-1000,7000);
		hexcg[i]->SetTitle(Form("Excitation Energy '%s' [with all gammas]; E_{exc} [keV];",opt));	  		
  	list[i]->Add(hexcg[i]);	   		 
  	  	
  	// pixel center theta plots		
  	hexcgamthetacm[i] = new TH3S(Form("ExcGamThetaCm_%s",opt),"",180,0,180,2000,0,4000,800,-1000,7000);
  	hexcgamthetacm[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Theta Cm For '%s'; #theta_{CM} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",opt));
  	list[i]->Add(hexcgamthetacm[i]);

		hexcthetacm[i] = new TH2F(Form("ExcThetaCm_%s",opt),"",180,0,180,800,-1000,7000);
		hexcthetacm[i]->SetTitle(Form("Excitation Energy vs. Theta Cm '%s'; #theta_{CM} [#circ]; E_{exc} [keV]",opt));	  		
  	list[i]->Add(hexcthetacm[i]);	  	
  	
  	hexcgamthetalab[i] = new TH3S(Form("ExcGamThetaLab_%s",opt),"",180,0,180,2000,0,4000,800,-1000,7000);
  	hexcgamthetalab[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Theta Lab For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",opt));
  	list[i]->Add(hexcgamthetalab[i]);	 	
  			
		// secondary hists to check everything is working
  	hthetatheta[i] = new TH2F(Form("ThetaLabThetaCm_%s",opt),"",180,0,180,180,0,180);
  	hthetatheta[i]->SetTitle(Form("Theta Lab vs. Theta Cm For '%s'; #theta_{CM} [#circ]; #theta_{LAB} [#circ]",opt));
		list[i]->Add(hthetatheta[i]);

		hphasetadd[i] = new TH2F(Form("PhaseTAdd_%s",opt),"",400,0,400,700,0,7);
		hphasetadd[i]->SetTitle(Form("Phase vs. TAdd For '%s'; TIGRESS time [ns]; RF Phase [ ~ SHARC time]",opt));
		list[i]->Add(hphasetadd[i]); 
	
		htaddtadd[i] = new TH2F(Form("TAddTAdd_%s",opt),"",400,0,400,400,0,400);
		htaddtadd[i]->SetTitle(Form("TAdd vs. TAdd For '%s'; TIGRESS time [ns]; TIGRESS time [ns]",opt));
		list[i]->Add(htaddtadd[i]);
		
		fflush(stdout);
		
	}
	printf(" Making empty histograms .. %4.1f%%     -- COMPLETE \n\n",100.0);	
	
	printf("\n List of Histograms :-");
	TObject *obj;
	for(int i=0; i<imax; i++){

  	if(i==0)
      printf("\n -> REACTION = ' DP '");
  	if(i==1)
      printf("\n -> REACTION = ' PP '");
  	if(i==2)
      printf("\n -> REACTION = ' DD '");
	
	  TIter iter(list[i]);
	  while((obj = iter.Next()))
	    printf("\n\t%s*    %s",obj->ClassName(),obj->GetName());

		for(int j=0; j<3; j++){
		  if(i>0)
		    break;

			if(j==0)
        printf("\n\t -> SECTION = ' UQ '");
			if(j==1)
        printf("\n\t -> SECTION = ' UB '");
			if(j==2)
        printf("\n\t -> SECTION = ' DB '");
				      
      TIter iter2(list_sec[i][j]);
      while((obj = iter2.Next()))
        printf("\n\t\t%s*    %s",obj->ClassName(),obj->GetName());	
    }    
  }
			
	Double_t tmin = 230.0, tmax = 265.0; // large time window
	Double_t thmin = 0;//135.0; // option to use only largest angle data : disabled

  printf("\n\n Using the following control parameters:-");
  printf("\n TIGTIME  :   %.2f <   tigtime [ns] < %.2f",tmin,tmax);
  printf("\n TIGTHETA :   tigtheta [deg] > %.2f\n",thmin);  

	printf("\n\n Filling histograms:-\n\n");
	Int_t nentries=chain->GetEntries(), sec_indx;	

	for(int n=0; n<nentries; n++){
	
		chain->GetEntry(n);
		
		if(n%250000==0){
			printf("\t processing event %8i/%i\t[%5.2f%%]\r",n,nentries,(Double_t)n/nentries*100.0);
			fflush(stdout);
		}
		
		if(!all_types && type>0)
		  continue;			
		
		// fill hist excluding dead strips and other options
		if(TSharcAnalysis::BadStrip(det,front,-1) || TSharcAnalysis::BadStrip(det,-1,back))
			continue;							
		
		hkin->Fill(thetalab,ekin);
	
		if(type==0 && exc>7000)//|| thetalab<60)) // gets rid of inelastic 95Sr stuff in (d,p)
			continue;		
			
		if(det<9) // detector section identifier
			sec_indx = 2;	
		else if(det<13)
			sec_indx = 1;			
		else if(det>=13)
			sec_indx = 0;
						
		hexc[type]->Fill(exc);	
		// sectional excitaion spectra [only dp]
		if(type==0) hexc_sec[0][sec_indx]->Fill(exc);	
				
		// fill the gamma energy = 0 bin with everything [including events not coincident with TIGRESS]				
		hexcthetacm[type]->Fill(thetacm,exc);	

		hexcgamthetalab[type]->Fill(thetalab,0.0,exc);
		hexcgamthetacm[type]->Fill(thetacm,0.0,exc);	

		hdengcm->Fill(thetacm,denergy);
		
		hthetatheta[type]->Fill(thetacm,thetalab);

		if(addsize>4) // only use first 3 gammas
			addsize=3;
			
		// we require a gamma energy for these hists
		for(int i=0; i<(int)addsize; i++){
			for(int j=i+1; j<(int)addsize; j++){
				htaddtadd[type]->Fill(tadd[i],tadd[j]);
				
				if(tadd[i]<tmin || tadd[i]>tmax || eadd[i]<50.0) 
					continue;								
								
				if(tadd[j]<tmin || tadd[j]>tmax || eadd[j]<50.0) 
					continue;				
					
				if(thetadd[i]<thmin || thetadd[j]<thmin) 
				  continue;						
					
				hexcgamgam[type]->Fill(eadd[i],eadd[j],exc);
				hexcgamgam[type]->Fill(eadd[j],eadd[i],exc);
													
				hgamgam[type]->Fill(eadd[i],eadd[j]);
				hgamgam[type]->Fill(eadd[j],eadd[i]);
				
				// fill in detector section plots
				if(type==0) hexcgamgam_sec[0][sec_indx]->Fill(eadd[i],eadd[j],exc);				
				if(type==0) hexcgamgam_sec[0][sec_indx]->Fill(eadd[j],eadd[i],exc);		
				
				if(type==0) hgamgam_sec[0][sec_indx]->Fill(eadd[i],eadd[j]);				
				if(type==0) hgamgam_sec[0][sec_indx]->Fill(eadd[j],eadd[i]);				

			}
			hphasetadd[type]->Fill(tadd[i],phase);
			if(type==0) hphasetadd_sec[0][sec_indx]->Fill(tadd[i],phase);
			
			if(tadd[i]<tmin || tadd[i]>tmax || eadd[i]<50.0) 
				continue;							
				
      if(thetadd[i]<thmin) 
        continue;					

			hexcgam[type]->Fill(eadd[i],exc);		
			hexcg[type]->Fill(exc);									
			hgam[type]->Fill(eadd[i]);
			
			if(type==0) hexcgam_sec[0][sec_indx]->Fill(eadd[i],exc);		
			if(type==0) hgam_sec[0][sec_indx]->Fill(eadd[i]);			
			
			hgamthetatig[type]->Fill(thetadd[i],eadd[i]);
			if(type==0) hgamthetatig_sec[0][sec_indx]->Fill(thetadd[i],eadd[i]);
			
			hexcgamthetatig[type]->Fill(thetadd[i],eadd[i],exc);
			if(type==0) hexcgamthetatig_sec[0][sec_indx]->Fill(thetadd[i],eadd[i],exc);			
			
			hexcgamthetalab[type]->Fill(thetalab,eadd[i],exc);
			hexcgamthetacm[type]->Fill(thetacm,eadd[i],exc);		
		}
		
	}
	
	hexcgamthetalab[0]->GetYaxis()->SetRangeUser(0,3);
	TH2F *h2 = (TH2F*)hexcgamthetalab[0]->Project3D("zx");
	TH1D *h = (TH1D*)h2->ProjectionX("px",30,70);
	list[0]->Add(h2);
	list[0]->Add(h);
	

  printf("\n\t Complete!\n\n\n Now writing results to file:-\n");
	fflush(stdout);
	
	TFile *file;
	file = new TFile(outfilename.c_str(),"RECREATE");	
	
	file->SetCompressionLevel(4);
	
	file->mkdir("dp");
	file->cd("dp");	
	list[0]->Write();
	printf("\n -> Wrote list to ' dp '");
	fflush(stdout);
	
  file->cd("/");		
	file->mkdir("dp/UQ");
	file->cd("dp/UQ");		
	list_sec[0][0]->Write();
	printf("\n -> Wrote list to ' dp/UQ '");
	fflush(stdout);	
	
  file->cd("/");		
	file->mkdir("dp/UB");
	file->cd("dp/UB");	
	list_sec[0][1]->Write();	
	printf("\n -> Wrote list to ' dp/UB '");
	fflush(stdout);	
	
  file->cd("/");		
	file->mkdir("dp/DB");
	file->cd("dp/DB");		
	list_sec[0][2]->Write();	
	printf("\n -> Wrote list to ' dp/DB '");
	fflush(stdout);	
	
	if(all_types){	
		file->cd("/");	
		file->mkdir("pp");
		file->cd("pp");	
		list[1]->Write();	
    printf("\n -> Wrote list to ' pp '");
	
		file->cd("/");	
		file->mkdir("dd");
		file->cd("dd");		
		list[2]->Write();	
		printf("\n -> Wrote list to ' dd '");
	
	}	

	file->cd("/");		
	file->mkdir("Acceptances");
	file->cd("Acceptances");
	acclist->Write();
	printf("\n -> Wrote acceptances");
	
	file->cd("/");			
	hkin->Write();
	hdengcm->Write();
	printf("\n -> Wrote auxiliary objects");
	
	printf("\n\n Made File ' %s '.. Closing\n\n",file->GetName());
	
	file->Close();
	
  return 0;		
}

