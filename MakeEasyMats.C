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
#include"TTigress.h"

/*
  ____MakeEasyMats____
  Prepares multi-dimensional histograms for TTigressAnalysis
  - A reaction can be selected using eg "dp" or "dp+pp+dt" or "dp dt"
  - If one reaction is selected, it is possible to make sectional (DB UB UQ) histograms
  - If more than one reaction is selected, sectional histograms cannot be made
*/

Double_t CalculateDoppler(Double_t eng, UInt_t tigdet, UInt_t tigcry, UInt_t tigseg, TVector3 pos_offs, Double_t rbeta){

  TVector3 segpos = TTigress::GetPosition(tigdet,tigcry,tigseg);
  Double_t theta = (segpos+pos_offs).Theta();
  Double_t doppler = eng/sqrt(1-rbeta*rbeta)*(1-rbeta*cos(theta));

  return doppler;
}

Bool_t MakeEasyMats(UInt_t A, std::string treename = "$DATADIR/redwood_Sr95_*", std::string outfilename = "", std::string reac="dp+dt", Bool_t low_res=false){
  
  Double_t beame; // beam enegry in middle of target
  Double_t targ;//target thickness
  std::string badstripsname;
 
  // I tested these on all three data sets and they work fine
  Double_t tmin=220.0, tmax=265.0;//  time window 
  Double_t minchgdiff=-500.0, maxchgdiff=200.0; // dcharge-dchargeb
  Double_t excmin=-2000, excmax=8000; // excitation energy range in histograms
  Double_t thmin=0.0; // option to use only largest angle data : disabled 
  Double_t pos[] = {0.0,0.0,0.0};
  TVector3 tig_offs(0.0,0.0,0.0);

  if(A == 94){ // beam mass number
    beame = 499.5; 
    targ = 5.0; 
    badstripsname = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips_sr94.txt";	
  } else if (A == 95){
    beame = 510.9; 
    targ = 4.5; 
    badstripsname = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";	   
 //   tig_offs.SetXYZ(0.0,0.0,-0.4);
  } else if (A == 96){
    beame = 513.5; 
    targ = 5.0; 
    badstripsname = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";	   
  //  tig_offs.SetXYZ(0.0,0.0,-0.4);
  }
  
	TChain *chain = new TChain("EasyTree");
  chain->Add(treename.c_str());		
  printf("\n* Loaded:	%i trees *\n",chain->GetNtrees());
	
  static const unsigned long npos = std::string::npos;	
  std::vector<std::string> opt(5);
  UInt_t nreac=0;
	TReaction *r[5];	
	if(reac.find("dp")!=npos){
    r[0] = new TReaction(Form("sr%i",A),"d","p",Form("sr%i",A+1),beame,0,true);  	
    printf("\n**  Made Reaction: %s\n",r[0]->GetNameFull());	
    opt.at(0).assign("dp"); nreac++;
  }
  if(reac.find("pp")!=npos){
    r[1] = new TReaction(Form("sr%i",A),"p","p",Form("sr%i",A),beame,0,true); 
    printf("\n**  Made Reaction: %s\n",r[1]->GetNameFull());
    opt.at(1).assign("pp"); nreac++;
  }
  if(reac.find("dd")!=npos){
    r[2] = new TReaction(Form("sr%i",A),"d","d",Form("sr%i",A),beame,0,true); 
    printf("\n**  Made Reaction: %s\n",r[2]->GetNameFull());
    opt.at(2).assign("dd"); nreac++;
  }
  if(reac.find("cc")!=npos){
    r[3] = new TReaction(Form("sr%i",A),"c12","c12",Form("sr%i",A),beame,0,true); 
    printf("\n**  Made Reaction: %s\n\n",r[3]->GetNameFull());
    opt.at(3).assign("cc"); nreac++;
  }
  if(reac.find("dt")!=npos){
    r[4] = new TReaction(Form("sr%i",A),"d","t",Form("sr%i",A-1),beame,0,true); 
    printf("\n**  Made Reaction: %s\n\n",r[4]->GetNameFull());
    opt.at(4).assign("dt"); nreac++;
  }  
  if(!nreac){
    printf("\n\n Error :  No reactions have been selected. Try ""dp+dt""\n\n");
    return false;
  }
  
	TSharcAnalysis::SetTarget(pos[0],pos[1],pos[2],targ,"cd2",0.5,0.5,0.5);
	TList *acclist = TSharcAnalysis::GetAcceptanceList(r[0],badstripsname.c_str());
	
	TSharcAnalysis::Print();

	// set tree branch addresses
	UInt_t det=0,front=0,back=0,type=0,addsize=0;
	Double_t thetalab=0,thetacm=0,phicorr=0,denergy=0,dcharge=0,dchargeb=0,penergy=0;
	Double_t ekin=0,exc=0,phase=0,rbeta=0;
	Double_t tadd[10], eadd[10], thetadd[10], tigeng[10];
  UShort_t tigdet[10], tigcry[10], tigseg[10];
  	
	chain->SetBranchAddress("det",&det);
	chain->SetBranchAddress("front",&front);
	chain->SetBranchAddress("back",&back);
	
	chain->SetBranchAddress("type",&type);	
	chain->SetBranchAddress("phase",&phase);	
				
	chain->SetBranchAddress("thetalab",&thetalab);
	chain->SetBranchAddress("thetacm",&thetacm);

	chain->SetBranchAddress("denergy",&denergy);
	chain->SetBranchAddress("dcharge",&dcharge);
	chain->SetBranchAddress("dchargeb",&dchargeb);
	chain->SetBranchAddress("penergy",&penergy);	

	chain->SetBranchAddress("ekin",&ekin);
	chain->SetBranchAddress("exc",&exc);

	chain->SetBranchAddress("rbeta",&rbeta);

	chain->SetBranchAddress("addsize",&addsize);		
	chain->SetBranchAddress("eadd",&eadd);
	chain->SetBranchAddress("tadd",&tadd);
	chain->SetBranchAddress("thetadd",&thetadd);
	chain->SetBranchAddress("tigeng",&tigeng);
	chain->SetBranchAddress("tigdet",&tigdet);
	chain->SetBranchAddress("tigcry",&tigcry);
	chain->SetBranchAddress("tigseg",&tigseg);

	// for different reactions
	TList *list[5];	 // list for each reaction, list[0]="dp", list[1]="pp", list[2]="dd", list[3]="cc", list[4]="dt"		
	
  TH3S *hexcgamgam[5], *hexcgamthetacm[5], *hexcgamthetalab[5];
  TH3S *hexcgamthetatig[5], *hexcthetatheta[5], *hphasetaddeadd[5];//, *hexcgamthetacm2[3], *hexcgamthetalab2[3];
  TH2F *hgamgam[5], *hexcgam[5], *hexcthetacm[5];
  TH2F *hgamthetatig[5], *hphasetadd[5], *htaddtadd[5];
	TH1D *hgam[5], *hexc[5], *hexcg[5];
	
	Bool_t no_sec=false;	
	if(nreac>1)
	  low_res=true;
	if(low_res){
	  printf("\n\n Sectional histograms [dbox,ubox,uqqq] will not be made.\n\n");
	  no_sec=true; // don't make sectional matrices	  
	}
	// for different sharc sections
	TList *list_sec[3];	 // list for each reaction and SHARC section combination		
	TH3S *hexcgamgam_sec[3], *hexcgamthetatig_sec[3];	
	TH2F *hexcgam_sec[3], *hgamgam_sec[3], *hgamthetatig_sec[3], *hphasetadd_sec[3];
	TH1D *hexc_sec[3], *hgam_sec[3];
	
  TH2F *hkin = new TH2F("KinVsThetaLab",Form("Kinematics For ^{%i}Sr; #theta_{LAB} [#circ]; Energy [keV]",A),720,0,180,3000,0,30000);
  TH2F *hdengcm = new TH2F("DelVsThetaCm","ElasticKinematicsCm; #theta_{CM} [#circ]; Delta Energy [keV]",720,0,180,1000,0,10000);    
    
	printf("\n\n Making empty histograms .. %4.1f%%\r",0.0);	
  fflush(stdout);  
  
  for(int i=0; i<5; i++){
  	
  	if(opt.at(i).length()<2)
  	  continue;
  	list[i] = new TList;	
  	const char *rname = opt.at(i).c_str();
  	
  	// useful for gamma gamma analysis	-- binning is large and gamma energy range is small
		hexcgamgam[i] = new TH3S(Form("ExcGamGam_%s",rname),"",1000,0,4000,1000,0,4000,500,excmin,excmax); 
		hexcgamgam[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Gamma Energy '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]; E_{exc} [keV]",rname));
		list[i]->Add(hexcgamgam[i]);				

  	hexcgamthetatig[i] = new TH3S(Form("ExcGamTigressTheta_%s",rname),"",180,0,180,2000,0,4000,1000,excmin,excmax);
  	hexcgamthetatig[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. TIGRESS Theta For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",rname));
  	list[i]->Add(hexcgamthetatig[i]);    	

  	hexcgamthetacm[i] = new TH3S(Form("ExcGamThetaCm_%s",rname),"",180,0,180,2000,0,4000,1000,excmin,excmax);
  	hexcgamthetacm[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Theta Cm For '%s'; #theta_{CM} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",rname));
  	list[i]->Add(hexcgamthetacm[i]);  	
  	
  	hgamthetatig[i] = new TH2F(Form("GamTigressTheta_%s",rname),"",180,0,180,4000,0,4000);
  	hgamthetatig[i]->SetTitle(Form("Gamma Energy vs. TIGRESS Theta For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]",rname));
  	list[i]->Add(hgamthetatig[i]);     	

		hexcthetacm[i] = new TH2F(Form("ExcThetaCm_%s",rname),"",180,0,180,1000,excmin,excmax);
		hexcthetacm[i]->SetTitle(Form("Excitation Energy vs. Theta Cm '%s'; #theta_{CM} [#circ]; E_{exc} [keV]",rname));	  		
  	list[i]->Add(hexcthetacm[i]);	  	

		hexcgam[i] = new TH2F(Form("ExcGam_%s",rname),"",4000,0,4000,1000,excmin,excmax);
		hexcgam[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy '%s'; E_{#gamma} [keV]; E_{exc} [keV]",rname));	  		
  	list[i]->Add(hexcgam[i]);	
		
		hgamgam[i] = new TH2F(Form("GamGam_%s",rname),"",4000,0,4000,4000,0,4000);
		hgamgam[i]->SetTitle(Form("Gamma Energy vs. Gamma Energy '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]",rname));	  		
  	list[i]->Add(hgamgam[i]);
  	
		hexc[i] = new TH1D(Form("Exc_%s",rname),"",1000,excmin,excmax);
		hexc[i]->SetTitle(Form("Excitation Energy '%s'; E_{exc} [keV];",rname));	  		
  	list[i]->Add(hexc[i]);
  	  	
		hgam[i] = new TH1D(Form("Gam_%s",rname),"",4000,0,4000);
		hgam[i]->SetTitle(Form("Gamma Energy '%s'; E_{#gamma} [keV];",rname));	  		
  	list[i]->Add(hgam[i]);	  	
	
  	// sectional matrices for DBOX, UBOX & UQQQ [dp only]  	
		if(nreac==1){
		  for(int j=0; j<3; j++){
		
      list_sec[j] = new TList;		
		  const char *sec;

			if(j==0)
				sec = "UQ";
			if(j==1)
				sec = "UB";	
			if(j==2)
				sec = "DB";		
				
			printf(" Making empty histograms .. %4.1f%%\r",(double)j/3.0*100.0);	
			fflush(stdout);
			
			hexcgamgam_sec[j] = new TH3S(Form("ExcGamGam%s_%s",sec,rname),"",1000,0,4000,1000,0,4000,500,excmin,excmax); 
			hexcgamgam_sec[j]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Gamma Energy [%s] '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]; E_{exc} [keV]",sec,rname));
			list_sec[j]->Add(hexcgamgam_sec[j]);	
			
      hexcgamthetatig_sec[j] = new TH3S(Form("ExcGamTigressTheta%s_%s",sec,rname),"",180,0,180,2000,0,4000,1000,excmin,excmax);
      hexcgamthetatig_sec[j]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. TIGRESS Theta [%s] For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",sec,rname));
      list_sec[j]->Add(hexcgamthetatig_sec[j]); 				

      hgamthetatig_sec[j] = new TH2F(Form("GamTigressTheta%s_%s",sec,rname),"",180,0,180,4000,0,4000);
      hgamthetatig_sec[j]->SetTitle(Form("Gamma Energy vs. TIGRESS Theta For [%s] '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]",sec,rname));
      list_sec[j]->Add(hgamthetatig_sec[j]);    
      			
			hexcgam_sec[j] = new TH2F(Form("ExcGam%s_%s",sec,rname),"",4000,0,4000,1000,excmin,excmax);
			hexcgam_sec[j]->SetTitle(Form("Excitation Energy vs. Gamma Energy [%s] '%s'; E_{#gamma} [keV]; E_{exc} [keV]",sec,rname));	  		
			list_sec[j]->Add(hexcgam_sec[j]);	
		
			hgamgam_sec[j] = new TH2F(Form("GamGam%s_%s",sec,rname),"",4000,0,4000,4000,0,4000);
			hgamgam_sec[j]->SetTitle(Form("Gamma Energy vs. Gamma Energy [%s] '%s'; E_{#gamma} [keV]; E_{#gamma} [keV]",sec,rname));	  		
			list_sec[j]->Add(hgamgam_sec[j]);
		
			hexc_sec[j] = new TH1D(Form("Exc%s_%s",sec,rname),"",1000,excmin,excmax);
			hexc_sec[j]->SetTitle(Form("Excitation Energy [%s] '%s'; E_{exc} [keV];",sec,rname));	  		
			list_sec[j]->Add(hexc_sec[j]);
				
			hgam_sec[j] = new TH1D(Form("Gam%s_%s",sec,rname),"",4000,0,4000);
			hgam_sec[j]->SetTitle(Form("Gamma Energy [%s] '%s'; E_{#gamma} [keV];",sec,rname));	  		
			list_sec[j]->Add(hgam_sec[j]);	  	 		
			
      hphasetadd_sec[j] = new TH2F(Form("PhaseTAdd%s_%s",sec,rname),"",400,0,400,700,0,7);
      hphasetadd_sec[j]->SetTitle(Form("Phase vs. TAdd [%s] For '%s'; TIGRESS time [ns]; RF Phase [ ~ SHARC time]",sec,rname));
      list_sec[j]->Add(hphasetadd_sec[j]); 		
      }		
  	}			
    
    if(low_res) // no need to build secondary histograms
      continue;
      
		hexcg[i] = new TH1D(Form("Exc_AllGam_%s",rname),"",1000,excmin,excmax);
		hexcg[i]->SetTitle(Form("Excitation Energy '%s' [with all gammas]; E_{exc} [keV];",rname));	  		
  	list[i]->Add(hexcg[i]);	   		   	
  	
  	hexcgamthetalab[i] = new TH3S(Form("ExcGamThetaLab_%s",rname),"",180,0,180,2000,0,4000,1000,excmin,excmax);
  	hexcgamthetalab[i]->SetTitle(Form("Excitation Energy vs. Gamma Energy vs. Theta Lab For '%s'; #theta_{LAB} [#circ]; E_{#gamma} [keV]; E_{exc} [keV]",rname));
  	list[i]->Add(hexcgamthetalab[i]);	 	
  			
		// secondary hists to check everything is working
  	hexcthetatheta[i] = new TH3S(Form("ExcThetaLabThetaCm_%s",rname),"",180,0,180,180,0,180,1000,excmin,excmax);
  	hexcthetatheta[i]->SetTitle(Form("Theta Lab vs. Theta Cm vs. Excitation energy For '%s'; #theta_{CM} [#circ]; #theta_{LAB} [#circ]; E_{exc} [keV]",rname));
		list[i]->Add(hexcthetatheta[i]);  

		hphasetadd[i] = new TH2F(Form("PhaseTAdd_%s",rname),"",400,0,400,700,0,7);
		hphasetadd[i]->SetTitle(Form("Phase vs. TAdd For '%s'; TIGRESS time [ns]; RF Phase [ ~ SHARC time]",rname));
		list[i]->Add(hphasetadd[i]); 
		
		hphasetaddeadd[i] = new TH3S(Form("PhaseTAddEadd_%s",rname),"",400,0,400,700,0,7,1000,0,4000);
		hphasetaddeadd[i]->SetTitle(Form("Phase vs. TAdd For '%s'; TIGRESS time [ns]; RF Phase [ ~ SHARC time]; E_{#gamma} [keV]",rname));
		list[i]->Add(hphasetaddeadd[i]); 
	
		htaddtadd[i] = new TH2F(Form("TAddTAdd_%s",rname),"",400,0,400,400,0,400);
		htaddtadd[i]->SetTitle(Form("TAdd vs. TAdd For '%s'; TIGRESS time [ns]; TIGRESS time [ns]",rname));
		list[i]->Add(htaddtadd[i]);
		
		fflush(stdout);		
	}
	printf(" Making empty histograms .. %4.1f%%     -- COMPLETE \n\n",100.0);	
	
	printf("\n List of Histograms :-");
	TObject *obj;
	for(int i=0; i<5; i++){

  	if(opt.at(i).length()<2)
  	  continue;
    printf("\n -> REACTION = ' %s '",opt.at(i).c_str());
	  TIter iter(list[i]);
	  while((obj = iter.Next()))
	    printf("\n\t%s*    %s",obj->ClassName(),obj->GetName());

    if(nreac==1){
      for(int j=0; j<3; j++){

        if(j==0)
          printf("\n\t -> SECTION = ' UQ '");
        if(j==1)
          printf("\n\t -> SECTION = ' UB '");
        if(j==2)
          printf("\n\t -> SECTION = ' DB '");
              
        TIter iter2(list_sec[j]);
        while((obj = iter2.Next()))
          printf("\n\t\t%s*    %s",obj->ClassName(),obj->GetName());	
      }    
    }
  }
			
  printf("\n\n Using the following control parameters:-");
  printf("\n TIGTIME  :   %.2f <   tigtime [ns] < %.2f",tmin,tmax);
  printf("\n CHGDIFF  :   %.2f <   fchg - bchg  < %.2f",minchgdiff,maxchgdiff);
  printf("\n TIGTHETA :   tigtheta [deg] > %.2f\n",thmin);  

	printf("\n\n Filling histograms:-\n\n");
	Int_t nentries=chain->GetEntries(), sec_indx;	
	Double_t excr;
  double d2r = 0.017453;
  double r2d = 57.2958;			
  UInt_t ncount[5];

	for(int n=0; n<nentries; n++){
	
		chain->GetEntry(n);
		
		if(n%250000==0){
			printf("\t processing event %8i/%i\t[%5.2f%%]\r",n,nentries,(Double_t)n/nentries*100.0);
			fflush(stdout);
		}
		
  	if(opt.at(type).length()<2)
  	  continue;		// if 'type' was not included in specified reactions, continue
		
		// fill hist excluding dead strips and other options
		if(TSharcAnalysis::BadStrip(det,front,-1) || TSharcAnalysis::BadStrip(det,-1,back))
			continue;							
		
		if(dcharge-dchargeb>maxchgdiff || dcharge-dchargeb<minchgdiff)
  		continue;
		
		hkin->Fill(thetalab,ekin);
		if(type==1 || type==2 || type==3){ // elastics must assume no excitation energy
        r[type]->SetExc(0.0);		
		    thetacm = r[type]->ConvertThetaLabToCm(thetalab*d2r,2)*r2d;			
		    hdengcm->Fill(thetacm,denergy);
		//    printf("\n %i type == %i thetalab = %.1f   thetacm = %.1f  180-2*thetalab = %.1f",n,type,thetalab,thetacm,180.0-2*thetalab);		    
    //    fflush(stdout);
    } else if(type==4 && penergy<100)
      continue;
      
//     if(!all_types && type>0) // only fill dp
//       continue;		
		  		
		// (p,p) can be reconstructed as (d,p) here too. 
	  //  if(type==1) // && penergy>50.0) // only with PID?
	  //    type=0;
    // reconstruct excitation energy using fixed TReaction
    // excr = r[type]->GetExcEnergy(ekin*1e-3,thetalab*d2r,2)*1e3;
    // r[type]->SetExc(excr*1e3);
    // exc = excr;    
	
		if(exc<excmin || exc>excmax)
			continue;		
			
    ncount[type]++;  
			
    // now set excitation energy to get angular conversions and stuff right
    r[type]->SetExc(exc*1e-3);

    // set thetacm using the excitation energy			
    thetacm = r[type]->ConvertThetaLabToCm(thetalab*d2r,2)*r2d;			
			
		if(det<9) // detector section identifier
			sec_indx = 2;	
		else if(det<13)
			sec_indx = 1;			
		else if(det>=13)
			sec_indx = 0;
						
		hexc[type]->Fill(exc);	
		// sectional excitation spectra [only one reaction]
		if(!no_sec) hexc_sec[sec_indx]->Fill(exc);	
				
		// fill the gamma energy = 0 bin with everything [including events not coincident with TIGRESS]				
		hexcthetacm[type]->Fill(thetacm,exc);	
		hexcgamthetacm[type]->Fill(thetacm,0.0,exc);	
		
    if(!low_res){
      hexcgamthetalab[type]->Fill(thetalab,0.0,exc);
      hexcthetatheta[type]->Fill(thetacm,thetalab,exc);
    }
//		if(addsize>3) // only use first 3 gammas
//			addsize=3;

		// we require a gamma energy for these hists
		for(int i=0; i<(int)addsize; i++){
      
      if(tig_offs.Mag())
        eadd[i] = CalculateDoppler(tigeng[i],tigdet[i],tigcry[i],tigseg[i],tig_offs,rbeta);
		
			for(int j=i+1; j<(int)addsize; j++){
				
				if(!low_res) htaddtadd[type]->Fill(tadd[i],tadd[j]);
				
				if(tadd[i]<tmin || tadd[i]>tmax || eadd[i]<50.0) 
					continue;								
								
				if(tadd[j]<tmin || tadd[j]>tmax || eadd[j]<50.0) 
					continue;				
					
				if(thetadd[i]<thmin || thetadd[j]<thmin) 
				  continue;						
				
        if(tig_offs.Mag())
          eadd[j] = CalculateDoppler(tigeng[j],tigdet[j],tigcry[j],tigseg[j],tig_offs,rbeta);
				
        hexcgamgam[type]->Fill(eadd[i],eadd[j],exc);
        hexcgamgam[type]->Fill(eadd[j],eadd[i],exc);
        
        hgamgam[type]->Fill(eadd[i],eadd[j]);
        hgamgam[type]->Fill(eadd[j],eadd[i]);
      
				// fill in detector section plots if required
				if(!no_sec){
				  hexcgamgam_sec[sec_indx]->Fill(eadd[i],eadd[j],exc);				
				  hexcgamgam_sec[sec_indx]->Fill(eadd[j],eadd[i],exc);		
				
				  hgamgam_sec[sec_indx]->Fill(eadd[i],eadd[j]);				
				  hgamgam_sec[sec_indx]->Fill(eadd[j],eadd[i]);				
        }
			}
			
			if(!low_res){
			  hphasetaddeadd[type]->Fill(tadd[i],phase,eadd[i]);
			  hphasetadd[type]->Fill(tadd[i],phase);
			}
			
			if(!no_sec) hphasetadd_sec[sec_indx]->Fill(tadd[i],phase);
			
			if(tadd[i]<tmin || tadd[i]>tmax || eadd[i]<50.0) 
				continue;							
				
      if(thetadd[i]<thmin) 
        continue;					

			hexcgamthetacm[type]->Fill(thetacm,eadd[i],exc);					

			hexcgam[type]->Fill(eadd[i],exc);		
			hgam[type]->Fill(eadd[i]);
			
      hgamthetatig[type]->Fill(thetadd[i],eadd[i]);
			hexcgamthetatig[type]->Fill(thetadd[i],eadd[i],exc);
			
			if(!low_res){
        hexcg[type]->Fill(exc);			
			  hexcgamthetalab[type]->Fill(thetalab,eadd[i],exc);			
      }			
			if(!no_sec){
        hexcgam_sec[sec_indx]->Fill(eadd[i],exc);		
        hgam_sec[sec_indx]->Fill(eadd[i]);			
        hgamthetatig_sec[sec_indx]->Fill(thetadd[i],eadd[i]);
        hexcgamthetatig_sec[sec_indx]->Fill(thetadd[i],eadd[i],exc);			
			}
		}
		
	}
	
	if(!low_res){
    hexcgamthetalab[0]->GetYaxis()->SetRangeUser(0,3);
    TH2F *h2 = (TH2F*)hexcgamthetalab[0]->Project3D("zx");
    TH1D *h = (TH1D*)h2->ProjectionX("px",30,70);
    list[0]->Add(h2);
    list[0]->Add(h);
	}

  printf("\n\t Complete!\n\n\n Now writing results to file:-\n");
	fflush(stdout);
	
	TFile *file;
	if(outfilename.size()==0)
	  outfilename.assign(Form("Results_%s.root",treename.c_str()));
	
	file = new TFile(outfilename.c_str(),"RECREATE");	
	file->SetCompressionLevel(1);
	
  for(int i=0; i<5; i++){
    
    printf("\n Number of ' %s ' events : %i",opt.at(i).c_str(),ncount[i]);
  
    if(opt.at(i).length()<2)
      continue;
		file->cd("/");	
    file->mkdir(opt.at(i).c_str());
    file->cd(opt.at(i).c_str());	
    list[i]->Write();
    printf("\n -> Wrote list to ' %s '",opt.at(i).c_str());
    fflush(stdout);
  }
	
	if(!no_sec){
    file->cd("/");		
    file->mkdir(Form("%s/UQ",opt.at(0).c_str()));
    file->cd(Form("%s/UQ",opt.at(0).c_str()));
    list_sec[0]->Write();
    printf("\n -> Wrote list to ' %s/UQ '",opt.at(0).c_str());
    fflush(stdout);	
  
    file->cd("/");		
    file->mkdir(Form("%s/UB",opt.at(0).c_str()));
    file->cd(Form("%s/UB",opt.at(0).c_str()));
    list_sec[1]->Write();
    printf("\n -> Wrote list to ' %s/UB '",opt.at(0).c_str());
    fflush(stdout);	
    
    file->cd("/");		
    file->mkdir(Form("%s/DB",opt.at(0).c_str()));
    file->cd(Form("%s/DB",opt.at(0).c_str()));
    list_sec[2]->Write();
    printf("\n -> Wrote list to ' %s/DB '",opt.at(0).c_str());
    fflush(stdout);	
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
	
  return true;		
}
