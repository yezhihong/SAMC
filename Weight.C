/*
  Weight some variables by cross section.
  -- Zhihong Ye, 08/08/2011
*/
//C++ headers
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cmath>
#include "unistd.h"
#include "math.h"
#include "stdio.h"
using namespace std;
//ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TStorage.h"
#include "TChain.h"

void Weight()
{
  const int Bin = 60; 
  gSystem->Exec("ls -l *.root");
  TString filename; 
  cerr <<"Which RootFiles? "; cin >> filename;

  // / gSystem->Exec("cp -vf "+filename+" "+filename+".old");

  TFile *f1 = new TFile(filename,"update");
  assert(f1);

  TTree *t1 =(TTree*)f1->Get("SAMC");  assert(t1);

  /*Define Variables{{{*/
  int Id; //Event Id
  int IsPassed; //if pass through the magnet, 1=true, 0=false(not pass)
  int IsQualified; //if rec var is in the range of gen, 1=true, 0=false(not

  /*only for RadCor{{{*/
  double mp_cor; //Multiple Photon Correction
  double tail; //radiative tail
  double tail_err; //error for radiative tail
  double R; //Phys.Rev.D 12 1884 A83

  double cs_qfs;//quasi-free scattering cross section
  double cs_del;//Delta
  double cs_x;  //DIS
  double cs_r1500; //Resonace 1500 MeV
  double cs_r1700; //Resonace 1700 MeV
  double cs_2N; //Dip region
  double cs_q;//sum of all the cross section above

  double cs_p; //sigma_p (peak approximation for internal radation)
  double cs_b_bt; //the part before target in sigma_b
  double cs_b_at; //the part after target in sigma_b
  double cs_b; //sigma_b (external bremsstrahlung ionization)
  double cs_ex; //sigma_ex (extract internal radation)
  double cs_qt;
  /*}}}*/

  /*Member Data from random{{{*/
  double beam_x; //cm
  double beam_y; //cm
  double reactz_gen; // interaction z position cm
  double E_s;//=incident beam energy (MeV)
  double E_p;//=incident beam energy (MeV)
  //gen=generator
  double th_tg_gen; //theta target
  double ph_tg_gen; //phi target
  double dp_gen; //dp at target
  /*}}}*/

  /*Member Data derived from variables above{{{*/
  double x_tg_gen; //
  double y_tg_gen; //
  /*}}}*/

  /*Member Data for refinement{{{*/
  double dp_ref; //dp for John.LeRose matrix
  double x_tg_ref; //x tg for John.LeRose matrix
  double y_tg_ref; //y tg for John.LeRose matrix
  double th_tg_ref; //th tg for John.LeRose matrix
  double ph_tg_ref; //ph tg for John.LeRose matrix
  /*}}}*/

  /*Member Data for focal plane {{{*/
  double x_fp; //
  double y_fp; //
  double th_fp; //
  double ph_fp; //
  double th_fp_no_ort; //th_fp without TXFIT orthogonalization
  double q1ex_x,q1ex_y;//0:x 1:y
  double dent_x,dent_y;
  double dext_x,dext_y;
  double q3en_x,q3en_y;
  double q3ex_x,q3ex_y;

  int IsPassedQ1Ex;
  int IsPassedDipoleEn;
  int IsPassedDipoleEx;
  int IsPassedQ3En;
  int IsPassedQ3Ex;

  /* bool IsPassedQ1Ex; */
  /* bool IsPassedDipoleEn; */
  /* bool IsPassedDipoleEx; */
  /* bool IsPassedQ3En; */
  /* bool IsPassedQ3Ex; */
  /*}}}*/

  /*Member Data for reconstruct target varibles {{{*/
  double x_tg_rec; //cm
  double y_tg_rec; //cm
  double th_tg_rec; //
  double ph_tg_rec; //
  double dp_rec;
  double Angle_rec;
  double reactz_rec;
  double rvalue;
  double cs_M;//mott cross section
  double cs_Final;//final cross section
  /*}}}*/
  /*}}}*/

  /***********************************************************/
  /**               Set BranchAddress                       **/
  /***********************************************************/
  t1->SetBranchAddress("Id",&Id); //Event Id
  t1->SetBranchAddress("IsPassed",&IsPassed); //if pass through the magnet, 1=true, 0=false(not pass)
  t1->SetBranchAddress("IsQualified",&IsQualified); //if rec var is in the range of gen, 1=true, 0=false(not
  
  /*only for RadCor{{{*/
  t1->SetBranchAddress("mp_cor",&mp_cor   ); //Multiple Photon Correction
  t1->SetBranchAddress("tail",&tail   ); //radiative tail
  t1->SetBranchAddress("R",&R   ); //Phys.Rev.D 12 1884 A83
 
  t1->SetBranchAddress("cs_qfs",&cs_qfs   );//quasi-free scattering cross section
  t1->SetBranchAddress("cs_del",&cs_del   );//Delta
  t1->SetBranchAddress("cs_x",&cs_x   );  //DIS
  t1->SetBranchAddress("cs_r1500",&cs_r1500   ); //Resonace 1500 MeV
  t1->SetBranchAddress("cs_r1700",&cs_r1700   ); //Resonace 1700 MeV
  t1->SetBranchAddress("cs_2N",&cs_2N   ); //Dip region
  t1->SetBranchAddress("cs_q",&cs_q   );//sum of all the cross section above
 
  t1->SetBranchAddress("cs_p",&cs_p   ); //sigma_p (peak approximation for internal radation)
  t1->SetBranchAddress("cs_b_bt",&cs_b_bt   ); //the part before target in sigma_b
  t1->SetBranchAddress("cs_b_at",&cs_b_at   ); //the part after target in sigma_b
  t1->SetBranchAddress("cs_b",&cs_b   ); //sigma_b (external bremsstrahlung ionization)
  t1->SetBranchAddress("cs_ex",&cs_ex   ); //sigma_ex (extract internal radation)
  t1->SetBranchAddress("cs_qt",&cs_qt   );
  /*}}}*/
  
  /*Member Data from random{{{*/
  t1->SetBranchAddress("beam_x",&beam_x    ); //cm
  t1->SetBranchAddress("beam_y",&beam_y    ); //cm
  t1->SetBranchAddress("reactz_gen",&reactz_gen    ); // interaction z position cm
  t1->SetBranchAddress("E_s",&E_s    );//=incident beam energy (MeV)
  t1->SetBranchAddress("E_p",&E_p    );//=incident beam energy (MeV)
  //gen=generator
  t1->SetBranchAddress("th_tg_gen",&th_tg_gen    ); //theta target
  t1->SetBranchAddress("ph_tg_gen",&ph_tg_gen    ); //phi target
  t1->SetBranchAddress("dp_gen",&dp_gen    ); //dp at target
  /*}}}*/
  
  /*Member Data derived from variables above{{{*/
  t1->SetBranchAddress("x_tg_gen",&x_tg_gen    ); //
  t1->SetBranchAddress("y_tg_gen",&y_tg_gen    ); //
  /*}}}*/
  
  /*Member Data for refinement{{{*/
  t1->SetBranchAddress("dp_ref",&dp_ref    ); //dp for John.LeRose matrix
  t1->SetBranchAddress("x_tg_ref",&x_tg_ref    ); //x tg for John.LeRose matrix
  t1->SetBranchAddress("y_tg_ref",&y_tg_ref    ); //y tg for John.LeRose matrix
  t1->SetBranchAddress("th_tg_ref",&th_tg_ref    ); //th tg for John.LeRose matrix
  t1->SetBranchAddress("ph_tg_ref",&ph_tg_ref    ); //ph tg for John.LeRose matrix
  /*}}}*/
  
  /*Member Data for focal plane {{{*/
  t1->SetBranchAddress("x_fp",&x_fp    ); //
  t1->SetBranchAddress("y_fp",&y_fp    ); //
  t1->SetBranchAddress("th_fp",&th_fp    ); //
  t1->SetBranchAddress("ph_fp",&ph_fp    ); //
  t1->SetBranchAddress("th_fp_no_ort",&th_fp_no_ort    ); //th_fp without TXFIT orthogonalization
 
  t1->SetBranchAddress("q1ex_x",&q1ex_x    );//0:x 1:y
  t1->SetBranchAddress("dent_x",&dent_x    );
  t1->SetBranchAddress("dext_x",&dext_x    );
  t1->SetBranchAddress("q3en_x",&q3en_x    );
  t1->SetBranchAddress("q3ex_x",&q3ex_x    );
 
  t1->SetBranchAddress("q1ex_y",&q1ex_y    );//0:x 1:y
  t1->SetBranchAddress("dent_y",&dent_y    );
  t1->SetBranchAddress("dext_y",&dext_y    );
  t1->SetBranchAddress("q3en_y",&q3en_y    );
  t1->SetBranchAddress("q3ex_y",&q3ex_y    );
  
  t1->SetBranchAddress("IsPassedQ1Ex",&IsPassedQ1Ex    );
  t1->SetBranchAddress("IsPassedDipoleEn",&IsPassedDipoleEn    );
  t1->SetBranchAddress("IsPassedDipoleEx",&IsPassedDipoleEx    );
  t1->SetBranchAddress("IsPassedQ3En",&IsPassedQ3En    );
  t1->SetBranchAddress("IsPassedQ3Ex",&IsPassedQ3Ex    );
  /*}}}*/
  
  /*Member Data for reconstruct target varibles {{{*/
  t1->SetBranchAddress("x_tg_rec",&x_tg_rec    ); //cm
  t1->SetBranchAddress("y_tg_rec",&y_tg_rec    ); //cm
  t1->SetBranchAddress("th_tg_rec",&th_tg_rec    ); //
  t1->SetBranchAddress("ph_tg_rec",&ph_tg_rec    ); //
  t1->SetBranchAddress("dp_rec",&dp_rec    );
  t1->SetBranchAddress("reactz_rec",&reactz_rec    );
  t1->SetBranchAddress("rvalue",&rvalue    );

  t1->SetBranchAddress("cs_M",&cs_M    );//mott cross section
  t1->SetBranchAddress("cs_Final",&cs_Final    );//final cross section
  
  double Theta,Xbj;
  t1->SetBranchAddress("E_s",&E_s);
  t1->SetBranchAddress("E_p",&E_p);
  t1->SetBranchAddress("Theta",&Theta);
  t1->SetBranchAddress("Xbj",&Xbj);
/*}}}*/

/***********************************************************/
/*       End of  Setting BranchAddress                     */
/***********************************************************/
  /*Define Histrograms{{{*/
  //TH1F *hxfp = new TH1F("hxfp","Weighted X_fp",Bin,-0.5,0.5);
  TH1F *hxfp = new TH1F("hxfp","",Bin,-0.5,0.5);
  hxfp->SetXTitle("x_{fp}");
  hxfp->GetXaxis()->CenterTitle(1);
  hxfp->GetXaxis()->SetTitleFont(32);
  hxfp->GetXaxis()->SetTitleSize(0.06);

  //TH1F *hyfp = new TH1F("hyfp","Weighted Y_fp",Bin,-0.05,0.05);
  TH1F *hyfp = new TH1F("hyfp","",Bin,-0.05,0.05);
  hyfp->SetXTitle("y_{fp}");
  hyfp->GetXaxis()->CenterTitle(1);
  hyfp->GetXaxis()->SetTitleFont(32);
  hyfp->GetXaxis()->SetTitleSize(0.06);

  //TH1F *hthfp = new TH1F("hthfp","Weighted Theta_fp",Bin,-0.02,0.02);
  TH1F *hthfp = new TH1F("hthfp","",Bin,-0.02,0.02);
  hthfp->SetXTitle("#theta_{fp}");
  hthfp->GetXaxis()->CenterTitle(1);
  hthfp->GetXaxis()->SetTitleFont(32);
  hthfp->GetXaxis()->SetTitleSize(0.06);

  TH1F *hphfp = new TH1F("hphfp","",Bin,-0.05,0.05);
  hphfp->SetXTitle("#phi_{fp}");
  hphfp->GetXaxis()->CenterTitle(1);
  hphfp->GetXaxis()->SetTitleFont(32);
  hphfp->GetXaxis()->SetTitleSize(0.06);


  TH1F *hxtg = new TH1F("hxtg","",Bin,-0.01,0.01);
  hxtg->SetXTitle("x_{tg}");
  hxtg->GetXaxis()->CenterTitle(1);
  hxtg->GetXaxis()->SetTitleFont(32);
  hxtg->GetXaxis()->SetTitleSize(0.06);

 TH1F *hytg;TH1F *hrpz; 
   	 if(filename.Contains("he3")||filename.Contains("h2")||filename.Contains("he4")){
    hytg = new TH1F("hytg","",Bin,-0.07,0.07);
    hrpz = new TH1F("hrpz","",Bin,-0.12,0.12);
  }
  else{
    hytg = new TH1F("hytg","",Bin,-0.01,0.01);
    hrpz = new TH1F("hrpz","",Bin,-0.01,0.01);
  }
  hytg->SetXTitle("y_{tg}");
  hytg->GetXaxis()->CenterTitle(1);
  hytg->GetXaxis()->SetTitleFont(32);
  hytg->GetXaxis()->SetTitleSize(0.06);

  hrpz->SetXTitle("z_{react}");
  hrpz->GetXaxis()->CenterTitle(1);
  hrpz->GetXaxis()->SetTitleFont(32);
  hrpz->GetXaxis()->SetTitleSize(0.06);

  TH1F *hthtg = new TH1F("hthtg","",Bin,-0.05,0.05);
  hthtg->SetXTitle("#theta_{tg}");
  hthtg->GetXaxis()->CenterTitle(1);
  hthtg->GetXaxis()->SetTitleFont(32);
  hthtg->GetXaxis()->SetTitleSize(0.06);

  TH1F *hphtg = new TH1F("hphtg","",Bin,-0.03,0.03);
  hphtg->SetXTitle("#phi_{tg}");
  hphtg->GetXaxis()->CenterTitle(1);
  hphtg->GetXaxis()->SetTitleFont(32);
  hphtg->GetXaxis()->SetTitleSize(0.06);

  TH1F *hdptg = new TH1F("hdptg","",Bin,-0.05,0.05);
  hdptg->SetXTitle("#delta p_{tg}");
  hdptg->GetXaxis()->CenterTitle(1);
  hdptg->GetXaxis()->SetTitleFont(32);
  hdptg->GetXaxis()->SetTitleSize(0.06);

  TH2F *hdpphtg = new TH2F("hdpphtg","",Bin,-0.07,0.07,Bin,-0.06,0.06);
  hdpphtg->SetXTitle("#delta p_{tg}");
  hdpphtg->GetXaxis()->CenterTitle(1);
  hdpphtg->GetXaxis()->SetTitleFont(32);
  hdpphtg->GetXaxis()->SetTitleSize(0.06);
  hdpphtg->SetYTitle("#phi_{tg}");
  hdpphtg->GetYaxis()->CenterTitle(1);
  hdpphtg->GetYaxis()->SetTitleFont(32);
  hdpphtg->GetYaxis()->SetTitleSize(0.06);

  TH2F *hxphtg = new TH2F("hxphtg","",Bin,0.5,8.0,Bin,-0.06,0.06);
  hxphtg->SetXTitle("x_{bj}");
  hxphtg->GetXaxis()->CenterTitle(1);
  hxphtg->GetXaxis()->SetTitleFont(32);
  hxphtg->GetXaxis()->SetTitleSize(0.06);
  hxphtg->SetYTitle("#phi_{tg}");
  hxphtg->GetYaxis()->CenterTitle(1);
  hxphtg->GetYaxis()->SetTitleFont(32);
  hxphtg->GetYaxis()->SetTitleSize(0.06);
 /*}}}*/
  Int_t Evt = t1->GetEntries();

  /*Fill Histrograms{{{*/
  double cs_avg = 0.0;
  for(int i=0;i<Evt;i++){
	  t1->GetEntry(i);
	  if(IsPassed
			  &&abs(y_tg_rec+0.00)<0.080
			  &&abs(th_tg_rec)<0.03
			  &&abs(ph_tg_rec)<0.02
			  &&abs(reactz_rec)<0.075
			  &&abs(dp_rec)<0.03)	
		  cs_avg += cs_Final;
  }
  cs_avg /= Evt;
  
  double weight=1.0;
  for(int j=0;j<Evt;j++){
	  t1->GetEntry(j);
	  //    if(IsPassed&&cs_Final>0.0){
	  if(IsPassed){
		  if(abs(y_tg_rec+0.00)<0.080
				  &&abs(th_tg_rec)<0.03
				  &&abs(ph_tg_rec)<0.02
				  &&abs(reactz_rec)<0.075
				  &&abs(dp_rec)<0.03){	
			  weight = cs_Final/cs_avg;

			  hxfp->Fill(x_fp,weight);
			  hyfp->Fill(y_fp,weight);
			  hthfp->Fill(th_fp,weight);
			  hphfp->Fill(ph_fp,weight);

			  hxtg->Fill(x_tg_rec,weight);
			  hytg->Fill(y_tg_rec,weight);
			  hthtg->Fill(th_tg_rec,weight);
			  hphtg->Fill(ph_tg_rec,weight);
			  hdptg->Fill(dp_rec,weight);
			  hrpz ->Fill(reactz_rec,weight);

			  if(!(j%1000))
				  cerr<<Form(" %d events have been proceeded! ", j)<<endl;

		  }
	  }

	  hdpphtg->Fill(dp_rec, ph_tg_rec, weight);
	  hxphtg->Fill(Xbj, ph_tg_rec, weight);
  }
  
  hxfp->Write(); hyfp->Write(); hthfp->Write(); hphfp->Write();
  hxtg->Write(); hytg->Write(); hthtg->Write(); hphtg->Write(); hdptg->Write();
  hrpz->Write(); hdpphtg->Write(); hxphtg->Write();

  /*}}}*/
  f1->Close();
}

