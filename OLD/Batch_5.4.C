#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#include "CommonParameters.h"
#include "ADC_TARGET_Pedestals.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr2.h"
#include "Pedestals.h"
#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"
#include "Plot_Event_Display.C"
#include "Event_Display_functions.cxx"
#include "Batch_Variables_5.4.h"
#endif

using namespace std;

//double SFT_Test_CR_Event(int Run_Number, int evt, double phi, bool to_print, double C2X_centroid, double length_in_target);
//vector<double> SFT_Test_CR_Event2(int Run_Number, int evt, double phi, double C2X_centroid, double length_in_target);
int kaon_fiber(double x_bar, double y_bar);

//void Batch_5_4(Int_t Run_Number=5, Int_t ievt=0, Int_t enable_cout=0){
void Batch_5_4(int Run_Number=5, int ievt_min=0, int ievt_max=20, int Switch_Display=1){

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);

  //bool flag = false;

  //if(Switch_Display==0 || Switch_Display==1 || Switch_Display==2){
  //  flag = true;
  //}

  //if(!flag){
  //  cout << endl;
  //  cout << "Flag Error! " << endl;
  //  cout << "Please, choose another flag" << endl;
  //  cout << "0: Screen Printout (Reduced Statistics)" << endl;
  //  cout << "1: File Printout (Reduced Statistics) " << endl;
  //  cout << "2: File Printout (Full Statistics)" << endl;
  //  cout << endl;
  //}

  /////////////////////   Dave's Time Walk Correction File  ////////////////////////
  double par_in[256][3] = {0.};
  double par_err[356][3] = {0.};
  Int_t ADCmax = 3450;
  double Yfi = 0.;
  double Ani = 0.;
  double Yni = 0.;
  double Tpi = 0.;
  //bool NoFile = false;

  char ParsTarg1[100];
  sprintf(ParsTarg1,"TimeWalk%d.dat",Run_Number);

  //if(!ifstream(ParsTarg1)) NoFile = true;

  /////////////////////////////////////////////////////////////////////////

  sprintf(path_input,"%s",path_merged);
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  char output[100];
  sprintf(output,"RUN_%d_DATA.txt",Run_Number);
  ofstream fout;
  fout.open(output);

  TChain *fChain=new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("TDC_Trig",tdc_trigger); // tdc_trigger for TDC_diff calculation

  fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
  fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
  fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D); 
  
  fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);    
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

  fChain->SetBranchAddress("MWPCADC",MwpcADC);

  fChain->SetBranchAddress("ADC_C2X_R",adc_c2x_r);
  fChain->SetBranchAddress("ADC_C2X_L",adc_c2x_l);
  fChain->SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
  fChain->SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
  fChain->SetBranchAddress("ADC_C3X_R",adc_c3x_r);
  fChain->SetBranchAddress("ADC_C3X_L",adc_c3x_l);
  fChain->SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
  fChain->SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
  fChain->SetBranchAddress("ADC_C4X_R",adc_c4x_r);
  fChain->SetBranchAddress("ADC_C4X_L",adc_c4x_l);
  fChain->SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
  fChain->SetBranchAddress("ADC_C4Y_L",adc_c4y_l);  

  fChain->SetBranchAddress("TDC_Ck", tdc_ck);
  fChain->SetBranchAddress("TDC_Cpi", tdc_cpi);

  fChain->SetBranchAddress("VT48_TDC",tdc_vt48);

  fChain->SetBranchAddress("EvFlag", Event_flag);
  fChain->SetBranchAddress("EvTag", Event_tag);
 
  int nentries;
  nentries = fChain->GetEntries();

  if(Switch_Display==2){
  	ievt_min=0;
  	ievt_max=nentries;
  }

  for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_HG[i]) + TARGET_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;


  ///////////////////////////////////////////////////////////////////////////////////////
  /// LOOP OVER EVENTS

	for(int ievt=ievt_min; ievt<ievt_max; ievt++){
		if(ievt%10000==0 && (Switch_Display==0 || Switch_Display==2)){
			cout << "***  " << ievt << " events done!" << endl;
		}

  		fChain->GetEntry(ievt);

    	for(int j_TARGET=0; j_TARGET<256; j_TARGET++){
    		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
    		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
    		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
    	}

    	////////////////////////////////////////////////////////////////////////
    	/// Time walk (by Dave)

      	ifstream parTARGdat(ParsTarg1,ios::in);
      	Int_t ij = 0;
      	Int_t ik = 0;

      	// Read in parameters and their errors. (errors not used)
      	for(Int_t ii = 0; ii < nBars; ii++){
        	parTARGdat >> ij >> ik >> par_in[ii][0] >> par_err[ii][0];
        	parTARGdat >> ij >> ik >> par_in[ii][1] >> par_err[ii][1];
        	parTARGdat >> ij >> ik >> par_in[ii][2] >> par_err[ii][2];
      	}

      	for(Int_t ii = 0; ii<256; ii++){

        	Yfi = par_in[ii][0] - par_in[ii][1]/sqrt(ADCmax - par_in[ii][2]);
        	Ani = adc_high_target[ii]-TARGET_ADC_Ped_HG[ii]; //SEB

        	if((Ani >= TARGET_ADC_Thr_HG_Offset) && (Ani < ADCmax)){ // SEB
          	Yni = par_in[ii][0] - par_in[ii][1]/sqrt(Ani - par_in[ii][2]);
          	Tpi = Yfi - Yni;
          	for(int jj=0; jj<16; jj++){
            	if(tdc_le_target[ii][jj]>0) tdc_le_target[ii][jj] = tdc_le_target[ii][jj] + Tpi;
          	}
        	}
      	}
      	///////////////////////////////////////////////////////////////////////////

    	///// SFT
    	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
      		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
     		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
    		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];

    		if(ADC_High_SFT[j_SFT]<0)	ADC_High_SFT_corr[j_SFT]=0;
    		if(ADC_High_SFT[j_SFT]>=0)	ADC_High_SFT_corr[j_SFT]=ADC_High_SFT[j_SFT];
    	}


    	///// TOF1 ADCs 
    	for(int i=0; i<12; i++){
      	ADC_TOF1[i] = ADC_tof1U[i]-TOF1U_ADC_Thr[i];
      	ADC_TOF1[i+12] = ADC_tof1D[i]-TOF1D_ADC_Thr[i];
    	}

    	///// TOF1 & TOF2 TDCs
    	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
    		TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
    		TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      		TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      		TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      		TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      		TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
    	}

    	///// C2X
    	for(Int_t j_C2=0; j_C2<56; j_C2++){
      		ADC_C2X_R[j_C2] = double(adc_c2x_r[j_C2])-(ADC_C2X_Thr_R[j_C2]+MWPC_Thr_Offset_C2X);
      		ADC_C2X_L[j_C2] = double(adc_c2x_l[j_C2])-(ADC_C2X_Thr_L[j_C2]+MWPC_Thr_Offset_C2X);
    	}

    	///// C3X
    	for(Int_t j_C3=0; j_C3<64; j_C3++){
    		ADC_C3X_R[j_C3] = double(adc_c3x_r[j_C3])-(ADC_C3X_Thr_R[j_C3]+MWPC_Thr_Offset_C3X);
      		ADC_C3X_L[j_C3] = double(adc_c3x_l[j_C3])-(ADC_C3X_Thr_L[j_C3]+MWPC_Thr_Offset_C3X);
    	}

    	///// C4X
    	for(Int_t j_C4=0; j_C4<72; j_C4++){
    		ADC_C4X_R[j_C4] = double(adc_c4x_r[j_C4])-(ADC_C4X_Thr_R[j_C4]+MWPC_Thr_Offset_C4X);
      		ADC_C4X_L[j_C4] = double(adc_c4x_l[j_C4])-(ADC_C4X_Thr_L[j_C4]+MWPC_Thr_Offset_C4X);
    	}

    	///// C2Y, C3Y & C4Y
    	for(Int_t j_CY=0; j_CY<16; j_CY++){
    		ADC_C2Y_R[j_CY] = double(adc_c2y_r[j_CY])-(ADC_C2Y_Thr_R[j_CY]+MWPC_Thr_Offset_C2Y);
      		ADC_C2Y_L[j_CY] = double(adc_c2y_l[j_CY])-(ADC_C2Y_Thr_L[j_CY]+MWPC_Thr_Offset_C2Y);
      		ADC_C3Y_R[j_CY] = double(adc_c3y_r[j_CY])-(ADC_C3Y_Thr_R[j_CY]+MWPC_Thr_Offset_C3Y);
      		ADC_C3Y_L[j_CY] = double(adc_c3y_l[j_CY])-(ADC_C3Y_Thr_L[j_CY]+MWPC_Thr_Offset_C3Y);
      		ADC_C4Y_R[j_CY] = double(adc_c4y_r[j_CY])-(ADC_C4Y_Thr_R[j_CY]+MWPC_Thr_Offset_C4Y);
      		ADC_C4Y_L[j_CY] = double(adc_c4y_l[j_CY])-(ADC_C4Y_Thr_L[j_CY]+MWPC_Thr_Offset_C4Y);
    	}


    	//******* GOOD TARGET EVENTS
    	Good_TARGET_Event = false;
    	count_TARGET_evts = 0;

    	for(int i=0; i<256; i++){
    		if((ADC_High_TARGET[i]>0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
           	   (ADC_High_TARGET[i]>0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
               (ADC_High_TARGET[i]>0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
               (ADC_High_TARGET[i]>0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
               (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0)) count_TARGET_evts++;
    	}

    	if(count_TARGET_evts >= n_hit) Good_TARGET_Event = true;
    	//if(enable_cout == 9) Good_TARGET_Event = true;


    	//******* GOOD TOF EVENTS
    	for(int ii=0; ii<12; ii++){   
    		Good_TOF1_ADC[ii]=false;
      		Good_TOF2_ADC[ii]=false;
      		Good_TOF1_TDC[ii]=false;
      		Good_TOF2_TDC[ii]=false;
      		Good_TOF1[ii]=false;
      		Good_TOF2[ii]=false;
      		Good_TOFs[ii]=false;
    	}
      
    	Good_TOF_Event = false;
    	Event_On_Blacklist = false;

    	//if(blacklist.fail()){
    	//  cout << "Error: Could not read blacklist file." << endl;  // TO CHECK !
    	//}
    	//else{
    	//  while(getline(blacklist,Current_Event)){
    	//    sscanf(Current_Event.c_str(), "%d", &current_event);
    	//    if(current_event == ievt){
    	//      Event_On_Blacklist = true;
   	 	//      break;
    	//    }
    	//  }
    	//} 
    	//blacklist.close();

    	for(int i =0; i<12; i++){
    		ADC_TOF1U[i] = ADC_TOF1[i];
     		ADC_TOF1D[i] = ADC_TOF1[i+12];
    	}

    	for(int i=0; i<12; i++){
    		if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;

    		if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
               (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;

      		if(Good_TOF1_TDC[i] || Good_TOF1_ADC[i]) Good_TOF1[i] = true;

      		if((ADC_TOF2AO[i]>=0 && ADC_TOF2AI[i]>=0) || (ADC_TOF2BO[i]>=0 && ADC_TOF2BI[i]>=0))  Good_TOF2_ADC[i] = true;

      		if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i])  &&
                (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
               ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  &&
                (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC[i] = true;

      		if(Good_TOF2_ADC[i] || Good_TOF2_TDC[i]) Good_TOF2[i] = true;
    	}

    	for(int k=0; k<12; k++){
      		if(k!=0 && k!=11){
        		if((Good_TOF2[k] && Good_TOF1[k-1]) || (Good_TOF2[k] && Good_TOF1[k]) || (Good_TOF2[k] && Good_TOF1[k+1])) Good_TOFs[k] = true;
      		}
    	}

    	if((Good_TOF2[0] && Good_TOF1[11]) || (Good_TOF2[0] && Good_TOF1[0]) || (Good_TOF2[0] && Good_TOF1[1]))  Good_TOFs[0] = true;

    	if((Good_TOF2[11] && Good_TOF1[10]) || (Good_TOF2[11] && Good_TOF1[11]) || (Good_TOF2[11] && Good_TOF1[0]))  Good_TOFs[11] = true;

    	for(int kk=0; kk<12; kk++){
    		if(Good_TOFs[kk]) Good_TOF_Event = true;
    	}
    	//if(enable_cout == 9) Good_TOF_Event = true;


    	//******* GOOD MWPC EVENTS
    	count_C2X = 0;    count_C2Y = 0;
    	count_C3X = 0;    count_C3Y = 0;
    	count_C4X = 0;    count_C4Y = 0;
    	Good_MWPC_Event = false;

    	for(int i=0; i<56; i++){
    		if(ADC_C2X_R[i]>0. || ADC_C2X_L[i]>0.) count_C2X++;
    	}

    	for(int ii=0; ii<16; ii++){
    		if(ADC_C2Y_R[ii]>0. || ADC_C2Y_L[ii]>0.) count_C2Y++;
    	}

    	for(int j=0; j<64; j++){
    		if(ADC_C3X_R[j]>0. || ADC_C3X_L[j]>0.) count_C3X++;
    	}

    	for(int jj=0; jj<16; jj++){
    		if(ADC_C3Y_R[jj]>0. || ADC_C3Y_L[jj]>0.) count_C3Y++;
    	}

    	for(int k=0; k<72; k++){
    		if(ADC_C4X_R[k]>0. || ADC_C4X_L[k]>0.) count_C4X++;
    	}

    	for(int kk=0; kk<16; kk++){
    		if(ADC_C4Y_R[kk]>0. || ADC_C4Y_L[kk]>0.) count_C4Y++;
    	}

    	if(count_C2X>0 && count_C2Y>0 && count_C3X>0 && count_C3Y>0 && count_C4X>0 && count_C4Y>0) Good_MWPC_Event = true;
    	//if(enable_cout == 9) Good_MWPC_Event = true;


    	//******* GOOD EVENTS
    	if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist){
    		Good_Event = true;
    		//cout << ievt << "  " << "GOOD !" << endl;
    	}
    	//else cout << ievt << "  " << "BAD !" << endl;


    	///// SKIP IF NOT A GOOD EVENT
    	if(Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event) continue;
    	if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event) continue;
    	if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event) continue;
    	if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event) continue;
    	if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event) continue;
    	if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event) continue;
    	if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event) continue;
    	
    	if(!Good_Event) return;


    	///// DETERMINE K-STOP TIME
    	TDC_average = -1;

  		for(int ii=0; ii<256; ii++){
    		max_index_all[ii]=0;
    		max_ADC_all[ii]=0;
  		}
  		max_index_flag=0;

  		for(int i = 0; i<256; i++){
    		max_ADC_all[i] = -100000000;
    		max_index_all[i] = -1;
  		}

  		for(int j = 0; j<256; j++){
    		for(int i=0; i<256; i++){
      			max_index_flag = 0;

      			for(int k = 0; k<256; k++){
        			if(i == max_index_all[k]){
        				max_index_flag = 1;
        				break;
        			}
      			}
      			if (max_index_flag == 1) continue;
      			else{
        			if(ADC_Low_TARGET[i]>max_ADC_all[j]){
          				max_index_all[j] = i;
          				max_ADC_all[j] = ADC_Low_TARGET[i];
        			}
      			}
    		}
  		}

  		//Calculate TDC Average
  		TDC_LG_max = -1;
  		TDC_LG_max2 = -1;
  		index_max1=0;
  		index_max2=0;

  		for(int i = 0; i<256; i++){
   			for(int j=0; j<6; j++){
      			if(TDC_LG_max == -1){
        			if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
          			TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        			index_max1 = max_index_all[i];
      			}
      			else if(TDC_LG_max2 == -1){
        			if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
          			TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
        			index_max2 = max_index_all[i];
      			}
    		}
 		}

  		if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
  		else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
  		else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;


  		kaon_TDC_min = TDC_average + TDC_Avg_Offset_min;
  		kaon_TDC_max = TDC_average + TDC_Avg_Offset_max;
  		TDC_min_Kstop = TDC_average + TDC_Kstop_Avg_Offset_min;
  		TDC_max_Kstop = TDC_average + TDC_Kstop_Avg_Offset_max;

  		TDC_min_TARGET = kaon_TDC_min;


  		for(int i=0; i<256; i++){
    		has_TDC_hit[i] = false;
    		has_TDC_hit_Kstop[i] = false;
  		}


  		for(Int_t i=0; i<256; i++){
    		for (Int_t k=0; k<4; k++) {
      			if((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
      			if((tdc_le_target[i][k]>=kaon_TDC_min) && (tdc_le_target[i][k]<=kaon_TDC_max)) has_TDC_hit[i] = true;
    		}
  		}


  		///// ???
  		Angle_ADC_cut = 0;
  		x_inc = 0;
  		y_inc = 0;
  		hit_count = 0;
  		int count = 0;

  		max_index = 0;
  		max_index2 = 0;
  		max_index3 = 0;
  		max_index4 = 0;

  		max_ADC = -100000000;
  		max_ADC2 = -100000000;
  		max_ADC3 = -100000000;
  		max_ADC4 = -100000000;

  		for(Int_t q=0; q<256; q++){
    		if(ADC_High_TARGET[q]>max_ADC){
      			max_index = q;
      			max_ADC = ADC_High_TARGET[q];
    		}
  		}

  		for(Int_t q=0; q<256; q++){
    		if(q == max_index) continue;
    		else{
      			if(ADC_High_TARGET[q]>max_ADC2){
        			max_index2 = q;
        			max_ADC2 = ADC_High_TARGET[q];
      			}
    		}
  		}

  		for(Int_t q=0; q<256; q++){
    		if((q == max_index) || (q == max_index2)) continue;
    		else{
      			if(ADC_High_TARGET[q]>max_ADC3) {
        			max_index3 = q;
        			max_ADC3 = ADC_High_TARGET[q];
      			}
    		}
  		}

  		for(Int_t q=0; q<256; q++){
    		if((q == max_index) || (q == max_index2) || (q == max_index3)) continue;
    		else{
      			if(ADC_High_TARGET[q]>max_ADC4) {
        			max_index4 = q;
        			max_ADC4 = ADC_High_TARGET[q];
      			}
    		}
  		}

  		x_cent = Xloc[max_index];
  		y_cent = Yloc[max_index];
  
  		for(int i=0; i<256; i++) hyp[i] = -1;

  		for(Int_t j=0; j<256; j++){
   	 		hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
  		}


  		///// DETERMINE TOF1 GAP
  		for(int j=0; j<12; j++){
    		has_ADC_TOF1_hit[j] = false;
    		has_TDC_TOF1_hit[j] = false;
    		has_ADC_TOF2_hit[j] = false;
    		has_TDC_TOF2_hit[j] = false;
    		has_both_ADC_TOF1_hit[j] = false;
    		has_both_TDC_TOF1_hit[j] = false;
  		}

  		/// GOOD TOF2
  		for(int i = 0; i < 12; i++){
    		if((ADC_TOF2AO[i]>0 && ADC_TOF2AI[i]>0) || (ADC_TOF2BO[i]>0 && ADC_TOF2BI[i]>0)) {has_ADC_TOF2_hit[i]=true;}
    		if((((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) && (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])))
    		||  (((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) && (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])))) {has_TDC_TOF2_hit[i]=true;}
  		}

  		/// GOOD TOF1
  		for(int i = 0; i < 12; i++){
    		if(ADC_TOF1U[i]>0 && ADC_TOF1D[i]>0) {has_both_ADC_TOF1_hit[i] = true;}
    		if(ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
    		if((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) && (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_both_TDC_TOF1_hit[i] = true;}
    		if((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
  		}

  		for(int jj=0; jj<12; jj++){
    		ADC_TOF1_hit[jj] = 0;
    		ADCTDC_TOF1_hit[jj] = 0;
    		ADC_TOF2_hit[jj] = 0;
    		ADCTDC_TOF2_hit[jj] = 0;
  		}

  		for(int k=0; k<12; k++) {
    		if(has_ADC_TOF1_hit[k]){
      			if(has_TDC_TOF1_hit[k]) ADCTDC_TOF1_hit[k]++;
      			else ADC_TOF1_hit[k]++;
    		}
    		if(has_ADC_TOF2_hit[k]){
      			if(has_TDC_TOF2_hit[k]) ADCTDC_TOF2_hit[k]++;
      			else ADC_TOF2_hit[k]++;
    		}
  		}

  		/// Determine which TOF2 is hit
   		selected_TOF2 = 0;

 		for(int i = 0; i<12; i++){
    		if(has_TDC_TOF2_hit[i] && has_ADC_TOF2_hit[i]) selected_TOF2 = i + 1;
  			}

  		if(selected_TOF2 == 0){
    		for(int i = 0; i<12; i++){
      			if(has_TDC_TOF2_hit[i] || has_ADC_TOF2_hit[i]) selected_TOF2 = i+1;
    		}
  		}

   		for(int kk=0; kk<12; kk++) gap_counter[kk] = 0;


  		/// GAP SCORING
		//if(scoring_type==2){
    	for(int i=0; i<12; i++){
      		if(ADC_TOF1U[i]>=0) gap_counter[i]++;
    	}

    	for(int i=0; i<12; i++){
      		if(ADC_TOF1D[i]>=0) gap_counter[i]++;
    	}

    	for(int i=0; i<12; i++){
      		if(TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) gap_counter[i]++;
    	}

    	for(int i=0; i<12; i++){
      		if(TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]) gap_counter[i]++;
    	}

    	for(int i=0; i<12; i++){
      		if(i!=0 && i!=11 && gap_counter[i]>0){
        		if(ADC_TOF2AO[i]>=0) gap_counter[i]++;
        		if(ADC_TOF2AO[i-1]>=0) gap_counter[i]++;
        		if(ADC_TOF2AO[i+1]>=0) gap_counter[i]++;
        
        		if(ADC_TOF2AI[i]>=0) gap_counter[i]++;
        		if(ADC_TOF2AI[i-1]>=0) gap_counter[i]++;
        		if(ADC_TOF2AI[i+1]>=0) gap_counter[i]++;
        
        		if(ADC_TOF2BO[i]>=0) gap_counter[i]++;
        		if(ADC_TOF2BO[i-1]>=0) gap_counter[i]++;
        		if(ADC_TOF2BO[i+1]>=0) gap_counter[i]++;
        
        		if(ADC_TOF2BI[i]>=0) gap_counter[i]++;
        		if(ADC_TOF2BI[i-1]>=0) gap_counter[i]++;
        		if(ADC_TOF2BI[i+1]>=0) gap_counter[i]++;
	    
	    		if(TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) gap_counter[i]++;
        		if(TDC_TOF2AO[i-1]>=TOF2AO_TDC_min[i-1] && TDC_TOF2AO[i-1]<=TOF2AO_TDC_max[i-1]) gap_counter[i]++;
        		if(TDC_TOF2AO[i+1]>=TOF2AO_TDC_min[i+1] && TDC_TOF2AO[i+1]<=TOF2AO_TDC_max[i+1]) gap_counter[i]++;

        		if(TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) gap_counter[i]++;
        		if(TDC_TOF2AI[i-1]>=TOF2AI_TDC_min[i-1] && TDC_TOF2AI[i-1]<=TOF2AI_TDC_max[i-1]) gap_counter[i]++;
        		if(TDC_TOF2AI[i+1]>=TOF2AI_TDC_min[i+1] && TDC_TOF2AI[i+1]<=TOF2AI_TDC_max[i+1]) gap_counter[i]++;

        		if(TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) gap_counter[i]++;
        		if(TDC_TOF2BO[i-1]>=TOF2BO_TDC_min[i-1] && TDC_TOF2BO[i-1]<=TOF2BO_TDC_max[i-1]) gap_counter[i]++;
        		if(TDC_TOF2BO[i+1]>=TOF2BO_TDC_min[i+1] && TDC_TOF2BO[i+1]<=TOF2BO_TDC_max[i+1]) gap_counter[i]++;

        		if(TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i]) gap_counter[i]++;
        		if(TDC_TOF2BI[i-1]>=TOF2BI_TDC_min[i-1] && TDC_TOF2BI[i-1]<=TOF2BI_TDC_max[i-1]) gap_counter[i]++;
        		if(TDC_TOF2BI[i+1]>=TOF2BI_TDC_min[i+1] && TDC_TOF2BI[i+1]<=TOF2BI_TDC_max[i+1]) gap_counter[i]++;
      		}
    	}

    	if(gap_counter[0]>0){
      		if(ADC_TOF2AO[0]>=0) gap_counter[0]++;
      		if(ADC_TOF2AO[1]>=0) gap_counter[0]++;
      		if(ADC_TOF2AO[11]>=0) gap_counter[0]++;

      		if(ADC_TOF2AI[0]>=0) gap_counter[0]++;
      		if(ADC_TOF2AI[1]>=0) gap_counter[0]++;
      		if(ADC_TOF2AI[11]>=0) gap_counter[0]++;

      		if(ADC_TOF2BO[0]>=0) gap_counter[0]++;
      		if(ADC_TOF2BO[1]>=0) gap_counter[0]++;
      		if(ADC_TOF2BO[11]>=0) gap_counter[0]++;

      		if(ADC_TOF2BI[0]>=0) gap_counter[0]++;
      		if(ADC_TOF2BI[1]>=0) gap_counter[0]++;
      		if(ADC_TOF2BI[11]>=0) gap_counter[0]++;

      		if(TDC_TOF2AO[0]>=TOF2AO_TDC_min[0] && TDC_TOF2AO[0]<=TOF2AO_TDC_max[0]) gap_counter[0]++;
      		if(TDC_TOF2AO[1]>=TOF2AO_TDC_min[1] && TDC_TOF2AO[1]<=TOF2AO_TDC_max[1]) gap_counter[0]++;
      		if(TDC_TOF2AO[11]>=TOF2AO_TDC_min[11] && TDC_TOF2AO[11]<=TOF2AO_TDC_max[11]) gap_counter[0]++;

      		if(TDC_TOF2AI[0]>=TOF2AI_TDC_min[0] && TDC_TOF2AI[0]<=TOF2AI_TDC_max[0]) gap_counter[0]++;
      		if(TDC_TOF2AI[1]>=TOF2AI_TDC_min[1] && TDC_TOF2AI[1]<=TOF2AI_TDC_max[1]) gap_counter[0]++;
      		if(TDC_TOF2AI[11]>=TOF2AI_TDC_min[11] && TDC_TOF2AI[11]<=TOF2AI_TDC_max[11]) gap_counter[0]++;

      		if(TDC_TOF2BO[0]>=TOF2BO_TDC_min[0] && TDC_TOF2BO[0]<=TOF2BO_TDC_max[0]) gap_counter[0]++;
      		if(TDC_TOF2BO[1]>=TOF2BO_TDC_min[1] && TDC_TOF2BO[1]<=TOF2BO_TDC_max[1]) gap_counter[0]++;
      		if(TDC_TOF2BO[11]>=TOF2BO_TDC_min[11] && TDC_TOF2BO[11]<=TOF2BO_TDC_max[11]) gap_counter[0]++;

      		if(TDC_TOF2BI[0]>=TOF2BI_TDC_min[0] && TDC_TOF2BI[0]<=TOF2BI_TDC_max[0]) gap_counter[0]++;
      		if(TDC_TOF2BI[1]>=TOF2BI_TDC_min[1] && TDC_TOF2BI[1]<=TOF2BI_TDC_max[1]) gap_counter[0]++;
      		if(TDC_TOF2BI[11]>=TOF2BI_TDC_min[11] && TDC_TOF2BI[11]<=TOF2BI_TDC_max[11]) gap_counter[0]++;
    	}	

    	if(gap_counter[11]>0){
      		if(ADC_TOF2AO[11]>=0) gap_counter[11]++;
      		if(ADC_TOF2AO[10]>=0) gap_counter[11]++;
      		if(ADC_TOF2AO[0]>=0) gap_counter[11]++;

      		if(ADC_TOF2AI[11]>=0) gap_counter[11]++;
      		if(ADC_TOF2AI[10]>=0) gap_counter[11]++;
      		if(ADC_TOF2AI[0]>=0) gap_counter[11]++;

      		if(ADC_TOF2BO[11]>=0) gap_counter[11]++;
      		if(ADC_TOF2BO[10]>=0) gap_counter[11]++;
      		if(ADC_TOF2BO[0]>=0) gap_counter[11]++;

      		if(ADC_TOF2BI[11]>=0) gap_counter[11]++;
      		if(ADC_TOF2BI[10]>=0) gap_counter[11]++;
      		if(ADC_TOF2BI[0]>=0) gap_counter[11]++;

      		if(TDC_TOF2AO[11]>=TOF2AO_TDC_min[11] && TDC_TOF2AO[11]<=TOF2AO_TDC_max[11]) gap_counter[11]++;
      		if(TDC_TOF2AO[10]>=TOF2AO_TDC_min[10] && TDC_TOF2AO[10]<=TOF2AO_TDC_max[10]) gap_counter[11]++;
      		if(TDC_TOF2AO[0]>=TOF2AO_TDC_min[0] && TDC_TOF2AO[0]<=TOF2AO_TDC_max[0]) gap_counter[11]++;

      		if(TDC_TOF2AI[11]>=TOF2AI_TDC_min[11] && TDC_TOF2AI[11]<=TOF2AI_TDC_max[11]) gap_counter[11]++;
      		if(TDC_TOF2AI[10]>=TOF2AI_TDC_min[10] && TDC_TOF2AI[10]<=TOF2AI_TDC_max[10]) gap_counter[11]++;
      		if(TDC_TOF2AI[0]>=TOF2AI_TDC_min[0] && TDC_TOF2AI[0]<=TOF2AI_TDC_max[0]) gap_counter[11]++;

      		if(TDC_TOF2BO[11]>=TOF2BO_TDC_min[11] && TDC_TOF2BO[11]<=TOF2BO_TDC_max[11]) gap_counter[11]++;
      		if(TDC_TOF2BO[10]>=TOF2BO_TDC_min[10] && TDC_TOF2BO[10]<=TOF2BO_TDC_max[10]) gap_counter[11]++;
      		if(TDC_TOF2BO[0]>=TOF2BO_TDC_min[0] && TDC_TOF2BO[0]<=TOF2BO_TDC_max[0]) gap_counter[11]++;

      		if(TDC_TOF2BI[11]>=TOF2BI_TDC_min[11] && TDC_TOF2BI[11]<=TOF2BI_TDC_max[11]) gap_counter[11]++;
      		if(TDC_TOF2BI[10]>=TOF2BI_TDC_min[10] && TDC_TOF2BI[10]<=TOF2BI_TDC_max[10]) gap_counter[11]++;
      		if(TDC_TOF2BI[0]>=TOF2BI_TDC_min[0] && TDC_TOF2BI[0]<=TOF2BI_TDC_max[0]) gap_counter[11]++;
    	}
  		//}


    	high_gap_hit = 0;
    	gap_to_fit = 0;
    	score_max = 0;
    	tof1_ties.clear();

    	for (int k=0; k<12; k++) {
      		//if((gap_counter[k]>=high_gap_hit) && (has_TDC_TOF1_hit[k] || has_ADC_TOF1_hit[k])){
      		if(gap_counter[k]>=high_gap_hit){
        		if(gap_counter[k] == high_gap_hit) tof1_ties.push_back(k); 
        		else{
          			tof1_ties.clear();
          			tof1_ties.push_back(k);
        		}
        		high_gap_hit = gap_counter[k];
		    }
    	}
    
    	if(tof1_ties.size() > 1){  // only when there are more than 2 elements in tof1_ties!
      		for(unsigned int j=0; j<tof1_ties.size(); j++){
        		for(int jj=0; jj<8; jj++){
          			if(ADC_High_TARGET[channel[tof1_ties[j]][jj]] > 0) gap_counter[tof1_ties[j]]++;
          			if(ADC_Low_TARGET[channel[tof1_ties[j]][jj]] > 0) gap_counter[tof1_ties[j]]++;
        		}
      		}
    	}

    	for(int k=0; k<12; k++){
      		if(gap_counter[k] > score_max){
        		score_max = gap_counter[k];
        		gap_to_fit = k+1;
      		}
    	} 
	    gap_to_fit_rotate = gap_to_fit;

    	//cout << ievt << "   " << gap_to_fit << endl; 

    	///// Determine K-STOP BARS
    	for(int jj=0; jj<256; jj++){
    		k_stop_bar[jj] = false;
    		k_stop_bar_initial[jj] = false;
  		}

  		good_k_stop_bars.clear();

  		for(int i = 0; i<256; i++){
    		if(ADC_High_TARGET[i]>0 || ADC_Low_TARGET[i]>0) vec_TARGET_bar_selected.push_back(i);
    		if(ADC_High_TARGET[i] > HG_KAON && ADC_Low_TARGET[i] > LG_KAON && has_TDC_hit_Kstop[i]){
      			good_k_stop_bars.push_back(i);
      			k_stop_bar[i] = true;
      			k_stop_bar_initial[i] = true;
    		}
  		}


  		///// KAON FIT
  		vec_xx_kaon.clear();	vec_yy_kaon.clear();
  		vec_ex_kaon.clear();	vec_ey_kaon.clear();
  		vec_kaon_bars.clear();

  		a_fit_kaon = 0.;
  		b_fit_kaon = 0.;
  		Chis_kaon = 0.;
  		ndf_kaon = 99;
  		kaon_bk = false;

  		for(int i = 0; i<256; i++){
    		if(k_stop_bar[i]){
      			vec_xx_kaon.push_back(Xloc[i]);
      			vec_ex_kaon.push_back(TARGET_Errors_X);

      			vec_yy_kaon.push_back(Yloc[i]);
      			vec_ey_kaon.push_back(TARGET_Errors_Y);

      			vec_kaon_bars.push_back(i);
    		}
  		}

  		gr_kaon = new TGraphErrors(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0],
                                     &vec_ey_kaon[0],&vec_ey_kaon[0]);

  		gr_kaon_bk = new TGraph(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0]);
  		gr_kaon_fit = new TF1("kaon_fit", "pol1");
  		gr_kaon->GetXaxis()->SetLimits(-50.,50.);

  		if(vec_xx_kaon.size()>=5){
    		gr_kaon->Fit("kaon_fit","Q");
    		gr_kaon_fit = gr_kaon->GetFunction("kaon_fit");
    		a_fit_kaon = gr_kaon_fit->GetParameter(1);
    		b_fit_kaon = gr_kaon_fit->GetParameter(0);
    		Chis_kaon = gr_kaon_fit->GetChisquare();
    		ndf_kaon = gr_kaon_fit->GetNDF();
  		}

  		if(a_fit_kaon > 1000){
    		kaon_bk = true;

	   		gr_kaon_bk->GetXaxis()->SetLimits(-50.,50.);
    		gr_kaon_bk->GetYaxis()->SetRangeUser(-50.,50.);

    		gr_kaon_bk->Fit("kaon_fit","Q");
    		gr_kaon_fit = gr_kaon_bk->GetFunction("kaon_fit");
    		a_fit_kaon = gr_kaon_fit->GetParameter(1);
    		b_fit_kaon = gr_kaon_fit->GetParameter(0);
  		}


  		///// CALCULATE C2X_CENTROID
  		C2X_centroid = 0.;
 
  		for(int i=0; i<16; i++){   
    		C2Y_L[i] = 0.;
    		C2Y_R[i] = 0.;
    		C3Y_L[i] = 0.;
    		C3Y_R[i] = 0.;
    		C4Y_L[i] = 0.;
    		C4Y_R[i] = 0.;
  		}

  		for(int j=0; j<56; j++){
    		C2X_L[j] = 0.;
    		C2X_R[j] = 0.;
  		}

  		for(int k=0; k<64; k++){
    		C3X_L[k] = 0.;
    		C3X_R[k] = 0.;
  		}

  		for(int l=0; l<72; l++){
    		C4X_L[l] = 0.;
    		C4X_R[l] = 0.;
  		}

  		//C2 Counters
  		for(int i=0; i<56; i++){
    		if(Good_Event && ADC_C2X_L[i]>0) C2X_L[i] = ADC_C2X_L[i];
    		if(Good_Event && ADC_C2X_R[i]>0) C2X_R[i] = ADC_C2X_R[i];
  		}

  		for(int i=0; i<64; i++){
    		if(Good_Event && ADC_C3X_L[i]>0) C3X_L[i] = ADC_C3X_L[i];
    		if(Good_Event && ADC_C3X_R[i]>0) C3X_R[i] = ADC_C3X_R[i];
  		}

  		for(int i=0; i<72; i++){
    		if(Good_Event && ADC_C4X_L[i]>0) C4X_L[i] = ADC_C4X_L[i];
    		if(Good_Event && ADC_C4X_R[i]>0) C4X_R[i] = ADC_C4X_R[i];
  		}

  		for(int i=0; i<16; i++){
    		if(Good_Event && ADC_C2Y_L[i]>0) C2Y_L[i] = ADC_C2Y_L[i];
    		if(Good_Event && ADC_C2Y_R[i]>0) C2Y_R[i] = ADC_C2Y_R[i];
    		if(Good_Event && ADC_C3Y_L[i]>0) C3Y_L[i] = ADC_C3Y_L[i];
    		if(Good_Event && ADC_C3Y_R[i]>0) C3Y_R[i] = ADC_C3Y_R[i];
    		if(Good_Event && ADC_C4Y_L[i]>0) C4Y_L[i] = ADC_C4Y_L[i];
    		if(Good_Event && ADC_C4Y_R[i]>0) C4Y_R[i] = ADC_C4Y_R[i];
    	}

  		first_cluster = true;
  		cluster_spacing = 0;
  		cluster_length_count = 0;

  		C2X_clusters = 0;
  		C2Y_clusters = 0;
  		C3X_clusters = 0;
  		C3Y_clusters = 0;
  		C4X_clusters = 0;
  		C4Y_clusters = 0;

  		C2X_cluster_index.clear(); 
  		C2X_cluster_length.clear();
  		C2Y_cluster_index.clear();
  		C2Y_cluster_length.clear();

  		C3X_cluster_index.clear();
  		C3X_cluster_length.clear();
  		C3Y_cluster_index.clear();
  		C3Y_cluster_length.clear();

  		C4X_cluster_index.clear();
  		C4X_cluster_length.clear();
  		C4Y_cluster_index.clear();
  		C4Y_cluster_length.clear();

  		if(selected_TOF2 > 6){  // LEFT
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<56; i++){
      			if(C2X_L[i] > 0. && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C2X_L[i] > 0.){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C2X_clusters++;
          				C2X_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 56 && C2X_L[i+1] <= 0. && C2X_L[i+2] <= 0.) C2X_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 56 && C2X_L[i+1] <= 0.) C2X_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 56) C2X_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C2X_L[i-1] <= 0.) cluster_length_count = 0;
        			else if (i != 0) cluster_length_count++;
      			}
    		}

    		// count C2Y clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<16; i++){
      			if(C2Y_L[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C2Y_L[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
        				C2Y_clusters++;
          				C2Y_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 16 && C2Y_L[i+1] <= 0 && C2Y_L[i+2] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 16 && C2Y_L[i+1] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 16) C2Y_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C2Y_L[i-1] <= 0) cluster_length_count = 0;
        			else if (i != 0) cluster_length_count++;
      			}
    		}

    		// count C3X clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<64; i++){
    	  		if(C3X_L[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C3X_L[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C3X_clusters++;
          				C3X_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 64 && C3X_L[i+1] <= 0 && C3X_L[i+2] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 64 && C3X_L[i+1] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 64) C3X_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C3X_L[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}

    		// count C3Y clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<16; i++){
      			if(C3Y_L[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C3Y_L[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C3Y_clusters++;
          				C3Y_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 16 && C3Y_L[i+1] <= 0 && C3Y_L[i+2] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 16 && C3Y_L[i+1] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 16) C3Y_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C3Y_L[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}

    		// count C4X clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<72; i++){
    	  		if(C4X_L[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C4X_L[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C4X_clusters++;
          				C4X_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 72 && C4X_L[i+1] <= 0 && C4X_L[i+2] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        			else if(i == 70 && C4X_L[71] <= 0) C4X_cluster_length.push_back(cluster_length_count);
       			 	else if(i == 71) C4X_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C4X_L[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}

    		// count C4Y clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<16; i++){
      			if(C4Y_L[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C4Y_L[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
        	  			C4Y_clusters++;
          				C4Y_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 16 && C4Y_L[i+1] <= 0 && C4Y_L[i+2] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 16 && C4Y_L[i+1] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 16) C4Y_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C4Y_L[i-1] <= 0) cluster_length_count = 0;
        			else if(i != 0) cluster_length_count++;
      			}
    		}
  		}
  		else{  // RIGHT
    		// count C2X clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<56; i++){
      			if(C2X_R[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C2X_R[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C2X_clusters++;
          				C2X_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 56 && C2X_R[i+1] <= 0 && C2X_R[i+2] <= 0) C2X_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 56 && C2X_R[i+1] <= 0) C2X_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 56)C2X_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C2X_R[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}

    		// count C2Y clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<16; i++){
      			if(C2Y_R[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C2Y_R[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C2Y_clusters++;
          				C2Y_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 16 && C2Y_R[i+1] <= 0 && C2Y_R[i+2] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 16 && C2Y_R[i+1] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 16) C2Y_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C2Y_R[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}

    		// count C3X clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<64; i++){
    			if(C3X_R[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C3X_R[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C3X_clusters++;
          				C3X_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 64 && C3X_R[i+1] <= 0 && C3X_R[i+2] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 64 && C3X_R[i+1] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 64) C3X_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C3X_R[i-1] <= 0) cluster_length_count = 0;
        			else if (i != 0) cluster_length_count++;
      			}
    		}

    		// count C3Y clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<16; i++){
      			if(C3Y_R[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C3Y_R[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C3Y_clusters++;
          				C3Y_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 16 && C3Y_R[i+1] <= 0 && C3Y_R[i+2] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 16 && C3Y_R[i+1] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 16) C3Y_cluster_length.push_back(cluster_length_count);
      			}
      			else{
       				cluster_spacing++;
        			if(i != 0 && C3Y_R[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}

    		// count C4X clusters
    		first_cluster = true;
    		cluster_spacing = 0;
    		cluster_length_count = 0;

    		for(int i = 0; i<72; i++){
      			if(C4X_R[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C4X_R[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C4X_clusters++;
          				C4X_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 72 && C4X_R[i+1] <= 0 && C4X_R[i+2] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        			else if(i == 70 && C4X_R[71] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        			else if(i == 71) C4X_cluster_length.push_back(cluster_length_count);
      			}	
      			else{
        			cluster_spacing++;
        			if(i != 0 && C4X_R[i-1] <= 0) cluster_length_count = 0;
        			else if ( i != 0) cluster_length_count++;
      			}
    		}


    		// count C4Y clusters
    		first_cluster = true;
    		cluster_spacing = 0;
   			cluster_length_count = 0;

    		for(int i = 0; i<16; i++){
      			if(C4Y_R[i] > 0 && first_cluster){
        			cluster_spacing = MWPC_cluster_separation + 1;
        			first_cluster = false;
      			}

      			if(C4Y_R[i] > 0){
        			if(cluster_spacing > MWPC_cluster_separation){
          				C4Y_clusters++;
          				C4Y_cluster_index.push_back(i);
        			}
        			cluster_length_count++;
        			cluster_spacing = 0;

        			if(i+2 < 16 && C4Y_R[i+1] <= 0 && C4Y_R[i+2] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 2 == 16 && C4Y_R[i+1] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        			else if(i + 1 == 16) C4Y_cluster_length.push_back(cluster_length_count);
      			}
      			else{
        			cluster_spacing++;
        			if(i != 0 && C4Y_R[i-1] <= 0) cluster_length_count = 0;
        			else if (i != 0) cluster_length_count++;
      			}
    		}
  		}

  		C2_centroid_num = 0;
  		C2_centroid_den = 0;

  		if(C2X_clusters == 1){
    		if(selected_TOF2>=7){
      			for(int i=0; i<56; i++){
      				if(C2X_L[i] > 0.){
          				C2_centroid_num += C2X_L[i]*C2_ZLoc[i];
          				C2_centroid_den += C2X_L[i];
        			}
      			}
    		}
    		else{
      			for(int i=0; i<56; i++){
        			if(C2X_R[i] > 0){
          				C2_centroid_num += C2X_R[i]*C2_ZLoc[i];
          				C2_centroid_den += C2X_R[i];
        			}
      			}
    		}

     		C2X_centroid = double(C2_centroid_num)/C2_centroid_den;
  		}


 		///// Ck AND Cpi
  		TDC_ck_sum = 0;       TDC_ck_avg = 0.;     TDC_ck_sigma = 0.;
  		TDC_cpi_sum = 0;      TDC_cpi_avg = 0.;    TDC_cpi_sigma = 0.;

  		TDC_ck_sigma2 = 0.;
  		TDC_cpi_sigma2 = 0.;

  		TDC_ck_sum2 = 0;    TDC_ck_avg2=0.;    TDC_ck_counter = 0;
  		TDC_cpi_sum2 = 0;   TDC_cpi_avg2=0.;   TDC_cpi_counter = 0;

  		 for(int i=0; i<14; i++){
    		TDC_ck_selected[i] = 0;
    		TDC_cpi_selected[i] = 0;
  		}
  		vec_Ck.clear();		vec_Cpi.clear();

  		for(int ic = 0; ic < 14; ic++){

    		for(int jc=0; jc<6; jc++){
      			if(tdc_ck[ic][jc]>=TDC_Ck_min && tdc_ck[ic][jc]<=TDC_Ck_max) TDC_ck_selected[ic]=tdc_ck[ic][jc];
      			if(tdc_cpi[ic][jc]>=TDC_Cpi_min && tdc_cpi[ic][jc]<=TDC_Cpi_max) TDC_cpi_selected[ic]=tdc_cpi[ic][jc];
    		}
    		if(TDC_ck_selected[ic]>0) vec_Ck.push_back(TDC_ck_selected[ic]);
    		if(TDC_cpi_selected[ic]>0) vec_Cpi.push_back(TDC_cpi_selected[ic]);
    		TDC_ck_selected[ic]=0;
    		TDC_cpi_selected[ic]=0;
  		}

  		for(unsigned int ik=0; ik<vec_Ck.size(); ik++) TDC_ck_sum += vec_Ck[ik];
  		for(unsigned int ip=0; ip<vec_Cpi.size(); ip++) TDC_cpi_sum += vec_Cpi[ip];

  		if(vec_Ck.size()>0) TDC_ck_avg = double(TDC_ck_sum)/double(vec_Ck.size());
  		else TDC_ck_avg = -1;

  		if(vec_Cpi.size()>0) TDC_cpi_avg = double(TDC_cpi_sum)/double(vec_Cpi.size());
  		else TDC_cpi_avg = -1;

  		for(unsigned int i=0; i<vec_Ck.size(); i++) TDC_ck_sigma2 += pow((vec_Ck[i]-TDC_ck_avg),2);
  		for(unsigned int j=0; j<vec_Cpi.size(); j++) TDC_cpi_sigma2 += pow((vec_Cpi[j]-TDC_cpi_avg),2);

  		if(vec_Ck.size()>0) TDC_ck_sigma = sqrt(TDC_ck_sigma2/vec_Ck.size());
  		else TDC_ck_sigma = -1;

  		if(vec_Cpi.size()>0) TDC_cpi_sigma = sqrt(TDC_cpi_sigma2/vec_Cpi.size());
  		else TDC_cpi_sigma = -1;

  		for(unsigned int i=0; i<vec_Ck.size(); i++){
    		if(abs(vec_Ck[i]-TDC_ck_avg) <= 1.4*TDC_ck_sigma){
      			TDC_ck_sum2 += vec_Ck[i];
      			TDC_ck_counter++;
    		}
  		}

  		for(unsigned int j=0; j<vec_Cpi.size(); j++){
    		if(abs(vec_Cpi[j]-TDC_cpi_avg) <= 1.4*TDC_cpi_sigma){
      			TDC_cpi_sum2 += vec_Cpi[j];
      			TDC_cpi_counter++;
    		}
  		}

	  	if(TDC_ck_counter>0) TDC_ck_avg2 = double(TDC_ck_sum2)/double(TDC_ck_counter);
  		else TDC_ck_avg2 = -1;

  		if(TDC_cpi_counter>0) TDC_cpi_avg2 = double(TDC_cpi_sum2)/double(TDC_cpi_counter);
  		else TDC_cpi_avg2 = -1;


  		///// DETERMINE IF A HIT TARGET HAS NEIGHBOURS
   		for(int i=0; i<256; i++) TARGET_High_has_neighbours[i] = false;

  		for(int i = 0; i<256; i++){
    		for(int j=0; j<8; j++){
      			if((TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] > Angle_ADC_cut && has_TDC_hit[TARGET_neighbours[i][j]] && !k_stop_bar[TARGET_neighbours[i][j]]) ||
         		   (TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] <0 && ADC_Low_TARGET[TARGET_neighbours[i][j]]>0 && !k_stop_bar[TARGET_neighbours[i][j]])){
          			TARGET_High_has_neighbours[i] = true;
          			break;
      			}
    		}
  		}


  		///// FILL LEPTON FIT 1
  		vec_xx_lepton.clear();		vec_ex_lepton.clear();
  		vec_yy_lepton.clear();		vec_ey_lepton.clear();
  		vec_lepton_bars.clear();	vec_lepton_size.clear();
  		count=0;	
  		hit_count=0;

  		for(Int_t i=0; i<256; i++){
    		if(TARGET_High_has_neighbours[i]){
      			if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
        			if(!k_stop_bar[i]){

          				if(gap_to_fit==6 || gap_to_fit==12){

            				vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
            				vec_ex_lepton.push_back(TARGET_Errors_X);

            				vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
            				vec_ey_lepton.push_back(TARGET_Errors_Y);

            				vec_lepton_bars.push_back(i);
          				}
          				else{
            				vec_xx_lepton.push_back(Xloc[i]);
            				vec_ex_lepton.push_back(TARGET_Errors_X);

            				vec_yy_lepton.push_back(Yloc[i]);
            				vec_ey_lepton.push_back(TARGET_Errors_Y);

            				vec_lepton_bars.push_back(i);
          				}

         				vec_lepton_size.push_back(Xloc[i]);
        			}
        			count++;

        			if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                  			channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                  			channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                  			channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
          				if(!k_stop_bar[i]){

            				if(gap_to_fit==6 || gap_to_fit==12){
              					vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
              					vec_ex_lepton.push_back(TARGET_Errors_X);
              					vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
              					vec_ex_lepton.push_back(TARGET_Errors_X);

              					vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
              					vec_ey_lepton.push_back(TARGET_Errors_Y);
              					vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
              					vec_ey_lepton.push_back(TARGET_Errors_Y);
            				}
            				else{
             				 	vec_xx_lepton.push_back(Xloc[i]);
              					vec_ex_lepton.push_back(TARGET_Errors_X);
              					vec_xx_lepton.push_back(Xloc[i]);
              					vec_ex_lepton.push_back(TARGET_Errors_X);

              					vec_yy_lepton.push_back(Yloc[i]);
              					vec_ey_lepton.push_back(TARGET_Errors_Y);
              					vec_yy_lepton.push_back(Yloc[i]);
              					vec_ey_lepton.push_back(TARGET_Errors_Y);
            				}

            				vec_lepton_size.push_back(Xloc[i]);
          				}
        			}
      			}

      			if(ADC_High_TARGET[i]>Angle_ADC_cut){
        			x_inc = x_inc + Xloc[i];
        			y_inc = y_inc + Yloc[i];
        			hit_count++;
      			}

      			if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0){
        			if(!k_stop_bar[i]){

          				if(gap_to_fit==6 || gap_to_fit==12){
            				vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
            				vec_ex_lepton.push_back(TARGET_Errors_X);

            				vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
            				vec_ey_lepton.push_back(TARGET_Errors_Y);

            				vec_lepton_bars.push_back(i);
          				}
          				else{
            				vec_xx_lepton.push_back(Xloc[i]);
            				vec_ex_lepton.push_back(TARGET_Errors_X);

            				vec_yy_lepton.push_back(Yloc[i]);
            				vec_ey_lepton.push_back(TARGET_Errors_Y);

            				vec_lepton_bars.push_back(i);
          				}

          				vec_lepton_size.push_back(Xloc[i]);
        			}
        			count++;
      			}
    		}
  		}


  		for(int g=0; g<12; g++){
    		Gap[g][0][0] = TOF_Xloc[3*g];
    		Gap[g][1][0] = TOF_Xloc[3*g+1];
    		Gap[g][2][0] = TOF_Xloc[3*g+2];

    		Gap[g][0][1] = TOF_Yloc[3*g];
    		Gap[g][1][1] = TOF_Yloc[3*g+1];
    		Gap[g][2][1] = TOF_Yloc[3*g+2];
  		}

  		for(int i = 0; i<3; i++){

    		if(gap_to_fit == 6){
      			vec_xx_lepton.push_back(Gap[2][1][0]);
      			vec_ex_lepton.push_back(TOF1_Errors_X[2][1]);

      			vec_yy_lepton.push_back(Gap[2][1][1]);
      			vec_ey_lepton.push_back(TOF1_Errors_Y[2][1]);
    		}	
    		else if(gap_to_fit == 12){
      			vec_xx_lepton.push_back(Gap[8][1][0]);
      			vec_ex_lepton.push_back(TOF1_Errors_X[2][1]);

      			vec_yy_lepton.push_back(Gap[8][1][1]);
      			vec_ey_lepton.push_back(TOF1_Errors_Y[2][1]);
    		}
    		else{
     		 	vec_xx_lepton.push_back(Gap[gap_to_fit-1][1][0]);
     		 	vec_ex_lepton.push_back(TOF1_Errors_X[gap_to_fit-1][1]);

      			vec_yy_lepton.push_back(Gap[gap_to_fit-1][1][1]);
      			vec_ey_lepton.push_back(TOF1_Errors_Y[gap_to_fit-1][1]);
    		}
  		}	

  		///// FIT LEPTON 1
  		a_lepton_fit_1 = 0.;
  		b_lepton_fit_1 = 0.;
  		Chis_lepton_fit_1 = 0.;
  		ndf_lepton_fit_1 = 0;

  		if(vec_xx_lepton.size()>=5){
  			gr_lepton_1 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
    								 	  &vec_ex_lepton[0],&vec_ey_lepton[0]);

  			func_lepton_fit_1 = new TF1("Lepton_fit_1", "pol1");

  			gr_lepton_1->GetXaxis()->SetLimits(-50.,50.);
  			gr_lepton_1->GetYaxis()->SetRangeUser(-50.,50.);

  			if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
    			func_lepton_fit_1->SetParameter(0,0);
    			func_lepton_fit_1->SetParameter(1,1);
  			}
  			if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
    			func_lepton_fit_1->SetParameter(0,0);
    			func_lepton_fit_1->SetParameter(1,-1);
  			}
  			else{
    			func_lepton_fit_1->SetParameter(0,0);
    			func_lepton_fit_1->SetParameter(1,1);
  			}

  			func_lepton_fit_1->SetParLimits(0,-50,50);
  			func_lepton_fit_1->SetParLimits(1,-50,50);

  			gr_lepton_1->Fit("Lepton_fit_1","Q");
  			func_lepton_fit_1 = gr_lepton_1->GetFunction("Lepton_fit_1");
  			a_lepton_fit_1 = func_lepton_fit_1->GetParameter(1);
  			b_lepton_fit_1 = func_lepton_fit_1->GetParameter(0);
  			Chis_lepton_fit_1 = func_lepton_fit_1->GetChisquare();
  			ndf_lepton_fit_1 = func_lepton_fit_1->GetNDF();
  		}
  		else continue;

  		///// FILL LEPTON FIT 2 
  		vec_xx_lepton.clear();            vec_ex_lepton.clear();
  		vec_yy_lepton.clear();            vec_ey_lepton.clear();
  		vec_xx_lepton_rotate.clear();     vec_ex_lepton_rotate.clear();
  		vec_yy_lepton_rotate.clear();     vec_ey_lepton_rotate.clear();
  		vec_lepton_size.clear();
  		lepton_counter = 0;
  		vec_lepton_bars_final.clear();

  		if(gap_to_fit == 6 || gap_to_fit == 12){
    		for(unsigned int i=0; i<vec_lepton_bars.size(); i++){
      			if(distance_to_line(Xloc[TARGET_Rotated_index[vec_lepton_bars[i]]],Yloc[TARGET_Rotated_index[vec_lepton_bars[i]]],a_lepton_fit_1,b_lepton_fit_1) <= max_dist
        		   && TARGET_High_has_neighbours[vec_lepton_bars[i]] && !k_stop_bar[vec_lepton_bars[i]]){
        			if(ADC_High_TARGET[vec_lepton_bars[i]]>Angle_ADC_cut && has_TDC_hit[vec_lepton_bars[i]]){

          				vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
          				vec_ex_lepton.push_back(TARGET_Errors_X);

          				vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
          				vec_ey_lepton.push_back(TARGET_Errors_Y);

          				vec_lepton_bars_rotate.push_back(vec_lepton_bars[i]);
          				vec_lepton_bars_final.push_back(vec_lepton_bars[i]);
          				lepton_counter++;

          				if(gap_to_fit == 6){
           	 				if(IsIn(TARGET_Rotated_index[vec_lepton_bars[i]],channel[3-1][0], channel[3-1][1],
                      			channel[3-1][2], channel[3-1][3],
                      			channel[3-1][4], channel[3-1][5],
                      			channel[3-1][6], channel[3-1][7])){

              					vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
              					vec_ex_lepton.push_back(TARGET_Errors_X);

              					vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
              					vec_ey_lepton.push_back(TARGET_Errors_Y);

              					//vec_lepton_bars_rotate.push_back(vec_lepton_bars[i]);
            				}
          				}

          				if(gap_to_fit == 12){
            				if(IsIn(TARGET_Rotated_index[vec_lepton_bars[i]],channel[9-1][0], channel[9-1][1],
                      			channel[9-1][2], channel[9-1][3],
                      			channel[9-1][4], channel[9-1][5],
                      			channel[9-1][6], channel[9-1][7])){

              					vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
              					vec_ex_lepton.push_back(TARGET_Errors_X);

              					vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
              					vec_ey_lepton.push_back(TARGET_Errors_Y);

              					//vec_lepton_bars_rotate.push_back(vec_lepton_bars[i]);
            				}
          				}
        			}

        			if(ADC_High_TARGET[vec_lepton_bars[i]]<0 && ADC_Low_TARGET[vec_lepton_bars[i]]>0){
          				vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
          				vec_ex_lepton.push_back(TARGET_Errors_X);

          				vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[vec_lepton_bars[i]]]);
          				vec_ey_lepton.push_back(TARGET_Errors_Y);

          				vec_lepton_bars_rotate.push_back(vec_lepton_bars[i]);
          				vec_lepton_bars_final.push_back(vec_lepton_bars[i]);
          				lepton_counter++;
        			}
      			}
    		}
  		}
  		else{
    		for(unsigned int i=0; i<vec_lepton_bars.size(); i++){
      			if(distance_to_line(Xloc[vec_lepton_bars[i]],Yloc[vec_lepton_bars[i]],a_lepton_fit_1,b_lepton_fit_1) <= max_dist
        			&& TARGET_High_has_neighbours[vec_lepton_bars[i]] && !k_stop_bar[vec_lepton_bars[i]]){
        			if(ADC_High_TARGET[vec_lepton_bars[i]]>Angle_ADC_cut && has_TDC_hit[vec_lepton_bars[i]]){

          				vec_xx_lepton.push_back(Xloc[vec_lepton_bars[i]]);
          				vec_ex_lepton.push_back(TARGET_Errors_X);

          				vec_yy_lepton.push_back(Yloc[vec_lepton_bars[i]]);
          				vec_ey_lepton.push_back(TARGET_Errors_Y);

          				vec_lepton_bars_final.push_back(vec_lepton_bars[i]);
          				lepton_counter++;

          				if(IsIn(vec_lepton_bars[i],channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],
                    		channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                    		channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                    		channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){

            				vec_xx_lepton.push_back(Xloc[vec_lepton_bars[i]]);
            				vec_ex_lepton.push_back(TARGET_Errors_X);

            				vec_yy_lepton.push_back(Yloc[vec_lepton_bars[i]]);
            				vec_ey_lepton.push_back(TARGET_Errors_Y);
          				}
        			}

        			if(ADC_High_TARGET[vec_lepton_bars[i]]<0 && ADC_Low_TARGET[vec_lepton_bars[i]]>0){
          				vec_xx_lepton.push_back(Xloc[vec_lepton_bars[i]]);
          				vec_ex_lepton.push_back(TARGET_Errors_X);

          				vec_yy_lepton.push_back(Yloc[vec_lepton_bars[i]]);
          				vec_ey_lepton.push_back(TARGET_Errors_Y);

          				vec_lepton_bars_final.push_back(vec_lepton_bars[i]);
          				lepton_counter++;
        			}
      			}
    		}
  		}

  		///// FIT LEPTON 2
  		a_lepton_fit_2 = 0.;
  		b_lepton_fit_2 = 0.;
  		Chis_lepton_fit_2 = 0.;
  		ndf_lepton_fit_2 = 0;

		gr_lepton_2 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
                                 &vec_ex_lepton[0],&vec_ey_lepton[0]);

  		func_lepton_fit_2 = new TF1("lepton_fit_2", "pol1");

  		if(vec_xx_lepton.size()>0){

    		gr_lepton_2->GetXaxis()->SetLimits(-50.,50.);
    		gr_lepton_2->GetYaxis()->SetRangeUser(-50.,50.);

    		if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
      			func_lepton_fit_2->SetParameter(0,0);
      			func_lepton_fit_2->SetParameter(1,1);
    		}
    		if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
      			func_lepton_fit_2->SetParameter(0,0);
      			func_lepton_fit_2->SetParameter(1,-1);
    		}
    		else{
      			func_lepton_fit_2->SetParameter(0,0);
      			func_lepton_fit_2->SetParameter(1,1);
    		}

    		func_lepton_fit_2->SetParLimits(0,-50,50);
    		func_lepton_fit_2->SetParLimits(1,-50,50);

    		gr_lepton_2->Fit("lepton_fit_2","Q");
    		func_lepton_fit_2 = gr_lepton_2->GetFunction("lepton_fit_2");

    		a_lepton_fit_2 = func_lepton_fit_2->GetParameter(1);
    		b_lepton_fit_2 = func_lepton_fit_2->GetParameter(0);
    		Chis_lepton_fit_2 = func_lepton_fit_2->GetChisquare();
    		ndf_lepton_fit_2 = func_lepton_fit_2->GetNDF();
  		}
  		//else cout << "Empty Fit 2" << endl;
  		else continue;


  		///// ???
  		for(int i=0; i<2; i++){
    		//x_int_TDC_selected_weighted[i] = 999.; y_int_TDC_selected_weighted[i] = 999.;
    		//x_int_GoodLG[i] = 999.;                y_int_GoodLG[i] = 999.;
    		//x_int_GoodLG_weighted[i] = 999.;       y_int_GoodLG_weighted[i] = 999.;
   			x_int_TDC_Gap_Fibers[i] = 999.;        y_int_TDC_Gap_Fibers[i] = 999.;
    		x_int_TDC_Gap_Fibers_SFT[i] = 999.;    y_int_TDC_Gap_Fibers_SFT[i] = 999.;
    		x_int_TARGET[i] = 999.;                y_int_TARGET[i] = 999.;
  		}

  		//x_int_GoodLG[0] = intersectx1(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);  // UNUSED
  		//x_int_GoodLG[1] = intersectx2(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
  		//y_int_GoodLG[0] = y1_int(x_int_GoodLG[0], a_fit_GoodLG, b_fit_GoodLG);
  		//y_int_GoodLG[1] = y2_int(x_int_GoodLG[1], a_fit_GoodLG, b_fit_GoodLG);

  		x_int_TDC_Gap_Fibers[0] = intersectx1(a_lepton_fit_2, b_lepton_fit_2, R_TOF1);
  		x_int_TDC_Gap_Fibers[1] = intersectx2(a_lepton_fit_2, b_lepton_fit_2, R_TOF1);
  		y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_lepton_fit_2, b_lepton_fit_2);
  		y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_lepton_fit_2, b_lepton_fit_2);

  		x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_lepton_fit_2, b_lepton_fit_2, R_SFT_L1);
  		x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_lepton_fit_2, b_lepton_fit_2, R_SFT_L1);
  		y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_lepton_fit_2, b_lepton_fit_2);
  		y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_lepton_fit_2, b_lepton_fit_2);

  		x_GoodLG_intersect1=0.;    y_GoodLG_intersect1=0.;
  		x_TDC_Gap_Fibers=0.;       y_TDC_Gap_Fibers=0.;

  		for(int i=0; i<2; i++){
    		dist1_GoodLG[i] = 999.;
    		dist1_TDC_Gap_Fibers[i] = 999.;
    		dist1_TDC_Gap_Fibers_SFT[i] = 999.;
  		}

  		for(int i=0; i<2; i++){
    		if(gap_to_fit == 6){
      			dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      			dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      			dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[3-1][1][0], Gap[3-1][1][1]);
    		}
    		else if(gap_to_fit == 12){
      			dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      			dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      			dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[9-1][1][0], Gap[9-1][1][1]);
    		}
    		else{
     			dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      			dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      			dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    		}
  		}

  		if(dist1_GoodLG[0] < dist1_GoodLG[1]){
    		x_GoodLG_intersect1 = x_int_GoodLG[0];
    		y_GoodLG_intersect1 = y_int_GoodLG[0];
  		}
  		else if(dist1_GoodLG[1] < dist1_GoodLG[0]){
    		x_GoodLG_intersect1 = x_int_GoodLG[1];
    		y_GoodLG_intersect1 = y_int_GoodLG[1];
  		}
  		//else cout << "ERROR !" << endl;


  		if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1]){
    		x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[0];
    		y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[0];
  		}
  		else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0]){
    		x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[1];
    		y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[1];
  		}
  		else cout << "ERROR !" << endl;

  		for(int j=0; j<3; j++){
  			dist2_TDC_selected[j] = 999.;
  			dist2_GoodLG[j] = 999.;
  		} 

  		dist2_TDC_selected_min = 999.;
  		dist2_GoodLG_min = 999.;
  		selected_TDC_selected = 0;

  		for(int ii=0; ii<3; ii++){
    		if(gap_to_fit == 6){
      			dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[3-1][ii][0], Gap[3-1][ii][1]);
      			dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[3-1][ii][0], Gap[3-1][ii][1]);
    		}
    		else if(gap_to_fit == 12){
      			dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[9-1][ii][0], Gap[9-1][ii][1]);
      			dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[9-1][ii][0], Gap[9-1][ii][1]);
    		}
    		else{
      			dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
      			dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
    		}

    		if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min){
      			dist2_TDC_selected_min = dist2_TDC_selected[ii];
      			selected_TDC_selected = ii;
    		}

    		if(dist2_GoodLG[ii] <= dist2_GoodLG_min) dist2_GoodLG_min = dist2_GoodLG[ii];
  		}


  		///// FILL LEPTON FIT 3		
  		for(int i = 0; i<3; i++){

    		if(gap_to_fit == 6){
      			vec_xx_lepton.push_back(Gap[3-1][selected_TDC_selected][0]);
      			vec_ex_lepton.push_back(TOF1_Errors_X[3-1][selected_TDC_selected]);

      			vec_yy_lepton.push_back(Gap[3-1][selected_TDC_selected][1]);
      			vec_ey_lepton.push_back(TOF1_Errors_Y[3-1][selected_TDC_selected]);
    		}
    		else if(gap_to_fit == 12){
      			vec_xx_lepton.push_back(Gap[9-1][selected_TDC_selected][0]);
      			vec_ex_lepton.push_back(TOF1_Errors_X[9-1][selected_TDC_selected]);

      			vec_yy_lepton.push_back(Gap[9-1][selected_TDC_selected][1]);
      			vec_ey_lepton.push_back(TOF1_Errors_Y[9-1][selected_TDC_selected]);
    		}
    		else{
    			vec_xx_lepton.push_back(Gap[gap_to_fit-1][selected_TDC_selected][0]);
    			vec_ex_lepton.push_back(TOF1_Errors_X[gap_to_fit-1][selected_TDC_selected]);

    			vec_yy_lepton.push_back(Gap[gap_to_fit-1][selected_TDC_selected][1]);
    			vec_ey_lepton.push_back(TOF1_Errors_Y[gap_to_fit-1][selected_TDC_selected]);
    		}
  		}


  		///// FIT LEPTON 3
  		a_lepton_fit_3 = 0.;
  		b_lepton_fit_3 = 0.;
  		Chis_lepton_fit_3 = 0.;
  		ndf_lepton_fit_3 = 0;
  		
  		ParError = 999.;
  		ChiS = 0.;
  		ndf = 0;

  		gr_lepton_3 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
                                         &vec_ex_lepton[0],&vec_ey_lepton[0]);
  
  		func_lepton_fit_3 = new TF1("lepton_fit_3", "pol1");
		gr_lepton_3->GetXaxis()->SetLimits(-50.,50.);
  		gr_lepton_3->GetYaxis()->SetRangeUser(-50.,50.);


  		if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
    		func_lepton_fit_3->SetParameter(0,0);
    		func_lepton_fit_3->SetParameter(1,1);
  		}
  		if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
    		func_lepton_fit_3->SetParameter(0,0);
    		func_lepton_fit_3->SetParameter(1,-1);
  		}
  		else{
      		func_lepton_fit_3->SetParameter(0,0);
      		func_lepton_fit_3->SetParameter(1,1);
    	}

  		func_lepton_fit_3->SetParLimits(0,-50,50);
  		func_lepton_fit_3->SetParLimits(1,-50,50);

  		gr_lepton_3->Fit("lepton_fit_3","Q");
  		func_lepton_fit_3 = gr_lepton_3->GetFunction("lepton_fit_3");

  		a_lepton_fit_3 = func_lepton_fit_3->GetParameter(1);
  		b_lepton_fit_3 = func_lepton_fit_3->GetParameter(0);

  		ParError = func_lepton_fit_3->GetParError(1);
  		Chis_lepton_fit_3 = func_lepton_fit_3->GetChisquare();
  		ndf_lepton_fit_3 = func_lepton_fit_3->GetNDF();


  		///// ???
		vec_xx_int_TDC_TARGET.clear();
		vec_yy_int_TDC_TARGET.clear();


  		x_int_TDC_Gap_Fibers[0] = intersectx1(a_lepton_fit_3, b_lepton_fit_3, R_TOF1);
  		x_int_TDC_Gap_Fibers[1] = intersectx2(a_lepton_fit_3, b_lepton_fit_3, R_TOF1);
  		x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_lepton_fit_3, b_lepton_fit_3, R_SFT_L1);
  		x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_lepton_fit_3, b_lepton_fit_3, R_SFT_L1);

  		y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_lepton_fit_3, b_lepton_fit_3);
  		y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_lepton_fit_3, b_lepton_fit_3);
  		y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_lepton_fit_3, b_lepton_fit_3);
  		y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_lepton_fit_3, b_lepton_fit_3);

  		x_int_TARGET[0] = intersectx1(a_lepton_fit_3, b_lepton_fit_3, R_TARGET);
  		x_int_TARGET[1] = intersectx2(a_lepton_fit_3, b_lepton_fit_3, R_TARGET);
  		y_int_TARGET[0] = y1_int(x_int_TARGET[0], a_lepton_fit_3, b_lepton_fit_3);
  		y_int_TARGET[1] = y2_int(x_int_TARGET[1], a_lepton_fit_3, b_lepton_fit_3);


  		x_TDC_Gap_Fibers_intersect1=0.;          y_TDC_Gap_Fibers_intersect1=0.;
  		x_TDC_Gap_Fibers_SFT_intersect1=0.;      y_TDC_Gap_Fibers_SFT_intersect1=0.;
  		x_TARGET_intersect=0;                    y_TARGET_intersect=0;
  		x_Arrows=0;                              y_Arrows=0;

  		dist1_TARGET_intersect[0] = 999.;
  		dist1_TARGET_intersect[1] = 999.;

  		for(int i=0; i<2; i++){
    		if(gap_to_fit==6){
      			dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      			dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      			dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[3-1][1][0], Gap[3-1][1][1]);
    		}
    		else if(gap_to_fit==12){
      			dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      			dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      			dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[9-1][1][0], Gap[9-1][1][1]);
    		}
    		else{
      			dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      			dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      			dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    		}
  		}

  		if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1]){
    		x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[0];
    		y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[0];
  		}
  		else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0]){
    		x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[1];
    		y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[1];
  		}
  		else cout << "ERROR !" << endl;

  		if(dist1_TDC_Gap_Fibers_SFT[0] < dist1_TDC_Gap_Fibers_SFT[1]){
    		x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[0];
    		y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[0];
  		}
  		else if(dist1_TDC_Gap_Fibers_SFT[1] < dist1_TDC_Gap_Fibers_SFT[0]){
    		x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[1];
    		y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[1];
  		}
  		else cout << "ERROR !" << endl;

  		if(dist1_TARGET_intersect[0] < dist1_TARGET_intersect[1]){  // UNUSED
    		x_TARGET_intersect = x_int_TARGET[0];
    		y_TARGET_intersect = y_int_TARGET[0];
    		x_Arrows = x_int_TARGET[1];
    		y_Arrows = y_int_TARGET[1];
  		}
  		else if(dist1_TARGET_intersect[1] < dist1_TARGET_intersect[0]){  // UNUSED
    		x_TARGET_intersect = x_int_TARGET[1];
    		y_TARGET_intersect = y_int_TARGET[1];
    		x_Arrows = x_int_TARGET[0];
    		y_Arrows = y_int_TARGET[0];
  		}
  		else cout << "ERROR !" << endl;


  		if(gap_to_fit==6 || gap_to_fit==12){  // TO MOVE

    		vec_xx_int_TDC_Gap_Fibers.push_back(y_TDC_Gap_Fibers_intersect1);
    		vec_yy_int_TDC_Gap_Fibers.push_back(-x_TDC_Gap_Fibers_intersect1);

    		vec_xx_int_TDC_Gap_Fibers_SFT.push_back(y_TDC_Gap_Fibers_SFT_intersect1);
    		vec_yy_int_TDC_Gap_Fibers_SFT.push_back(-x_TDC_Gap_Fibers_SFT_intersect1);

    		vec_xx_int_TDC_TARGET.push_back(y_TARGET_intersect);
    		vec_yy_int_TDC_TARGET.push_back(-x_TARGET_intersect);
  		}
  		else{   // TO MOVE
    		vec_xx_int_TDC_Gap_Fibers.push_back(x_TDC_Gap_Fibers_intersect1);
    		vec_yy_int_TDC_Gap_Fibers.push_back(y_TDC_Gap_Fibers_intersect1);

    		vec_xx_int_TDC_Gap_Fibers_SFT.push_back(x_TDC_Gap_Fibers_SFT_intersect1);
    		vec_yy_int_TDC_Gap_Fibers_SFT.push_back(y_TDC_Gap_Fibers_SFT_intersect1);

    		vec_xx_int_TDC_TARGET.push_back(x_TARGET_intersect);
    		vec_yy_int_TDC_TARGET.push_back(y_TARGET_intersect);
  		}


  		///// ANGLE CALCULATION
  		a_final_guide = 0.;
  		alpha_guide = 0.;
  		tanalpha_guide = 0.;
  		angle_final_guide = 0.;
  		Delta_phi = 999.99;  
  		Delta_phi_deg = 999.99;

  		if(gap_to_fit==6 || gap_to_fit==12) a_final_guide = a_lepton_fit_3;
  		else a_final_guide = a_lepton_fit_3;

  		tanalpha_guide = a_final_guide;
  		alpha_guide = atan(tanalpha_guide);

  		/// Determination of Delta Phi
  		Delta_phi = sin(alpha_guide)*cos(alpha_guide)*(ParError/a_final_guide);
  		Delta_phi_deg = (180/PI)*Delta_phi;


  		if(gap_to_fit==6 || gap_to_fit==12){
    		alpha_guide = alpha_guide*(180/PI) + 90;

    		if((-x_TDC_Gap_Fibers_intersect1 + x_TARGET_intersect) < 0.) angle_final_guide = 180. + alpha_guide;
    		else{

      			if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0) angle_final_guide = alpha_guide;
      			else angle_final_guide = alpha_guide + 360.;
    		}
  		}
  		else{

    		if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) < 0) angle_final_guide = 180. + alpha_guide * (180./PI);
    		else{

      			if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0) angle_final_guide = alpha_guide * (180./PI);
      			else angle_final_guide = alpha_guide * (180./PI) + 360.;
    		}
  		}

  		if(angle_final_guide > 360.) angle_final_guide -= 360.;


  		///// ???
  		Yth_TOF1=999.;
  		New_ChiS=999.;
  		ChiS_cos=999.;
  		a_fit_lepton_rotate = 999.;
  		b_fit_lepton_rotate = 999.;

  		a_fit_lepton_rotate = -1/a_lepton_fit_3;
  		b_fit_lepton_rotate = b_lepton_fit_3/a_lepton_fit_3;

  		// K Stop line / length in target
  		TLine *k_stop_line;

  		if(X_BAR != -10000){
    		if(angle_final_guide >= 90 && angle_final_guide <= 270){ // second x point in the left
      			k_stop_line = new TLine(X_BAR, Y_BAR, -50, tan(angle_final_guide*M_PI/180.0)*(-50-X_BAR) + Y_BAR);
    		}
    		else{ // second x point in the right
      			k_stop_line = new TLine(X_BAR, Y_BAR, 50, tan(angle_final_guide*M_PI/180.0)*(50-X_BAR) + Y_BAR);
    		}
  		}
  		else{
    		k_stop_line = new TLine(0,0,0,0);
  		}

  		k_stop_line->SetLineColor(3);
  		k_stop_line->SetLineWidth(2);

  		if(X_BAR != -10000){
    		x_tof1_intersect_1 = 999.;
    		y_tof1_intersect_1 = 999.;
    		x_tof1_intersect_2 = 999.;
    		y_tof1_intersect_2 = 999.;
    		x_tof1_intersect = 999.;
    		y_tof1_intersect = 999.;

    		alpha = angle_final_guide - 90.0;


    		//if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    		if(gap_to_fit_rotate==6 || gap_to_fit_rotate==12){
      			x_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
      			x_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
      			y_tof1_intersect_1 = -x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
      			y_tof1_intersect_2 = -x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);

      			if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
        			x_tof1_intersect = x_tof1_intersect_1;
        			y_tof1_intersect = y_tof1_intersect_1;
      			}
      			else{
        			x_tof1_intersect = x_tof1_intersect_2;
       	 			y_tof1_intersect = y_tof1_intersect_2;
      			}
    		}
    		else{
      			x_tof1_intersect_1 = x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
      			x_tof1_intersect_2 = x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
      			y_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
      			y_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);

      			if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
        			x_tof1_intersect = x_tof1_intersect_1;
        			y_tof1_intersect = y_tof1_intersect_1;
      			}
      			else{
        			x_tof1_intersect = x_tof1_intersect_2;
        			y_tof1_intersect = y_tof1_intersect_2;
      			}
    		}
  		}
  		else{
    		cout << "No K-Stop for track length." << endl;
  		}


  		///// CALCULATE XBAR and YBAR
  		X_weights = 999.;
  		Y_weights = 999.;
  		total_energy = 999.;

  		vec_kaon_centroid_coordinates.clear();
  		vec_fit_lines_intersect.clear();
  		vec_k_stop_coordinates.clear();
  		vec_xx_kaon_stop.clear();
  		vec_yy_kaon_stop.clear();

  		for(unsigned int k=0; k<vec_kaon_bars.size(); k++){
    		X_weights += ADC_Low_TARGET[vec_kaon_bars[k]]*Xloc[vec_kaon_bars[k]];
    		Y_weights += ADC_Low_TARGET[vec_kaon_bars[k]]*Yloc[vec_kaon_bars[k]];
    		total_energy += ADC_Low_TARGET[vec_kaon_bars[k]];
  		}

  		if(vec_kaon_bars.size()>0 || total_energy!=0){
    		vec_kaon_centroid_coordinates.push_back(X_weights/total_energy);
    		vec_kaon_centroid_coordinates.push_back(Y_weights/total_energy);
  		}
  		else{
    		vec_kaon_centroid_coordinates.push_back(999.99);
    		vec_kaon_centroid_coordinates.push_back(999.99);
  		}

  		if(vec_kaon_bars.size()>0){
    		if(gap_to_fit==6 || gap_to_fit==12){
      			vec_fit_lines_intersect = _2lines_intersect(a_fit_lepton_rotate, b_fit_lepton_rotate, a_fit_kaon, b_fit_kaon);
    		}
    		else{
      			vec_fit_lines_intersect = _2lines_intersect(a_lepton_fit_3, b_lepton_fit_3, a_fit_kaon, b_fit_kaon);
    		}
  		}
  		else{
    		vec_fit_lines_intersect.push_back(999.99);
    		vec_fit_lines_intersect.push_back(999.99);
  		}

  		if(vec_kaon_bars.size() < 5){
    		vec_k_stop_coordinates.push_back(vec_kaon_centroid_coordinates[0]);
    		vec_k_stop_coordinates.push_back(vec_kaon_centroid_coordinates[1]);
  		}
  		else{
    		vec_k_stop_coordinates.push_back(vec_fit_lines_intersect[0]);
    		vec_k_stop_coordinates.push_back(vec_fit_lines_intersect[1]);
  		}

  		vec_xx_kaon_stop.push_back(vec_k_stop_coordinates[0]);
  		vec_yy_kaon_stop.push_back(vec_k_stop_coordinates[1]);


  		i_kaon_bar = 999;
  		i_kaon_bar = kaon_fiber(vec_k_stop_coordinates[0],vec_k_stop_coordinates[1]);

  		length_in_target = 999.;
		length_in_target = distance(vec_xx_kaon_stop[0],vec_yy_kaon_stop[0],
                           vec_xx_int_TDC_TARGET[0],vec_yy_int_TDC_TARGET[0]); //CORRECT


  		///// CALCULATE TDC DIFF
  		tdc_ck_corr = tdc_trigger[0][0] - TDC_ck_avg2;
  		TDC_diff = 550 - (0.91*TDC_average - 0.625*tdc_ck_corr);


  		///// CALCULATE SUM_ADC_HG_LEPTON
  		sum_ADC_HG_lepton = 0;
 		for(unsigned int i=0; i<vec_lepton_bars_final.size(); i++){
    		sum_ADC_HG_lepton += ADC_High_TARGET[vec_lepton_bars_final[i]]+TARGET_ADC_Thr_HG_Offset;
  		}


  		///// CALCULATE AVERAGE TDCs
		sum_TDC_lepton = 0.;
  		sum_TDC_kaon = 0.;

  		counter_TDC_lepton = 0.;
  		counter_TDC_kaon = 0.;

  		Average_TDC_lepton = -1.;
  		Average_TDC_kaon = -1.;


  		for(unsigned int i=0; i<vec_lepton_bars_final.size(); i++){
    		for(int j=0; j<16; j++){
      			if(tdc_le_target[vec_lepton_bars_final[i]][j]>kaon_TDC_min && tdc_le_target[vec_lepton_bars_final[i]][j]<kaon_TDC_max){
        			sum_TDC_lepton += 550-(0.91*tdc_le_target[vec_lepton_bars_final[i]][j] - 0.625*(tdc_vt48[0][0]-TDC_ck_avg2));
        			counter_TDC_lepton++;
      			}
    		}
  		}

  		for(unsigned int i=0; i<vec_kaon_bars.size(); i++){
   			for(int j=0; j<16; j++){
      			if(tdc_le_target[vec_kaon_bars[i]][j]>TDC_min_Kstop && tdc_le_target[vec_kaon_bars[i]][j]<TDC_max_Kstop){
        			sum_TDC_kaon += 550-(0.91*tdc_le_target[vec_kaon_bars[i]][j] - 0.625*(tdc_vt48[0][0]-TDC_ck_avg2));
        			counter_TDC_kaon++;
      			}
    		}
  		}


  		if(counter_TDC_lepton>0) Average_TDC_lepton = sum_TDC_lepton/counter_TDC_lepton;
  		else Average_TDC_lepton = -1;

  		if(counter_TDC_kaon >0) Average_TDC_kaon = sum_TDC_kaon/counter_TDC_kaon;
  		else Average_TDC_kaon = -1;

  		///////////////////////////////////////////////////////////////////////////////////////////


  		/// DEBUG
  		if(Switch_Display==1){
 			cout << fixed;
  			cout << setw(4) << Run_Number << "  ";
  			cout << setw(7) << ievt << "  ";
  			cout << setw(2) << Event_tag[0] << "  ";
  			cout << setw(2) << gap_to_fit_rotate << "  ";
  			cout << setw(2) << selected_TOF2 << "  ";	
 			cout << setw(7) << setprecision(3) << angle_final_guide << "  ";
  			cout << setw(5) << setprecision(2) << Delta_phi_deg << "  ";
  			cout << setw(8) << setprecision(2) << Chis_lepton_fit_3 << "  ";
  			cout << setw(3) << setprecision(0) << ndf_lepton_fit_3 << "  ";
  			cout << setw(7) << setprecision(2) << Chis_lepton_fit_3/ndf_lepton_fit_3 << "  ";
  			cout << setw(3) << lepton_counter << "  ";
  			cout << setw(7) << setprecision(2) << vec_xx_int_TDC_TARGET[0] << "  ";
  			cout << setw(7) << setprecision(2) << vec_yy_int_TDC_TARGET[0] << "  ";
			cout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[0] << "  ";
  			cout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[1] << "  ";
  			cout << setw(6) << setprecision(2) << vec_fit_lines_intersect[0] << "  ";
  			cout << setw(6) << setprecision(2) << vec_fit_lines_intersect[1] << "  ";
  			cout << setw(6) << setprecision(2) << vec_k_stop_coordinates[0] << "  ";
  			cout << setw(6) << setprecision(2) << vec_k_stop_coordinates[1] << "  ";
  			cout << setw(3) << i_kaon_bar << "  ";
  			cout << setw(3) << vec_kaon_bars.size() << "  ";
  			cout << setw(6) << setprecision(2) << Chis_kaon/ndf_kaon << "  ";
  			cout << setw(3) << vec_Ck.size() << "  ";
  			cout << setw(3) << vec_Cpi.size() << "  ";
  			cout << setw(8) << setprecision(3) << length_in_target << "  ";
  			cout << setw(8) << setprecision(3) << C2X_centroid << "  ";
  			cout << setw(6) << setprecision(1) << TDC_diff << "  ";
  			cout << setw(7) << sum_ADC_HG_lepton << "  "; 
  			cout << setw(7) << setprecision(2) << Average_TDC_lepton << "  ";
  			cout << setw(7) << setprecision(2) << Average_TDC_kaon << "  ";
  			cout << endl;
  		}


  		if(Switch_Display==0 || Switch_Display==1 || Switch_Display==2){  		
 			fout << fixed;
  			fout << setw(4) << Run_Number << "  ";
  			fout << setw(7) << ievt << "  ";
  			fout << setw(2) << gap_to_fit_rotate << "  ";
  			fout << setw(2) << selected_TOF2 << "  "; 		
 			fout << setw(7) << setprecision(3) << angle_final_guide << "  ";
  			fout << setw(5) << setprecision(2) << Delta_phi_deg << "  ";
  			fout << setw(8) << setprecision(2) << Chis_lepton_fit_3 << "  ";
  			fout << setw(3) << setprecision(0) << ndf_lepton_fit_3 << "  ";
  			fout << setw(7) << setprecision(2) << Chis_lepton_fit_3/ndf_lepton_fit_3 << "  ";
  			fout << setw(3) << lepton_counter << "  ";
  			fout << setw(7) << setprecision(2) << vec_xx_int_TDC_TARGET[0] << "  ";
  			fout << setw(7) << setprecision(2) << vec_yy_int_TDC_TARGET[0] << "  ";
			fout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[0] << "  ";
  			fout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[1] << "  ";
  			fout << setw(6) << setprecision(2) << vec_fit_lines_intersect[0] << "  ";
  			fout << setw(6) << setprecision(2) << vec_fit_lines_intersect[1] << "  ";
  			fout << setw(6) << setprecision(2) << vec_k_stop_coordinates[0] << "  ";
  			fout << setw(6) << setprecision(2) << vec_k_stop_coordinates[1] << "  ";
  			fout << setw(3) << i_kaon_bar << "  ";
  			fout << setw(3) << vec_kaon_bars.size() << "  ";
  			fout << setw(6) << setprecision(2) << Chis_kaon/ndf_kaon << "  ";
  			fout << setw(3) << vec_Ck.size() << "  ";
  			fout << setw(3) << vec_Cpi.size() << "  ";
  			fout << setw(8) << setprecision(3) << length_in_target << "  ";
  			fout << setw(8) << setprecision(3) << C2X_centroid << "  ";
  			fout << setw(6) << setprecision(1) << TDC_diff << "  ";
  			fout << setw(7) << sum_ADC_HG_lepton << "  "; 
  			fout << setw(7) << setprecision(2) << Average_TDC_lepton << "  ";
  			fout << setw(7) << setprecision(2) << Average_TDC_kaon << "  ";
  			fout << endl;
  		}

  		delete gr_kaon;
		delete gr_kaon_bk;
		//delete gr_kaon_fit;
		delete gr_lepton_1;
		//delete func_lepton_fit_1;
		delete gr_lepton_2;
		//delete func_lepton_fit_2;
		delete gr2_Leptons_rotate;
		delete gr_lepton_3;
		//delete func_lepton_fit_3;
		delete gr3_Leptons_rotate;

	} // EndLoop ever Events
  ////////////////////////////////////////////////////////////////////////////////////////



}






/*
double SFT_Test_CR_Event(int Run_Number, int evt, double phi, bool to_print, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];
  Int_t adc_low_sft[128];
  Int_t tdc_le_sft[128][16];
  Int_t tdc_te_sft[128][16];

  Int_t HG_SFT_ADC_Thr[128] = {0};

  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

  Double_t ADC_High_SFT_corr[128];
  Int_t has_TDC_SFT_hit[128] = {0};

  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char path_input[200];
  sprintf(path_input,"%s",path_merged);

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  TChain *fChain= new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->GetEntry(evt);

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0) ADC_High_SFT_corr[j]=0;
    else ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  double sft_z_selected = 0.;
  sft_z_selected = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, to_print, 0, false, C2X_centroid, length_in_target);

  return sft_z_selected;
}

vector<double> SFT_Test_CR_Event2(int Run_Number, int evt, double phi, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];
  Int_t adc_low_sft[128];
  Int_t tdc_le_sft[128][16];
  Int_t tdc_te_sft[128][16];

  Int_t HG_SFT_ADC_Thr[128] = {0};


  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

  Double_t ADC_High_SFT_corr[128];
  Int_t has_TDC_SFT_hit[128] = {0};

  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char path_input[200];
  sprintf(path_input,"%s",path_merged);

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);


  TChain *fChain= new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->GetEntry(evt);

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0) ADC_High_SFT_corr[j]=0;
    else ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }

  for(Int_t ii=0; ii<128; ii++){
    for(Int_t qq=0; qq<6; qq++) {
      if(tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  //double sft_z_selected = 0.;
  vector<double> ZZ;
  ZZ = Z_Avg(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, phi, true, 0, false, C2X_centroid, length_in_target);

  //return sft_z_selected;
  return ZZ;
}
*/
