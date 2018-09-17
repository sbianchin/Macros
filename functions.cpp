#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#endif

bool areneighbours(int bar1, int bar2){

  bool bar_has_neighbours = false;

  for(int i=0; i<8; i++){
    if(TARGET_neighbours[bar1][i]==bar2) bar_has_neighbours = true;
  }
  
  return bar_has_neighbours;
}

int fun1(Int_t ADC_TOF1U[12], Int_t ADC_TOF1D[12], Int_t TDC_TOF1U[12], Int_t TDC_TOF1D[12],
		  Int_t ADC_TOF2AO[12], Int_t ADC_TOF2AI[12], Int_t ADC_TOF2BO[12], Int_t ADC_TOF2BI[12],
		  Int_t TDC_TOF2AO[12], Int_t TDC_TOF2AI[12], Int_t TDC_TOF2BO[12], Int_t TDC_TOF2BI[12]){

    int gap_counter[12] = {0};
  	int high_gap_hit = 0;
  	int gap_to_fit = 0;
  	int score_max = 0;

  	vector<int> tof1_ties;

    for(int i=0; i<12; i++){
      if(ADC_TOF1U[i]>=0) gap_counter[i]++;
      if(ADC_TOF1D[i]>=0) gap_counter[i]++;
      if(TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) gap_counter[i]++;
      if(TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]) gap_counter[i]++;
      if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) && 
        (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i])) gap_counter[i]++;
      if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) && 
        (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]) && 
        (ADC_TOF1U[i]>=0 && ADC_TOF1D[i]>=0)) gap_counter[i]++;
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


  tof1_ties.clear();

  for (int k=0; k<12; k++) {
    if(gap_counter[k]>=high_gap_hit){
      if(gap_counter[k] == high_gap_hit) tof1_ties.push_back(k);
      else{
        tof1_ties.clear();
        tof1_ties.push_back(k);
      }
      high_gap_hit = gap_counter[k];
    }
  }

  bool tie_breaker = false;




  if(tof1_ties.size() > 1){  // only when there are more than 2 elements in tof1_ties!
  	//cout << "THERE IS A TIE HERE !! " << endl;
    tie_breaker = true;
    for(unsigned int j=0; j<tof1_ties.size(); j++){
      for(int jj=0; jj<12; jj++){
        if(ADC_High_TARGET[channel[tof1_ties[j]][jj]]>=0){
          if(has_TDC_hit[channel[tof1_ties[j]][jj]]) gap_counter[tof1_ties[j]]++;
        }
        else{
          if(ADC_Low_TARGET[channel[tof1_ties[j]][jj]]>0) gap_counter[tof1_ties[j]]++;
        }
      }
    }
  }

  for(int k=0; k<12; k++){
  	//cout << "SCORES FUNCTION : " << k+1 << "  " << gap_counter[k] << endl;
  	if(gap_counter[k] > score_max){
  		score_max = gap_counter[k];
  		gap_to_fit = k+1;
  	}
  }

  //cout << endl;
  //return gap_to_fit;

  //cout << endl;
  //cout << "########################################################" << endl;
  //cout << endl;
  //cout << "TEST FUNCTION ALBANE : " << gap_to_fit << endl;
  //cout << endl;
  //cout << "########################################################" << endl;
  //cout << endl;

  return gap_to_fit;

}

vector<int> select_kaon(bool display=false){

  vector<int> vec_candidates;
  vector<int> vec_candidates_LG;
  vector<int> vec_output;

  int max1 = -100;  int kaon_max1 = -1;
  int max2 = -100;  int kaon_max2 = -1;
  int kaon_sub = -1; 

  int TDC_output = -1;
  int TDC_max1 = -1;
  int TDC_max2 = -1;


  for(int i=0; i<256; i++){
    if(ADC_High_TARGET[i]>HG_KAON){
      for(int j=0; j<6; j++){
        if(tdc_le_target[i][j]>=TDC_Thr_min && tdc_le_target[i][j]<=TDC_Thr_max){
          vec_candidates.push_back(i); 
        }         
      }
    }
  }

  if(vec_candidates.size()==0){
    for(int i=0; i<256; i++){
      //if(ADC_Low_TARGET[i]>max1){
      if(ADC_High_TARGET[i]>max1){
        for(int j=0; j<6; j++){
          if(tdc_le_target[i][j]>=TDC_Thr_min && tdc_le_target[i][j]<=TDC_Thr_max){
            //cout << "tdc_le_target[187] : " << tdc_le_target[187][j] << endl;
            //max1 = ADC_Low_TARGET[i];
            max1 = ADC_High_TARGET[i];
            kaon_max1 = i;
            kaon_sub = i;
          }
        }
      }
    }

    for(int j=0; j<6; j++){
      if(tdc_le_target[kaon_max1][j]>=TDC_Thr_min && tdc_le_target[kaon_max1][j]<=TDC_Thr_max){
        TDC_output = tdc_le_target[kaon_max1][j];
        TDC_max1 = tdc_le_target[kaon_max1][j];
      }
    }
    if(display){
      cout << "NO KAON !" << endl;
      cout << endl;
      printf("                  Fiber  ADC(HG)  ADC(LG)   TDC\n");
      printf("Kaon 1st Max. :   %4d    %4d     %4d    %4d\n",kaon_max1, ADC_High_TARGET[kaon_max1], ADC_Low_TARGET[kaon_max1], TDC_max1);
      printf("Kaon 2nd Max. :    ---    ----     ----    ----\n");
      printf("TDC Average   :  %4d\n",TDC_output);
      cout << "KAON SUBSTITUTE : " << kaon_sub << endl;
      cout << endl;
    }
  }

  else if(vec_candidates.size()==1){
    max1 = ADC_Low_TARGET[vec_candidates[0]];
    kaon_max1 = vec_candidates[0];

    for(int j=0; j<6; j++){
      if(tdc_le_target[kaon_max1][j]>=TDC_Thr_min && tdc_le_target[kaon_max1][j]<=TDC_Thr_max){
        TDC_output = tdc_le_target[kaon_max1][j];
        TDC_max1 = tdc_le_target[kaon_max1][j];
      }
    }

    if(max1<LG_KAON || (TDC_max1<TDC_Thr_min || TDC_max1>TDC_Thr_max)) kaon_sub = vec_candidates[0];

    if(display){
      cout << "ONLY 1 KAON !" << endl;
      cout << endl;
      printf("                  Fiber  ADC(HG)  ADC(LG)   TDC\n");
      printf("Kaon 1st Max. :   %4d    %4d     %4d    %4d\n",kaon_max1, ADC_High_TARGET[kaon_max1], ADC_Low_TARGET[kaon_max1], TDC_max1);
      printf("Kaon 2nd Max. :    ---    ----     ----    ----\n");
      printf("TDC Average   :  %4d\n",TDC_output);
      cout << "KAON SUBSTITUTE : " << kaon_sub << endl;
      cout << endl;
    }
  }
  else if(vec_candidates.size()>1){

    for(unsigned int i=0; i<vec_candidates.size(); i++){
      if(ADC_Low_TARGET[vec_candidates[i]]>max1){
        max1 = ADC_Low_TARGET[vec_candidates[i]];
        kaon_max1 = vec_candidates[i];  
      }
    }  

    vec_candidates.erase(remove(vec_candidates.begin(), vec_candidates.end(), kaon_max1), vec_candidates.end());

    for(unsigned int i=0; i<vec_candidates.size(); i++){
      if(ADC_Low_TARGET[vec_candidates[i]]>max2){
        max2 = ADC_Low_TARGET[vec_candidates[i]];
        kaon_max2 = vec_candidates[i];  
      }
    }

    for(int j=0; j<6; j++){
      if(tdc_le_target[kaon_max1][j]>=TDC_Thr_min && tdc_le_target[kaon_max1][j]<=TDC_Thr_max){
        TDC_max1 = tdc_le_target[kaon_max1][j];
      }
      if(tdc_le_target[kaon_max2][j]>=TDC_Thr_min && tdc_le_target[kaon_max2][j]<=TDC_Thr_max){
        TDC_max2 = tdc_le_target[kaon_max2][j];
      }
    }

    if(abs(TDC_max1 - TDC_max2) <= T_limit) TDC_output = (TDC_max1 + TDC_max2)/2;
    else if(max1 > max2) TDC_output = TDC_max1;
    else if(max1 < max2) TDC_output = TDC_max2;
    
    if(max1<LG_KAON || (TDC_max1<TDC_Thr_min || TDC_max1>TDC_Thr_max)) kaon_sub = kaon_max1;

    if(display){
      cout << "MORE THAN ONE KAON !" << endl;
      cout << endl;
      printf("                  Fiber  ADC(HG)  ADC(LG)   TDC\n");
      printf("Kaon 1st Max. :   %4d    %4d     %4d    %4d\n",kaon_max1, ADC_High_TARGET[kaon_max1], ADC_Low_TARGET[kaon_max1], TDC_max1);
      printf("Kaon 2nd Max. :   %4d    %4d     %4d    %4d\n",kaon_max2, ADC_High_TARGET[kaon_max2], ADC_Low_TARGET[kaon_max2], TDC_max2);
      printf("TDC Average   :  %4d\n",TDC_output);
      cout << "KAON SUBSTITUTE : " << kaon_sub << endl;
      cout << endl;
    }
  }  

  vec_output.push_back(kaon_sub);
  vec_output.push_back(TDC_output);

  return vec_output;
}

vector<double> get_slope(double x1, double y1, double x2, double y2){

  vector<double> v_output;
  double a;
  double b;
  
  a = (y2-y1)/(x2-x1);
  b = y1 - a*x1;

  v_output.push_back(a);
  v_output.push_back(b);

  return v_output;
}

double angle_calculation(double slope, double xstop, double ystop, double xtof1, double ytof1){

  double angle = -1.;
  double alpha = -1.;

  bool quad1 = false;
  bool quad2 = false;
  bool quad3 = false;
  bool quad4 = false;

  if(slope>=0 && xtof1-xstop>=0 && ytof1-ystop>=0) quad1 = true;
  if(slope<=0 && xtof1-xstop<=0 && ytof1-ystop>=0) quad2 = true;
  if(slope>=0 && xtof1-xstop<=0 && ytof1-ystop<=0) quad3 = true;
  if(slope<=0 && xtof1-xstop>=0 && ytof1-ystop<=0) quad4 = true;

  alpha = atan(slope)*(180/PI);

  if(quad1) angle = alpha;
  else if(quad2) angle = alpha + 180.;
  else if(quad3) angle = alpha + 180.;
  else if(quad4) angle = alpha + 360.;

  return angle;
}

vector<double> lineLeptonToTOF1(vector<int> vec_in, int TOF1_gap, double X_Kstop, double Y_Kstop){

  vector<double> vec_out;
  
  bool bars_are_neighbours = false;
  bool b_bar1 = false;
  bool b_bar2 = false;
  double dist1 = -1;
  double dist2 = -1;


  if(vec_in.size() == 0){
    vec_out.push_back(X_Kstop);
    vec_out.push_back(Y_Kstop);
    vec_out.push_back(TOF1_Xloc[TOF1_gap-1][1]);
    vec_out.push_back(TOF1_Yloc[TOF1_gap-1][1]);
  }

  if(vec_in.size() == 1){
    vec_out.push_back(Xloc[vec_in[0]]);
    vec_out.push_back(Yloc[vec_in[0]]);
    vec_out.push_back(TOF1_Xloc[TOF1_gap-1][1]);
    vec_out.push_back(TOF1_Yloc[TOF1_gap-1][1]);
  }

  if(vec_in.size() == 2){

    bars_are_neighbours = areneighbours(vec_in[0], vec_in[1]);

    if(bars_are_neighbours){
      vec_out.push_back(0.5*(Xloc[vec_in[0]]+Xloc[vec_in[1]]));
      vec_out.push_back(0.5*(Yloc[vec_in[0]]+Yloc[vec_in[1]]));
      vec_out.push_back(TOF1_Xloc[TOF1_gap-1][1]);   
      vec_out.push_back(TOF1_Yloc[TOF1_gap-1][1]);   
    }
    else{
      dist1 = sqrt(pow(Xloc[vec_in[0]]-X_Kstop,2)+pow(Yloc[vec_in[0]]-Y_Kstop,2));
      dist2 = sqrt(pow(Xloc[vec_in[1]]-X_Kstop,2)+pow(Yloc[vec_in[1]]-Y_Kstop,2));

      if(dist1 <= 10.) b_bar1 = true;
      if(dist2 <= 10.) b_bar2 = true;

      if(b_bar1 && !b_bar2){
        vec_out.push_back(Xloc[vec_in[0]]);
        vec_out.push_back(Yloc[vec_in[0]]);
      }
      else if(!b_bar1 && b_bar2){
        vec_out.push_back(Xloc[vec_in[1]]);
        vec_out.push_back(Yloc[vec_in[1]]);
      }
      else if(b_bar1 && b_bar2){
        if(dist1 < dist2){
          vec_out.push_back(Xloc[vec_in[0]]);
          vec_out.push_back(Yloc[vec_in[0]]);
        }
        else{
          vec_out.push_back(Xloc[vec_in[1]]);
          vec_out.push_back(Yloc[vec_in[1]]);        
        }
      }
      else if(!b_bar1 && !b_bar2){
        vec_out.push_back(X_Kstop);
        vec_out.push_back(Y_Kstop);
      }
    }

    vec_out.push_back(TOF1_Xloc[TOF1_gap-1][1]);
    vec_out.push_back(TOF1_Yloc[TOF1_gap-1][1]);
  }

  return vec_out;
}

double distance_2points(double x1, double y1, double x2, double y2){

  double dist;
  dist = sqrt(pow((x1-x2),2)+pow((y1-y2),2));
  return dist;
}

vector<double> intersect_circle(double a=0., double b=0., double R=44.5, double xtof1=-1., double ytof1=-1.){

  vector<double> intersect;
  double x1, y1;
  double x2, y2;
  double dist1, dist2;


  x1 = (-2*a*b - sqrt(4*a*a*b*b - 4*(a*a+1)*(b*b-R*R))) / (2*(a*a+1));
  y1 = a*x1+b;
  x2 = (-2*a*b + sqrt(4*a*a*b*b - 4*(a*a+1)*(b*b-R*R))) / (2*(a*a+1));
  y2 = a*x2+b;

  dist1 = distance_2points(x1, y1, xtof1, ytof1);
  dist2 = distance_2points(x2, y2, xtof1, ytof1);

  if(dist1 <= dist2){
    intersect.push_back(x1);
    intersect.push_back(y1);
  }
  else{
    intersect.push_back(x2);
    intersect.push_back(y2);
  }

  return intersect;
}


