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
#include "CommonParameters_7.2.h"
#include "ADC_TARGET_Pedestals.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr2.h"
#include "Pedestals.h"
#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"
#include "Plot_Event_Display_7.0.C"
#include "Event_Display_functions.cxx"
#include "Batch_Variables_7.2.h"
#include "functions_7.2.cpp"
//#include "functions_7.2.h"
#include "pruning_functions.cpp"
#include "drawing_pruning.cpp"
#include "drawing_fitting.cpp"
#include "drawing_Event_Display.cpp"
#include "printouts.cpp"

#include "functions_dev.h"

#endif

TChain fChain("Tree");

struct uncertainty{
  double dx;
  double dy;
};

double SFT_Test_CR_Event(int Run_Number, int evt, double phi, bool to_print, double C2X_centroid, double length_in_target);
vector<double> SFT_Test_CR_Event2(int Run_Number, int evt, double phi, double C2X_centroid, double length_in_target);
int kaon_fiber(double x_bar, double y_bar);

uncertainty calculate_uncertainty(double x, double y, double m, double b, double dm, double db) {
  double dy = pow(x * dm, 2) + pow(db, 2) + x * dm * db;
  double dx = pow((y-b)/(m*m) * dm, 2) + pow(db/m, 2) + (y-b)/(m*m*m)*dm*db;
  uncertainty d;
  d.dx = dx;
  d.dy = dy;
  return d;
}

bool check_vertical_kaon(vector<double> vec_xx_kaon) {
  bool vertical_kaon = true;
  for (unsigned int i = 1; i < vec_xx_kaon.size(); ++i) {
    if (vec_xx_kaon[i] != vec_xx_kaon[0]){
      vertical_kaon = false;
      break;
    }
  }
  return vertical_kaon;
}



void Event_Display(int Run_Number=3994, int ievt=0, int flag_pruning=0, int enable_cout=0, int Switch_Printout = 1,
  int Switch_Display = 1, int print_output = 1, int batch = 0, int batch_flag=0){

  // TH2F *h_empty = new TH2F("EMPTY", "EMPTY", 500, -50, 50, 500, -50, 50);
  // TCanvas *cc1;
  // cc1 = new TCanvas("EMPTY","EMPTY",50,50,800,800);
  // cc1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)"); 
  // cc1->Divide(2,2);


  bool flag_pruning_line;
  if(flag_pruning==0) flag_pruning_line = true;
  else flag_pruning_line = false;
  //bool flag_pruning_line = false;

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  has_edge_bars=false;
  bool intersect_outside_target = false;

  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int walk=1;
  bool Zlist_flag = false;
  char Version[100] = "Version X";


  if(enable_cout!=0 && enable_cout!=1 && enable_cout!=2 && enable_cout!=9){
    cout << "  " << endl;
    cout << "Flag Error !" << endl;
    cout << "  " << endl;
    return;
  }


  /////////////////////   Dave's Time Walk Correction File  ////////////////////////
  double par_in[256][3] = {0.};
  double par_err[356][3] = {0.};
  Int_t ADCmax = 3450;
  double Yfi = 0.;
  double Ani = 0.;
  double Yni = 0.;
  double Tpi = 0.;
  bool NoFile = false;

  char ParsTarg1[100];
  sprintf(ParsTarg1,"TimeWalk%d.dat",Run_Number);

  if(!ifstream(ParsTarg1)) NoFile = true;

  //////////////////////////////////////////////////////////////////////////////////


  const char* environment = getenv("path_merged");
  sprintf(Name_finput,"%s/Run%dMS.root",environment, Run_Number);

  ifstream ifile(Name_finput);
  if(!ifile.good()){ 
    cout << endl; 
    cout << endl;
    cout << "SORRY, THIS FILE DOESN'T EXIST !" << endl;
    cout << "(PLEASE, CHOOSE ANOTHER ONE)" << endl;
    cout << endl;
    cout << endl;
    return;
  }

  for(int ii=0; ii<256; ii++){
    HG_TARGET_ADC_Thr[ii] = 0;
    LG_TARGET_ADC_Thr[ii] = 0;
  }

  for(int jj=0; jj<128; jj++){
    HG_SFT_ADC_Thr[jj] = 0;
    LG_SFT_ADC_Thr[jj] = 0;
  }

  for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_HG[i]) + TARGET_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  TDC_min_SFT = SFT_TDC_min[0];
  TDC_max_SFT = SFT_TDC_max[0];

  for(int kk=0; kk<256; kk++){
    TARGET_High_has_neighbours[kk] = false;
    kaon_has_neighbours[kk] = false;
  }

  n_hit = 2;   // ## Number of hits required in the TARGET

  ADC_High_corr_max = 0;

  ADC_cut_TARGET2 = 850;

  TDC_diff = -1;
  tdc_ck_corr = 0.;

  for(int i=0; i<128; i++) has_TDC_SFT_hit[i] = 0;
  for(int j=0; j<40; j++) Event_flag[j] = false;
  for(int k=0; k<12; k++) Z_TOF1[k] = 999.99;

  sprintf(file_mapping,"../../Mapping");
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);
  sprintf(par_finput2,"%s/MWPC_map2.txt",file_mapping);

  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  ifstream fdat2(par_finput2,ios::in);
  for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
  fdat2.close();

  error_flag=1;

  vec_xx_lepton_pr.clear();
  vec_yy_lepton_pr.clear();
  vec_ex_lepton_pr.clear();
  vec_ey_lepton_pr.clear();

  vec_xx_lepton.clear();
  vec_yy_lepton.clear();
  vec_ex_lepton.clear();
  vec_ey_lepton.clear();

  vec_xx_lepton_test.clear();
  vec_yy_lepton_test.clear();

  vec_xx_lepton_rotate.clear();
  vec_yy_lepton_rotate.clear();
  vec_ex_lepton_rotate.clear();
  vec_ey_lepton_rotate.clear();

  vector<int> vec_lepton_bars_norepeat;
  vec_lepton_bars.clear();
  vec_lepton_bars_final.clear();
  vec_lepton_bars_rotate.clear();
  vec_lepton_bars_norepeat.clear();
  vec_kaon_bars.clear();

  vec_lepton_rotate_size.clear();
  vec_bar.clear();
  vec_bar_rotate.clear();
  vec_yprime.clear();
  vec_yprime_rotate.clear();
  vec_Dy.clear();
  vec_Dy_rotate.clear();

  sumS = 999.;
  sumS_rotate = 999.;

  vec_xx_kaon.clear();
  vec_yy_kaon.clear();
  vec_ex_kaon.clear();
  vec_ey_kaon.clear();

  vec_kaon_xx_rotate.clear();
  vec_kaon_yy_rotate.clear();
  vec_kaon_ex_rotate.clear();
  vec_kaon_ey_rotate.clear();

  vec_Ck.clear();
  vec_Cpi.clear();

  vec_xx_TOF1_Marker.clear();
  vec_yy_TOF1_Marker.clear();

  vec_xx_TOF1.clear();
  vec_yy_TOF1.clear();

  vec_xx_TOF1_closest.clear();
  vec_yy_TOF1_closest.clear();

  vec_xx_int_TDC_Gap_Fibers.clear();
  vec_yy_int_TDC_Gap_Fibers.clear();

  vec_xx_int_TDC_Gap_Fibers_SFT.clear();
  vec_yy_int_TDC_Gap_Fibers_SFT.clear();

  vec_xx_int_TDC_TARGET.clear();
  vec_yy_int_TDC_TARGET.clear();
  
  vec_xx_int_TOF1.clear();
  vec_yy_int_TOF1.clear();

  vec_xx_kaon_stop.clear();
  vec_yy_kaon_stop.clear();

  for(int i=0; i<14; i++){
    TDC_ck_selected[i] = 0;
    TDC_cpi_selected[i] = 0;
  }

  TDC_ck_sum = 0;       TDC_ck_avg = 0.;     TDC_ck_sigma = 0.;
  TDC_cpi_sum = 0;      TDC_cpi_avg = 0.;    TDC_cpi_sigma = 0.;

  TDC_ck_sigma2 = 0.;
  TDC_cpi_sigma2 = 0.;

  TDC_ck_sum2 = 0;    TDC_ck_avg2=0.;    TDC_ck_counter = 0;
  TDC_cpi_sum2 = 0;   TDC_cpi_avg2=0.;   TDC_cpi_counter = 0;

  vec_tdc_b0_6.clear();
  vec_tdc_b0_7.clear();
  vec_TARGET_bar_selected.clear();

  for(int i=0; i<12; i++){
    for(int j=0; j<5; j++){
      for(int k=0; k<2; k++){
        Gap[i][j][k] = 0;
      }
    }
  }

  if (batch == 0) {

  fChain.Reset();
  fChain.Add(Name_finput);
  fChain.SetMakeClass(1);
  fChain.SetBranchAddress("TDC_Trig",tdc_trigger); // tdc_trigger for TDC_diff calculation

  fChain.SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain.SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain.SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain.SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain.SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain.SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain.SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain.SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain.SetBranchAddress("ADC_TOF1U",ADC_tof1U);
  fChain.SetBranchAddress("ADC_TOF1D",ADC_tof1D);
  fChain.SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain.SetBranchAddress("TDC_TOF1D",TDC_tof1D);

  fChain.SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain.SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain.SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain.SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);
  fChain.SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain.SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain.SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain.SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);

  fChain.SetBranchAddress("MWPCADC",MwpcADC);

  fChain.SetBranchAddress("ADC_C2X_R",adc_c2x_r);
  fChain.SetBranchAddress("ADC_C2X_L",adc_c2x_l);
  fChain.SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
  fChain.SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
  fChain.SetBranchAddress("ADC_C3X_R",adc_c3x_r);
  fChain.SetBranchAddress("ADC_C3X_L",adc_c3x_l);
  fChain.SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
  fChain.SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
  fChain.SetBranchAddress("ADC_C4X_R",adc_c4x_r);
  fChain.SetBranchAddress("ADC_C4X_L",adc_c4x_l);
  fChain.SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
  fChain.SetBranchAddress("ADC_C4Y_L",adc_c4y_l);

  fChain.SetBranchAddress("TDC_Ck", tdc_ck);
  fChain.SetBranchAddress("TDC_Cpi", tdc_cpi);

  fChain.SetBranchAddress("VT48_TDC",tdc_vt48);

  fChain.SetBranchAddress("EvFlag", Event_flag);
  fChain.SetBranchAddress("EvTag", Event_tag);
  }

  Int_t nentries = 0;
  nentries = (Int_t)fChain.GetEntries();

  Good_Event=false;

  ofstream fout;
  char batch_output[50];

  ofstream fout_csv;
  char batch_output_csv[50];

  ofstream fout2;
  char batch_output2[50];


  if(batch==1 && (batch_flag== 0 || batch_flag==2 || batch_flag==3)){ // {print_output
    sprintf(batch_output, "Run_%d_Batch_output.txt", Run_Number);
    fout.open(batch_output, ios::app);

    sprintf(batch_output_csv, "Run_%d_Batch_output.csv", Run_Number);
    fout_csv.open(batch_output_csv, ios::app);   

    sprintf(batch_output2, "Run_%d_Zlist.txt", Run_Number);
    fout2.open(batch_output2, ios::app);  
  }
 

  for(int ivt=ievt; ivt<ievt+1; ivt++){
    fChain.GetEntry(ivt);
 
    ////////////////////////  Dave's Time Walk Correction /////////////////////// 

    if(walk==1){
      ifstream parTARGdat(ParsTarg1,ios::in);
      Int_t ij = 0;
      Int_t ik = 0;

      // Read in parameters and their errors (errors not used)
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
    }

    /////////////////////////////////////////////////////////////////////////////

    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
      TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
    }

    for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
      ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
      TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
      ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
    }

    for(int i=0; i<12; i++){
      ADC_TOF1[i] = ADC_tof1U[i]-TOF1U_ADC_Thr[i];
      ADC_TOF1[i+12] = ADC_tof1D[i]-TOF1D_ADC_Thr[i];
    }

    for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
      ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-TOF2AO_ADC_Thr[j_TOF2];
      ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-TOF2BO_ADC_Thr[j_TOF2];
      ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-TOF2AI_ADC_Thr[j_TOF2];
      ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-TOF2BI_ADC_Thr[j_TOF2];
   }

    for(Int_t j_C2=0; j_C2<56; j_C2++){
      ADC_C2X_R[j_C2] = double(adc_c2x_r[j_C2])-(ADC_C2X_Thr_R[j_C2]+MWPC_Thr_Offset_C2X);
      ADC_C2X_L[j_C2] = double(adc_c2x_l[j_C2])-(ADC_C2X_Thr_L[j_C2]+MWPC_Thr_Offset_C2X);
    }

    for(Int_t j_C3=0; j_C3<64; j_C3++){
      ADC_C3X_R[j_C3] = double(adc_c3x_r[j_C3])-(ADC_C3X_Thr_R[j_C3]+MWPC_Thr_Offset_C3X);
      ADC_C3X_L[j_C3] = double(adc_c3x_l[j_C3])-(ADC_C3X_Thr_L[j_C3]+MWPC_Thr_Offset_C3X);
    }

    for(Int_t j_C4=0; j_C4<72; j_C4++){
      ADC_C4X_R[j_C4] = double(adc_c4x_r[j_C4])-(ADC_C4X_Thr_R[j_C4]+MWPC_Thr_Offset_C4X);
      ADC_C4X_L[j_C4] = double(adc_c4x_l[j_C4])-(ADC_C4X_Thr_L[j_C4]+MWPC_Thr_Offset_C4X);
    }

    for(Int_t j_CY=0; j_CY<16; j_CY++){
      ADC_C2Y_R[j_CY] = double(adc_c2y_r[j_CY])-(ADC_C2Y_Thr_R[j_CY]+MWPC_Thr_Offset_C2Y);
      ADC_C2Y_L[j_CY] = double(adc_c2y_l[j_CY])-(ADC_C2Y_Thr_L[j_CY]+MWPC_Thr_Offset_C2Y);
      ADC_C3Y_R[j_CY] = double(adc_c3y_r[j_CY])-(ADC_C3Y_Thr_R[j_CY]+MWPC_Thr_Offset_C3Y);
      ADC_C3Y_L[j_CY] = double(adc_c3y_l[j_CY])-(ADC_C3Y_Thr_L[j_CY]+MWPC_Thr_Offset_C3Y);
      ADC_C4Y_R[j_CY] = double(adc_c4y_r[j_CY])-(ADC_C4Y_Thr_R[j_CY]+MWPC_Thr_Offset_C4Y);
      ADC_C4Y_L[j_CY] = double(adc_c4y_l[j_CY])-(ADC_C4Y_Thr_L[j_CY]+MWPC_Thr_Offset_C4Y);
    }

    for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
    }

    //********* GOOD TARGET EVENTS
    Good_TARGET_Event = false;
    count_TARGET_evts = 0;
    for(int i=0; i<256; i++)
    {
       if((ADC_High_TARGET[i]>0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0)) count_TARGET_evts++;
    }

    if(count_TARGET_evts >= n_hit) Good_TARGET_Event = true;

    if(enable_cout == 9) Good_TARGET_Event = true;

    //********* GOOD TOF EVENTS
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

    if(enable_cout == 9) Good_TOF_Event = true;


    //********* GOOD MWPC EVENTS
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

    if(count_C2X>0 && count_C2Y>0 && count_C3X>0 && count_C3Y>0 && count_C4X>0 && count_C4Y>0)  Good_MWPC_Event = true;

    if(enable_cout == 9) Good_MWPC_Event = true;

    if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist)  Good_Event = true;

    if(Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO MWPC!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      } 
      break;
    }

    if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO TOF!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }

    if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO TARGET!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }

    if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO TOF and NO MWPC!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }

    if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO TARGET and NO MWPC!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }

    if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO TARGET and NO TOF!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }

    if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << "   --  NO TARGET and NO TOF and NO MWPC!" << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }

    if(Event_On_Blacklist){    // TO CHECK
      if(Switch_Printout != 0){
        cout << endl;
        cout << "Event: "<< ievt << " is on the blacklist." << endl;
        cout << " " << endl;
        cout << ">>>>>  Please, choose another event" << endl;
        cout << " " << endl;
      }
      break;
    }
  }

  if(!Good_Event){
    //cout << "1111" << endl;
    if(batch==1){
      if(batch_flag==0 || batch_flag==2 || batch_flag==3){
        fout << endl;
        fout_csv << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," << endl;
      }
      if(batch_flag==1) cout << endl;
    }
  return;
  }
 
  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0;
    if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }

  ////////////////////  GEOMETRY !!!  ////////////////////////////////////////////

  //////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE
  TDC_average = -1;
  vector<int> k_select;   k_select.clear();
  int TDC_average_NEW = -1;
  int kaon_substitute = -1;
  bool has_kaon_sub = false;
  vector<double> vec_xx_kaon_substitute;            vec_xx_kaon_substitute.clear();
  vector<double> vec_yy_kaon_substitute;            vec_yy_kaon_substitute.clear();
  vector<double> vec_xx_kaon_substitute_rotated;    vec_xx_kaon_substitute_rotated.clear(); 
  vector<double> vec_yy_kaon_substitute_rotated;    vec_yy_kaon_substitute_rotated.clear(); 

  /// NEW WAY
  k_select = select_kaon(false);
  kaon_substitute = k_select[0];
  TDC_average_NEW = k_select[1];
  if(kaon_substitute > -1) has_kaon_sub = true;

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
      if(max_index_flag == 1) continue;

      else {
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
        if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max){
          TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        }
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

  int kaon_TDC_min;
  int kaon_TDC_max;
  int TDC_min_Kstop;
  int TDC_max_Kstop;

  if(TDC_average_NEW > 0){
    kaon_TDC_min = TDC_average_NEW + TDC_Avg_Offset_min;
    kaon_TDC_max = TDC_average_NEW + TDC_Avg_Offset_max;
    TDC_min_Kstop = TDC_average_NEW + TDC_Kstop_Avg_Offset_min;
    TDC_max_Kstop = TDC_average_NEW + TDC_Kstop_Avg_Offset_max;
  }
  else{
    kaon_TDC_min = TDC_Thr_min;
    kaon_TDC_max = TDC_Thr_max;
    TDC_min_Kstop = TDC_Thr_min;
    TDC_max_Kstop = TDC_Thr_max;   
  }

  TDC_min_TARGET = kaon_TDC_min;

  for(int i=0; i<256; i++){
    has_TDC_hit[i] = false;
    has_TDC_hit_Kstop[i] = false;
  }

  for(Int_t i=0; i<256; i++){
    for (Int_t k=0; k<6; k++) {
      if ((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
      if ((tdc_le_target[i][k]>=kaon_TDC_min) && (tdc_le_target[i][k]<=kaon_TDC_max)) has_TDC_hit[i] = true;
    }
  }

  for(Int_t ii=0; ii<128; ii++){
    if(ADC_High_SFT_corr[ii]>ADC_High_corr_max) ADC_High_corr_max=ADC_High_SFT_corr[ii];
  }

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++){
      if(tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }

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
    if (q == max_index) continue;
    else{
      if(ADC_High_TARGET[q]>max_ADC2){
        max_index2 = q;
        max_ADC2 = ADC_High_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2)) continue;
    else{
      if(ADC_High_TARGET[q]>max_ADC3) {
        max_index3 = q;
        max_ADC3 = ADC_High_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;
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

  for(int j=0; j<12; j++){
    has_ADC_TOF1_hit[j] = false;
    has_TDC_TOF1_hit[j] = false;
    has_ADC_TOF2_hit[j] = false;
    has_TDC_TOF2_hit[j] = false;
    has_both_ADC_TOF1_hit[j] = false;
    has_both_TDC_TOF1_hit[j] = false;
  }

  ///Set TOF2 Lines
  for(int i = 0; i < 12; i++){
    if ((ADC_TOF2AO[i]>0 && ADC_TOF2AI[i]>0) || (ADC_TOF2BO[i]>0 && ADC_TOF2BI[i]>0)) {has_ADC_TOF2_hit[i]=true;}
    if ((((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) && (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])))
    ||  (((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) && (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])))) {has_TDC_TOF2_hit[i]=true;}
  }

  ///Set TOF1 Lines
  for(int i = 0; i < 12; i++){
    if (ADC_TOF1U[i]>0 && ADC_TOF1D[i]>0) {has_both_ADC_TOF1_hit[i] = true;}
    if (ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) && (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_both_TDC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
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

  selected_TOF2 = 0;

  // Determine which TOF2 is hit
  for(int i = 0; i<12; i++){
    if(has_TDC_TOF2_hit[i] && has_ADC_TOF2_hit[i]) selected_TOF2 = i + 1;
  }

  if(selected_TOF2 == 0){
    for(int i = 0; i<12; i++){
      if(has_TDC_TOF2_hit[i] || has_ADC_TOF2_hit[i]) selected_TOF2 = i+1;
    }
  }

  // for(int kk=0; kk<12; kk++) gap_counter[kk] = 0;


  //// GAP SCORING !

  vector<int> vec_TOF1_Gap;   vec_TOF1_Gap.clear();
  int gap_to_fit;


  vec_TOF1_Gap = _get_TOF1_Gap(ADC_TOF1U, ADC_TOF1D, TDC_TOF1U, TDC_TOF1D,
                               ADC_TOF2AO, ADC_TOF2AI, ADC_TOF2BO, ADC_TOF2BI,
                               TDC_TOF2AO, TDC_TOF2AI, TDC_TOF2BO, TDC_TOF2BI);

  gap_to_fit = vec_TOF1_Gap[0];
  gap_to_fit_rotated = TOF1_rotated[gap_to_fit-1]+1;


  for(int jj=0; jj<256; jj++){
    k_stop_bar[jj] = false;
    k_stop_bar_initial[jj] = false;
  }

  good_k_stop_bars.clear();

  for(int i = 0; i<256; i++){
    if(ADC_High_TARGET[i]>=0 || ADC_Low_TARGET[i]>=0) vec_TARGET_bar_selected.push_back(i);
    if(ADC_High_TARGET[i] >= HG_KAON && ADC_Low_TARGET[i] >= LG_KAON && has_TDC_hit_Kstop[i]){
      good_k_stop_bars.push_back(i);
      k_stop_bar[i] = true;
      k_stop_bar_initial[i] = true;
    }
  }


  /// PRUNING 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  /// DETERMINE KAON BARS //////////////////////////////////////////////////////////////////////////////
  
  vector<int> V_kaon_bars_0;  V_kaon_bars_0.clear();
  vector<int> V_kaon_bars;    V_kaon_bars.clear();
  bool is_kaon[256];
  
  // Fill vector with EVERY possible kaon candidates
  for(int i = 0; i<256; i++){
    is_kaon[i] = false;
    if(ADC_High_TARGET[i] >= HG_KAON && ADC_Low_TARGET[i] >= LG_KAON && has_TDC_hit_Kstop[i]){
      is_kaon[i] = true;
      V_kaon_bars_0.push_back(i);  // keep tracks of every intial candidates
      V_kaon_bars.push_back(i);
    }
  }

  // Remove candidates that are single IF enough kaon candidates
  if(V_kaon_bars.size() > kaon_cluster_size){
    V_kaon_bars = _remove_singles(V_kaon_bars, is_kaon, is_kaon, kaon_substitute);
  }

  // Add substitute kaon if V_kaon_bars is EMPTY
  if(V_kaon_bars.size()==0 && has_kaon_sub){
    V_kaon_bars.push_back(kaon_substitute);
  }

  /// Fill Kaon Graph
  vector<double> vec_x_kaons;     vec_x_kaons.clear();
  vector<double> vec_y_kaons;     vec_y_kaons.clear();
  vector<double> vec_ex_kaons;    vec_ex_kaons.clear();
  vector<double> vec_ey_kaons;    vec_ey_kaons.clear();

  vector<double> vec_x_kaons_rot;     vec_x_kaons_rot.clear();
  vector<double> vec_y_kaons_rot;     vec_y_kaons_rot.clear();
  vector<double> vec_ex_kaons_rot;    vec_ex_kaons_rot.clear();
  vector<double> vec_ey_kaons_rot;    vec_ey_kaons_rot.clear();

  for(unsigned int i=0; i<V_kaon_bars.size(); i++){
    vec_x_kaons.push_back(Xloc[V_kaon_bars[i]]);
    vec_y_kaons.push_back(Yloc[V_kaon_bars[i]]);
    vec_ex_kaons.push_back(TARGET_Errors_X);
    vec_ey_kaons.push_back(TARGET_Errors_Y);
  }

  TGraph *gra_kaon_0;
  TGraph *gr_kaons_unfitted;
  TGraph *gr_kaons_final;

  
  gra_kaon_0 = new TGraphErrors(V_kaon_bars.size(), &vec_x_kaons[0], &vec_y_kaons[0],
                                &vec_ex_kaons[0], &vec_ey_kaons[0]);
  
  gr_kaons_unfitted = new TGraphErrors(V_kaon_bars.size(), &vec_x_kaons[0], &vec_y_kaons[0],
                                &vec_ex_kaons[0], &vec_ey_kaons[0]);

  gra_kaon_0->SetMarkerColor(4);
  gra_kaon_0->SetMarkerStyle(21);
  gra_kaon_0->SetMarkerSize(0.8);
  gra_kaon_0->GetYaxis()->SetTitle("Y title");
  gra_kaon_0->GetXaxis()->SetLimits(-50.,50.);
  gra_kaon_0->GetYaxis()->SetRangeUser(-50.,50.);

  gr_kaons_unfitted->SetMarkerColor(4);
  gr_kaons_unfitted->SetMarkerStyle(21);
  gr_kaons_unfitted->SetMarkerSize(0.8);
  gr_kaons_unfitted->GetYaxis()->SetTitle("Y title");
  gr_kaons_unfitted->GetXaxis()->SetLimits(-50.,50.);
  gr_kaons_unfitted->GetYaxis()->SetRangeUser(-50.,50.);
 
  double a_fit_kaon_0 = 999.99;
  double b_fit_kaon_0 = 999.99;
  double ChiS_kaon_0 = 999.99;
  int ndf_kaon_0 = 999;

  double a_fit_kaon_final = 999.99;
  double b_fit_kaon_final = 999.99;
  double ChiS_kaon_final = 999.99;
  int ndf_kaon_final = 999;

  bool to_rotate_kaons = false;

  TF1 *func_kaon_fit_0;
  func_kaon_fit_0 = new TF1("gr_kaon_0", "pol1");

  if(V_kaon_bars.size()>kaon_cluster_size){
    gra_kaon_0->Fit("gr_kaon_0","QW");  // <--- OPTION "W" TO BE REMOVED IF TOO MANY ERROR MESSAGES !
    gra_kaon_0->Fit("gr_kaon_0","Q");

    func_kaon_fit_0 = gra_kaon_0->GetFunction("gr_kaon_0");
    func_kaon_fit_0->SetLineWidth(2);
    func_kaon_fit_0->SetLineColor(4);

    a_fit_kaon_0 = func_kaon_fit_0->GetParameter(1);
    b_fit_kaon_0 = func_kaon_fit_0->GetParameter(0);
    ChiS_kaon_0 = func_kaon_fit_0->GetChisquare();
    ndf_kaon_0 = func_kaon_fit_0->GetNDF();

    if(abs(a_fit_kaon_0) > 10) to_rotate_kaons = true;
  }

  TGraph *gr_kaons_rot;
  TLine *line_fit_kaons;

  //to_rotate_kaons = true; /// <--- TO REMOVE 
  //to_rotate_kaons = false; /// <--- TO REMOVE 

  if(!to_rotate_kaons){
    gr_kaons_final = gra_kaon_0;
    line_fit_kaons = new TLine(99., 99., 99., 99.);
    
    a_fit_kaon_final = a_fit_kaon_0;
    b_fit_kaon_final = b_fit_kaon_0;
    ChiS_kaon_final = ChiS_kaon_0;
    ndf_kaon_final = ndf_kaon_0;
  }  
  else if(to_rotate_kaons){
    for(unsigned int j=0; j<V_kaon_bars.size(); j++){
      vec_x_kaons_rot.push_back(Xloc[TARGET_Rotated_index[V_kaon_bars[j]]]);
      vec_y_kaons_rot.push_back(Yloc[TARGET_Rotated_index[V_kaon_bars[j]]]);
      vec_ex_kaons_rot.push_back(TARGET_Errors_X);
      vec_ey_kaons_rot.push_back(TARGET_Errors_Y);
    }

    gr_kaons_rot = new TGraphErrors(V_kaon_bars.size(), &vec_x_kaons_rot[0], &vec_y_kaons_rot[0],
                                         &vec_ex_kaons_rot[0], &vec_ey_kaons_rot[0]);

    TF1 *func_kaon_fit_rot;
    func_kaon_fit_rot = new TF1("gr_kaon_rot", "pol1");

    gr_kaons_rot->Fit("gr_kaon_rot","QW");  // <--- OPTION "W" TO BE REMOVED IF TOO MANY ERROR MESSAGES !
    gr_kaons_rot->Fit("gr_kaon_rot","Q");

    func_kaon_fit_rot = gr_kaons_rot->GetFunction("gr_kaon_rot");
    func_kaon_fit_rot->SetLineWidth(2);
    func_kaon_fit_rot->SetLineColor(4);

    a_fit_kaon_0 = func_kaon_fit_rot->GetParameter(1);
    b_fit_kaon_0 = func_kaon_fit_rot->GetParameter(0);
    ChiS_kaon_0 = func_kaon_fit_rot->GetChisquare();
    ndf_kaon_0 = func_kaon_fit_rot->GetNDF();

    a_fit_kaon_final = -1/a_fit_kaon_0;  
    b_fit_kaon_final = b_fit_kaon_0/a_fit_kaon_0;   
    ChiS_kaon_final = ChiS_kaon_0 ;   
    ndf_kaon_final = ndf_kaon_0;   

    line_fit_kaons = new TLine(50.*a_fit_kaon_0 + b_fit_kaon_0, -50., 
                              -50.*a_fit_kaon_0 + b_fit_kaon_0, 50.);

    gr_kaons_final = gr_kaons_unfitted; 
  }



  //////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  /// DETERMINE LEPTON BARS ////////////////////////////////////////////////////////////////////////////

  bool is_lepton[256];
  vector<int> V_lepton_bars_0;  V_lepton_bars_0.clear();
  vector<int> V_lepton_bars;    V_lepton_bars.clear();


  for(int i=0; i<256; i++){
    is_lepton[i] = false;
    if(!is_kaon[i]){
      if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
        is_lepton[i] = true;
        V_lepton_bars_0.push_back(i);
        V_lepton_bars.push_back(i);
      }
      else if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0){
        is_lepton[i] = true;
        V_lepton_bars_0.push_back(i);
        V_lepton_bars.push_back(i);
      }
    }
  }


  // Remove single bars (bars with no neighbours)
  V_lepton_bars = _remove_singles(V_lepton_bars, is_lepton, is_kaon, kaon_substitute);


  //Remove kaon substitute if present
  if(has_kaon_sub){
    V_lepton_bars.erase(find(V_lepton_bars.begin(), V_lepton_bars.end(), kaon_substitute));
  }


  ///////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// CHECK IF NEED ROTATION  /////////////////////////////////////////////////////////////////////////////////////

  // LEPTONS
  vector<int> V_lepton_0;         V_lepton_0.clear();
  vector<double> V_x_lepton_0;    V_x_lepton_0.clear();
  vector<double> V_y_lepton_0;    V_y_lepton_0.clear();
  vector<double> V_ex_lepton_0;   V_ex_lepton_0.clear();
  vector<double> V_ey_lepton_0;   V_ey_lepton_0.clear();

  for(unsigned int i=0; i<V_lepton_bars.size(); i++){
    V_lepton_0.push_back(V_lepton_bars[i]);
    V_x_lepton_0.push_back(Xloc[V_lepton_bars[i]]);
    V_ex_lepton_0.push_back(TARGET_Errors_X);
    V_y_lepton_0.push_back(Yloc[V_lepton_bars[i]]);
    V_ey_lepton_0.push_back(TARGET_Errors_Y);   
  }

  // Add Edge Bars
  for(unsigned int j=0; j<V_lepton_bars.size(); j++){
    if(IsIn(V_lepton_bars[j],channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
              channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
              channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
              channel[gap_to_fit-1][6], channel[gap_to_fit-1][7],
              channel[gap_to_fit-1][8], channel[gap_to_fit-1][9],
              channel[gap_to_fit-1][10], channel[gap_to_fit-1][11])){

      has_edge_bars = true;
      
      for(int iw=0; iw<weight_edge_bars_fit1; iw++){  
        V_lepton_0.push_back(V_lepton_bars[j]);

        V_x_lepton_0.push_back(Xloc[V_lepton_bars[j]]);
        V_ex_lepton_0.push_back(TARGET_Errors_X);

        V_y_lepton_0.push_back(Yloc[V_lepton_bars[j]]);
        V_ey_lepton_0.push_back(TARGET_Errors_Y);
      }
    }
  }

  // Add TOF! with weight
  for(int k=0; k<weight_TOF1_fit1; k++){
    V_x_lepton_0.push_back(TOF1_Xloc[gap_to_fit-1][2]);
    V_ex_lepton_0.push_back(TOF1_Errors_X[gap_to_fit-1][2]);
  
    V_y_lepton_0.push_back(TOF1_Yloc[gap_to_fit-1][2]);
    V_ey_lepton_0.push_back(TOF1_Errors_Y[gap_to_fit-1][2]);
  }

  TGraph *graph_lepton_0;
  TF1 *func_lepton_fit_0;

  TGraph *N_graph_lepton_0;
  TF1 *N_func_lepton_fit_0;

  graph_lepton_0 = new TGraphErrors(V_x_lepton_0.size(),&V_x_lepton_0[0],&V_y_lepton_0[0],
                                   &V_ex_lepton_0[0],&V_ey_lepton_0[0]);

  double a_fit_lepton_0 = 999.;
  double b_fit_lepton_0 = 999.;
  double ChiS_lepton_0 = 999.;
  double ndf_lepton_0 = 999.;

  func_lepton_fit_0 = new TF1("gr_lepton_0", "pol1");

  graph_lepton_0->SetMarkerColor(2);
  graph_lepton_0->SetMarkerStyle(21);
  graph_lepton_0->SetMarkerSize(0.8);
  graph_lepton_0->GetYaxis()->SetTitle("Y title");
  graph_lepton_0->GetXaxis()->SetLimits(-50.,50.);
  graph_lepton_0->GetYaxis()->SetRangeUser(-50.,50.);

  //if(vec_lepton_bars_norepeat.size()>lepton_cluster_size){ 
  if(V_lepton_bars.size()>lepton_cluster_size){ 
    graph_lepton_0->Fit("gr_lepton_0","QW");  // <--- OPTION "W" TO BE REMOVED IF TOO MANY ERROR MESSAGES !
    graph_lepton_0->Fit("gr_lepton_0","Q");

    func_lepton_fit_0 = graph_lepton_0->GetFunction("gr_lepton_0");
    func_lepton_fit_0->SetLineWidth(2);
    func_lepton_fit_0->SetLineColor(2);
    
    a_fit_lepton_0 = func_lepton_fit_0->GetParameter(1);
    b_fit_lepton_0 = func_lepton_fit_0->GetParameter(0);
    ChiS_lepton_0 = func_lepton_fit_0->GetChisquare();
    ndf_lepton_0 = func_lepton_fit_0->GetNDF();
  }
  
  bool to_rotate_leptons = false;
  if(V_lepton_bars.size()>lepton_cluster_size && abs(a_fit_lepton_0) > 10) to_rotate_leptons = true;

  vector<double> vec_x_leptons_rot;     vec_x_leptons_rot.clear();
  vector<double> vec_y_leptons_rot;     vec_y_leptons_rot.clear();
  vector<double> vec_ex_leptons_rot;    vec_ex_leptons_rot.clear();
  vector<double> vec_ey_leptons_rot;    vec_ey_leptons_rot.clear();

  TGraph *gr_leptons_unfitted;
  gr_leptons_unfitted = new TGraphErrors(V_x_lepton_0.size(), &V_x_lepton_0[0], &V_y_lepton_0[0],
                                        &V_ex_lepton_0[0], &V_ey_lepton_0[0]);

  double a_fit_lepton_final = 999.99;
  double b_fit_lepton_final = 999.99;
  double ChiS_lepton_final = 999.99;
  int ndf_lepton_final = 999;

  TGraph *gr_leptons_0_NEW;
  TGraph *gr_leptons_rot;
  TLine *line_fit_leptons;
  int TOF1_gap = -1; 

  if(!to_rotate_leptons){
    gr_leptons_0_NEW = graph_lepton_0;
    line_fit_leptons = new TLine(99., 99., 99., 99.);
    
    a_fit_lepton_final = a_fit_lepton_0;
    b_fit_lepton_final = b_fit_lepton_0;
    ChiS_lepton_final = ChiS_lepton_0;
    ndf_lepton_final = ndf_lepton_0;
  }  
  else if(to_rotate_leptons){

    TOF1_gap = TOF1_rotated[gap_to_fit-1]+1;

    for(unsigned int i=0; i<V_lepton_bars.size(); i++){
      vec_x_leptons_rot.push_back(Xloc[TARGET_Rotated_index[V_lepton_bars[i]]]);
      vec_ex_leptons_rot.push_back(TARGET_Errors_X);
      vec_y_leptons_rot.push_back(Yloc[TARGET_Rotated_index[V_lepton_bars[i]]]);
      vec_ey_leptons_rot.push_back(TARGET_Errors_Y);   
    }

  // Add Edge Bars
  for(unsigned int j=0; j<V_lepton_bars.size(); j++){
    if(IsIn(V_lepton_bars[j],channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
              channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
              channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
              channel[gap_to_fit-1][6], channel[gap_to_fit-1][7],
              channel[gap_to_fit-1][8], channel[gap_to_fit-1][9],
              channel[gap_to_fit-1][10], channel[gap_to_fit-1][11])){

      has_edge_bars = true;
      
      for(int iw=0; iw<weight_edge_bars_fit1; iw++){  
        vec_x_leptons_rot.push_back(Xloc[TARGET_Rotated_index[V_lepton_bars[j]]]);
        vec_ex_leptons_rot.push_back(TARGET_Errors_X);

        vec_y_leptons_rot.push_back(Yloc[TARGET_Rotated_index[V_lepton_bars[j]]]);
        vec_ey_leptons_rot.push_back(TARGET_Errors_Y);
      }
    }
  }

  // Add TOF1 with weight
  for(int k=0; k<weight_TOF1_fit1; k++){
    vec_x_leptons_rot.push_back(TOF1_Xloc[TOF1_gap-1][2]);
    vec_ex_leptons_rot.push_back(TOF1_Errors_X[TOF1_gap-1][2]);
  
    vec_y_leptons_rot.push_back(TOF1_Yloc[TOF1_gap-1][2]);
    vec_ey_leptons_rot.push_back(TOF1_Errors_Y[TOF1_gap-1][2]);
  }

    gr_leptons_rot = new TGraphErrors(vec_x_leptons_rot.size(), &vec_x_leptons_rot[0], &vec_y_leptons_rot[0],
                                    &vec_ex_leptons_rot[0], &vec_ey_leptons_rot[0]);

    TF1 *func_lepton_fit_rot;
    func_lepton_fit_rot = new TF1("gr_lepton_rot", "pol1");

    //gr_leptons_rot->Fit("gr_lepton_rot","QW");  // <--- OPTION "W" TO BE REMOVED IF TOO MANY ERROR MESSAGES !
    gr_leptons_rot->Fit("gr_lepton_rot","Q");

    func_lepton_fit_rot = gr_leptons_rot->GetFunction("gr_lepton_rot");
    func_lepton_fit_rot->SetLineWidth(2);
    func_lepton_fit_rot->SetLineColor(4);

    a_fit_lepton_0 = func_lepton_fit_rot->GetParameter(1);
    b_fit_lepton_0 = func_lepton_fit_rot->GetParameter(0);
    ChiS_lepton_0 = func_lepton_fit_rot->GetChisquare();
    ndf_lepton_0 = func_lepton_fit_rot->GetNDF();

    if(a_fit_lepton_0 !=0 && a_fit_lepton_0 !=0){
      a_fit_lepton_final = -1/a_fit_lepton_0;  
      b_fit_lepton_final = a_fit_lepton_0/b_fit_kaon_0;   
      ChiS_lepton_final = ChiS_lepton_0 ;   
      ndf_lepton_final = ndf_lepton_0;
    }
    else{
      cout << "ERROR nan !" << endl;
      return;
    }


     line_fit_leptons = new TLine(50.*a_fit_lepton_0 + b_fit_lepton_0, -50., 
                               -50.*a_fit_lepton_0 + b_fit_lepton_0, 50.);

     gr_leptons_0_NEW = gr_leptons_unfitted; 
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  /// DETERMINE K_STOP_0
  /// KAON_PR
  vector<int> V_kaon_pr;         V_kaon_pr.clear();
  vector<double> V_x_kaon_pr;    V_x_kaon_pr.clear();
  vector<double> V_y_kaon_pr;    V_y_kaon_pr.clear();
  vector<double> V_ex_kaon_pr;   V_ex_kaon_pr.clear();
  vector<double> V_ey_kaon_pr;   V_ey_kaon_pr.clear();
  vector<double> vec_X_kaon_pruned_NEW;   vec_X_kaon_pruned_NEW.size();
  vector<double> vec_Y_kaon_pruned_NEW;   vec_Y_kaon_pruned_NEW.size();


  if(to_rotate_kaons || to_rotate_leptons){
    for(unsigned int i=0; i<V_kaon_bars.size(); i++){
      V_kaon_pr.push_back(TARGET_Rotated_index[V_kaon_bars[i]]);
      V_x_kaon_pr.push_back(Xloc[TARGET_Rotated_index[V_kaon_bars[i]]]);
      V_y_kaon_pr.push_back(Yloc[TARGET_Rotated_index[V_kaon_bars[i]]]);
      V_ex_kaon_pr.push_back(TARGET_Errors_X);
      V_ey_kaon_pr.push_back(TARGET_Errors_Y);
      vec_X_kaon_pruned_NEW.push_back(Xloc[V_kaon_bars[i]]);
      vec_Y_kaon_pruned_NEW.push_back(Yloc[V_kaon_bars[i]]);
    }
  }
  else{
    for(unsigned int i=0; i<V_kaon_bars.size(); i++){
      V_kaon_pr.push_back(V_kaon_bars[i]);     
      V_x_kaon_pr.push_back(Xloc[V_kaon_bars[i]]);
      V_y_kaon_pr.push_back(Yloc[V_kaon_bars[i]]);
      V_ex_kaon_pr.push_back(TARGET_Errors_X);
      V_ey_kaon_pr.push_back(TARGET_Errors_Y);
      vec_X_kaon_pruned_NEW.push_back(Xloc[V_kaon_bars[i]]);
      vec_Y_kaon_pruned_NEW.push_back(Yloc[V_kaon_bars[i]]);
    }   
  }

  // LEPTON_PR 
  vector<int> V_lepton_pr;         V_lepton_pr.clear();
  vector<double> V_x_lepton_pr;    V_x_lepton_pr.clear();
  vector<double> V_y_lepton_pr;    V_y_lepton_pr.clear();
  vector<double> V_ex_lepton_pr;   V_ex_lepton_pr.clear();
  vector<double> V_ey_lepton_pr;   V_ey_lepton_pr.clear();

  if(to_rotate_leptons){
    for(unsigned int i=0; i<V_lepton_0.size(); i++){
      V_lepton_pr.push_back(TARGET_Rotated_index[V_lepton_0[i]]);
      V_x_lepton_pr.push_back(Xloc[TARGET_Rotated_index[V_lepton_0[i]]]);
      V_y_lepton_pr.push_back(Yloc[TARGET_Rotated_index[V_lepton_0[i]]]);
      V_ex_lepton_pr.push_back(TARGET_Errors_X);
      V_ey_lepton_pr.push_back(TARGET_Errors_Y);
    }

    //Add TOF1
    for(int k=0; k<weight_TOF1_fit1; k++){
      V_x_lepton_pr.push_back(TOF1_Xloc[TOF1_rotated[gap_to_fit-1]][2]);
      V_ex_lepton_pr.push_back(TOF1_Errors_X[TOF1_rotated[gap_to_fit-1]][2]);
  
      V_y_lepton_pr.push_back(TOF1_Yloc[TOF1_rotated[gap_to_fit-1]][2]);
      V_ey_lepton_pr.push_back(TOF1_Errors_Y[TOF1_rotated[gap_to_fit-1]][2]);
    }
  }
  else{
    for(unsigned int i=0; i<V_lepton_0.size(); i++){
      V_lepton_pr.push_back(V_lepton_0[i]);     
      V_x_lepton_pr.push_back(Xloc[V_lepton_0[i]]);
      V_y_lepton_pr.push_back(Yloc[V_lepton_0[i]]);
      V_ex_lepton_pr.push_back(TARGET_Errors_X);
      V_ey_lepton_pr.push_back(TARGET_Errors_Y);
    } 
    
    //Add TOF1
    for(int k=0; k<weight_TOF1_fit1; k++){
      V_x_lepton_pr.push_back(TOF1_Xloc[gap_to_fit-1][2]);
      V_ex_lepton_pr.push_back(TOF1_Errors_X[gap_to_fit-1][2]);
  
      V_y_lepton_pr.push_back(TOF1_Yloc[gap_to_fit-1][2]);
      V_ey_lepton_pr.push_back(TOF1_Errors_Y[gap_to_fit-1][2]);
    }
  }

  TGraph *N_graph_kaon_pr;
  N_graph_kaon_pr = new TGraphErrors(V_x_kaon_pr.size(),&V_x_kaon_pr[0],&V_y_kaon_pr[0],
                                      &V_ex_kaon_pr[0],&V_ey_kaon_pr[0]);
  TGraph *N_graph_lepton_pr;
  N_graph_lepton_pr = new TGraphErrors(V_x_lepton_pr.size(),&V_x_lepton_pr[0],&V_y_lepton_pr[0],
                                      &V_ex_lepton_pr[0],&V_ey_lepton_pr[0]);


  ////////////////////////////////////////////////////////////////////////////////////////////

  TF1 *N_func_lepton_fit_pr;
  N_func_lepton_fit_pr = new TF1("N_gr_lepton_pr", "pol1");
  double N_a_fit_lepton_pr = 999.;
  double N_b_fit_lepton_pr = 999.;
  double N_ChiS_lepton_pr = 999.;
  double N_ndf_lepton_pr = 999.;

  if(V_lepton_pr.size()>lepton_cluster_size){
    //graph_lepton_pr->Fit("gr_lepton_pr","QW");  // <--- OPTION "W" TO BE REMOVED IF TOO MANY ERROR MESSAGES !
    N_graph_lepton_pr->Fit("N_gr_lepton_pr","Q");

    N_func_lepton_fit_pr = N_graph_lepton_pr->GetFunction("N_gr_lepton_pr");
    N_func_lepton_fit_pr->SetLineWidth(2);
    N_func_lepton_fit_pr->SetLineColor(2);
  
    N_a_fit_lepton_pr = N_func_lepton_fit_pr->GetParameter(1);
    N_b_fit_lepton_pr = N_func_lepton_fit_pr->GetParameter(0);
    N_ChiS_lepton_pr = N_func_lepton_fit_pr->GetChisquare();
    N_ndf_lepton_pr = N_func_lepton_fit_pr->GetNDF();
  }

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// DETERMINE KAON CENTROID

  vector<double> V_x_kaon_centroid;   V_x_kaon_centroid.clear();
  vector<double> V_y_kaon_centroid;   V_y_kaon_centroid.clear();


  // X_weights = 0.;
  // Y_weights = 0.;
  // total_energy = 0.;


  vec_kaon_centroid_coordinates.clear();
  //vec_fit_lines_intersect.clear();
  vec_k_stop_coordinates.clear();

  vector<double> vec_xx_kaon_centroid;      vec_xx_kaon_centroid.clear();
  vector<double> vec_yy_kaon_centroid;      vec_yy_kaon_centroid.clear();

  vector<double> V_centroid;
  V_centroid.clear();
  double N_centroid_X = 0.;
  double N_centroid_Y = 0.;

  V_centroid = _get_kaon_centroid(V_kaon_bars);
  N_centroid_X = V_centroid[0];
  N_centroid_Y = V_centroid[1];


  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// DETERMINE KSTOP 0

  vector<double> N_test;      N_test.clear();
  vector<double> N_Kstop_0;

  N_Kstop_0 = _get_kstop_0(V_lepton_bars, V_kaon_bars, V_centroid,
                          a_fit_kaon_final, b_fit_kaon_final, a_fit_lepton_final, b_fit_lepton_final);


  // N_x_Kstop_0.push_back(N_Kstop_0[0]);
  // N_y_Kstop_0.push_back(N_Kstop_0[1]);

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// TRIANGULAR PRUNING

  vector<int> TOF1_neighbours;
  TOF1_neighbours = _neighbour_TOF1(gap_to_fit);

  vector<double> vec_xx_pruning_area;   vec_xx_pruning_area.clear();
  vector<double> vec_yy_pruning_area;   vec_yy_pruning_area.clear();

  vec_xx_pruning_area = _pruning_area(N_Kstop_0[0], TOF1_Xloc[TOF1_neighbours[0]][2], TOF1_Xloc[TOF1_neighbours[1]][2]);
  vec_yy_pruning_area = _pruning_area(N_Kstop_0[1], TOF1_Yloc[TOF1_neighbours[0]][2], TOF1_Yloc[TOF1_neighbours[1]][2]);

  TGraph *gr_pruning_area;
  gr_pruning_area = new TGraph(vec_xx_pruning_area.size(), &vec_xx_pruning_area[0], &vec_yy_pruning_area[0]);
  gr_pruning_area->SetLineWidth(2);
  gr_pruning_area->SetLineColor(3);
  gr_pruning_area->SetMarkerStyle(34);
  gr_pruning_area->SetMarkerSize(1.3);
  gr_pruning_area->SetMarkerColor(2);

  // cc1->cd(1);
  // h_empty->Draw();
  // gr_pruning_area->Draw("same");

  vector<int> vec_pruned_triangle_NEW;
  vec_pruned_triangle_NEW.clear();
 
  vec_pruned_triangle_NEW = _Triangular_pruning2(gap_to_fit, N_Kstop_0[0], N_Kstop_0[1], 
                                        V_lepton_bars, V_kaon_bars);
  


  vec_X_kaon_pruned_NEW.clear();
  vec_Y_kaon_pruned_NEW.clear();


  for(unsigned int i=0; i<vec_pruned_triangle_NEW.size(); i++){
    vec_X_kaon_pruned_NEW.push_back(Xloc[vec_pruned_triangle_NEW[i]]);
    vec_Y_kaon_pruned_NEW.push_back(Yloc[vec_pruned_triangle_NEW[i]]);
  }
 
  TGraph *gr_kaon_pruned_NEW;
  gr_kaon_pruned_NEW = new TGraph(vec_pruned_triangle_NEW.size(), &vec_X_kaon_pruned_NEW[0], &vec_Y_kaon_pruned_NEW[0]);
  gr_kaon_pruned_NEW->SetMarkerStyle(21);
  gr_kaon_pruned_NEW->SetMarkerSize(0.8);
  gr_kaon_pruned_NEW->SetMarkerColor(3);


  //////////////////////////////////////////////////////////////////////////////////////////////
  /// LINE PRUNING

  vector<double> vec_pruning_line;  vec_pruning_line.clear();
  vector<int> vec_pruned_line_NEW;  vec_pruned_line_NEW.clear();

  vec_pruned_line_NEW = _Linear_pruning(to_rotate_leptons, N_Kstop_0[0], N_Kstop_0[1], gap_to_fit, 7.5,V_lepton_bars);

  ///////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  for(int g=0; g<12; g++){
    // Gap[g][0][0] = TOF_Xloc[3*g];
    // Gap[g][1][0] = TOF_Xloc[3*g+1];
    // Gap[g][2][0] = TOF_Xloc[3*g+2];

    // Gap[g][0][1] = TOF_Yloc[3*g];
    // Gap[g][1][1] = TOF_Yloc[3*g+1];
    // Gap[g][2][1] = TOF_Yloc[3*g+2];

    Gap[g][0][0] = TOF1_Xloc[g][0];
    Gap[g][1][0] = TOF1_Xloc[g][1];
    Gap[g][2][0] = TOF1_Xloc[g][2];
    Gap[g][3][0] = TOF1_Xloc[g][3];
    Gap[g][4][0] = TOF1_Xloc[g][4];

    Gap[g][0][1] = TOF1_Yloc[g][0];
    Gap[g][1][1] = TOF1_Yloc[g][1];
    Gap[g][2][1] = TOF1_Yloc[g][2];
    Gap[g][3][1] = TOF1_Yloc[g][3];
    Gap[g][4][1] = TOF1_Yloc[g][4];
   }


  ///////////////////////////////////////////////////////////////////////////////////////
  /// WORKING WITH PRUNED DATA

  vector<int> vec_leptons_pruned_NEW;       vec_leptons_pruned_NEW.clear();
  vector<double> vec_x_leptons_pruned_NEW;  vec_x_leptons_pruned_NEW.clear();
  vector<double> vec_y_leptons_pruned_NEW;  vec_y_leptons_pruned_NEW.clear();

  if(flag_pruning_line){
    vec_leptons_pruned_NEW = vec_pruned_line_NEW;
  }
  else{
    vec_leptons_pruned_NEW = vec_pruned_triangle_NEW;
  }

  if(to_rotate_leptons){
    for(unsigned int i=0; i<vec_leptons_pruned_NEW.size(); i++){
      vec_x_leptons_pruned_NEW.push_back(Xloc[TARGET_Rotated_index[vec_leptons_pruned_NEW[i]]]);
      vec_y_leptons_pruned_NEW.push_back(Yloc[TARGET_Rotated_index[vec_leptons_pruned_NEW[i]]]);
    }
  }
  else{
    for(unsigned int i=0; i<vec_leptons_pruned_NEW.size(); i++){
      vec_x_leptons_pruned_NEW.push_back(Xloc[vec_leptons_pruned_NEW[i]]);
      vec_y_leptons_pruned_NEW.push_back(Yloc[vec_leptons_pruned_NEW[i]]);
    }
  }

  TGraph *gr_test_pruned;
  gr_test_pruned = new TGraph(vec_leptons_pruned_NEW.size(),&vec_x_leptons_pruned_NEW[0],&vec_y_leptons_pruned_NEW[0]);

  gr_test_pruned->SetMarkerColor(2);
  gr_test_pruned->SetMarkerStyle(21);
  gr_test_pruned->SetMarkerSize(0.8);

  vector<double> vec_line_parameters;
  vec_line_parameters.clear();


  //////////////////////////////////////////////////////////////////////////////////////////////
  /// LEPTONS FIT 1

  if(to_rotate_leptons) TOF1_gap = TOF1_rotated[gap_to_fit-1]+1;
  else TOF1_gap = gap_to_fit;

  int EB_addition_weight = 0;

  vector<int> vec_leptons_fit1;
  vec_leptons_fit1.clear();

  /// GET SELECTED LEPTONS FOR FIT 1
  /// CAREFUL !!! : These lepon bars can be rotated !
  vec_leptons_fit1 = _select_leptons_to_fit(to_rotate_leptons, gap_to_fit, vec_leptons_pruned_NEW,
                                            vec_line_parameters, EB_addition_weight);

  bool to_restore = false;
  if(to_rotate_leptons) to_restore = true;

  /// PREPARE TGRAPHERRORS gr_leptons_fit1 
  vector<double> vec_x_leptons_fit1;   vec_x_leptons_fit1.clear();
  vector<double> vec_y_leptons_fit1;   vec_y_leptons_fit1.clear();
  vector<double> vec_ex_leptons_fit1;  vec_ex_leptons_fit1.clear();
  vector<double> vec_ey_leptons_fit1;  vec_ey_leptons_fit1.clear();


  //*** Fill in the lepton bars (including "edge bars")
  for(unsigned int i=0; i<vec_leptons_fit1.size(); i++){
    vec_x_leptons_fit1.push_back(Xloc[vec_leptons_fit1[i]]);
    vec_y_leptons_fit1.push_back(Yloc[vec_leptons_fit1[i]]);
    vec_ex_leptons_fit1.push_back(TARGET_Errors_X);
    vec_ey_leptons_fit1.push_back(TARGET_Errors_Y);
  }

  //*** Add TOF1 gap (weighted)
  /// ADD TOF1 GAP (WEIGHTED) TO THE FIT
  for(int i=0; i<weight_TOF1_fit1; i++){
    vec_x_leptons_fit1.push_back(TOF1_Xloc[TOF1_gap-1][2]);
    vec_y_leptons_fit1.push_back(TOF1_Yloc[TOF1_gap-1][2]);
    vec_ex_leptons_fit1.push_back(TOF1_Errors_X[TOF1_gap-1][2]);
    vec_ey_leptons_fit1.push_back(TOF1_Errors_Y[TOF1_gap-1][2]);
  }

  TGraph *gr_leptons_fit1;
  gr_leptons_fit1 = new TGraphErrors(vec_x_leptons_fit1.size(), &vec_x_leptons_fit1[0], &vec_y_leptons_fit1[0],
                                    &vec_ex_leptons_fit1[0], &vec_ey_leptons_fit1[0]);

  gr_leptons_fit1->SetMarkerStyle(21);
  gr_leptons_fit1->SetMarkerColor(2);
  gr_leptons_fit1->SetMarkerSize(0.8);
  gr_leptons_fit1->GetXaxis()->SetLimits(-50.,50.);
  gr_leptons_fit1->GetYaxis()->SetRangeUser(-50.,50.);

  double a_fit1;
  double b_fit1;
  double Chis_fit1;
  double ndf_fit1;

  TF1 *f_leptons_fit_1 = new TF1("Leptons_fit1", "pol1");
  //gr_leptons_fit1->Fit("Leptons_fit1","QW");
  gr_leptons_fit1->Fit("Leptons_fit1","Q");
  delete func_lepton_fit_1;
  
  TF1 *f_leptons_fit1 = gr_leptons_fit1->GetFunction("Leptons_fit1");
  a_fit1 = f_leptons_fit1->GetParameter(1);
  b_fit1 = f_leptons_fit1->GetParameter(0);
  Chis_fit1 = f_leptons_fit1->GetChisquare();
  ndf_fit1 = f_leptons_fit1->GetNDF();
  f_leptons_fit1->SetLineWidth(2);
  f_leptons_fit1->SetLineColor(2);

  /// CHECK IF A ROTATION IS NEEDED 
  // *** Or simply use the double fit method with first fit using option QW
  bool to_rotate_leptons_post_pruning = false;
  if(abs(a_fit1)>15) to_rotate_leptons_post_pruning = true;

  if(!to_rotate_leptons && to_rotate_leptons_post_pruning){
    
    to_restore = true;

    vec_leptons_fit1 = _rotate_vector(vec_leptons_fit1);

    vec_x_leptons_fit1.clear();
    vec_y_leptons_fit1.clear();
    vec_ex_leptons_fit1.clear();
    vec_ey_leptons_fit1.clear();

    for(unsigned int i=0; i<vec_leptons_fit1.size(); i++){
      vec_x_leptons_fit1.push_back(Xloc[vec_leptons_fit1[i]]);
      vec_y_leptons_fit1.push_back(Yloc[vec_leptons_fit1[i]]);
      vec_ex_leptons_fit1.push_back(TARGET_Errors_X);
      vec_ey_leptons_fit1.push_back(TARGET_Errors_Y);
    }

    TOF1_gap = TOF1_rotated[gap_to_fit-1]+1;

    for(int i=0; i<weight_TOF1_fit1; i++){
      vec_x_leptons_fit1.push_back(TOF1_Xloc[TOF1_gap-1][2]);
      vec_y_leptons_fit1.push_back(TOF1_Yloc[TOF1_gap-1][2]);
      vec_ex_leptons_fit1.push_back(TOF1_Errors_X[TOF1_gap-1][2]);
      vec_ey_leptons_fit1.push_back(TOF1_Errors_Y[TOF1_gap-1][2]);
    }

    gr_leptons_fit1->Set(0);
    gr_leptons_fit1 = new TGraphErrors(vec_x_leptons_fit1.size(),&vec_x_leptons_fit1[0],&vec_y_leptons_fit1[0],
                                   &vec_ex_leptons_fit1[0],&vec_ey_leptons_fit1[0]);

    gr_leptons_fit1->SetMarkerStyle(21);
    gr_leptons_fit1->SetMarkerColor(2);
    gr_leptons_fit1->SetMarkerSize(0.8);
    gr_leptons_fit1->GetXaxis()->SetLimits(-50.,50.);
    gr_leptons_fit1->GetYaxis()->SetRangeUser(-50.,50.);
  
    TF1 *f_leptons_fit_1 = new TF1("Leptons_fit1", "pol1");
    gr_leptons_fit1->Fit("Leptons_fit1","Q");
    delete func_lepton_fit_1;
  
    TF1 *f_leptons_fit1 = gr_leptons_fit1->GetFunction("Leptons_fit1");
    a_fit1 = f_leptons_fit1->GetParameter(1);
    b_fit1 = f_leptons_fit1->GetParameter(0);
    Chis_fit1 = f_leptons_fit1->GetChisquare();
    ndf_fit1 = f_leptons_fit1->GetNDF();
    f_leptons_fit1->SetLineWidth(2);
    f_leptons_fit1->SetLineColor(2);
  }

  vector<double> vec_line_parameters_fit1;  vec_line_parameters_fit1.clear();
  vec_line_parameters_fit1.push_back(a_fit1);
  vec_line_parameters_fit1.push_back(b_fit1);

  bool To_rotate_TEST = false;
  if(to_rotate_leptons || to_rotate_leptons_post_pruning) To_rotate_TEST = true;


  ////////////////////////////////////////////////////////////////////////////////////////////
  /// LEPTONS FIT 2
  /// GET SELECTED LEPTONS FOR FIT 2
  /// CAREFUL !!! : These lepon bars can be rotated !
  vector<int> vec_leptons_fit2;   vec_leptons_fit2.clear();
  EB_addition_weight = 1; 

  vector<double> vec_x_leptons_fit2;    vec_x_leptons_fit2.clear();
  vector<double> vec_y_leptons_fit2;    vec_y_leptons_fit2.clear();
  vector<double> vec_ex_leptons_fit2;   vec_ex_leptons_fit2.clear();
  vector<double> vec_ey_leptons_fit2;   vec_ey_leptons_fit2.clear();

  vec_leptons_fit2 = _select_leptons_to_fit(false, TOF1_gap, vec_leptons_fit1,
                                            vec_line_parameters_fit1, EB_addition_weight);



  for(unsigned int i=0; i<vec_leptons_fit2.size(); i++){
    vec_x_leptons_fit2.push_back(Xloc[vec_leptons_fit2[i]]);
    vec_y_leptons_fit2.push_back(Yloc[vec_leptons_fit2[i]]);
    vec_ex_leptons_fit2.push_back(TARGET_Errors_X);
    vec_ey_leptons_fit2.push_back(TARGET_Errors_Y);
  }

  TGraph *gr_leptons_fit2;
  gr_leptons_fit2 = new TGraphErrors(vec_leptons_fit2.size(), &vec_x_leptons_fit2[0], &vec_y_leptons_fit2[0], 
                                    &vec_ex_leptons_fit2[0], &vec_ey_leptons_fit2[0]);

  double a_fit2 = 999.99;
  double b_fit2 = 999.99;
  double Chis_fit2 = 999.99;
  double ndf_fit2 = 999.99;
  vector<int> vec_leptons_pruned_unique;
  vec_leptons_pruned_unique.clear();
  int Lepton_counter = -1; 

  gr_leptons_fit2->SetMarkerStyle(21);
  gr_leptons_fit2->SetMarkerColor(2);
  gr_leptons_fit2->SetMarkerSize(0.8);
  gr_leptons_fit2->GetXaxis()->SetLimits(-50.,50.);
  gr_leptons_fit2->GetYaxis()->SetRangeUser(-50.,50.);

  TF1 *f_leptons_fit_2 = new TF1("Leptons_fit2", "pol1");
  if(vec_leptons_fit2.size()>0){
    gr_leptons_fit2->Fit("Leptons_fit2","QW");
    gr_leptons_fit2->Fit("Leptons_fit2","Q");
    delete func_lepton_fit_2;
  
    TF1 *f_leptons_fit2 = gr_leptons_fit2->GetFunction("Leptons_fit2");
    a_fit2 = f_leptons_fit2->GetParameter(1);
    b_fit2 = f_leptons_fit2->GetParameter(0);
    Chis_fit2 = f_leptons_fit2->GetChisquare();
    ndf_fit2 = f_leptons_fit2->GetNDF();
    f_leptons_fit2->SetLineWidth(2);
    f_leptons_fit2->SetLineColor(2);
  }

  vec_leptons_pruned_unique = _remove_doubles(vec_leptons_fit2);
  Lepton_counter = vec_leptons_pruned_unique.size();



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// ???
  for(int k=0; k<256; k++){
    for(int l=0; l<6; l++){
        TDC_LE_TARGET_corrected[k][l] = 9999;
        //if(tdc_le_target[k][l]>-1)  TDC_LE_TARGET_corrected[k][l] = tdc_le_target[k][l] - TDC_average;
        if(tdc_le_target[k][l]>-1)  TDC_LE_TARGET_corrected[k][l] = tdc_le_target[k][l] - TDC_average_NEW;
    }
  }

  for(int k=0; k<256; k++){
    for(int l=0; l<6; l++){
        if(TDC_LE_TARGET_corrected[k][l] == 9999) sprintf(TDC_LE_TARGET_corr[k][l],"----");
        if(TDC_LE_TARGET_corrected[k][l] != 9999) sprintf(TDC_LE_TARGET_corr[k][l],"%4d",TDC_LE_TARGET_corrected[k][l]);
    }
  }

  //// Don't do anything if the event has less than n_hits hits in the TARGET
  if (enable_cout == 9) count = n_hit+1;
  if(count<-10){
    if (Switch_Printout != 0) {
      cout << " >>>>  Event "<< ievt << " has fewer than "<< n_hit << " high gain hits in the TARGET (outliers and K stop removed)!" << endl;
      cout << " >>>>  Please, choose another event" << endl;
      cout << " " << endl;
    }
    return;
  }

  // Arrays for finding mwpc clustering

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

  C2X_centroid = 0.;

  //C2 Counters
  for(int i=0; i<56; i++){
    if (Good_Event && ADC_C2X_L[i]>0) C2X_L[i] = ADC_C2X_L[i];
    if (Good_Event && ADC_C2X_R[i]>0) C2X_R[i] = ADC_C2X_R[i];
  }

  for(int i=0; i<64; i++){
    if (Good_Event && ADC_C3X_L[i]>0) C3X_L[i] = ADC_C3X_L[i];
    if (Good_Event && ADC_C3X_R[i]>0) C3X_R[i] = ADC_C3X_R[i];
  }

  for(int i=0; i<72; i++){
    if (Good_Event && ADC_C4X_L[i]>0) C4X_L[i] = ADC_C4X_L[i];
    if (Good_Event && ADC_C4X_R[i]>0) C4X_R[i] = ADC_C4X_R[i];
  }

  for(int i=0; i<16; i++){
    if (Good_Event && ADC_C2Y_L[i]>0) C2Y_L[i] = ADC_C2Y_L[i];
    if (Good_Event && ADC_C2Y_R[i]>0) C2Y_R[i] = ADC_C2Y_R[i];
    if (Good_Event && ADC_C3Y_L[i]>0) C3Y_L[i] = ADC_C3Y_L[i];
    if (Good_Event && ADC_C3Y_R[i]>0) C3Y_R[i] = ADC_C3Y_R[i];
    if (Good_Event && ADC_C4Y_L[i]>0) C4Y_L[i] = ADC_C4Y_L[i];
    if (Good_Event && ADC_C4Y_R[i]>0) C4Y_R[i] = ADC_C4Y_R[i];
  }

  // Find clustering of mwpcs
  ////////////////////////////
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
        else if ( i != 0) cluster_length_count++;
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
        else if ( i != 0) cluster_length_count++;
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
        else if ( i != 0) cluster_length_count++;
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

  // B0 Counter
  for(int i=0; i<16; i++){
    if(tdc_vt48[6][i]>TDC_B0_min && tdc_vt48[6][i]<TDC_B0_max) vec_tdc_b0_6.push_back(tdc_vt48[6][i]);
    else vec_tdc_b0_6.push_back(-1);

    if(tdc_vt48[7][i]>TDC_B0_min && tdc_vt48[7][i]<TDC_B0_max) vec_tdc_b0_7.push_back(tdc_vt48[7][i]);
    else vec_tdc_b0_7.push_back(-1);
  }

  sort(vec_tdc_b0_6.begin(), vec_tdc_b0_6.end());
  sort(vec_tdc_b0_7.begin(), vec_tdc_b0_7.end());


  // ### Display Ck and Cpi infomation

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

  if(vec_Ck.size()>1) TDC_ck_sigma = sqrt(TDC_ck_sigma2/(vec_Ck.size()-1));
  else TDC_ck_sigma = 0;

  if(vec_Cpi.size()>1) TDC_cpi_sigma = sqrt(TDC_cpi_sigma2/(vec_Cpi.size()-1));
  else TDC_cpi_sigma = 0;

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

  n_bar_selected=999;
  n_bar_selected = vec_TARGET_bar_selected.size();
  
  //char TDC_LE_TARGET_sel_string[n_bar_selected][16][20];
  char TDC_LE_TARGET_sel_string[256][16][20];

  for(unsigned int i=0; i<vec_TARGET_bar_selected.size(); i++){
    for(int j=0; j<16; j++){
      if(tdc_le_target[vec_TARGET_bar_selected[i]][j] > 0){
        sprintf(TDC_LE_TARGET_sel_string[i][j],"%5.1f",550 - (0.91*tdc_le_target[vec_TARGET_bar_selected[i]][j] - 0.625*(tdc_vt48[0][0]-TDC_ck_avg2)));
      }
      else sprintf(TDC_LE_TARGET_sel_string[i][j],"-----");
    }
  }

  par0_ADC = 0;  par0_TDC = 0.;
  par1_ADC = 0;  par1_TDC = 0.;

  X_BAR = 999.;
  Y_BAR = 999.;

  /// Data counters
  has_data_TDC2 = 0;
  has_data_ADC2 = 0;
  has_data_ADC3 = 0;
  has_data_ADC4 = 0;
  has_data_ADCA = 0;

  ///// Select edge fiber for track fitting
  Xloc_gap = 0;
  Yloc_gap = 0;

  if (gap_to_fit == 0){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[0][j]]>0) && (has_TDC_hit[channel[0][j]])){
        Xloc_gap = Xloc[channel[0][j]];
        Yloc_gap = Yloc[channel[0][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[0][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[0][3]];
  }
  if (gap_to_fit == 1){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[1][j]]>0) && (has_TDC_hit[channel[1][j]])){
        Xloc_gap = Xloc[channel[1][j]];
        Yloc_gap = Yloc[channel[1][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[1][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[1][3]];
  }
  if (gap_to_fit == 2){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[2][j]]>0) && (has_TDC_hit[channel[2][j]])){
        Xloc_gap = Xloc[channel[2][j]];
        Yloc_gap = Yloc[channel[2][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[2][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[2][3]];
  }
  if (gap_to_fit == 3){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[3][j]]>0) && (has_TDC_hit[channel[3][j]])){
        Xloc_gap = Xloc[channel[3][j]];
        Yloc_gap = Yloc[channel[3][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[3][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[3][3]];
  }
  if (gap_to_fit == 4){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[4][j]]>0) && (has_TDC_hit[channel[4][j]])){
        Xloc_gap = Xloc[channel[4][j]];
        Yloc_gap = Yloc[channel[4][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[4][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[4][3]];
  }
  if (gap_to_fit == 5){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[5][j]]>0) && (has_TDC_hit[channel[5][j]])){
        Xloc_gap = Xloc[channel[5][j]];
        Yloc_gap = Yloc[channel[5][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[5][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[5][3]];
  }
  if (gap_to_fit == 6){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[6][j]]>0) && (has_TDC_hit[channel[6][j]])){
        Xloc_gap = Xloc[channel[6][j]];
        Yloc_gap = Yloc[channel[6][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[6][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[6][3]];
  }
  if (gap_to_fit == 7){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[7][j]]>0) && (has_TDC_hit[channel[7][j]])){
        Xloc_gap = Xloc[channel[7][j]];
        Yloc_gap = Yloc[channel[7][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[7][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[7][3]];
  }
  if (gap_to_fit == 8){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[8][j]]>0) && (has_TDC_hit[channel[8][j]])){
        Xloc_gap = Xloc[channel[8][j]];
        Yloc_gap = Yloc[channel[8][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[8][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[8][3]];
  }
  if (gap_to_fit == 9){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[9][j]]>0) && (has_TDC_hit[channel[9][j]])){
        Xloc_gap = Xloc[channel[9][j]];
        Yloc_gap = Yloc[channel[9][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[9][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[9][3]];
  }
  if (gap_to_fit == 10){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[10][j]]>0) && (has_TDC_hit[channel[10][j]])){
        Xloc_gap = Xloc[channel[10][j]];
        Yloc_gap = Yloc[channel[10][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[10][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[10][3]];
  }
  if (gap_to_fit == 11){
    for (int j=0; j<12; j++){
      if ((ADC_High_TARGET[channel[11][j]]>0) && (has_TDC_hit[channel[11][j]])){
        Xloc_gap = Xloc[channel[11][j]];
        Yloc_gap = Yloc[channel[11][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[11][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[11][3]];
  }

  ////Fill TARGET tracking histograms with data
  for(Int_t i=0; i<256; i++){
    if ((i == max_index) || (i == max_index2) || (i == max_index3) || (i == max_index4)) continue;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (hyp[i] > 12.5)) has_data_ADC2++;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) has_data_TDC2++;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) has_data_ADC3++;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (hyp[i] < 12.5)) has_data_ADC4++;
  }

  if (has_data_ADC3 > 1){

    xcoord = 0.;
    unique_x = 0;

    if (has_data_ADCA > 2){

      /// Ring intercept coordinates
      determinant = 4*(pow(par0_ADC,2))*(pow(par1_ADC,2)) - 4*(pow(par1_ADC,2) + 1)*(pow(par0_ADC,2)-1600);
      x_circle_int1 = (-2*(par0_ADC)*(par1_ADC) + sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
      y_circle_int1 = (par1_ADC)*(x_circle_int1) + par0_ADC;
      x_circle_int2 = (-2*(par0_ADC)*(par1_ADC) - sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
      y_circle_int2 = (par1_ADC)*(x_circle_int2) + par0_ADC;

      if (unique_x == 1){
        x_circle_int1 = xcoord;
        y_circle_int1 = sqrt(1600-(pow(xcoord,2)));
        x_circle_int2 = xcoord;
        y_circle_int2 = sqrt(1600-(pow(xcoord,2)))*-1;
      }

      SFTxdistance1 = pow((x_circle_int1-Xloc_gap),2);
      SFTydistance1 = pow((y_circle_int1-Yloc_gap),2);

      SFTxdistance2 = pow((x_circle_int2-Xloc_gap),2);
      SFTydistance2 = pow((y_circle_int2-Yloc_gap),2);

      SFTxhyp1 = double(sqrt(double(SFTxdistance1) + double(SFTydistance1)));
      SFTxhyp2 = double(sqrt(double(SFTxdistance2) + double(SFTydistance2)));

      SFT_x_intercept = 999.;
      SFT_y_intercept = 999.;

      if (SFTxhyp1 < SFTxhyp2){
        SFT_x_intercept = x_circle_int1;
        SFT_y_intercept = y_circle_int1;
      }
      else {
        SFT_x_intercept = x_circle_int2;
        SFT_y_intercept = y_circle_int2;
      }

      cout << "" << endl;
      cout << "SFT Circle Intercept" << endl;
      cout << "X coordinate: " << SFT_x_intercept << " -- Y coordinate: " << SFT_y_intercept << endl;

      /// Determine angle phi between centre of TARGET and SFT ring intercept
      SFT_phi = 999.;

      if (SFT_x_intercept >0) SFT_phi = (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;
      else SFT_phi = 360 + (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;

      cout << "" << endl;
      cout << "Angle between SFT Layer Ring Intercept and Centre of TARGET: " << SFT_phi << " deg." << endl;
      cout << "" << endl;

      if (has_data_ADC4 > 1){

        x_intercept = double(double(par0_ADC - par0_TDC)/double(par1_TDC - par1_ADC));
        y_intercept = double(par1_ADC)*double(x_intercept) + par0_ADC;

        for(int i=0; i<256; i++){
          x_distance4[i] = 0.;
          y_distance4[i] = 0.;
          distances4[i] = 0;
        }

        for (int q=0; q<256; q++) {
          x_distance4[q] = pow((Xloc[q]-x_intercept),2);
          y_distance4[q] = pow((Yloc[q]-y_intercept),2);
          distances4[q] = double(sqrt(double(x_distance4[q]) + double(y_distance4[q])));
        }

        //double min_distance = 10000.0;
        min_distance = 10000.0;

        for (int q=0; q<256; q++) {
          x_distance4[q] = pow((Xloc[q]-x_intercept),2);
          y_distance4[q] = pow((Yloc[q]-y_intercept),2);
          distances4[q] = double(sqrt(double(x_distance4[q]) + double(y_distance4[q])));
        }

        for (int q=0; q<256; q++) {
          if (distances4[q] < min_distance) {
            min_distance = distances4[q];
          }
        }
      }
    }
  }

  a_fit_GoodLG=0.;                     b_fit_GoodLG=0.;
  a_fit_GoodLG_weighted=0.;            b_fit_GoodLG_weighted=0.;
  a_fit_TDC_selected_weighted=0.;      b_fit_TDC_selected_weighted=0.;

  for(int i=0; i<12; i++)
  {
    vec_xx_TOF1_Marker.push_back(Gap[i][0][0]);  // TO MOVE
    vec_xx_TOF1_Marker.push_back(Gap[i][1][0]);
    vec_xx_TOF1_Marker.push_back(Gap[i][2][0]);

    vec_yy_TOF1_Marker.push_back(Gap[i][0][1]);
    vec_yy_TOF1_Marker.push_back(Gap[i][1][1]);
    vec_yy_TOF1_Marker.push_back(Gap[i][2][1]);
  }

  for(int g=0; g<5; g++){
    vec_xx_TOF1.push_back(Gap[gap_to_fit-1][g][0]);  // TO MOVE
    vec_yy_TOF1.push_back(Gap[gap_to_fit-1][g][1]);
  }

  vec_xx_lepton.clear();            vec_ex_lepton.clear();
  vec_yy_lepton.clear();            vec_ey_lepton.clear();
  vec_xx_lepton_rotate.clear();     vec_ex_lepton_rotate.clear();
  vec_yy_lepton_rotate.clear();     vec_ey_lepton_rotate.clear();
  //vec_lepton_size.clear();
  lepton_counter = 0;


  //////////////////////////////////////////////////////////////////////////////////////////////
  /// LEPTONS FIT 3

  vector<double> vec_intersect_TOF1;  vec_intersect_TOF1.clear();
  int selected_TOF1_section = -1;

  if(vec_leptons_fit2.size()>0){
    vec_intersect_TOF1 = _intersect_TOF1(TOF1_gap-1, a_fit2, b_fit2);
  }

  vector<double> int_with_circle;
  int_with_circle.clear();

  if(vec_leptons_fit2.size()>0){
    int_with_circle = _good_intersect_with_circle(a_fit2, b_fit2, R_TOF1, TOF1_gap);
  }

  if(vec_leptons_fit2.size()>0){
    selected_TOF1_section = _which_TOF1_section(TOF1_gap, int_with_circle);
  }
  else selected_TOF1_section = -1;

  /// CAREFUL !!! : These lepon bars can be rotated !
  vector<int> vec_leptons_fit3;   vec_leptons_fit3.clear();
  EB_addition_weight = 1; 

  vector<double> vec_x_leptons_fit3;    vec_x_leptons_fit3.clear();
  vector<double> vec_y_leptons_fit3;    vec_y_leptons_fit3.clear();
  vector<double> vec_ex_leptons_fit3;   vec_ex_leptons_fit3.clear();
  vector<double> vec_ey_leptons_fit3;   vec_ey_leptons_fit3.clear();

  vec_x_leptons_fit3 = vec_x_leptons_fit2;
  vec_y_leptons_fit3 = vec_y_leptons_fit2;
  vec_ex_leptons_fit3 = vec_ex_leptons_fit2;
  vec_ey_leptons_fit3 = vec_ey_leptons_fit2;

  if(selected_TOF1_section != -1){
    for(int i=0; i<weight_TOF1_fit3; i++){
      vec_x_leptons_fit3.push_back(TOF1_Xloc[TOF1_gap-1][selected_TOF1_section]);
      vec_y_leptons_fit3.push_back(TOF1_Yloc[TOF1_gap-1][selected_TOF1_section]);
      vec_ex_leptons_fit3.push_back(TOF1_Errors_X[TOF1_gap-1][selected_TOF1_section]);
      vec_ey_leptons_fit3.push_back(TOF1_Errors_Y[TOF1_gap-1][selected_TOF1_section]);
    }
  }

  TGraph *gr_leptons_fit3;
  gr_leptons_fit3 = new TGraphErrors(vec_x_leptons_fit3.size(), &vec_x_leptons_fit3[0], &vec_y_leptons_fit3[0], 
                                    &vec_ex_leptons_fit3[0], &vec_ey_leptons_fit3[0]);

  double a_fit3 = 999.99;
  double b_fit3 = 999.99;
  double Chis_fit3 = 999.99;
  double ndf_fit3 = 999.99;
  double a_fit3_error = 999.99;
  double b_fit3_error = 999.99;


  gr_leptons_fit3->SetMarkerStyle(21);
  gr_leptons_fit3->SetMarkerColor(2);
  gr_leptons_fit3->SetMarkerSize(0.8);
  gr_leptons_fit3->GetXaxis()->SetLimits(-50.,50.);
  gr_leptons_fit3->GetYaxis()->SetRangeUser(-50.,50.);

  TF1 *f_leptons_fit_3 = new TF1("Leptons_fit3", "pol1");

  if(vec_x_leptons_fit3.size()>0){
    //gr_leptons_fit3->Fit("Leptons_fit3","QW");
    gr_leptons_fit3->Fit("Leptons_fit3","Q");
    delete func_lepton_fit_3;
  
    TF1 *f_leptons_fit3 = gr_leptons_fit3->GetFunction("Leptons_fit3");
    a_fit3 = f_leptons_fit3->GetParameter(1);
    b_fit3 = f_leptons_fit3->GetParameter(0);
    Chis_fit3 = f_leptons_fit3->GetChisquare();
    ndf_fit3 = f_leptons_fit3->GetNDF();
    a_fit3_error = f_leptons_fit3->GetParError(1);
    b_fit3_error = f_leptons_fit3->GetParError(0);
    f_leptons_fit3->SetLineWidth(2);
    f_leptons_fit3->SetLineColor(2);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////
  //// LEPTON GRAPH FINAL !

  vector<double> vec_x_leptons_final;    vec_x_leptons_final.clear();
  vector<double> vec_y_leptons_final;    vec_y_leptons_final.clear();
  vector<double> vec_ex_leptons_final;   vec_ex_leptons_final.clear();
  vector<double> vec_ey_leptons_final;   vec_ey_leptons_final.clear();

  TGraph *gr_leptons_final;

  if(to_restore){
    for(unsigned int i=0; i<vec_leptons_fit2.size(); i++){
      vec_x_leptons_final.push_back(Xloc[TARGET_Rotated_index_inverse[vec_leptons_fit2[i]]]);
      vec_y_leptons_final.push_back(Yloc[TARGET_Rotated_index_inverse[vec_leptons_fit2[i]]]);
      vec_ex_leptons_final.push_back(TARGET_Errors_X);
      vec_ey_leptons_final.push_back(TARGET_Errors_Y);
    }

    for(int i=0; i<weight_TOF1_fit3; i++){
      vec_x_leptons_final.push_back(TOF1_Xloc[gap_to_fit-1][selected_TOF1_section]);
      vec_y_leptons_final.push_back(TOF1_Yloc[gap_to_fit-1][selected_TOF1_section]);
      vec_ex_leptons_final.push_back(TOF1_Errors_X[gap_to_fit-1][selected_TOF1_section]);
      vec_ey_leptons_final.push_back(TOF1_Errors_Y[gap_to_fit-1][selected_TOF1_section]);
    }   

    gr_leptons_final = new TGraphErrors(vec_x_leptons_final.size(), &vec_x_leptons_final[0], &vec_y_leptons_final[0],
                                      &vec_ex_leptons_final[0], &vec_ey_leptons_final[0]);
  }
  else gr_leptons_final = gr_leptons_fit3;

  vector<double> final_fit_XY;
  TLine *restored_fitl_line;
  if(to_restore){
    final_fit_XY = _restore_fit_line(a_fit3, b_fit3);

    restored_fitl_line = new TLine(final_fit_XY[0], final_fit_XY[1],final_fit_XY[2],final_fit_XY[3]);

    restored_fitl_line->SetLineWidth(2);
    restored_fitl_line->SetLineColor(2);
  }



  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// DETERMINE FINAL KSTOP  <---- TO FINISH !!!

  vector<double> vec_kstop_FINAL; vec_kstop_FINAL.clear();
  int i_kaon_bar = -1;
  double R_Kstop = 999.99;
  
  //V_kaon_bars
  vec_kstop_FINAL = _get_kstop_final(to_restore,
                                     vec_leptons_fit2, V_kaon_bars, V_centroid,
                                     a_fit_kaon_final, b_fit_kaon_final, a_fit3, b_fit3,
                                     has_kaon_sub, kaon_substitute);

  i_kaon_bar = kaon_fiber(vec_kstop_FINAL[0],vec_kstop_FINAL[1]);

  R_Kstop = sqrt(vec_kstop_FINAL[0]*vec_kstop_FINAL[0] + vec_kstop_FINAL[1]*vec_kstop_FINAL[1]);

 ////////////////////////////////////////////////////////////////////////////////////////////////
/*  /// DRAW FITTING CANVAS

  _Drawing_fitting(Run_Number, ievt, To_rotate_TEST,
                   gr_leptons_fit1, gr_leptons_fit2, gr_leptons_fit3,
                   gr_kaons_final, line_fit_kaons,
                   int_with_circle);
*/
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// IF NO LEPTON LEFT 
  bool b_bbar = false;
  vector<double> vec_bbar;          vec_bbar.clear();
  vector<double> vec_slope_bbar;    vec_slope_bbar.clear();

  if(vec_leptons_fit2.size()<2){
  //if(vec_leptons_fit2.size()<2 || (vec_leptons_fit2.size()==2 && angle_final_guide==999.)){
      vec_bbar = lineLeptonToTOF1(vec_leptons_fit2, gap_to_fit, vec_kstop_FINAL[0], vec_kstop_FINAL[1]);
      b_bbar = true;
    }

    if(b_bbar){
      vec_slope_bbar = get_slope(vec_bbar[0], vec_bbar[1], vec_bbar[2], vec_bbar[3]);    
    }


  ////////////////////////////////////////////////////////////////////////////////////////////////


  double a_final = 999.99;
  double b_final = 999.99;
  double ChiS_final = 999.99;
  double ndf_final = 999.99;
  double a_final_error = 999.99;
  double b_final_error = 999.99;

  if(to_restore){
    if(a_fit3 !=0){
      a_final = -1/a_fit3;  
      b_final = b_fit3/a_fit3;   
      ChiS_final = Chis_fit3 ;   
      ndf_final = ndf_fit3; 
      a_final_error = a_fit3_error;
      b_final_error = b_fit3_error;
    }  
    else cout << "NAN ALERT ! -- Evt : " << ievt << endl;
  }
  else{
    a_final = a_fit3;  
    b_final = b_fit3;   
    ChiS_final = Chis_fit3 ;     
    ndf_final = ndf_fit3;
    a_final_error = a_fit3_error;
    b_final_error = b_fit3_error;
  }

  //////////////////////////////////////////////////////////////////////////////
  /// GET THE INTERSECT BETWEEN THE TWO FITTING LINES
  vector<double >vec_fit_lines_intersect;   vec_fit_lines_intersect.clear();
  vec_fit_lines_intersect = _2lines_intersect(a_fit_kaon_final, b_fit_kaon_final, a_fit3, b_fit3);

  /// CALCULATE INTERSECT WITH TARGET, TOF1 and SFT
  // Intersect with TARGET
  vector<double> vec_intersect_TARGET_final;   vec_intersect_TARGET_final.clear();
  vector<double> vec_X_int_TARGET;             vec_X_int_TARGET.clear();
  vector<double> vec_Y_int_TARGET;             vec_Y_int_TARGET.clear();

  if(b_bbar){    
    vec_intersect_TARGET_final = intersect_circle(vec_slope_bbar[0], vec_slope_bbar[1], R_TARGET, TOF1_Xloc[gap_to_fit-1][2], TOF1_Yloc[gap_to_fit-1][2]);
  }
  else{
    vec_intersect_TARGET_final = intersect_circle(a_final, b_final, R_TARGET, TOF1_Xloc[gap_to_fit-1][2], TOF1_Yloc[gap_to_fit-1][2]);    
  }
  vec_X_int_TARGET.push_back(vec_intersect_TARGET_final[0]);
  vec_Y_int_TARGET.push_back(vec_intersect_TARGET_final[1]);

  // Intersect with TOF1
  vector<double> vec_intersect_TOF1_final;    vec_intersect_TOF1_final.clear();
  vector<double> vec_X_int_TOF1;              vec_X_int_TOF1.clear();
  vector<double> vec_Y_int_TOF1;              vec_Y_int_TOF1.clear();

  if(b_bbar){
    vec_intersect_TOF1_final = intersect_circle(vec_slope_bbar[0], vec_slope_bbar[1], R_TOF1, TOF1_Xloc[gap_to_fit-1][2], TOF1_Yloc[gap_to_fit-1][2]);
  }
  else{
    vec_intersect_TOF1_final = intersect_circle(a_final, b_final, R_TOF1, TOF1_Xloc[gap_to_fit-1][2], TOF1_Yloc[gap_to_fit-1][2]);   
  }
  vec_X_int_TOF1.push_back(vec_intersect_TOF1_final[0]);
  vec_Y_int_TOF1.push_back(vec_intersect_TOF1_final[1]);

  // Intersect with SFT L1
  vector<double> vec_intersect_SFT_final;     vec_intersect_SFT_final.clear();
  vector<double> vec_X_int_SFT;               vec_X_int_SFT.clear();
  vector<double> vec_Y_int_SFT;               vec_Y_int_SFT.clear();

  if(b_bbar){
    vec_intersect_SFT_final = intersect_circle(vec_slope_bbar[0], vec_slope_bbar[1], R_SFT_L1, TOF1_Xloc[gap_to_fit-1][2], TOF1_Yloc[gap_to_fit-1][2]);
  }
  else{
    vec_intersect_SFT_final = intersect_circle(a_final, b_final, R_SFT_L1, TOF1_Xloc[gap_to_fit-1][2], TOF1_Yloc[gap_to_fit-1][2]);   
  }
  vec_X_int_SFT.push_back(vec_intersect_SFT_final[0]);
  vec_Y_int_SFT.push_back(vec_intersect_SFT_final[1]);

  /// LENGTH IN TARGET
  double length_in_target = -1.;
  length_in_target = distance(vec_kstop_FINAL[0],vec_kstop_FINAL[1],
                               vec_X_int_TARGET[0],vec_Y_int_TARGET[0]);


  /// TDC DIFF
  tdc_ck_corr = tdc_trigger[0][0] - TDC_ck_avg2;
  TDC_diff = 550 - (0.91*TDC_average_NEW - 0.625*tdc_ck_corr);


  /// SUM ADC HG
  sum_ADC_HG_lepton = 0;

  for(unsigned int i=0; i<vec_leptons_pruned_unique.size(); i++){
    sum_ADC_HG_lepton += ADC_High_TARGET[vec_leptons_pruned_unique[i]]+TARGET_ADC_Thr_HG_Offset;
  }




  /// AVERAGE TDCs
  // Leptons
  sum_TDC_lepton = 0.;
  counter_TDC_lepton = 0.;
  Average_TDC_lepton = 0.;

  for(unsigned int i=0; i<vec_leptons_pruned_unique.size(); i++){
    for(int j=0; j<16; j++){
      if(tdc_le_target[vec_leptons_pruned_unique[i]][j]>kaon_TDC_min && tdc_le_target[vec_leptons_pruned_unique[i]][j]<kaon_TDC_max){
        sum_TDC_lepton += 550-(0.91*tdc_le_target[vec_leptons_pruned_unique[i]][j] - 0.625*(tdc_vt48[0][0]-TDC_ck_avg2));
        counter_TDC_lepton++;
      }
    }
  }
  if (counter_TDC_lepton > 0) Average_TDC_lepton = sum_TDC_lepton/counter_TDC_lepton;
  else Average_TDC_lepton = -1;
  
  // Kaons
  sum_TDC_kaon = 0.;
  counter_TDC_kaon = 0.;
  Average_TDC_kaon = 0.;

  for(unsigned int i=0; i<V_kaon_bars.size(); i++){
    for(int j=0; j<16; j++){
      if(tdc_le_target[V_kaon_bars[i]][j]>TDC_min_Kstop && tdc_le_target[V_kaon_bars[i]][j]<TDC_max_Kstop){
        sum_TDC_kaon += 550-(0.91*tdc_le_target[V_kaon_bars[i]][j] - 0.625*(tdc_vt48[0][0]-TDC_ck_avg2));
        counter_TDC_kaon++;
      }
    }
  }

  if (counter_TDC_kaon > 0) Average_TDC_kaon = sum_TDC_kaon/counter_TDC_kaon;
  else Average_TDC_kaon = -1;




  /// CALCULATE FINAL TRACK ANGLE
  double alpha_guide = 999.99;
  double angle_final = 999.99;
  double delta_phi = 999.99;
  double delta_phi_deg = 999.99;

  alpha_guide = atan(a_final);
  angle_final = angle_calculation(a_final, vec_kstop_FINAL[0], vec_kstop_FINAL[1], vec_X_int_TOF1[0], vec_Y_int_TOF1[0]);

  delta_phi = sin(alpha_guide)*cos(alpha_guide)*(a_final_error/a_final);
  delta_phi_deg = (180./PI)*delta_phi;

  if(b_bbar){
    angle_final = angle_calculation(vec_slope_bbar[0], vec_bbar[0], vec_bbar[1], vec_bbar[2], vec_bbar[3]);
    delta_phi_deg = 99.99;
    ChiS_final = -1.;
    ndf_final = 1;
  }


  // UNCERTAINTY
  double delta_X = 999.99;
  double delta_Y = 999.99;

  delta_X = sqrt(pow(TARGET_Errors_X,2) + pow((sin(angle_final)*R_TARGET*delta_phi),2));
  delta_Y = sqrt(pow(TARGET_Errors_Y,2) + pow((cos(angle_final)*R_TARGET*delta_phi),2));

  /// HAS EDGE BARS  <---- TO REVIEW !
  int HasEdgeBars;
  if(has_edge_bars) HasEdgeBars = 1;
  else HasEdgeBars = 0;

  /// CHECK IF GOOD TOF1 EVENT  <--- TO REVIEW !
  bool IsGoodTOF1 = false;
  IsGoodTOF1 = _Good_TOF1(gap_to_fit-1,vec_intersect_TOF1_final);





















//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////












//   // UNUSEDPo
//   if(Zlist_flag){

//     ZZ.clear();
//     Z_selected.clear();

//     ZZ = SFT_Test_CR_Event2(Run_Number, ievt, angle_final_guide, C2X_centroid, length_in_target);

//     /// TOF1_Z_range
//     Z_TOF1[gap_to_fit-1] =
//     Slope_correction_TOF1[gap_to_fit-1]*((0.025*150*0.5*(TDC_TOF1U[gap_to_fit-1] - TDC_TOF1D[gap_to_fit-1])) + Offset_TOF1[gap_to_fit-1]);

//     if(Switch_Printout !=0){
//       cerr << "Z_TOF1  =  " << Z_TOF1[gap_to_fit-1] << " mm" << endl;
//     }

//     for(unsigned i=0; i<ZZ.size(); i++){
//         Z_selected.push_back(ZZ[i]);
//     }

//     /// SORT VECTOR ELEMENTS
//     sort(Z_selected.begin(), Z_selected.end());

//     if(Switch_Printout !=0){
//       cerr << "Z selected:  ";
//       for(unsigned k=0; k<Z_selected.size(); k++){
//         cerr << Z_selected[k] << "  ";
//       }
//     }
//   }



  /// PRINTOUTS
  int TDC_print[256];
  for(int i=0; i<256; i++) TDC_print[i] = -1;


  //cout << endl;
  
  if(Switch_Printout != 0){
    //enable_cout = 2;
    _ED_printouts(Name_finput, NoFile, ParsTarg1,
                  nentries, Version,
                  Good_Event, Switch, enable_cout, ievt,
                  kaon_TDC_min, kaon_TDC_max, TDC_min_Kstop, TDC_max_Kstop,
                  TDC_average_NEW, TDC_LE_TARGET_sel_string, gap_to_fit, 
                  angle_final, delta_phi_deg, ChiS_final, ndf_final, 
                  vec_X_int_SFT, vec_Y_int_SFT, 
                  vec_X_int_TOF1, vec_Y_int_TOF1,
                  vec_TOF1_Gap, length_in_target);

    //cout << endl;
  }







  /// OUTPUTS
  if(print_output != 0){
    cout << fixed;
    cout << setw(4) << Run_Number << " ";
    cout << setw(7) << ievt << " ";
    //cout << setw(1) << error_flag << " ";
    cout << setw(2) << Event_tag[0] << " ";
    cout << setw(2) << gap_to_fit << " ";
    cout << setw(2) << selected_TOF2 << " ";
    cout << setw(6) << setprecision(2) << angle_final << " ";
    cout << setw(6) << setprecision(2) << delta_phi_deg << " ";
    cout << setw(6) << setprecision(2) << delta_X << " ";
    cout << setw(6) << setprecision(2) << delta_Y << " ";
    cout << setw(7) << setprecision(2) << ChiS_final << " ";
    cout << setw(3) << setprecision(0) << ndf_final << " ";
    cout << setw(6) << setprecision(2) << ChiS_final/ndf_final << " ";
    cout << setw(3) << Lepton_counter << "  ";
    cout << setw(6) << setprecision(2) << vec_X_int_TARGET[0] << " ";
    cout << setw(6) << setprecision(2) << vec_Y_int_TARGET[0] << " " ;
    cout << setw(7) << setprecision(2) << N_centroid_X << " ";
    cout << setw(7) << setprecision(2) << N_centroid_Y << " ";
    cout << setw(7) << setprecision(2) << vec_fit_lines_intersect[0] << " ";
    cout << setw(7) << setprecision(2) << vec_fit_lines_intersect[1] << " ";
    cout << setw(7) << setprecision(2) << vec_kstop_FINAL[0] << " ";             // 21
    cout << setw(7) << setprecision(2) << vec_kstop_FINAL[1] << " ";
    cout << setw(7) << setprecision(2) << R_Kstop << " ";             // 22
    cout << setw(3) << i_kaon_bar << " ";
    cout << setw(3) << V_kaon_bars.size() << " ";
    cout << setw(5) << setprecision(2) << ChiS_kaon_final/ndf_kaon_final << " ";
    cout << setw(3) << vec_Ck.size() << " ";
    cout << setw(3) << vec_Cpi.size() << " ";
    cout << setw(7) << setprecision(2) << length_in_target << " ";
    cout << setw(7) << setprecision(2) << C2X_centroid << " ";
    cout << setw(6) << setprecision(1) << TDC_diff << " ";
    cout << setw(6) <<  sum_ADC_HG_lepton << " ";
    cout << setw(7) << setprecision(2) << Average_TDC_lepton << " ";
    cout << setw(7) << setprecision(2) << Average_TDC_kaon << " ";
    cout << setw(1) << setprecision(0) << HasEdgeBars << " ";
    cout << setw(1) << setprecision(0) << IsGoodTOF1 << " ";
    cout << endl;

    if(batch == 0){
      if(flag_pruning_line) cout << "LINE PRUNING USED" << endl;
      else cout << "TRIANGLE PRUNING USED " << endl;

      cout << endl;
      cout << endl;
    }
  }


  if(batch==1 && (batch_flag== 0 || batch_flag==2 || batch_flag==3)){ 
    fout << fixed;
    fout << setw(4) << Run_Number << " ";
    fout << setw(7) << ievt << " ";
    //fout << setw(1) << error_flag << " ";
    fout << setw(2) << Event_tag[0] << " ";
    fout << setw(2) << gap_to_fit << " ";
    fout << setw(2) << selected_TOF2 << " ";
    fout << setw(6) << setprecision(2) << angle_final << " ";
    fout << setw(6) << setprecision(2) << delta_phi_deg << " ";
    fout << setw(6) << setprecision(2) << delta_X << " ";
    fout << setw(6) << setprecision(2) << delta_Y << " ";
    fout << setw(7) << setprecision(2) << ChiS_final << " ";
    fout << setw(3) << setprecision(0) << ndf_final << " ";
    fout << setw(6) << setprecision(2) << ChiS_final/ndf_final << " ";
    fout << setw(3) << Lepton_counter << "  ";
    fout << setw(6) << setprecision(2) << vec_X_int_TARGET[0] << " ";
    fout << setw(6) << setprecision(2) << vec_Y_int_TARGET[0] << " " ;
    fout << setw(7) << setprecision(2) << N_centroid_X << " ";
    fout << setw(7) << setprecision(2) << N_centroid_Y << " ";
    fout << setw(7) << setprecision(2) << vec_fit_lines_intersect[0] << " ";
    fout << setw(7) << setprecision(2) << vec_fit_lines_intersect[1] << " ";
    fout << setw(7) << setprecision(2) << vec_kstop_FINAL[0] << " ";             // 21
    fout << setw(7) << setprecision(2) << vec_kstop_FINAL[1] << " ";
    fout << setw(7) << setprecision(2) << R_Kstop << " ";             // 22
    fout << setw(3) << i_kaon_bar << " ";
    fout << setw(3) << V_kaon_bars.size() << " ";
    fout << setw(5) << setprecision(2) << ChiS_kaon_final/ndf_kaon_final << " ";
    fout << setw(3) << vec_Ck.size() << " ";
    fout << setw(3) << vec_Cpi.size() << " ";
    fout << setw(7) << setprecision(2) << length_in_target << " ";
    fout << setw(7) << setprecision(2) << C2X_centroid << " ";
    fout << setw(6) << setprecision(1) << TDC_diff << " ";
    fout << setw(6) <<  sum_ADC_HG_lepton << " ";
    fout << setw(7) << setprecision(2) << Average_TDC_lepton << " ";
    fout << setw(7) << setprecision(2) << Average_TDC_kaon << " ";
    fout << setw(1) << setprecision(0) << HasEdgeBars << " ";
    fout << setw(1) << setprecision(0) << IsGoodTOF1 << " ";
    fout << endl;

    fout_csv << fixed;
    fout_csv << Run_Number << ",";                                               //  1
    fout_csv << ievt << ",";                                                     //  2
    fout_csv << error_flag << ",";                                               //  3
    fout_csv << Event_tag[0] << ",";                                             //  4
    fout_csv << gap_to_fit << ",";                                               //  5
    fout_csv << selected_TOF2 << ",";                                            //  6
    fout_csv << setprecision(2) << angle_final << ",";                     //  7
    fout_csv << setprecision(2) << delta_phi_deg << ",";                         //  8
    fout_csv << setprecision(2) << delta_X << ",";                               //  9
    fout_csv << setprecision(2) << delta_Y << ",";                               // 10
    fout_csv << setprecision(2) << ChiS_final << ",";                     // 11
    fout_csv << setprecision(0) << ndf_final << ",";                      // 12
    fout_csv << setprecision(2) << ChiS_final/ndf_final << ",";    // 13
    fout_csv << Lepton_counter << ",";                                           // 14
    fout_csv << setprecision(2) << vec_X_int_TARGET[0] << ",";              // 15
    fout_csv << setprecision(2) << vec_Y_int_TARGET[0] << "," ;             // 16
    fout_csv << setprecision(2) << N_centroid_X << ",";      // 17
    fout_csv << setprecision(2) << N_centroid_Y << ",";      // 18
    fout_csv << setprecision(2) << vec_fit_lines_intersect[0] << ",";            // 19
    fout_csv << setprecision(2) << vec_fit_lines_intersect[1] << ",";            // 20
    fout_csv << setprecision(2) << vec_kstop_FINAL[0] << ",";                 // 21
    fout_csv << setprecision(2) << vec_kstop_FINAL[1] << ",";                 // 22
    fout_csv << setw(7) << setprecision(2) << R_Kstop << ",";             // 22
    fout_csv << i_kaon_bar << ",";                                               // 23
    fout_csv << V_kaon_bars.size() << ",";                                     // 24
    fout_csv << setprecision(2) << ChiS_kaon_final/ndf_kaon_final << ",";                    // 25
    fout_csv << vec_Ck.size() << ",";                                            // 26
    fout_csv << vec_Cpi.size() << ",";                                           // 27
    fout_csv << setprecision(2) << length_in_target << ",";                      // 28
    fout_csv << setprecision(2) << C2X_centroid << ",";                          // 29
    fout_csv << setprecision(1) << TDC_diff << ",";                              // 30
    fout_csv << sum_ADC_HG_lepton << ",";                                        // 31
    fout_csv << setprecision(2) << Average_TDC_lepton << ",";                    // 32
    fout_csv << setprecision(2) << Average_TDC_kaon << ",";                      // 33
    fout_csv << setprecision(0) << HasEdgeBars << ",";                         // 34
    fout_csv << setprecision(0) << IsGoodTOF1;                          // 35
    fout_csv << endl;
    fout_csv.close(); 
  }


    /// DRAW EVENT DISPLAY CANVAS
    if(Switch_Display!=0){

      /// DRAW FITTING CANVAS
      _Drawing_fitting(Run_Number, ievt, To_rotate_TEST,
                       gr_leptons_fit1, gr_leptons_fit2, gr_leptons_fit3,
                       gr_kaons_final, line_fit_kaons,
                       int_with_circle);


      /// DRAW PRUNING CANVAS
      _Drawing_pruning(Run_Number, ievt,
                       V_kaon_bars, N_Kstop_0,
                       vec_pruned_triangle_NEW, vec_pruned_line_NEW,
                       gr_kaons_final, gr_kaons_unfitted, line_fit_kaons,
                       gr_leptons_0_NEW, line_fit_leptons);
 
      /// DRAW EVENT DISPLAY
      _Drawing_Event_Display(Run_Number, ievt, to_restore,
                            gap_to_fit, selected_TOF1_section,  
                            gr_leptons_final,gr_kaons_final,line_fit_kaons,
                            restored_fitl_line,
                            TARGET_ADC_Thr_HG_Offset, kaon_TDC_min, kaon_TDC_max,
                            HG_KAON, LG_KAON, TDC_min_Kstop, TDC_max_Kstop, i_kaon_bar,
                            has_kaon_sub, has_TDC_hit_Kstop,
                            has_ADC_TOF1_hit, has_TDC_TOF1_hit,
                            has_ADC_TOF2_hit, has_TDC_TOF2_hit,
                            has_TDC_hit, Switch, 
                            ADC_High_TARGET, ADC_Low_TARGET,
                            angle_final, delta_phi_deg,
                            ChiS_final, ndf_final,
                            vec_X_int_SFT, vec_Y_int_SFT,
                            vec_X_int_TARGET, vec_Y_int_TARGET,
                            vec_X_int_TOF1, vec_Y_int_TOF1,
                            vec_kstop_FINAL,
                            vec_intersect_TARGET_final,
                            b_bbar, vec_bbar);
    }

    return;

  } // End void


// // FUNCTIONS
// double SFT_Test_CR_Event(int Run_Number, int evt, double phi, bool to_print, double C2X_centroid, double length_in_target){
//   Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];
//   Int_t adc_low_sft[128];
//   Int_t tdc_le_sft[128][16];
//   Int_t tdc_te_sft[128][16];

//   Int_t HG_SFT_ADC_Thr[128] = {0};

//   for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

//   Double_t ADC_High_SFT_corr[128];
//   Int_t has_TDC_SFT_hit[128] = {0};

//   char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

//   char file_mapping[200];
//   sprintf(file_mapping,"../Mapping");

//   char par_finput[200];
//   sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

//   Int_t par_temp[2][128];
//   ifstream fdat(par_finput,ios::in);
//   for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
//   fdat.close();

//   char path_input[200];
//   sprintf(path_input,"%s",path_merged);

//   char Name_finput[200];
//   sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

//   TChain *fChain= new TChain("Tree");
//   fChain->Add(Name_finput);
//   fChain->SetMakeClass(1);

//   fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
//   fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
//   fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
//   fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

//   fChain->GetEntry(evt);

//   for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
//     ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
//   }

//   for(int j=0 ; j<128 ; j++){
//     if(ADC_High_SFT[j]<0) ADC_High_SFT_corr[j]=0;
//     else ADC_High_SFT_corr[j]=ADC_High_SFT[j];
//   }

//   for(Int_t ii=0; ii<128; ii++){
//     for (Int_t qq=0; qq<6; qq++) {
//       if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
//     }
//   }

//   double sft_z_selected = 0.;
//   sft_z_selected = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, to_print, 0, false, C2X_centroid, length_in_target);

//   return sft_z_selected;
// }

// vector<double> SFT_Test_CR_Event2(int Run_Number, int evt, double phi, double C2X_centroid, double length_in_target){
//   Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];
//   Int_t adc_low_sft[128];
//   Int_t tdc_le_sft[128][16];
//   Int_t tdc_te_sft[128][16];

//   Int_t HG_SFT_ADC_Thr[128] = {0};

//   for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

//   Double_t ADC_High_SFT_corr[128];
//   Int_t has_TDC_SFT_hit[128] = {0};

//   char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

//   char file_mapping[200];
//   sprintf(file_mapping,"../Mapping");

//   char par_finput[200];
//   sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

//   Int_t par_temp[2][128];
//   ifstream fdat(par_finput,ios::in);
//   for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
//   fdat.close();

//   char path_input[200];
//   sprintf(path_input,"%s",path_merged);

//   char Name_finput[200];
//   sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

//   TChain *fChain= new TChain("Tree");
//   fChain->Add(Name_finput);
//   fChain->SetMakeClass(1);
//   fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
//   fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
//   fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
//   fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

//   fChain->GetEntry(evt);

//   for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
//     ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
//   }

//   for(int j=0 ; j<128 ; j++){
//     if(ADC_High_SFT[j]<0) ADC_High_SFT_corr[j]=0;
//     else ADC_High_SFT_corr[j]=ADC_High_SFT[j];
//   }

//   for(Int_t ii=0; ii<128; ii++){
//     for(Int_t qq=0; qq<6; qq++) {
//       if(tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
//     }
//   }

//   vector<double> ZZ;
//   ZZ = Z_Avg(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, phi, true, 0, false, C2X_centroid, length_in_target);

//   return ZZ;
// }

