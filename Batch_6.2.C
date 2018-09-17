#include <stdio.h>
#include <chrono>
// #include "Batch_Variables.h"
#include "Event_Display_6.4.C"


using namespace std;


Int_t StartChain(TChain &fChain, char Name_finput[200]) {
 
 	fChain.Reset();
    fChain.Add(Name_finput);
    fChain.SetMakeClass(1);
    fChain.SetBranchAddress("VT48_TDC",tdc_vt48);

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
    fChain.SetBranchAddress("EvFlag", Event_flag);
    fChain.SetBranchAddress("EvTag", Event_tag);

    Int_t nentries = 0;
    nentries = (Int_t)fChain.GetEntries();
    return nentries;
}


void Batch_Job(int Run_Number=3994, int min_ievt=0, int max_ievt=1, int flag=1, bool FROM_FILE=false) {
  gROOT->SetBatch(1);

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);

  char batch_output[50];
  char batch_output2[50];
  bool Zlist_flag=false;

  int Read_Run_Number;
  int counter = 0; 
  ifstream fdat_read(file_batch, ios::in); 

  
  if(flag==0 || flag==2 || flag==3){
    sprintf(batch_output, "Run_%d_Batch_output.txt", Run_Number);
    ofstream fout;
    fout.open(batch_output);
    fout << "";
    fout.close();
  }

  if(Zlist_flag){
    if(flag==0 || flag==2 || flag==3){
      sprintf(batch_output2, "Run_%d_Zlist.txt", Run_Number);
      ofstream fout2;
      fout2.open(batch_output2);
      fout2 << "";
      fout2.close();
    }
  }

  char Name_finput[200];
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




  cout << endl;
  cout << "File opened:  " << Name_finput << endl;
  cout << endl;

  Int_t nentries = StartChain(fChain, Name_finput);
  std::chrono::high_resolution_clock::time_point t1, t2;
  t1 = std::chrono::high_resolution_clock::now();
  
  if (flag == 0) {
    int i = min_ievt;
    if(nentries>0){
      if(FROM_FILE){
        while(i<max_ievt && fdat_read.good()){
          fdat_read >> Read_Run_Number;
          counter++;
          i = Read_Run_Number;
          //if(fdat_read.eof()) break;
            
          if(counter%evt_counter==0 && i>0) cout << counter << " events processed!" << endl;
          Event_Display(Run_Number, Read_Run_Number, 0, 0, 0, 0, 1, 0);
        }
      }
      else{
        while(i<=max_ievt){
          if(i%evt_counter==0 && i>0) cout << i << " events processed!" << endl;
          Event_Display(Run_Number, i, 0, 0, 0, 0, 1, 0);
          i++;
        }
      }
    }
    else{
      cout << endl;
      cout << "ERROR : THIS FILE DOES NOT EXIT !" << endl;
      cout << endl;
      cout << endl;
      return;
    }
  }
  else if (flag == 1) {
    int i = min_ievt;
    if(nentries>0){
      if(FROM_FILE){
        while(i<max_ievt && fdat_read.good()){
          fdat_read >> Read_Run_Number;
          counter++;
          i = Read_Run_Number;
          //if(fdat_read.eof()) break;
            
          if(counter%evt_counter==0 && i>0) cout << counter << " events processed!" << endl;
          Event_Display(Run_Number, Read_Run_Number, 0, 0, 0, 1, 1, 1);
        }
      }
      else{
        while(i<=max_ievt){
          if(i%evt_counter==0 && i>0) cout << i << " events processed!" << endl;
          Event_Display(Run_Number, i, 0, 0, 0, 1, 1, 1);
          i++;
        }
      }
    }
    else{
      cout << endl;
      cout << "ERROR : THIS FILE DOES NOT EXIT !" << endl;
      cout << endl;
      cout << endl;
      return;
    }
  }
  else if (flag == 2) {
    int i = 0;
    if(nentries>0){
      if(FROM_FILE){
        while(i<=nentries && fdat_read.good()){
          fdat_read >> Read_Run_Number;
          counter++;
          //if(fdat_read.eof()) break;
            
          if(counter%evt_counter==0 && i>0) cout << counter << " events processed!" << endl;
          Event_Display(Run_Number, Read_Run_Number, 0, 0, 0, 0, 1, 2);
          i = Read_Run_Number;
        }
      }
      else{
        while(i<nentries){
          if(i%evt_counter==0 && i>0) cout << i << " events processed!" << endl;
          Event_Display(Run_Number, i, 0, 0, 0, 0, 1, 2);
          i++;
        }
      }
    }
    else{
      cout << endl;
      cout << "ERROR : THIS FILE DOES NOT EXIT !" << endl;
      cout << endl;
      cout << endl;
      return;
    }
  }
 /*
  else if (flag == 3){
    char events[200];
    int read_counter = 0;
    sprintf(events,"test.txt");
    ifstream fdat_read(events, ios::in);
    int Read_ievt;
    if(nentries>0){
      while(fdat_read.good()){
        //if(read_counter%10000==0 && read_counter>0){
        if(read_counter%evt_counter==0 && read_counter>0){
          cout << read_counter << " events processed!" << endl;
        }
        fdat_read >> Read_ievt;
        Event_Display(Run_Number, Read_ievt, 0, 0, 0, 0, 1, 3);
        read_counter++;
      }
    }
    else{
      cout << endl;
      cout << "ERROR : THIS FILE DOES NOT EXIT !" << endl;
      cout << endl;
      cout << endl;
      return;    
    }
  }
  */

  
  else cout << "Error: Flag needs to be 0, 1, 2, or 3.";
  t2 = std::chrono::high_resolution_clock::now();

  cout << endl;
  cout << "Time taken for all events to be processed : "
                << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
                << " seconds\n";
  cout << endl;
}
