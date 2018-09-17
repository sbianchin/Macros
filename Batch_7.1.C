#include <stdio.h>
#include <chrono>
// #include "Batch_Variables.h"
#include "Event_Display_7.1.C"

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
  char batch_output_csv[50];
  char batch_output2[50];
  char Root_Name[50];
  bool Zlist_flag=false;

  int Read_Run_Number;
  int counter = 0; 
  ifstream fdat_read(file_batch, ios::in);   
  TFile *f;

  sprintf(Root_Name, "TGT_Run%d.root", Run_Number);   
  
  if(flag==0 || flag==2){
    f = new TFile(Root_Name,"RECREATE"); 
  }

  TTree tree("TGT_tree","TARGET Tree");

  tree.Branch("RunNumber",&RunNumber);
  tree.Branch("EvtNumber",&EvtNumber);
  tree.Branch("FlagErr",&error_flag);
  tree.Branch("EvTag",&Event_tag[0]);
  tree.Branch("TOF1Gap",&gap_to_fit);
  tree.Branch("TOF2Gap",&selected_TOF2);
  tree.Branch("Phi",&angle_final_guide);
  tree.Branch("dPhi",&Delta_phi_deg);
  tree.Branch("dX",&Delta_X); // !!!
  tree.Branch("dY",&Delta_Y); // !!!
  tree.Branch("LeptonChiS",&Chis_lepton_fit_3);
  tree.Branch("LeptonNDF",&ndf_lepton_fit_3); 
  tree.Branch("LeptonRedChiS",&reduced_Chis);
  tree.Branch("LeptonNb",&lepton_counter);  // !!!
  tree.Branch("xIntTgt",&xx_int_TDC_TARGET);
  tree.Branch("yIntTgt",&yy_int_TDC_TARGET);
  tree.Branch("xKaonCentroid",&xx_kaon_centroid_coordinates);
  tree.Branch("yKaonCentroid",&yy_kaon_centroid_coordinates);
  tree.Branch("xFitLineIntersect",&xx_fit_lines_intersect);
  tree.Branch("yFitLineIntersect",&yy_fit_lines_intersect);
  tree.Branch("xKstop",&xx_kstop_final);
  tree.Branch("yKstop",&yy_kstop_final);
  tree.Branch("R_Kstop",&R_Kstop);
  tree.Branch("KaonBar",&i_kaon_bar);
  tree.Branch("KaonSize",&kaon_bar_size);
  tree.Branch("KaonRedChiS",&reduced_kaon_Chis);
  tree.Branch("CkSize",&Ck_size);
  tree.Branch("CpiSize",&Cpi_size);
  tree.Branch("LengthInTarget",&length_in_target);
  tree.Branch("C2XCentroid",&C2X_centroid);
  tree.Branch("TDCDiff",&TDC_diff);
  tree.Branch("SumADCHGLepton",&sum_ADC_HG_lepton);
  tree.Branch("AverageTDCLepton",&Average_TDC_lepton); // D??
  tree.Branch("AverageTDCKaon",&Average_TDC_kaon);
  tree.Branch("HasEdgeBars",&HasEdgeBars);
  tree.Branch("GoodTOF1Intersect",&Good_TOF1_intersect);



  if(flag==0 || flag==2 || flag==3){
    sprintf(batch_output, "Run_%d_Batch_output.txt", Run_Number);
    ofstream fout;
    fout.open(batch_output);
    fout << "";
    fout.close();

    sprintf(batch_output_csv, "Run_%d_Batch_output.csv", Run_Number);
    ofstream fout_csv;
    fout_csv.open(batch_output_csv);
    fout_csv << "";
    fout_csv.close();      
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
  
  if(flag == 0) {
    int i = min_ievt;
    if(nentries>0){
      if(FROM_FILE){
        while(counter<max_ievt && fdat_read.good()){
        //while(i<max_ievt && fdat_read.good()){
          fdat_read >> Read_Run_Number;
          counter++;
          i = Read_Run_Number;
          //if(fdat_read.eof()) break;
          //cout << "TEST : " << counter << endl; 
          if(counter%evt_counter==0 && i>0) cout << counter << " events processed!" << endl;
          Event_Display(Run_Number, Read_Run_Number, 0, 0, 0, 0, 1, 0);
          if(Good_Event) tree.Fill();
        }
        tree.Write();
      }
      else{
        while(i<=max_ievt){
          //cout << "TEST : " << i << endl;
          if(i%evt_counter==0 && i>0) cout << i << " events processed!" << endl;
          Event_Display(Run_Number, i, 0, 0, 0, 0, 1, 0);
          if(Good_Event) tree.Fill();
          i++;
        }
        tree.Write();
      }
    }
    else{
      cout << endl;
      cout << "ERROR : THIS FILE DOES NOT EXIST !" << endl;
      cout << endl;
      cout << endl;
      return;
    }
  }
  else if(flag == 1){
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
          //tree.Fill();
        }
        //tree.Write();
      }
      else{
        while(i<=max_ievt){
          if(i%evt_counter==0 && i>0) cout << i << " events processed!" << endl;
            Event_Display(Run_Number, i, 0, 0, 0, 1, 1, 1);
            //if(Good_Event) tree.Fill();
          i++;
        }
        //tree.Write();
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
  else if(flag == 2) {
    int i = 0;
    if(nentries>0){
      if(FROM_FILE){
        while(i<=nentries && fdat_read.good()){
          fdat_read >> Read_Run_Number;
          counter++;
          //if(fdat_read.eof()) break;
            
          if(counter%evt_counter==0 && i>0) cout << counter << " events processed!" << endl;
          Event_Display(Run_Number, Read_Run_Number, 0, 0, 0, 0, 1, 2);
          if(Good_Event) tree.Fill();
          i = Read_Run_Number;
        }
        tree.Write();
      }
      else{
        while(i<nentries){
          if(i%evt_counter==0 && i>0) cout << i << " events processed!" << endl;
          Event_Display(Run_Number, i, 0, 0, 0, 0, 1, 2);
          if(Good_Event) tree.Fill();
          i++;
        }
        tree.Write();
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
