
void TOF1_Detector_Attributes() {
  Gap1l->SetLineWidth(10);
  Gap1l->SetLineColor(15);

  Gap2l->SetLineWidth(10);
  Gap2l->SetLineColor(15);

  Gap3l->SetLineWidth(10);
  Gap3l->SetLineColor(15);

  Gap4l->SetLineWidth(10);
  Gap4l->SetLineColor(15);

  Gap5l->SetLineWidth(10);
  Gap5l->SetLineColor(15);

  Gap6l->SetLineWidth(10);
  Gap6l->SetLineColor(15);

  Gap7l->SetLineWidth(10);
  Gap7l->SetLineColor(15);

  Gap8l->SetLineWidth(10);
  Gap8l->SetLineColor(15);

  Gap9l->SetLineWidth(10);
  Gap9l->SetLineColor(15);

  Gap10l->SetLineWidth(10);
  Gap10l->SetLineColor(15);

  Gap11l->SetLineWidth(10);
  Gap11l->SetLineColor(15);

  Gap12l->SetLineWidth(10);
  Gap12l->SetLineColor(15);
}

void _Drawing_fitting(int Run_Number, int ievt, bool to_rotate,
                      TGraph *gr_leptons_fit1, TGraph *gr_leptons_fit2, TGraph *gr_leptons_fit3,
                      TGraph *gr_kaons_final, TLine *line_fit_kaons,
                      vector<double> vec_intersect_TOF1)
{

  TOF1_Detector_Attributes();
  //void _Drawing_pruning(int Run_Number, int ievt, bool to_rotate,
  //            TGraph *gr_kaons_pr2, TGraph *gr_leptons_pr2,
  //            vector<int> vec_kaons, vector<int> vec_leptons, vector<double> vec_kstop_0,
  //            vector<int> vec_triangle_leptons, vector<int> vec_line_leptons){
	
  TCanvas *c_fitting;
  c_fitting = new TCanvas("Fitting","FITTING PROCEDURE",0,200,700,700);
  c_fitting->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
  c_fitting->Divide(2,2);

  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");

  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0); 
  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);


  char str_leptons_fit1[100];
  char str_leptons_fit2[100];
  char str_leptons_fit3[100];
  char str_kaons_fit[100];

  if(to_rotate){
    sprintf(str_leptons_fit1,"Lepton Fit 1 (ROTATED 90^{o} CCW)  |  Run %d  --  Event %d",Run_Number,ievt);
    sprintf(str_leptons_fit2,"Lepton Fit 2 (ROTATED 90^{o} CCW)  |  Run %d  --  Event %d",Run_Number,ievt);
    sprintf(str_leptons_fit3,"Lepton Fit 3 (ROTATED 90^{o} CCW)  |  Run %d  --  Event %d",Run_Number,ievt);
  }
  else{
    sprintf(str_leptons_fit1,"Lepton Fit 1  |  Run %d  --  Event %d",Run_Number,ievt);
    sprintf(str_leptons_fit2,"Lepton Fit 2  |  Run %d  --  Event %d",Run_Number,ievt);
    sprintf(str_leptons_fit3,"Lepton Fit 3  |  Run %d  --  Event %d",Run_Number,ievt);
  }

  sprintf(str_kaons_fit,"Kaon Fit  |  Run %d  --  Event %d",Run_Number,ievt);

  TH2F *h_leptons_fit1 = new TH2F("Leptons Fit1", str_leptons_fit1, 500, -50, 50, 500, -50, 50);
  TH2F *h_leptons_fit2 = new TH2F("Lepton Fit2", str_leptons_fit2, 500, -50, 50, 500, -50, 50);
  TH2F *h_leptons_fit3 = new TH2F("Lepton Fit3", str_leptons_fit3, 500, -50, 50, 500, -50, 50);
  TH2F *h_kaons_fit = new TH2F("Kaon Fit", str_kaons_fit, 500, -50, 50, 500, -50, 50);


  //TARGET CENTER
  vector<double> vec_xx_Target_Center;  vec_xx_Target_Center.clear();
  vector<double> vec_yy_Target_Center;  vec_yy_Target_Center.clear();
  vec_xx_Target_Center.push_back(0.);
  vec_yy_Target_Center.push_back(0.);

  TGraph *gr_Target_Center = new TGraph(vec_xx_Target_Center.size(),&vec_xx_Target_Center[0],&vec_yy_Target_Center[0]);
  gr_Target_Center->SetMarkerStyle(5);
  gr_Target_Center->SetMarkerColor(1);
  gr_Target_Center->SetMarkerSize(1);


  //KAONS
/*  vector<double> vec_x_kaons;   vec_x_kaons.clear();
  vector<double> vec_y_kaons;   vec_y_kaons.clear();
  vector<double> vec_ex_kaons;  vec_ex_kaons.clear();
  vector<double> vec_ey_kaons;  vec_ey_kaons.clear();

  for(unsigned int i=0; i<vec_kaons.size(); i++){
    vec_x_kaons.push_back(Xloc[vec_kaons[i]]);
    vec_y_kaons.push_back(Yloc[vec_kaons[i]]);
    vec_ex_kaons.push_back(TARGET_Errors_X);
    vec_ey_kaons.push_back(TARGET_Errors_Y);
  }

  TGraph *gr_kaons;
  gr_kaons = new TGraphErrors(vec_kaons.size(),&vec_x_kaons[0],&vec_y_kaons[0],
                  &vec_ex_kaons[0],&vec_ey_kaons[0]);
    
  gr_kaons->SetMarkerColor(4);
  gr_kaons->SetMarkerStyle(21);
  gr_kaons->SetMarkerSize(0.8);
*/

  //FIT2 INTERSECT WITH TOF1
  vector<double> vec_x_intersect_TOF1;  vec_x_intersect_TOF1.clear();
  vector<double> vec_y_intersect_TOF1;  vec_y_intersect_TOF1.clear();

  if(vec_intersect_TOF1.size()>0){
    vec_x_intersect_TOF1.push_back(vec_intersect_TOF1[0]);
    vec_y_intersect_TOF1.push_back(vec_intersect_TOF1[1]);
  }

  TGraph *gr_intersect_TOF1;
  gr_intersect_TOF1 = new TGraph(vec_x_intersect_TOF1.size(), 
                      &vec_x_intersect_TOF1[0], &vec_y_intersect_TOF1[0]);

  gr_intersect_TOF1->SetMarkerColor(4);
  gr_intersect_TOF1->SetMarkerStyle(20);
  gr_intersect_TOF1->SetMarkerSize(0.8);

  /// PANEL 1
  c_fitting->cd(1);
  h_leptons_fit1->Draw();
  A1->Draw();
  A2->Draw();
  gr_Target_Center->Draw("sameP");

  //TOF1 Detectors
  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");
  ell_L1->Draw("same");

  gr_leptons_fit1->Draw("sameP");

  //////////////////////////////////////

  /// PANEL 2
  c_fitting->cd(2);
  h_kaons_fit->Draw();
  A1->Draw();
  A2->Draw();
  gr_Target_Center->Draw("sameP");

  //TOF1 Detectors
  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");
  ell_L1->Draw("same");

  gr_kaons_final->Draw("sameP");
  line_fit_kaons->Draw("same");

  //////////////////////////////////////

  /// PANEL 3
  c_fitting->cd(3);
  h_leptons_fit2->Draw();
  A1->Draw();
  A2->Draw();
  gr_Target_Center->Draw("sameP");

  //TOF1 Detectors
  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");
  ell_L1->Draw("same");

  gr_leptons_fit2->Draw("sameP");
  if(vec_intersect_TOF1.size()>0) gr_intersect_TOF1->Draw("sameP");
 
  //////////////////////////////////////

  /// PANEL 4
  c_fitting->cd(4);
  h_leptons_fit3->Draw();
  A1->Draw();
  A2->Draw();
  gr_Target_Center->Draw("sameP");

  //TOF1 Detectors
  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");
  ell_L1->Draw("same");

  gr_leptons_fit3->Draw("sameP");
}

