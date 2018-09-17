#ifndef __CINT__
#include <stdio.h>
#include <math.h>
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
#include "CommonParameters.h"
#endif


using namespace std;

double calculate_chi2(vector<double> vec_xx, vector<double> vec_yy, double m, double b) {
  double chi2 = 0.;
  vector<double> residuals;
  for (int i = 0; i < vec_xx.size(); ++i) {
    double residual = pow((vec_yy[i] - m*vec_xx[i] - b),2);// / (m*vec_xx[i] + b));
    chi2 += residual;
    cout << "Residual is: " << residual << endl;
  }
  // chi2 = accumulate(residuals.begin(), residuals.end(), 0);
  return chi2;
}
void Chis_test2(){

   gROOT->Clear();
   gROOT->Reset();

   gStyle->Clear();
   TH1::AddDirectory(kFALSE);
   gStyle->SetOptStat(0);


   Float_t Xloc[256] = { -7.75,-4.65,-1.55,1.55,4.65,7.75,
    -13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,
    -17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,
    -20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,
    -23.25,-20.15, -17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,
    -23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,
    -26.35,-23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,26.35,
    -26.35,-23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,26.35,
    -26.35,-23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,26.35,
    -26.35,-23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,26.35,
    -26.35,-23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,26.35,
    -26.35,-23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,26.35,
    -23.25,-20.15, -17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,
    -23.25,-20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,23.25,
    -20.15,-17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,20.15,
    -17.05,-13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,17.05,
    -13.95,-10.85,-7.75,-4.65,-1.55,1.55,4.65,7.75,10.85,13.95,
    -7.75,-4.65,-1.55,1.55,4.65,7.75};

   Float_t Yloc[256] = {         26.35,26.35,26.35,26.35,26.35,26.35,
    23.25,23.25,23.25,23.25,23.25,23.25,23.25,23.25,23.25,23.25,
    20.15,20.15,20.15,20.15,20.15,20.15,20.15,20.15,20.15,20.15,20.15,20.15,
    17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,17.05,
    13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,13.95,
    10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,10.85,
    7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75, 7.75,
    4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65, 4.65,
    1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55,
    -1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,-1.55,
    -4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,-4.65,
    -7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,-7.75,
    -10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,-10.85,
    -13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,-13.95,
    -17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,-17.05,
    -20.15,-20.15,-20.15,-20.15,-20.15,-20.15,-20.15,-20.15,-20.15,-20.15,-20.15,-20.15,
    -23.25,-23.25,-23.25,-23.25,-23.25,-23.25,-23.25,-23.25,-23.25,-23.25,
    -26.35,-26.35,-26.35,-26.35,-26.35,-26.35};

   double TARGET_Errors_X = sqrt(pow(3.2,2)/12);

   double TARGET_Errors_Y = sqrt(pow(3.2,2)/12);


   int n_bar, TOF1_counter;


/*  int bar_hit[20] = {3,10,32,132,4,
                     11,20,21,33,46,
                     47,62,95,97,113,
                     115,130,147,148,165};
*/
  vector<int> bar_hit = {109, 125, 126, 127, 142, 143, 158, 159, 160, 174,
    175, 176, 190, 191, 193, 194, 205, 210};
    n_bar = bar_hit.size();
   float R_TARGET = 29.0;
   TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
   TH2F *h_NULL = new TH2F("NULL", "NULL", 100,-50,50,100,-50,50);

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,700);


   int bar_number;
   double ChiS;
   int ndf;
   double a_fit, b_fit;
   double sum = 0;
   vector<double> vec_xx;        vec_xx.clear();
   vector<double> vec_yy;        vec_yy.clear();
   vector<double> vec_ex;        vec_ex.clear();
   vector<double> vec_ey;        vec_ey.clear();
   vector<int> barn;             barn.clear();
   vector<double> Y_expected;    Y_expected.clear();


   for(int i=0; i<n_bar; i++){
      vec_xx.push_back(Xloc[bar_hit[i]]);
      vec_yy.push_back(Yloc[bar_hit[i]]);
      // vec_ey.push_back(ex);
      // vec_ex.push_back(ey);
      // barn.push_back(bar_hit[i]);
   }

   TOF1_counter = 3;

/*   for(int j=0; j<3; j++){
      vec_xx.push_back(TOF1_Xloc[TOF1_counter-1][1]);
      //vec_ex.push_back(0.);
      vec_ex.push_back(TOF1_Errors_X[TOF1_counter-1][1]);

      vec_yy.push_back(TOF1_Yloc[TOF1_counter-1][1]);
      //vec_ey.push_back(0.);
      vec_ey.push_back(TOF1_Errors_Y[TOF1_counter-1][1]);
   }
*/


  //  TGraph *graph = new TGraphErrors(vec_xx.size(),&vec_xx[0],&vec_yy[0],
  //     &vec_ex[0],&vec_ey[0]);
  TGraph *graph = new TGraph(vec_xx.size(),&vec_xx[0],&vec_yy[0]);

   TF1 *func_fit = new TF1("gr", "pol1");

   graph->SetMarkerColor(4);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.8);
   graph->GetYaxis()->SetTitle("Y title");
   graph->GetXaxis()->SetLimits(-50.,50.);

   //gr->Fit("pol1","Q");
   //gr_fit = gr->GetFunction("pol1");


   //func_fit->SetParameter(0,0);
   //func_fit->SetParameter(1,1);
   //func_fit->SetParLimits(0,-50,50);
   //func_fit->SetParLimits(1,-50,50);

   graph->Fit("gr","Q");
   func_fit = graph->GetFunction("gr");
   func_fit->SetLineWidth(2);
   func_fit->SetLineColor(2);

   a_fit = func_fit->GetParameter(1);
   b_fit = func_fit->GetParameter(0);
   ChiS = func_fit->GetChisquare();
   ndf = func_fit->GetNDF();


   cout << endl;
   cout << endl;

   cout << "BAR NUMBERS : " << endl;
   cout << fixed;
   for(unsigned int j=0; j<barn.size(); j++){
      Y_expected.push_back((a_fit*Xloc[barn[j]])+b_fit);
      sum += pow(Yloc[barn[j]]-Y_expected[j],2);
      cout << setw(4) << barn[j] << "  ";
      cout << setw(6) << setprecision(2) << Xloc[barn[j]] << "  ";
      cout << setw(6) << setprecision(2) << Yloc[barn[j]] << "  ";
      cout << setw(7) << setprecision(5) << TARGET_Errors_X << "  ";
      cout << setw(7) << setprecision(5) << TARGET_Errors_Y << "     ";
      cout << setw(6) << setprecision(2) << Y_expected[j] << "  ";
      cout << setw(6) << setprecision(2) << Yloc[barn[j]]-Y_expected[j] << "  ";
      cout << endl;
   }
   cout << "SUM = " << setw(6) << setprecision(2) << sum << "  ";
   cout << endl;
   cout << endl;


   cout << endl;
   cout << endl;

   cout << "ChiS = " << ChiS << endl;
   double chi2 = calculate_chi2(vec_xx, vec_yy, a_fit, b_fit);
   cout << "Calculated chi2: " << chi2 << endl;
   cout << "NdF = " << ndf << endl;
   cout << endl;
   cout << "a_fit = " << a_fit << endl;
   cout << "b_fit = " << b_fit << endl;
   cout << endl;
   cout << endl;

   char ChiS_string[30];
   char ndf_string[30];
   sprintf(ChiS_string,"#chi^{2} = %3.2f", ChiS);
   sprintf(ndf_string,"NdF = %3i", ndf);

   TLatex *tex_ChiS;
   tex_ChiS = new TLatex(-45.,43.,ChiS_string);
   tex_ChiS->SetTextSize(0.05);
   tex_ChiS->SetLineWidth(2);

   TLatex *tex_ndf;
   tex_ndf = new TLatex(-45.,35.,ndf_string);
   tex_ndf->SetTextSize(0.05);
   tex_ndf->SetLineWidth(2);


   c1->cd();
   h_NULL->Draw();
   ell_Target->Draw("same");
   graph->Draw("sameP");
   tex_ChiS->Draw("same");
   tex_ndf->Draw("same");
}
