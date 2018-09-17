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
#include "ANAPATH.h"
#include "CommonParameters.h"
#include "ADC_Thresholds.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
//#include "MWPC_Thr.h"
#include "MWPC_Thr2.h"
#include "Pedestals.h"
#endif

#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"

using namespace std;
 
void Chis_test(){

   gStyle->SetOptStat(0);

   float R_TARGET = 29.0; 
   TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
   TH2F *h_NULL = new TH2F("NULL", "NULL", 100,-50,50,100,-50,50); 

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,700);


   int bar_number;
   double ChiS, ndf;
   vector<double> vec_xx;
   vector<double> vec_yy;
   vector<double> vec_ex;
   vector<double> vec_ey;
   vector<int> barn;

   cout << endl;
   cout << "Enter Bar Number :   ";
   
   cin >> bar_number;

   while(bar_number<256){

      cout << "Enter Bar Number :   ";  
      cin >> bar_number;   
      //test.push_back(i);

      if(bar_number>255){
         cout << "Thank you!" << endl;
         break;
      }
      
      vec_xx.push_back(Xloc[bar_number]);
      vec_yy.push_back(Yloc[bar_number]);
      vec_ex.push_back(1.6);
      vec_ey.push_back(1.6);
      barn.push_back(bar_number);
   }

   cout << endl;
   cout << endl;

   cout << "BAR NUMBERS : ";
   for(unsigned int j=0; j<barn.size(); j++){
      cout << barn[j] << "  ";
   }
   cout << endl;

   TGraph *gr = new TGraphErrors(vec_xx.size(),&vec_xx[0],&vec_yy[0],&vec_ex[0], &vec_ey[0]);
   TF1 *gr_fit = new TF1("gr", "pol1");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->GetYaxis()->SetTitle("Y title");
   gr->GetXaxis()->SetLimits(-50.,50.);
   gr->Draw("AP");
   gr->Fit("pol1","Q");
  
   gr_fit = gr->GetFunction("pol1");
   gr_fit->SetLineWidth(2);
   gr_fit->SetLineColor(2);

   ChiS = gr_fit->GetChisquare();
   ndf = gr_fit->GetNDF();

   cout << endl;
   cout << endl;

   cout << "ChiS = " << ChiS << endl;
   cout << "NdF = " << ndf << endl;
   cout << endl;
   


   c1->cd();
   h_NULL->Draw();
   ell_Target->Draw("same");
   gr->Draw("sameP");



}

