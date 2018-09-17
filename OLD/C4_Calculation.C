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
//#include "TSpectrum.h"
#include "TMarker.h"
//#include "Event_Display_MS.h"
#include "ANAPATH.h"
#include "Thresholds.h"
//#include "ListGoodEvents.h"
#endif

void C4_Calculation(float phi=90, float Z1_0 = 0., float p1 = 1., float p2 = 1.){

	float C1[16] = {0.};	//float Z1_0 = ; 
	float C2[16] = {0.};	float Z2_0 = 0.;   Z2_0 = Z1_0;
	float C3 = 0.;	float Z3 = 99999.;
	float C4 = 0.;	float Z4 = 99999.;
	float n3 = 7.;	float p3 = 1.;
	float n4 = 7.; 	float p4 = 1.;
	float pi = 3.14159265359;
	const float R_SFT[4] = {40.25, 41.1, 42.1, 42.95};
	const float alpha[4] = {3.57, 3.50, 3.87, 3.79};
	const float w_SFT[4] = {15.75, 15.75, 17.85, 17.85};
	const float l_SFT[4] = {w_SFT[0]/cos(alpha[0]*pi/180.), 
							w_SFT[1]/cos(alpha[1]*pi/180.), 
							w_SFT[2]/cos(alpha[2]*pi/180.),
							w_SFT[3]/cos(alpha[3]*pi/180.)};
	//float l4 = w4/cos(alpha[3]*pi/180.);
	//float l3 = w3/cos(alpha[2]*pi/180.);
	//float l2 = w2/cos(alpha[1]*pi/180.);
	//float l1 = w1/cos(alpha[0]*pi/180.);

	//float p1 = 8.;	float p2 = 13.;


	C4 = -(R_SFT[3] *(360-phi)*(pi/180)*tan(alpha[3]*pi/180) + (n4-1.)*l_SFT[3] + (1/17.)*(17.-p4)*l_SFT[3] + 0.5*l_SFT[3]/17.);

	Z4 = C4 + R_SFT[3]*(360-phi)*(pi/180)*tan(alpha[3]*pi/180) + (n4-1.)*l_SFT[3] + (1/17.)*(17.-p4)*l_SFT[3] + 0.5*l_SFT[3]/17.;

	cout << " " << endl;
	cout << "C4 = " << C4 << endl; 
	cout << "Z4 = " << Z4 << endl;

	cout << " " << endl;

	C3 = -(R_SFT[2] *(360-phi)*(pi/180)*tan(alpha[2]*pi/180) + (n3-1.)*l_SFT[2] + (1/17.)*(17.-p3)*l_SFT[2] + 0.5*l_SFT[2]/17.);
	Z3 = C3 + R_SFT[2]*(360-phi)*(pi/180)*tan(alpha[2]*pi/180) + (n3-1.)*l_SFT[2] + (1/17.)*(17.-p3)*l_SFT[2] + 0.5*l_SFT[2]/17.;

	cout << " " << endl;
	cout << "C3 = " << C3 << endl; 
	cout << "Z3 = " << Z3 << endl;

	cout << " " << endl;


	for(int n1=1; n1<16; n1++){
		C1[n1] = Z1_0 - (R_SFT[0]*phi*(pi/180)*tan(alpha[0]*pi/180) + (n1-1.)*l_SFT[0] + (1/15.)*(p1-1.)*l_SFT[0] + 0.5*l_SFT[0]/15.);
	}

	for(int n2=1; n2<16; n2++){
		C2[n2] = Z2_0 - (R_SFT[1]*phi*(pi/180)*tan(alpha[1]*pi/180) + (n2-1.)*l_SFT[1] + (1/15.)*(p2-1.)*l_SFT[1] + 0.5*l_SFT[1]/15.);
	}

	cout << " " << endl;
	for(int i=0; i<15; i++){
		cout << "n = " << i+1 << "  " << C1[i+1] << "  " << C2[i+1] << endl; 
	}

	cout << " " << endl;

	float c1 = -125.298;
	float c2 = -130.764;
	float N1 = 8.;
	float N2 = 8.;

	p1 = 1-((15./l_SFT[0])*(c1 + R_SFT[0]*phi*(pi/180)*tan(alpha[0]*pi/180) + (N1-1.)*l_SFT[0] + 0.5*l_SFT[0]/15));
	p2 = 1-((15./l_SFT[1])*(c2 + R_SFT[1]*phi*(pi/180)*tan(alpha[1]*pi/180) + (N2-1.)*l_SFT[1] + 0.5*l_SFT[1]/15));

	cout << "############  Z = 0 for n1 = " << N1 << " and p1 = " << p1 << "  (if C1 = " << c1 << ")" << endl;
	cout << "############  Z = 0 for n2 = " << N2 << " and p2 = " << p2 << "  (if C2 = " << c2 << ")" << endl;
}