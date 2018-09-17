#ifndef COMMONPARAMETERS_H
#define COMMONPARAMETERS_H

unsigned int lepton_cluster_size = 2;
unsigned int kaon_cluster_size = 3;    // mod from 5 to 3 by MH  sept 29/17
int n_hit = 2;
int weight_edge_bars_fit1 = 0;
int weight_edge_bars_fit2 = 1;

int weight_TOF1_fit1 = 3;
int weight_TOF1_fit3 = 3;

int evt_counter = 1000;

const float R_TARGET = 29.0;
const float R_TOF1 = 47.1;
//const float R_SFT_L1 = 40.0;
const float R_SFT_L1 = 41.6;
const float R_Cut = 29.0;

bool dist_pruned = false;
const double DISTANCE_MAX_TO_KAON_PRUNED = 6.;
const double DISTANCE_MAX_TO_KAON_FINAL = 6.;

//Some Global Variables
    //char data_path[] = "/triumfcs/trshare/trek/E36/Data/Merged_Trees/";
    //char data_path[] = "../Data/root/";
    //char anal_save_path[] = "Results/";
    //char histo_save_path[] = "RootFiles/";
  const Int_t nBars = 256;
  const Int_t nFibrs = 128;
  const Int_t nEasiTDC = 16;
  const Int_t nTOF1 = 24;
  const Int_t nTOF2 = 48;
  const Int_t nGaps = 12;
  const Double_t PI = 3.14159265;
  const double TOF1_Z_cut = 15;

const int MWPC_cluster_separation = 1;

Double_t Xhit[nBars];
Double_t Yhit[nBars];
Double_t Ehit[nBars];
Double_t XKhit[nBars];
Double_t YKhit[nBars];
// Bad bars info
Int_t NoBar[nBars];
Int_t nNoBars = 0;
// thresholds etc.
    Double_t HGpedT[nBars];
    Double_t LGpedT[nBars];
    Double_t HGpedF[nFibrs];
    Double_t LGpedF[nFibrs];
    Double_t TpedL, TpedH;
    Double_t TcutLT[nBars], TcutHT[nBars];
    Double_t TcutLF[nFibrs], TcutHF[nFibrs];
// Event number
Int_t nEvent;
Int_t nhits = 0;  // count number of high gain above pedestal hits found
Int_t nKhits = 0;  // count number of Kaon hits found
Int_t nhitsorig = 0; //number of high gain pedestal hits before drop "bad" hits
Int_t nhitsL = 0; // count of number of low gain pedestal hits found
Int_t nhits_t = 0; // count of number of "real" high gain hits found for tdc
Int_t nhitsF = 0;  // count number fibre high gain pedestal hits found
Int_t nhitsLF = 0; // count fibre number of low gain pedestal hits found
Int_t nhits_Ft = 0; // count fibre number of "real" high gain hits found for tdc
//Bar positions
Float_t Xloc[nBars] = {
                                                    - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,
                                    -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,
                            -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,
                    -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,
            -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,
            -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,
    -26.35, -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,  26.35,
    -26.35, -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,  26.35,
    -26.35, -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,  26.35,
    -26.35, -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,  26.35,
    -26.35, -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,  26.35,
    -26.35, -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,  26.35,
            -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,
            -23.25, -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,  23.25,
                    -20.15, -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,  20.15,
                            -17.05, -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,  17.05,
                                    -13.95, -10.85, - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75,  10.85,  13.95,
                                                    - 7.75, - 4.65, - 1.55,   1.55,   4.65,   7.75
};
Float_t Yloc[nBars] = {
                                                     26.35,  26.35,  26.35,  26.35,  26.35,  26.35,
                                     23.25,  23.25,  23.25,  23.25,  23.25,  23.25,  23.25,  23.25,  23.25,  23.25,
                             20.15,  20.15,  20.15,  20.15,  20.15,  20.15,  20.15,  20.15,  20.15,  20.15,  20.15,  20.15,
                     17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,  17.05,
             13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,  13.95,
             10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,  10.85,
      7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,   7.75,
      4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,   4.65,
      1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,   1.55,
    - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55, - 1.55,
    - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65, - 4.65,
    - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75, - 7.75,
            -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85, -10.85,
            -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95, -13.95,
                    -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05, -17.05,
                            -20.15, -20.15, -20.15, -20.15, -20.15, -20.15, -20.15, -20.15, -20.15, -20.15, -20.15, -20.15,
                                    -23.25, -23.25, -23.25, -23.25, -23.25, -23.25, -23.25, -23.25, -23.25, -23.25,
                                                    -26.35, -26.35, -26.35, -26.35, -26.35, -26.35
};
Int_t BarRow[nBars] = {
                              0,  0,  0,  0,  0,  0,
                      1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
              3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
          4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
          5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
      6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
      7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
      8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,
      9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
     11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
         12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
         13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
             14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
                 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
                     16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
                             17, 17, 17, 17, 17, 17
};

Int_t BarCol[nBars] = {
                              6,  7,  8,  9, 10, 11,
                      4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
                  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
              2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
          1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
          1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
          1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
          1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
              2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
                  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
                      4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
                              6,  7,  8,  9, 10, 11
};
Int_t TOF1_Gap[nTOF1] = {1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12};
Int_t TOF1_UD[nTOF1] = {1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2};
Int_t TOF2_Gap[nTOF2] = {1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,
                         1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12};
Int_t TOF2_UD[nTOF2] = {1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,
                        3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4};
//Partners list???

/*
//OLD
Float_t C2_ZLoc[56] =
{-275.0, -265.0, -255.0, -245.0, -235.0, -225.0, -215.0, -205.0, -195.0, -185.0,
 -175.0, -165.0, -155.0, -145.0, -135.0, -125.0, -115.0, -105.0,  -95.0,  -85.0,
  -75.0,  -65.0,  -55.0,  -45.0,  -35.0,  -25.0,  -15.0,   -5.0,    5.0,   15.0,
   25.0,   35.0,   45.0,   55.0,   65.0,   75.0,   85.0,   95.0,  105.0,  115.0,
  125.0,  135.0,  145.0,  155.0,  165.0,  175.0,  185.0,  195.0,  205.0,  215.0,
  225.0,  235.0,  245.0,  255.0,  265.0,  275.0};
*/

//NEW
Float_t C2_ZLoc[56] =
{275.0, 265.0, 255.0, 245.0, 235.0, 225.0, 215.0, 205.0, 195.0, 185.0,
 175.0, 165.0, 155.0, 145.0, 135.0, 125.0, 115.0, 105.0,  95.0,  85.0,
  75.0,  65.0,  55.0,  45.0,  35.0,  25.0,  15.0,   5.0,    -5.0,   -15.0,
 -25.0,   -35.0,   -45.0,   -55.0,   -65.0,   -75.0,  -85.0,   -95.0,  -105.0,  -115.0,
-125.0,  -135.0,  -145.0,  -155.0,  -165.0,  -175.0,  -185.0,  -195.0,  -205.0,  -215.0,
-225.0,  -235.0,  -245.0,  -255.0,  -265.0,  -275.0};


Float_t C2_YLoc[16] =
{75.0, 65.0, 55.0, 45.0, 35.0, 25.0, 15.0, 5.0, -5.0, -15.0, -25.0, -35.0, -45.0, -55.0, -65.0, -75.0};


double TOF_Gap1XLoc[3] = {47.1*cos(350*PI/180), 47.1*cos(0*PI/180), 47.1*cos(10*PI/180)};
double TOF_Gap1YLoc[3] = {47.1*sin(350*PI/180), 47.1*sin(0*PI/180), 47.1*sin(10*PI/180)};

double TOF_Gap2XLoc[3] = {47.1*cos(30*PI/180), 47.1*cos(30*PI/180), 47.1*cos(30*PI/180)};
double TOF_Gap2YLoc[3] = {47.1*sin(30*PI/180), 47.1*sin(30*PI/180), 47.1*sin(30*PI/180)};

//double TOF_Gap2XLoc[3] = {47.1*cos(30*PI/180) - (25/3)*sin(30*PI/180), 47.1*cos(30*PI/180), 47.1*cos(30*PI/180) + (25/3)*sin(30*PI/180)};
//double TOF_Gap2YLoc[3] = {47.1*sin(30*PI/180) + (25/3)*cos(30*PI/180), 47.1*sin(30*PI/180), 47.1*sin(30*PI/180) - (25/3)*cos(30*PI/180)};

double TOF_Gap3XLoc[3] = {47.1*cos(50*PI/180), 47.1*cos(60*PI/180), 47.1*cos(70*PI/180)};
double TOF_Gap3YLoc[3] = {47.1*sin(50*PI/180), 47.1*sin(60*PI/180), 47.1*sin(70*PI/180)};

double TOF_Gap4XLoc[3] = {47.1*cos(80*PI/180), 47.1*cos(90*PI/180), 47.1*cos(100*PI/180)};
double TOF_Gap4YLoc[3] = {47.1*sin(80*PI/180), 47.1*sin(90*PI/180), 47.1*sin(100*PI/180)};

double TOF_Gap5XLoc[3] = {47.1*cos(110*PI/180), 47.1*cos(120*PI/180), 47.1*cos(130*PI/180)};
double TOF_Gap5YLoc[3] = {47.1*sin(110*PI/180), 47.1*sin(120*PI/180), 47.1*sin(130*PI/180)};

double TOF_Gap6XLoc[3] = {47.1*cos(140*PI/180), 47.1*cos(150*PI/180), 47.1*cos(160*PI/180)};
double TOF_Gap6YLoc[3] = {47.1*sin(140*PI/180), 47.1*sin(150*PI/180), 47.1*sin(160*PI/180)};

double TOF_Gap7XLoc[3] = {47.1*cos(170*PI/180), 47.1*cos(180*PI/180), 47.1*cos(190*PI/180)};
double TOF_Gap7YLoc[3] = {47.1*sin(170*PI/180), 47.1*sin(180*PI/180), 47.1*sin(190*PI/180)};

double TOF_Gap8XLoc[3] = {47.1*cos(200*PI/180), 47.1*cos(210*PI/180), 47.1*cos(220*PI/180)};
double TOF_Gap8YLoc[3] = {47.1*sin(200*PI/180), 47.1*sin(210*PI/180), 47.1*sin(220*PI/180)};

double TOF_Gap9XLoc[3] = {47.1*cos(230*PI/180), 47.1*cos(240*PI/180), 47.1*cos(250*PI/180)};
double TOF_Gap9YLoc[3] = {47.1*sin(230*PI/180), 47.1*sin(240*PI/180), 47.1*sin(250*PI/180)};

double TOF_Gap10XLoc[3] = {47.1*cos(260*PI/180), 47.1*cos(270*PI/180), 47.1*cos(280*PI/180)};
double TOF_Gap10YLoc[3] = {47.1*sin(260*PI/180), 47.1*sin(270*PI/180), 47.1*sin(280*PI/180)};

double TOF_Gap11XLoc[3] = {47.1*cos(290*PI/180), 47.1*cos(300*PI/180), 47.1*cos(310*PI/180)};
double TOF_Gap11YLoc[3] = {47.1*sin(290*PI/180), 47.1*sin(300*PI/180), 47.1*sin(310*PI/180)};

double TOF_Gap12XLoc[3] = {47.1*cos(320*PI/180), 47.1*cos(330*PI/180), 47.1*cos(340*PI/180)};
double TOF_Gap12YLoc[3] = {47.1*sin(320*PI/180), 47.1*sin(330*PI/180), 47.1*sin(340*PI/180)};


double sigma_par = sqrt(pow(8.22,2)/12);   // sigma_parallel
double sigma_perp = sqrt(pow(5.,2)/12);    // sigma perpendicular

// AFTER 28/06/2017
double TOF1_Errors_X[12][5] = {sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), 
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180),
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180),
                               sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180),
                               sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180),
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180),
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180),
                               sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180)};


double TOF1_Errors_Y[12][5] = {sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), 
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180),
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180),
                               sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180),
                               sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180),
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180),
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180),
                               sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180)};


// BEFORE 28/06/2017
/*double TOF1_Errors_X[12][3] = {sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180),
                               sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180),
                               0., 0., 0.,
                               sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180),
                               sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180),
                               sigma_par, sigma_par, sigma_par,
                               sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180),
                               sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180),
                               0., 0., 0.,
                               sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180), sigma_par*sin(30*PI/180),
                               sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180), sigma_par*sin(60*PI/180),
                               sigma_par, sigma_par, sigma_par};


double TOF1_Errors_Y[12][3] = {sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180),
                               sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180),
                               sigma_par, sigma_par, sigma_par,
                               sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180),
                               sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180),
                               0., 0., 0.,
                               sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180),
                               sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180),
                               sigma_par, sigma_par, sigma_par,
                               sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180), sigma_par*cos(30*PI/180),
                               sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180), sigma_par*cos(60*PI/180),
                               0., 0., 0.};
  */                   

double TARGET_Errors_X = sqrt(pow(3.2,2)/12); 

double TARGET_Errors_Y = sqrt(pow(3.2,2)/12);

// R = 44.5mm
//float TOF_Xloc[36] = {15.220,22.25,28.604,34.089,38.538,41.816,43.824,44.5,43.824,41.816,38.538,34.089,
//28.604,22.25,15.220,7.727,0,-7.727,-15.220,-22.25,-28.604,-34.089,-38.538,-41.816,
//-43.824,-44.5,-43.824,-41.816,-38.538,-34.089,-28.604,-22.25,-15.220,-7.727,0,7.727};

// R = 47.1mm
/*float TOF_Xloc[36] = {16.109,23.55,30.275,36.081,40.790,44.259,46.384,47.1,46.384,44.259,40.790,36.081,
30.275,23.55,16.109,8.179,0.,-8.179,-16.109,-23.55,-30.275,-36.081,-40.790,-44.259,-46.384,-47.1,-46.384,
-44.259,-40.790,-36.081,-30.275,-23.55,-16.109,-8.179,0,8.179};
*/

//const double TOF1_length = 25.;
const double TOF1_length = 25.5;

double TOF_Xloc[36] = {47.1*cos(60*PI/180) - (TOF1_length/3.)*sin(60*PI/180), 47.1*cos(60*PI/180), 47.1*cos(60*PI/180) + (TOF1_length/3.)*sin(60*PI/180),
                       47.1*cos(30*PI/180) - (TOF1_length/3.)*sin(30*PI/180), 47.1*cos(30*PI/180), 47.1*cos(30*PI/180) + (TOF1_length/3.)*sin(30*PI/180),
                       47.1*cos(0*PI/180) - (TOF1_length/3.)*sin(0*PI/180), 47.1*cos(0*PI/180), 47.1*cos(0*PI/180) + (TOF1_length/3.)*sin(0*PI/180),
                       47.1*cos(330*PI/180) - (TOF1_length/3.)*sin(330*PI/180), 47.1*cos(330*PI/180), 47.1*cos(330*PI/180) + (TOF1_length/3.)*sin(330*PI/180),
                       47.1*cos(300*PI/180) - (TOF1_length/3.)*sin(300*PI/180), 47.1*cos(300*PI/180), 47.1*cos(300*PI/180) + (TOF1_length/3.)*sin(300*PI/180),
                       47.1*cos(270*PI/180) - (TOF1_length/3.)*sin(270*PI/180), 47.1*cos(270*PI/180), 47.1*cos(270*PI/180) + (TOF1_length/3.)*sin(270*PI/180),
                       47.1*cos(240*PI/180) - (TOF1_length/3.)*sin(240*PI/180), 47.1*cos(240*PI/180), 47.1*cos(240*PI/180) + (TOF1_length/3.)*sin(240*PI/180),
                       47.1*cos(210*PI/180) - (TOF1_length/3.)*sin(210*PI/180), 47.1*cos(210*PI/180), 47.1*cos(210*PI/180) + (TOF1_length/3.)*sin(210*PI/180),
                       47.1*cos(180*PI/180) - (TOF1_length/3.)*sin(180*PI/180), 47.1*cos(180*PI/180), 47.1*cos(180*PI/180) + (TOF1_length/3.)*sin(180*PI/180),
                       47.1*cos(150*PI/180) - (TOF1_length/3.)*sin(150*PI/180), 47.1*cos(150*PI/180), 47.1*cos(150*PI/180) + (TOF1_length/3.)*sin(150*PI/180),
                       47.1*cos(120*PI/180) - (TOF1_length/3.)*sin(120*PI/180), 47.1*cos(120*PI/180), 47.1*cos(120*PI/180) + (TOF1_length/3.)*sin(120*PI/180),
                       47.1*cos(90*PI/180) - (TOF1_length/3.)*sin(90*PI/180), 47.1*cos(90*PI/180), 47.1*cos(90*PI/180) + (TOF1_length/3.)*sin(90*PI/180)};


double TOF_Yloc[36] = {47.1*sin(60*PI/180) + (TOF1_length/3.)*cos(60*PI/180), 47.1*sin(60*PI/180), 47.1*sin(60*PI/180) - (TOF1_length/3.)*cos(60*PI/180),
                       47.1*sin(30*PI/180) + (TOF1_length/3.)*cos(30*PI/180), 47.1*sin(30*PI/180), 47.1*sin(30*PI/180) - (TOF1_length/3.)*cos(30*PI/180),
                       47.1*sin(0*PI/180) + (TOF1_length/3.)*cos(0*PI/180), 47.1*sin(0*PI/180), 47.1*sin(0*PI/180) - (TOF1_length/3.)*cos(0*PI/180),
                       47.1*sin(330*PI/180) + (TOF1_length/3.)*cos(330*PI/180), 47.1*sin(330*PI/180), 47.1*sin(330*PI/180) - (TOF1_length/3.)*cos(330*PI/180),
                       47.1*sin(300*PI/180) + (TOF1_length/3.)*cos(300*PI/180), 47.1*sin(300*PI/180), 47.1*sin(300*PI/180) - (TOF1_length/3.)*cos(300*PI/180),
                       47.1*sin(270*PI/180) + (TOF1_length/3.)*cos(270*PI/180), 47.1*sin(270*PI/180), 47.1*sin(270*PI/180) - (TOF1_length/3.)*cos(270*PI/180),
                       47.1*sin(240*PI/180) + (TOF1_length/3.)*cos(240*PI/180), 47.1*sin(240*PI/180), 47.1*sin(240*PI/180) - (TOF1_length/3.)*cos(240*PI/180),
                       47.1*sin(210*PI/180) + (TOF1_length/3.)*cos(210*PI/180), 47.1*sin(210*PI/180), 47.1*sin(210*PI/180) - (TOF1_length/3.)*cos(210*PI/180),
                       47.1*sin(180*PI/180) + (TOF1_length/3.)*cos(180*PI/180), 47.1*sin(180*PI/180), 47.1*sin(180*PI/180) - (TOF1_length/3.)*cos(180*PI/180),
                       47.1*sin(150*PI/180) + (TOF1_length/3.)*cos(150*PI/180), 47.1*sin(150*PI/180), 47.1*sin(150*PI/180) - (TOF1_length/3.)*cos(150*PI/180),
                       47.1*sin(120*PI/180) + (TOF1_length/3.)*cos(120*PI/180), 47.1*sin(120*PI/180), 47.1*sin(120*PI/180) - (TOF1_length/3.)*cos(120*PI/180),
                       47.1*sin(90*PI/180) + (TOF1_length/3.)*cos(90*PI/180), 47.1*sin(90*PI/180), 47.1*sin(90*PI/180) - (TOF1_length/3.)*cos(90*PI/180)};



double TEST_X[12][3] = {47.1*cos(60*PI/180) - (TOF1_length/3.)*sin(60*PI/180), 47.1*cos(60*PI/180), 47.1*cos(60*PI/180) + (TOF1_length/3.)*sin(60*PI/180),
                           47.1*cos(30*PI/180) - (TOF1_length/3.)*sin(30*PI/180), 47.1*cos(30*PI/180), 47.1*cos(30*PI/180) + (TOF1_length/3.)*sin(30*PI/180),
                           47.1*cos(0*PI/180) - (TOF1_length/3.)*sin(0*PI/180), 47.1*cos(0*PI/180), 47.1*cos(0*PI/180) + (TOF1_length/3.)*sin(0*PI/180),
                           47.1*cos(330*PI/180) - (TOF1_length/3.)*sin(330*PI/180), 47.1*cos(330*PI/180), 47.1*cos(330*PI/180) + (TOF1_length/3.)*sin(330*PI/180),
                           47.1*cos(300*PI/180) - (TOF1_length/3.)*sin(300*PI/180), 47.1*cos(300*PI/180), 47.1*cos(300*PI/180) + (TOF1_length/3.)*sin(300*PI/180),
                           47.1*cos(270*PI/180) - (TOF1_length/3.)*sin(270*PI/180), 47.1*cos(270*PI/180), 47.1*cos(270*PI/180) + (TOF1_length/3.)*sin(270*PI/180),
                           47.1*cos(240*PI/180) - (TOF1_length/3.)*sin(240*PI/180), 47.1*cos(240*PI/180), 47.1*cos(240*PI/180) + (TOF1_length/3.)*sin(240*PI/180),
                           47.1*cos(210*PI/180) - (TOF1_length/3.)*sin(210*PI/180), 47.1*cos(210*PI/180), 47.1*cos(210*PI/180) + (TOF1_length/3.)*sin(210*PI/180),
                           47.1*cos(180*PI/180) - (TOF1_length/3.)*sin(180*PI/180), 47.1*cos(180*PI/180), 47.1*cos(180*PI/180) + (TOF1_length/3.)*sin(180*PI/180),
                           47.1*cos(150*PI/180) - (TOF1_length/3.)*sin(150*PI/180), 47.1*cos(150*PI/180), 47.1*cos(150*PI/180) + (TOF1_length/3.)*sin(150*PI/180),
                           47.1*cos(120*PI/180) - (TOF1_length/3.)*sin(120*PI/180), 47.1*cos(120*PI/180), 47.1*cos(120*PI/180) + (TOF1_length/3.)*sin(120*PI/180),
                           47.1*cos(90*PI/180) - (TOF1_length/3.)*sin(90*PI/180), 47.1*cos(90*PI/180), 47.1*cos(90*PI/180) + (TOF1_length/3.)*sin(90*PI/180)};

double TEST_Y[12][3] = {47.1*sin(60*PI/180) + (TOF1_length/3.)*cos(60*PI/180), 47.1*sin(60*PI/180), 47.1*sin(60*PI/180) - (TOF1_length/3.)*cos(60*PI/180),
                           47.1*sin(30*PI/180) + (TOF1_length/3.)*cos(30*PI/180), 47.1*sin(30*PI/180), 47.1*sin(30*PI/180) - (TOF1_length/3.)*cos(30*PI/180),
                           47.1*sin(0*PI/180) + (TOF1_length/3.)*cos(0*PI/180), 47.1*sin(0*PI/180), 47.1*sin(0*PI/180) - (TOF1_length/3.)*cos(0*PI/180),
                           47.1*sin(330*PI/180) + (TOF1_length/3.)*cos(330*PI/180), 47.1*sin(330*PI/180), 47.1*sin(330*PI/180) - (TOF1_length/3.)*cos(330*PI/180),
                           47.1*sin(300*PI/180) + (TOF1_length/3.)*cos(300*PI/180), 47.1*sin(300*PI/180), 47.1*sin(300*PI/180) - (TOF1_length/3.)*cos(300*PI/180),
                           47.1*sin(270*PI/180) + (TOF1_length/3.)*cos(270*PI/180), 47.1*sin(270*PI/180), 47.1*sin(270*PI/180) - (TOF1_length/3.)*cos(270*PI/180),
                           47.1*sin(240*PI/180) + (TOF1_length/3.)*cos(240*PI/180), 47.1*sin(240*PI/180), 47.1*sin(240*PI/180) - (TOF1_length/3.)*cos(240*PI/180),
                           47.1*sin(210*PI/180) + (TOF1_length/3.)*cos(210*PI/180), 47.1*sin(210*PI/180), 47.1*sin(210*PI/180) - (TOF1_length/3.)*cos(210*PI/180),
                           47.1*sin(180*PI/180) + (TOF1_length/3.)*cos(180*PI/180), 47.1*sin(180*PI/180), 47.1*sin(180*PI/180) - (TOF1_length/3.)*cos(180*PI/180),
                           47.1*sin(150*PI/180) + (TOF1_length/3.)*cos(150*PI/180), 47.1*sin(150*PI/180), 47.1*sin(150*PI/180) - (TOF1_length/3.)*cos(150*PI/180),
                           47.1*sin(120*PI/180) + (TOF1_length/3.)*cos(120*PI/180), 47.1*sin(120*PI/180), 47.1*sin(120*PI/180) - (TOF1_length/3.)*cos(120*PI/180),
                           47.1*sin(90*PI/180) + (TOF1_length/3.)*cos(90*PI/180), 47.1*sin(90*PI/180), 47.1*sin(90*PI/180) - (TOF1_length/3.)*cos(90*PI/180)};


double TOF1_Xloc[12][5] = {47.1*cos(60*PI/180)-(2*TOF1_length/5.)*sin(60*PI/180),   47.1*cos(60*PI/180)-(TOF1_length/5.)*sin(60*PI/180),   47.1*cos(60*PI/180),  47.1*cos(60*PI/180)+(TOF1_length/5.)*sin(60*PI/180),   47.1*cos(60*PI/180)+(2*TOF1_length/5.)*sin(60*PI/180),
                        47.1*cos(30*PI/180)-(2*TOF1_length/5.)*sin(30*PI/180),   47.1*cos(30*PI/180)-(TOF1_length/5.)*sin(30*PI/180),   47.1*cos(30*PI/180),  47.1*cos(30*PI/180)+(TOF1_length/5.)*sin(30*PI/180),   47.1*cos(30*PI/180)+(2*TOF1_length/5.)*sin(30*PI/180),
                        47.1*cos(0*PI/180)-(2*TOF1_length/5.)*sin(0*PI/180),     47.1*cos(0*PI/180)-(TOF1_length/5.)*sin(0*PI/180),     47.1*cos(0*PI/180),   47.1*cos(0*PI/180)+(TOF1_length/5.)*sin(0*PI/180),     47.1*cos(0*PI/180)+(2*TOF1_length/5.)*sin(0*PI/180),
                        47.1*cos(330*PI/180)-(2*TOF1_length/5.)*sin(330*PI/180), 47.1*cos(330*PI/180)-(TOF1_length/5.)*sin(330*PI/180), 47.1*cos(330*PI/180), 47.1*cos(330*PI/180)+(TOF1_length/5.)*sin(330*PI/180), 47.1*cos(330*PI/180)+(2*TOF1_length/5.)*sin(330*PI/180),
                        47.1*cos(300*PI/180)-(2*TOF1_length/5.)*sin(300*PI/180), 47.1*cos(300*PI/180)-(TOF1_length/5.)*sin(300*PI/180), 47.1*cos(300*PI/180), 47.1*cos(300*PI/180)+(TOF1_length/5.)*sin(300*PI/180), 47.1*cos(300*PI/180)+(2*TOF1_length/5.)*sin(300*PI/180),
                        47.1*cos(270*PI/180)-(2*TOF1_length/5.)*sin(270*PI/180), 47.1*cos(270*PI/180)-(TOF1_length/5.)*sin(270*PI/180), 47.1*cos(270*PI/180), 47.1*cos(270*PI/180)+(TOF1_length/5.)*sin(270*PI/180), 47.1*cos(270*PI/180)+(2*TOF1_length/5.)*sin(270*PI/180),
                        47.1*cos(240*PI/180)-(2*TOF1_length/5.)*sin(240*PI/180), 47.1*cos(240*PI/180)-(TOF1_length/5.)*sin(240*PI/180), 47.1*cos(240*PI/180), 47.1*cos(240*PI/180)+(TOF1_length/5.)*sin(240*PI/180), 47.1*cos(240*PI/180)+(2*TOF1_length/5.)*sin(240*PI/180),
                        47.1*cos(210*PI/180)-(2*TOF1_length/5.)*sin(210*PI/180), 47.1*cos(210*PI/180)-(TOF1_length/5.)*sin(210*PI/180), 47.1*cos(210*PI/180), 47.1*cos(210*PI/180)+(TOF1_length/5.)*sin(210*PI/180), 47.1*cos(210*PI/180)+(2*TOF1_length/5.)*sin(210*PI/180),
                        47.1*cos(180*PI/180)-(2*TOF1_length/5.)*sin(180*PI/180), 47.1*cos(180*PI/180)-(TOF1_length/5.)*sin(180*PI/180), 47.1*cos(180*PI/180), 47.1*cos(180*PI/180)+(TOF1_length/5.)*sin(180*PI/180), 47.1*cos(180*PI/180)+(2*TOF1_length/5.)*sin(180*PI/180),
                        47.1*cos(150*PI/180)-(2*TOF1_length/5.)*sin(150*PI/180), 47.1*cos(150*PI/180)-(TOF1_length/5.)*sin(150*PI/180), 47.1*cos(150*PI/180), 47.1*cos(150*PI/180)+(TOF1_length/5.)*sin(150*PI/180), 47.1*cos(150*PI/180)+(2*TOF1_length/5.)*sin(150*PI/180),
                        47.1*cos(120*PI/180)-(2*TOF1_length/5.)*sin(120*PI/180), 47.1*cos(120*PI/180)-(TOF1_length/5.)*sin(120*PI/180), 47.1*cos(120*PI/180), 47.1*cos(120*PI/180)+(TOF1_length/5.)*sin(120*PI/180), 47.1*cos(120*PI/180)+(2*TOF1_length/5.)*sin(120*PI/180),
                        47.1*cos(90*PI/180)-(2*TOF1_length/5.)*sin(90*PI/180),   47.1*cos(90*PI/180)-(TOF1_length/5.)*sin(90*PI/180),   47.1*cos(90*PI/180),  47.1*cos(90*PI/180)+(TOF1_length/5.)*sin(90*PI/180),   47.1*cos(90*PI/180)+(2*TOF1_length/5.)*sin(90*PI/180)};


double TOF1_Yloc[12][5] = {47.1*sin(60*PI/180)+(2*TOF1_length/5.)*cos(60*PI/180),   47.1*sin(60*PI/180)+(TOF1_length/5.)*cos(60*PI/180),   47.1*sin(60*PI/180),  47.1*sin(60*PI/180)-(TOF1_length/5.)*cos(60*PI/180),   47.1*sin(60*PI/180)-(2*TOF1_length/5.)*cos(60*PI/180),
                        47.1*sin(30*PI/180)+(2*TOF1_length/5.)*cos(30*PI/180),   47.1*sin(30*PI/180)+(TOF1_length/5.)*cos(30*PI/180),   47.1*sin(30*PI/180),  47.1*sin(30*PI/180)-(TOF1_length/5.)*cos(30*PI/180),   47.1*sin(30*PI/180)-(2*TOF1_length/5.)*cos(30*PI/180),
                        47.1*sin(0*PI/180)+(2*TOF1_length/5.)*cos(0*PI/180),     47.1*sin(0*PI/180)+(TOF1_length/5.)*cos(0*PI/180),     47.1*sin(0*PI/180),   47.1*sin(0*PI/180)-(TOF1_length/5.)*cos(0*PI/180),     47.1*sin(0*PI/180)-(2*TOF1_length/5.)*cos(0*PI/180),
                        47.1*sin(330*PI/180)+(2*TOF1_length/5.)*cos(330*PI/180), 47.1*sin(330*PI/180)+(TOF1_length/5.)*cos(330*PI/180), 47.1*sin(330*PI/180), 47.1*sin(330*PI/180)-(TOF1_length/5.)*cos(330*PI/180), 47.1*sin(330*PI/180)-(2*TOF1_length/5.)*cos(330*PI/180),
                        47.1*sin(300*PI/180)+(2*TOF1_length/5.)*cos(300*PI/180), 47.1*sin(300*PI/180)+(TOF1_length/5.)*cos(300*PI/180), 47.1*sin(300*PI/180), 47.1*sin(300*PI/180)-(TOF1_length/5.)*cos(300*PI/180), 47.1*sin(300*PI/180)-(2*TOF1_length/5.)*cos(300*PI/180),
                        47.1*sin(270*PI/180)+(2*TOF1_length/5.)*cos(270*PI/180), 47.1*sin(270*PI/180)+(TOF1_length/5.)*cos(270*PI/180), 47.1*sin(270*PI/180), 47.1*sin(270*PI/180)-(TOF1_length/5.)*cos(270*PI/180), 47.1*sin(270*PI/180)-(2*TOF1_length/5.)*cos(270*PI/180),
                        47.1*sin(240*PI/180)+(2*TOF1_length/5.)*cos(240*PI/180), 47.1*sin(240*PI/180)+(TOF1_length/5.)*cos(240*PI/180), 47.1*sin(240*PI/180), 47.1*sin(240*PI/180)-(TOF1_length/5.)*cos(240*PI/180), 47.1*sin(240*PI/180)-(2*TOF1_length/5.)*cos(240*PI/180),
                        47.1*sin(210*PI/180)+(2*TOF1_length/5.)*cos(210*PI/180), 47.1*sin(210*PI/180)+(TOF1_length/5.)*cos(210*PI/180), 47.1*sin(210*PI/180), 47.1*sin(210*PI/180)-(TOF1_length/5.)*cos(210*PI/180), 47.1*sin(210*PI/180)-(2*TOF1_length/5.)*cos(210*PI/180),
                        47.1*sin(180*PI/180)+(2*TOF1_length/5.)*cos(180*PI/180), 47.1*sin(180*PI/180)+(TOF1_length/5.)*cos(180*PI/180), 47.1*sin(180*PI/180), 47.1*sin(180*PI/180)-(TOF1_length/5.)*cos(180*PI/180), 47.1*sin(180*PI/180)-(2*TOF1_length/5.)*cos(180*PI/180),
                        47.1*sin(150*PI/180)+(2*TOF1_length/5.)*cos(150*PI/180), 47.1*sin(150*PI/180)+(TOF1_length/5.)*cos(150*PI/180), 47.1*sin(150*PI/180), 47.1*sin(150*PI/180)-(TOF1_length/5.)*cos(150*PI/180), 47.1*sin(150*PI/180)-(2*TOF1_length/5.)*cos(150*PI/180),
                        47.1*sin(120*PI/180)+(2*TOF1_length/5.)*cos(120*PI/180), 47.1*sin(120*PI/180)+(TOF1_length/5.)*cos(120*PI/180), 47.1*sin(120*PI/180), 47.1*sin(120*PI/180)-(TOF1_length/5.)*cos(120*PI/180), 47.1*sin(120*PI/180)-(2*TOF1_length/5.)*cos(120*PI/180),
                        47.1*sin(90*PI/180)+(2*TOF1_length/5.)*cos(90*PI/180),   47.1*sin(90*PI/180)+(TOF1_length/5.)*cos(90*PI/180),   47.1*sin(90*PI/180),  47.1*sin(90*PI/180)-(TOF1_length/5.)*cos(90*PI/180),   47.1*sin(90*PI/180)-(2*TOF1_length/5.)*cos(90*PI/180)};


double X_TOF1_line[12][2] = {47.1*cos(60*PI/180)  - 0.5*TOF1_length*sin(60*PI/180),   47.1*cos(60*PI/180)  + 0.5*TOF1_length*sin(60*PI/180),
                             47.1*cos(30*PI/180)  - 0.5*TOF1_length*sin(30*PI/180),   47.1*cos(30*PI/180)  + 0.5*TOF1_length*sin(30*PI/180),
                             47.1*cos(0*PI/180)   - 0.5*TOF1_length*sin(0*PI/180),    47.1*cos(0*PI/180)   + 0.5*TOF1_length*sin(0*PI/180),
                             47.1*cos(330*PI/180) - 0.5*TOF1_length*sin(330*PI/180),  47.1*cos(330*PI/180) + 0.5*TOF1_length*sin(330*PI/180),
                             47.1*cos(300*PI/180) - 0.5*TOF1_length*sin(300*PI/180),  47.1*cos(300*PI/180) + 0.5*TOF1_length*sin(300*PI/180),
                             47.1*cos(270*PI/180) - 0.5*TOF1_length*sin(270*PI/180),  47.1*cos(270*PI/180) + 0.5*TOF1_length*sin(270*PI/180),
                             47.1*cos(240*PI/180) - 0.5*TOF1_length*sin(240*PI/180),  47.1*cos(240*PI/180) + 0.5*TOF1_length*sin(240*PI/180),
                             47.1*cos(210*PI/180) - 0.5*TOF1_length*sin(210*PI/180),  47.1*cos(210*PI/180) + 0.5*TOF1_length*sin(210*PI/180),
                             47.1*cos(180*PI/180) - 0.5*TOF1_length*sin(180*PI/180),  47.1*cos(180*PI/180) + 0.5*TOF1_length*sin(180*PI/180),
                             47.1*cos(150*PI/180) - 0.5*TOF1_length*sin(150*PI/180),  47.1*cos(150*PI/180) + 0.5*TOF1_length*sin(150*PI/180),
                             47.1*cos(120*PI/180) - 0.5*TOF1_length*sin(120*PI/180),  47.1*cos(120*PI/180) + 0.5*TOF1_length*sin(120*PI/180),
                             47.1*cos(90*PI/180)  - 0.5*TOF1_length*sin(90*PI/180),   47.1*cos(90*PI/180)  + 0.5*TOF1_length*sin(90*PI/180)};

double Y_TOF1_line[12][2] = {47.1*sin(60*PI/180)  + 0.5*TOF1_length*cos(60*PI/180),   47.1*sin(60*PI/180)  - 0.5*TOF1_length*cos(60*PI/180),
                             47.1*sin(30*PI/180)  + 0.5*TOF1_length*cos(30*PI/180),   47.1*sin(30*PI/180)  - 0.5*TOF1_length*cos(30*PI/180),
                             47.1*sin(0*PI/180)   + 0.5*TOF1_length*cos(0*PI/180),    47.1*sin(0*PI/180)   - 0.5*TOF1_length*cos(0*PI/180),
                             47.1*sin(330*PI/180) + 0.5*TOF1_length*cos(330*PI/180),  47.1*sin(330*PI/180) - 0.5*TOF1_length*cos(330*PI/180),
                             47.1*sin(300*PI/180) + 0.5*TOF1_length*cos(300*PI/180),  47.1*sin(300*PI/180) - 0.5*TOF1_length*cos(300*PI/180),
                             47.1*sin(270*PI/180) + 0.5*TOF1_length*cos(270*PI/180),  47.1*sin(270*PI/180) - 0.5*TOF1_length*cos(270*PI/180),
                             47.1*sin(240*PI/180) + 0.5*TOF1_length*cos(240*PI/180),  47.1*sin(240*PI/180) - 0.5*TOF1_length*cos(240*PI/180),
                             47.1*sin(210*PI/180) + 0.5*TOF1_length*cos(210*PI/180),  47.1*sin(210*PI/180) - 0.5*TOF1_length*cos(210*PI/180),
                             47.1*sin(180*PI/180) + 0.5*TOF1_length*cos(180*PI/180),  47.1*sin(180*PI/180) - 0.5*TOF1_length*cos(180*PI/180),
                             47.1*sin(150*PI/180) + 0.5*TOF1_length*cos(150*PI/180),  47.1*sin(150*PI/180) - 0.5*TOF1_length*cos(150*PI/180),
                             47.1*sin(120*PI/180) + 0.5*TOF1_length*cos(120*PI/180),  47.1*sin(120*PI/180) - 0.5*TOF1_length*cos(120*PI/180),
                             47.1*sin(90*PI/180)  + 0.5*TOF1_length*cos(90*PI/180),   47.1*sin(90*PI/180)  - 0.5*TOF1_length*cos(90*PI/180)};


//double X_line_test[1][2] = {47.1*cos(0*PI/180)-0.5*TOF1_length*sin(0*PI/180), 47.1*cos(0*PI/180)+0.5*TOF1_length*sin(0*PI/180)};
//double Y_line_test[1][2] = {47.1*sin(0*PI/180)+0.5*TOF1_length*cos(0*PI/180), 47.1*sin(0*PI/180)-0.5*TOF1_length*cos(0*PI/180)};
double X_line_test[1][2] = {47.1*cos(0*PI/180)-0.5*25.*sin(0*PI/180), 47.1*cos(0*PI/180)+0.5*25.*sin(0*PI/180)};
double Y_line_test[1][2] = {47.1*sin(0*PI/180)+0.5*25.*cos(0*PI/180), 47.1*sin(0*PI/180)-0.5*25.*cos(0*PI/180)};
TLine *Gap_test = new TLine(X_line_test[0][0], Y_line_test[0][0], X_line_test[0][1], Y_line_test[0][1]);


TLine *Gap1l = new TLine(X_TOF1_line[0][0],Y_TOF1_line[0][0],X_TOF1_line[0][1],Y_TOF1_line[0][1]);         //  12.5
TLine *Gap2l = new TLine(X_TOF1_line[1][0],Y_TOF1_line[1][0],X_TOF1_line[1][1],Y_TOF1_line[1][1]);         //  12.5
TLine *Gap3l = new TLine(X_TOF1_line[2][0],Y_TOF1_line[2][0],X_TOF1_line[2][1],Y_TOF1_line[2][1]);         //  12.5
TLine *Gap4l = new TLine(X_TOF1_line[3][0],Y_TOF1_line[3][0],X_TOF1_line[3][1],Y_TOF1_line[3][1]);         //  12.5
TLine *Gap5l = new TLine(X_TOF1_line[4][0],Y_TOF1_line[4][0],X_TOF1_line[4][1],Y_TOF1_line[4][1]);         //  12.5
TLine *Gap6l = new TLine(X_TOF1_line[5][0],Y_TOF1_line[5][0],X_TOF1_line[5][1],Y_TOF1_line[5][1]);         //  12.5
TLine *Gap7l = new TLine(X_TOF1_line[6][0],Y_TOF1_line[6][0],X_TOF1_line[6][1],Y_TOF1_line[6][1]);         //  12.5
TLine *Gap8l = new TLine(X_TOF1_line[7][0],Y_TOF1_line[7][0],X_TOF1_line[7][1],Y_TOF1_line[7][1]);         //  12.5
TLine *Gap9l = new TLine(X_TOF1_line[8][0],Y_TOF1_line[8][0],X_TOF1_line[8][1],Y_TOF1_line[8][1]);         //  12.5
TLine *Gap10l = new TLine(X_TOF1_line[9][0],Y_TOF1_line[9][0],X_TOF1_line[9][1],Y_TOF1_line[9][1]);        //  12.5
TLine *Gap11l = new TLine(X_TOF1_line[10][0],Y_TOF1_line[10][0],X_TOF1_line[10][1],Y_TOF1_line[10][1]);    //  12.5
TLine *Gap12l = new TLine(X_TOF1_line[11][0],Y_TOF1_line[11][0],X_TOF1_line[11][1],Y_TOF1_line[11][1]);    //  12.5

//TLine *Gap1l = new TLine(12.725, 47.04, 34.375, 34.54);                //  12.5
//TLine *Gap2l = new TLine(34.45, 34.375, 47.04, 12.725);                //  12.5
//TLine *Gap3l = new TLine(47.1, -12.5, 47.1, 12.5);                     //  12.5
//TLine *Gap4l = new TLine(34.45, -34.375, 47.04, -12.725);              //  12.5
//TLine *Gap5l = new TLine(12.725, -47.04, 34.375, -34.54);              //  12.5
//TLine *Gap6l = new TLine(-12.5, -47.1, 12.5, -47.1);                   //  12.5
//TLine *Gap7l = new TLine(-12.725, -47.04, -34.375, -34.54);            //  12.5
//TLine *Gap8l = new TLine(-34.45, -34.375, -47.04, -12.725);            //  12.5
//TLine *Gap9l = new TLine(-47.1, -12.5, -47.1, 12.5);                   //  12.5
//TLine *Gap10l = new TLine(-34.45, 34.375, -47.04, 12.725);             //  12.5
//TLine *Gap11l = new TLine(-12.725, 47.04, -34.375, 34.54);             //  12.5
//TLine *Gap12l = new TLine(-12.5, 47.1, 12.5, 47.1);                    //  12.5

const double TOF1_line_slope[12]={(Y_TOF1_line[0][0]-Y_TOF1_line[0][1])/(X_TOF1_line[0][0]-X_TOF1_line[0][1]),
                                  (Y_TOF1_line[1][0]-Y_TOF1_line[1][1])/(X_TOF1_line[1][0]-X_TOF1_line[1][1]),
                                  999.,
                                  (Y_TOF1_line[3][0]-Y_TOF1_line[3][1])/(X_TOF1_line[3][0]-X_TOF1_line[3][1]),
                                  (Y_TOF1_line[4][0]-Y_TOF1_line[4][1])/(X_TOF1_line[4][0]-X_TOF1_line[4][1]),
                                  0.,
                                  (Y_TOF1_line[6][0]-Y_TOF1_line[6][1])/(X_TOF1_line[6][0]-X_TOF1_line[6][1]),
                                  (Y_TOF1_line[7][0]-Y_TOF1_line[7][1])/(X_TOF1_line[7][0]-X_TOF1_line[7][1]),
                                  999.,
                                  (Y_TOF1_line[9][0]-Y_TOF1_line[9][1])/(X_TOF1_line[9][0]-X_TOF1_line[9][1]),
                                  (Y_TOF1_line[10][0]-Y_TOF1_line[10][1])/(X_TOF1_line[10][0]-X_TOF1_line[10][1]),
                                  0.};

/*int channel[12][8]={5,13,14,15,26,27,40,41,
27,40,41,56,57,73,90,91,
90,91,109,127,145,163,180,181,
180,181,197,212,213,226,227,239,
226,227,238,239,247,248,249,255,
242,247,250,251,252,253,254,255,
214,215,228,229,240,241,242,250,
164,165,182,198,199,214,215,228,
74,75,92,110,128,146,164,165,
16,28,29,42,43,58,74,75,
0,6,7,8,16,17,28,29,
0,1,2,3,4,5,8,13};
*/

int channel[12][12]={  5, 13, 14, 15, 24, 25, 26, 27, 38, 39, 40, 41,
                      27, 40, 41, 55, 56, 57, 71, 72, 73, 89, 90, 91,
                      90, 91,108,109,126,127,144,145,162,163,180,181,
                     179,180,181,195,196,197,211,212,213,226,227,239,
                     224,225,226,227,236,237,238,239,247,248,249,255,
                     242,243,244,245,246,247,250,251,252,253,254,255,
                     214,215,216,217,228,229,230,231,240,241,242,250,
                     164,165,166,182,183,184,198,199,200,214,215,228,
                      74, 75, 92, 93,110,111,128,129,146,147,164,165,
                      16, 28, 29, 42, 43, 44, 58, 59, 60, 74, 75, 76,
                       0,  6,  7,  8, 16, 17, 18, 19, 28, 29, 30, 31,
                       0,  1,  2,  3,  4,  5,  8,  9, 10, 11, 12, 13};           


int TOF1_rotated[12] = {9,10,11,0,1,2,3,4,5,6,7,8};

// Rotates the target by -90 deg.
Int_t TARGET_Rotated_index[256] = {
            164,146,128,110, 92,74,
        198,182,165,147,129,111, 93,75,58,42,
      214,199,183,166,148,130,112, 94,76,59,43,28,
    228,215,200,184,167,149,131,113, 95,77,60,44,29,16,
  240,229,216,201,185,168,150,132,114, 96,78,61,45,30,17, 6,
  241,230,217,202,186,169,151,133,115, 97,79,62,46,31,18, 7,
250,242,231,218,203,187,170,152,134,116, 98,80,63,47,32,19, 8,0,
251,243,232,219,204,188,171,153,135,117, 99,81,64,48,33,20, 9,1,
252,244,233,220,205,189,172,154,136,118,100,82,65,49,34,21,10,2,
253,245,234,221,206,190,173,155,137,119,101,83,66,50,35,22,11,3,
254,246,235,222,207,191,174,156,138,120,102,84,67,51,36,23,12,4,
255,247,236,223,208,192,175,157,139,121,103,85,68,52,37,24,13,5,
  248,237,224,209,193,176,158,140,122,104,86,69,53,38,25,14,
  249,238,225,210,194,177,159,141,123,105,87,70,54,39,26,15,
    239,226,211,195,178,160,142,124,106,88,71,55,40,27,
      227,212,196,179,161,143,125,107,89,72,56,41,
        213,197,180,162,144,126,108,90,73,57,
            181,163,145,127,109,91};
            
// Rotates the target by +90 deg.           
Int_t TARGET_Rotated_index_inverse[256] = {  
           91,109,127,145,163,181,
       57,73,90,108,126,144,162,180,197,213,
    41,56,72,89,107,125,143,161,179,196,212,227,
   27,40,55,71,88,106,124,142,160,178,195,211,226,239,
  15,26,39,54,70,87,105,123,141,159,177,194,210,225,238,249,
  14,25,38,53,69,86,104,122,140,158,176,193,209,224,237,248,
5,13,24,37,52,68,85,103,121,139,157,175,192,208,223,236,247,255,
4,12,23,36,51,67,84,102,120,138,156,174,191,207,222,235,246,254,
3,11,22,35,50,66,83,101,119,137,155,173,190,206,221,234,245,253,
2,10,21,34,49,65,82,100,118,136,154,172,189,205,220,233,244,252,
1, 9,20,33,48,64,81, 99,117,135,153,171,188,204,219,232,243,251,
0, 8,19,32,47,63,80, 98,116,134,152,170,187,203,218,231,242,250,
   7,18,31,46,62,79, 97,115,133,151,169,186,202,217,230,241,
   6,17,30,45,61,78, 96,114,132,150,168,185,201,216,229,240,
   16,29,44,60,77, 95,113,131,149,167,184,200,215,228,
      28,43,59,76, 94,112,130,148,166,183,199,214,
       42,58,75, 93,111,129,147,165,182,198,
         74, 92,110,128,146,164};

Int_t SFT_channel_to_layer[128] =
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 3, 3, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3
};

Int_t SFT_channel_to_fiber[128] =
{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
  11, 12, 13, 14, 15, 1, 2, 3, 4,
  5, 6, 7, 8, 9, 10, 11, 12, 13,
  14, 15, 1, 2, 17, 16, 15, 14,
  13, 12, 11, 10, 9, 8, 7, 6, 5,
  4, 3, 2, 1, 17, 16, 15, 14, 13,
  12, 11, 10, 9, 8, 7, 6, 5, 4,
  3, 1, 2, 3, 4, 5, 6, 7, 8, 9,
  10, 11, 12, 13, 14, 15,  1, 2,
  3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
  13, 14, 15, 1, 2, 17, 16, 15,
  14, 13, 12, 11, 10, 9, 8, 7, 6,
  5, 4, 3, 2, 1, 17, 16, 15, 14, 13,
  12, 11, 10, 9, 8, 7, 6, 5, 4, 3
};

char SFT_channel_to_US_DS[128] =
{
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
   'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U',
};


Int_t TARGET_neighbours[256][8] =
{
  -1,-1,-1,1,9,8,7,-1, -1,-1,-1,2,10,9,8,0, -1,-1,-1,3,11,10,9,1, -1,-1,-1,4,12,11,10,2, -1,-1,-1,5,13,12,11,3, // 5
    -1,-1,-1,-1,14,13,12,4, -1,-1,-1,7,18,17,16,-1, -1,-1,0,8,19,18,17,6, -1,0,1,9,20,19,18,7, 0,1,2,10,21,20,19,8,// 10
    1,2,3,11,22,21,20,9, 2,3,4,12,23,22,21,10, 3,4,5,13,24,23,22,11, 4,5,-1,14,25,24,23,12, 5,-1,-1,15,26,25,24,13, //15
    -1,-1,-1,-1,27,26,25,14, -1,-1,6,17,30,29,28,-1, -1,6,7,18,31,30,29,16, 6,7,8,19,32,31,30,15, 7,8,9,20,33,32,31,16, //20
    8,9,10,21,34,33,32,19, 9,10,11,22,35,34,33,20, 10,11,12,23,36,35,34,21, 11,12,13,24,37,36,35,22, 12,13,14,25,38,37,36,23, // 25 
    13,14,15,26,39,38,37,24, 14,15,-1,27,40,39,38,25, 15,-1,-1,-1,41,40,39,26, -1,-1,16,29,44,43,42,-1, -1,16,17,30,45,44,43,28, //30
    16,17,18,31,46,45,44,29, 17,18,19,32,47,46,45,30, 18,19,20,33,48,47,46,31, 19,20,21,34,49,48,47,32, 20,21,22,35,50,49,48,33, // 35
    21,22,23,36,51,50,49,34, 22,23,24,37,52,51,50,35, 23,24,25,38,53,52,51,36, 24,25,26,39,54,53,52,37, 25,26,27,40,55,54,53,38, // 40
    26,27,-1,41,56,55,54,39, 27,-1,-1,-1,57,56,55,40, -1,-1,28,43,59,58,-1,-1, -1,28,29,44,60,59,58,42, 28,29,30,45,61,60,59,43, //45
    29,30,31,46,62,61,60,44, 30,31,32,47,63,62,61,45, 31,32,33,48,64,63,62,46, 32,33,34,49,65,64,63,47, 33,34,35,50,66,65,64,48, //50
    34,35,36,51,67,66,65,49, 35,36,37,52,68,67,66,50, 36,37,38,53,69,68,67,51, 37,38,39,54,70,69,68,52, 38,39,40,55,71,70,69,53, //55
    39,40,41,56,72,71,70,54, 40,41,-1,57,73,72,71,55, 41,-1,-1,-1,-1,73,72,56, -1,42,43,59,76,75,74,-1, 42,43,44,60,77,76,75,58, //60
    43,44,45,61,78,77,76,59, 44,45,46,62,79,78,77,60, 45,46,47,63,80,79,78,61, 46,47,48,64,81,80,79,62, 47,48,49,65,82,81,80,63,//65
    48,49,50,66,83,82,81,64, 49,50,51,67,84,83,82,65, 50,51,52,68,85,84,83,66, 51,52,53,69,86,85,84,67, 52,53,54,70,87,86,85,68,//70
    53,54,55,71,88,87,86,69, 54,55,56,72,89,88,87,70, 55,56,57,73,90,89,88,71, 56,57,-1,-1,91,90,89,72, -1,-1,58,75,93,92,-1,-1,//75
    -1,58,59,76,94,93,92,74, 58,59,60,77,95,94,93,75, 59,60,61,78,96,95,94,76, 60,61,62,79,97,96,95,77, 61,62,63,80,98,97,96,78,//80
    62,63,64,81,99,98,97,79, 63,64,65,82,100,99,98,80, 64,65,66,83,101,100,99,81, 65,66,67,84,102,101,100,82, 66,67,68,85,103,102,101,83, //85
    67,68,69,86,104,103,102,84, 68,69,70,87,105,104,103,85, 69,70,71,88,106,105,104,86, 70,71,72,89,107,106,105,87, 71,72,73,90,108,107,106,88, //90 
    72,73,-1,91,109,108,107,89, 73,-1,-1,-1,-1,109,108,90, -1,74,75,93,111,110,-1,-1, 74,75,76,94,112,111,110,92, 75,76,77,95,113,112,11,93, //95
    76,77,78,96,114,113,112,94, 77,78,79,97,115,114,113,95, 78,79,80,98,116,115,114,96, 79,80,81,99,117,116,115,97, 80,81,82,100,118,117,116,98, //100
    81,82,83,101,119,118,117,99, 82,83,84,102,120,119,118,100, 83,84,85,103,121,120,119,101, 84,85,86,104,122,121,120,102, 85,86,87,105,123,122,121,103, //105
    86,87,88,106,124,123,122,104, 87,88,89,107,125,124,123,105, 88,89,90,108,126,125,124,106, 89,90,91,109,127,126,125,107, 90,91,-1,-1,-1,127,126,108, //110
    -1,92,93,111,129,128,-1,-1, 92,93,94,112,130,129,128,110, 93,94,95,113,131,130,129,111, 94,95,96,114,132,131,130,112, 95,96,97,115,133,132,131,113, //115
    96,97,98,116,134,133,132,114, 97,98,99,117,135,134,133,115, 98,99,100,118,136,135,134,116, 99,100,101,119,137,136,135,117, 100,101,102,120,138,137,136,118, //120 
    101,102,103,121,139,138,137,119, 102,103,104,122,140,139,138,120, 103,104,105,123,141,140,139,121, 104,105,106,124,142,141,140,122, 105,106,107,125,143,142,141,123, //125
    106,107,108,126,144,143,142,124, 107,108,109,127,145,144,143,125, 108,109,-1,-1,-1,145,144,126, -1,110,111,129,147,146,-1,-1, 110,111,112,130,148,147,146,128, //130
    111,112,113,131,149,148,147,129, 112,113,114,132,150,149,148,130, 113,114,115,133,151,150,149,131, 114,115,116,134,152,151,150,132, 115,116,117,135,153,152,151,133, //135
    116,117,118,136,154,153,152,134, 117,118,119,137,155,154,153,135, 118,119,120,138,156,155,154,136, 119,120,121,139,157,156,155,137, 120,121,122,140,158,157,156,138, //140
    121,122,123,141,159,158,157,139, 122,123,124,142,160,159,158,140, 123,124,125,143,161,160,159,141, 124,125,126,144,162,161,160,142, 125,126,127,145,163,162,161,143, //145
    126,127,-1,-1,-1,163,162,144,  -1,128,129,147,165,164,-1,-1, 128,129,130,148,166,165,164,146, 129,130,131,149,167,166,165,147, 130,131,132,150,168,167,166,148, //150
    131,132,133,151,169,168,167,149, 132,133,134,152,170,169,168,150, 133,134,135,153,171,170,169,151, 134,135,136,154,172,171,170,152, 135,136,137,155,173,172,171,153, //155
    136,137,138,156,174,173,172,154, 137,138,139,157,175,174,173,155, 138,139,140,158,176,175,174,156, 139,140,141,159,177,176,175,157, 140,141,142,160,178,177,176,158, //160
    141,142,143,161,179,178,177,159, 142,143,144,162,180,179,178,160, 143,144,145,163,181,180,179,161, 144,145,-1,-1,-1,181,180,162, -1,146,147,165,182,-1,-1,-1, //165
    146,147,148,166,183,182,-1,164, 147,148,149,167,184,183,182,165, 148,149,150,168,185,184,183,166, 149,150,151,169,186,185,184,167, 150,151,152,170,187,186,185,168, //170 
    151,152,153,171,188,187,186,169, 152,153,154,172,189,188,187,170, 153,154,155,173,190,189,188,171, 154,155,156,174,191,190,189,172, 155,156,157,175,192,191,190,173, //175
    156,157,158,176,193,192,191,174, 157,158,159,177,194,193,192,175, 158,169,160,178,195,194,193,176, 159,160,161,179,196,195,194,177, 160,161,162,180,197,196,195,178, //180
    161,162,163,181,-1,197,196,179, 162,163,-1,-1,-1,-1,197,180, 164,165,166,183,199,198,-1,-1, 165,166,167,184,200,199,198,182, 166,167,168,185,201,200,199,183, //185
    167,168,169,186,202,201,200,184, 168,169,170,187,203,202,201,185, 169,170,171,188,204,203,202,186, 170,171,172,189,205,204,203,187, 171,172,173,190,206,205,204,188,  //190 
    172,173,174,191,207,206,205,189, 173,174,175,192,208,207,206,190, 174,175,176,193,209,208,207,191, 175,176,177,194,210,209,208,192, 176,177,178,195,211,210,209,193, //195
    177,178,177,196,212,211,210,194, 178,179,180,197,213,212,211,195, 179,180,181,-1,-1,213,212,196, -1,182,183,199,214,-1,-1,-1, 182,183,184,200,215,214,-1,198,  //200
    183,184,185,201,216,215,214,199, 184,185,186,202,217,216,215,200, 185,186,187,203,218,217,216,201, 186,187,188,204,219,218,217,202, 187,188,189,205,220,219,218,203, //205
    188,189,190,206,221,220,219,204, 189,190,191,207,222,221,220,205, 190,191,192,208,223,222,221,206, 191,192,193,209,224,223,222,207, 192,193,194,210,225,224,223,208, //210
    193,194,195,211,226,225,224,209, 194,195,196,212,227,226,225,210, 195,196,197,213,-1,227,226,211, 196,197,-1,-1,-1,-1,227,212, 198,199,200,215,228,-1,-1,-1, //215
    199,200,201,210,229,228,-1,214, 200,201,202,217,230,229,228,215, 201,202,203,218,231,230,229,210, 202,203,204,219,232,231,230,217, 203,204,205,220,233,232,231,218, //220 
    204,205,206,221,234,233,232,219, 205,206,207,222,235,234,233,220, 206,207,208,223,236,235,234,221, 207,208,209,224,237,236,235,222, 208,209,210,225,238,237,236,223,  //225
    209,210,211,226,239,238,237,224, 210,211,212,227,-1,239,238,225, 211,212,213,-1,-1,-1,239,226, 214,215,216,229,240,-1,-1,-1, 215,216,217,230,241,240,-1,228, //230
    216,217,218,231,242,241,240,229, 217,218,219,232,243,242,241,230, 218,219,220,233,244,243,242,231, 219,220,221,234,245,244,243,232, 220,221,222,235,246,245,244,233,  //235 
    221,222,223,236,247,246,245,234, 222,223,224,237,248,247,246,235, 223,224,225,238,249,248,247,236, 224,225,226,239,-1,249,248,237, 225,226,227,-1,-1,-1,249,238, //240
    228,229,230,241,-1,-1,-1,-1, 229,230,231,242,250,-1,-1,240, 230,231,232,243,251,250,-1,241, 231,232,233,244,252,251,250,242, 232,233,234,245,253,252,251,243,  //245
    233,234,235,246,254,253,252,244, 234,235,246,247,255,254,253,245, 235,236,237,248,-1,255,254,246, 236,237,238,249,-1,-1,255,247, 237,238,239,-1,-1,-1,-1,248, //250
    241,242,243,251,-1,-1,-1,-1, 242,243,244,252,-1,-1,-1,250, 243,244,245,253,-1,-1,-1,251, 244,245,246,254,-1,-1,-1,252, 245,246,247,255,-1,-1,-1,253, //255
    246,247,248,-1,-1,-1,-1,254
};

Int_t C2XL_Channels[56] = {237, 128, 173, 192, 236, 129, 172, 193, 235, 130,
                           171, 194, 234, 131, 170, 195, 233, 132, 169, 196,
                           232, 133, 168, 197, 231, 134, 167, 198, 230, 135, 
                           166, 199, 229, 136, 165, 200, 228, 137, 164, 201, 
                           227, 138, 163, 202, 226, 139, 162, 203, 225, 140,
                           161, 204, 224, 141, 160, 205};  

Int_t C2XR_Channels[56] = {253, 144, 189, 208, 252, 145, 188, 209, 251, 146,
                           187, 210, 250, 147, 186, 211, 249, 148, 185, 212, 
                           248, 149, 184, 213, 247, 150, 183, 214, 246, 151, 
                           182, 215, 245, 152, 181, 216, 244, 153, 180, 217, 
                           243, 154, 179, 218, 242, 155, 178, 219, 241, 156,
                           177, 220, 240, 157, 176, 221};

Int_t C2YL_Channels[16] = {111, 110, 109, 108, 107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96};

Int_t C2YR_Channels[16] = {127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117, 116, 115, 114, 113, 112};



Int_t C3XL_Channels[64] = {64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495};

Int_t C3XR_Channels[64] = {80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511};

Int_t C3YL_Channels[16] = {463, 462, 461, 460, 459, 458, 457, 456, 455, 454, 453, 452, 451, 450, 449, 448};

Int_t C3YR_Channels[16] = {479, 478, 477, 476, 475, 474, 473, 472, 471, 470, 469, 468, 467, 466, 465, 464};



Int_t C4XL_Channels[72] = {416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431,
                           384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399,
                           352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367,
                           320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335,
                           288, 289, 290, 291, 292, 293, 294, 295};

Int_t C4XR_Channels[72] = {432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447,
                           400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415,
                           368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383,
                           336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351,
                           304, 305, 306, 307, 308, 309, 310, 311};

Int_t C4YL_Channels[16] = {271, 270, 269, 268, 267, 266, 265, 264, 263, 262, 261, 260, 259, 258, 257, 256};

Int_t C4YR_Channels[16] = {287, 286, 285, 284, 283, 282, 281, 280, 279, 278, 277, 276, 275, 274, 273, 272};



#endif // COMMONPARAMETERS_H

