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

vector<int> _neighbour_TOF1(int TOF1_gap){
	
	vector<int> TOF1_neighbours;
	int TOF1_minus;
	int TOF1_plus;

	if(TOF1_gap-1 == 0){
		TOF1_minus = 11;
		TOF1_plus = 1;
	}
	else if(TOF1_gap-1 == 11){
		TOF1_minus = 10;
		TOF1_plus = 0;
	}
	else{
		TOF1_minus = (TOF1_gap-1)-1;
		TOF1_plus = (TOF1_gap-1)+1;
	}

	TOF1_neighbours.push_back(TOF1_minus);
	TOF1_neighbours.push_back(TOF1_plus);

	return TOF1_neighbours;
}

vector<double> _pruning_area(double x1, double x2, double x3){

	vector<double> vec_xx_pruning_area;
	vec_xx_pruning_area.clear();

	vec_xx_pruning_area.push_back(x1);
	vec_xx_pruning_area.push_back(x2);
	vec_xx_pruning_area.push_back(x3);
	vec_xx_pruning_area.push_back(x1);

	/// DEBUG
	// cout << "pruning area : " << vec_xx_pruning_area[0] << "  ";
	// cout << vec_xx_pruning_area[1] << "  ";
	// cout << vec_xx_pruning_area[2] << "  ";
	// cout << vec_xx_pruning_area[3] << "  " << endl;

	return vec_xx_pruning_area;
}

vector<int> _Triangular_pruning(bool to_rotate, TGraph *gr_pruning_area, vector<int> leptons_unpruned, vector<int> kaons){

	vector<int> vec_pruned_temp;	vec_pruned_temp.clear();
	vector<int> vec_pruned;			vec_pruned.clear();


		if(kaons.size()==0){
			if(to_rotate){
				for(unsigned int i=0; i<leptons_unpruned.size(); i++){
					vec_pruned.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
				}
			}
			else{
				for(unsigned int i=0; i<leptons_unpruned.size(); i++){
					vec_pruned.push_back(leptons_unpruned[i]);
				}			
			}
		}
		else{
			if(to_rotate){
				for(unsigned int i=0; i<leptons_unpruned.size(); i++){
					if(gr_pruning_area->IsInside(Xloc[TARGET_Rotated_index[leptons_unpruned[i]]], Yloc[TARGET_Rotated_index[leptons_unpruned[i]]])==1){
						vec_pruned_temp.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
					}
				}
			}
			else{
				for(unsigned int i=0; i<leptons_unpruned.size(); i++){
					if(gr_pruning_area->IsInside(Xloc[leptons_unpruned[i]], Yloc[leptons_unpruned[i]])==1){
						vec_pruned_temp.push_back(leptons_unpruned[i]);
					}
				}
			}
		}

		if(vec_pruned_temp.size() >= 2){
			vec_pruned = vec_pruned_temp;
		}
		else{
			if(kaons.size()>0){
				if(to_rotate){
					for(unsigned int i=0; i<leptons_unpruned.size(); i++){
						vec_pruned.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
					}
				}
				else{
					for(unsigned int i=0; i<leptons_unpruned.size(); i++){
						vec_pruned.push_back(leptons_unpruned[i]);
					}			
				}
			}
		}

		return vec_pruned;
}

vector<double> V_X_pruning_area;
vector<double> V_Y_pruning_area;

vector<int> _Triangular_pruning2(int TOF1_Gap, double x_kstop, double y_kstop, 
	                             vector<int> leptons_unpruned, vector<int> kaons){

	vector<int> vec_pruned_temp;	vec_pruned_temp.clear();
	vector<int> vec_pruned;			vec_pruned.clear();
	vector<int> vec_temp;			vec_temp.clear();


	///DETERMINE TOF1 NEIGHBOURS USING _neighbour_TOF1 function
	vector<int> vec_TOF1_neighbours;	vec_TOF1_neighbours.clear();

	// if(to_rotate) vec_TOF1_neighbours = _neighbour_TOF1(TOF1_rotated[gap_to_fit-1]+1);
	// else vec_TOF1_neighbours = _neighbour_TOF1(gap_to_fit);

	//vec_TOF1_neighbours = _neighbour_TOF1(gap_to_fit);
	vec_TOF1_neighbours = _neighbour_TOF1(TOF1_Gap);

	/// DETERMINE PRUNING AREA USING _Ppruning_area function
	//vector<double> V_X_pruning_area; V_X_pruning_area.clear();
	//vector<double> V_Y_pruning_area; V_Y_pruning_area.clear();
	V_X_pruning_area.clear();
	V_Y_pruning_area.clear();

	V_X_pruning_area = _pruning_area(x_kstop, TOF1_Xloc[vec_TOF1_neighbours[0]][2], TOF1_Xloc[vec_TOF1_neighbours[1]][2]);
  	V_Y_pruning_area = _pruning_area(y_kstop, TOF1_Yloc[vec_TOF1_neighbours[0]][2], TOF1_Yloc[vec_TOF1_neighbours[1]][2]);


  	TGraph *gr_pruning_area;
  	gr_pruning_area = new TGraph(V_X_pruning_area.size(), &V_X_pruning_area[0], &V_Y_pruning_area[0]);
  	//   gr_pruning_area->SetLineWidth(2);
  	//   gr_pruning_area->SetLineColor(4);
  	//   gr_pruning_area->SetMarkerStyle(34);
  	//   gr_pruning_area->SetMarkerSize(1.3);
  	//   gr_pruning_area->SetMarkerColor(2);
  	//   gr_pruning_area->Draw("same");

  	/// PRUNING DATA
/*	if(kaons.size()==0){
		if(to_rotate){
			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
				vec_pruned.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
			}
		}
		else{
			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
				vec_pruned.push_back(leptons_unpruned[i]);
			}			
		}
	}
	else{
		if(to_rotate){
			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
				if(gr_pruning_area->IsInside(Xloc[TARGET_Rotated_index[leptons_unpruned[i]]], Yloc[TARGET_Rotated_index[leptons_unpruned[i]]])==1){
					vec_pruned_temp.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
				}
			}
		}
		else{
			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
				if(gr_pruning_area->IsInside(Xloc[leptons_unpruned[i]], Yloc[leptons_unpruned[i]])==1){
					vec_pruned_temp.push_back(leptons_unpruned[i]);
				}
			}
		}
	}
*/
	if(kaons.size()==0){
		// if(to_rotate){
		// 	for(unsigned int i=0; i<leptons_unpruned.size(); i++){
		// 		vec_pruned.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
		// 	}
		// }
		//else{
			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
				vec_pruned.push_back(leptons_unpruned[i]);
			}			
		//}
	}
	else{
		// if(to_rotate){
		// 	for(unsigned int i=0; i<leptons_unpruned.size(); i++){
		// 		if(gr_pruning_area->IsInside(Xloc[TARGET_Rotated_index[leptons_unpruned[i]]], Yloc[TARGET_Rotated_index[leptons_unpruned[i]]])==1){
		// 			vec_pruned_temp.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
		// 		}
		// 	}
		// }
		//else{
			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
				if(gr_pruning_area->IsInside(Xloc[leptons_unpruned[i]], Yloc[leptons_unpruned[i]])==1){
					vec_pruned_temp.push_back(leptons_unpruned[i]);
				}
			}
		//}
	}


	// if(vec_pruned_temp.size() >= 2){
	// 	vec_pruned = vec_pruned_temp;
	// }
	// else{
	// 	if(kaons.size()>0){
	// 		if(to_rotate){
	// 			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
	// 				vec_pruned.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
	// 			}
	// 		}
	// 		else{
	// 			for(unsigned int i=0; i<leptons_unpruned.size(); i++){
	// 				vec_pruned.push_back(leptons_unpruned[i]);
	// 			}			
	// 		}
	// 	}
	// }

	if(vec_pruned_temp.size() >= 2){
		vec_pruned = vec_pruned_temp;
	}
	else{
		if(kaons.size()>0){
			// if(to_rotate){
			// 	for(unsigned int i=0; i<leptons_unpruned.size(); i++){
			// 		vec_pruned.push_back(TARGET_Rotated_index[leptons_unpruned[i]]);
			// 	}
			// }
			//else{
				for(unsigned int i=0; i<leptons_unpruned.size(); i++){
					vec_pruned.push_back(leptons_unpruned[i]);
				}			
			//}
		}
	}

	// if(to_rotate){
	// 	for(unsigned int i=0; i<vec_pruned.size(); i++){
	// 		vec_temp.push_back(TARGET_Rotated_index_inverse[vec_pruned[i]]);
	// 	}
	// 	vec_pruned.clear();
	// 	vec_pruned = vec_temp;
	// }

	vector<double> V_X_pruned;	V_X_pruned.clear();
	vector<double> V_Y_pruned;	V_Y_pruned.clear();

	for(unsigned int i=0; i<vec_pruned.size(); i++){
		// if(to_rotate){
		// 	V_X_pruned.push_back(Xloc[TARGET_Rotated_index[vec_pruned[i]]]);
		// 	V_Y_pruned.push_back(Yloc[TARGET_Rotated_index[vec_pruned[i]]]);
		// }
		//else{
			V_X_pruned.push_back(Xloc[vec_pruned[i]]);
			V_Y_pruned.push_back(Yloc[vec_pruned[i]]);
		//}
	}

  	TGraph *gr_pruned_lepton;
  	gr_pruned_lepton = new TGraph(vec_pruned.size(), &V_X_pruned[0], &V_Y_pruned[0]);
  	// gr_pruned_lepton->SetMarkerStyle(21);
  	// gr_pruned_lepton->SetMarkerSize(0.8);
  	// gr_pruned_lepton->SetMarkerColor(4);
  	// gr_pruned_lepton->Draw("sameP");


	return vec_pruned;
}


vector<double> vec_pruning_line_coordinates;
vector<double> vec_pruning_line;

vector<double> _pruning_line_OLD(bool to_rotate, double x_kstop, double y_kstop, int TOF1_gap){

	double x1, y1;
	double x2, y2;
	double b;
	double a_out = -1.;
	double b_out = -1.;

	//vector<double> vec_pruning_line;
	vec_pruning_line.clear();

	vec_pruning_line_coordinates.clear();

	if(to_rotate) TOF1_gap = TOF1_rotated[TOF1_gap-1]+1;

	//Define min and max
	if(TOF1_gap==3 || TOF1_gap==9){

		x1 = x_kstop;
		y1 = -50.;
		
		x2 = x_kstop;
		y2 = 50.;

		a_out = 999.99;
		b_out = 999.99;
	}
	else{

		b = y_kstop - TOF1_line_slope[TOF1_gap-1]*x_kstop;

		x1 = -50.;
		y1 = TOF1_line_slope[TOF1_gap-1]*x1 + b;

		x2 = 50.;
		y2 = TOF1_line_slope[TOF1_gap-1]*x2 + b;

		a_out = TOF1_line_slope[TOF1_gap-1];
		b_out = b;
	}

	vec_pruning_line_coordinates.push_back(x1);
	vec_pruning_line_coordinates.push_back(y1);
	vec_pruning_line_coordinates.push_back(x2);
	vec_pruning_line_coordinates.push_back(y2);

	vec_pruning_line.push_back(x1);
	vec_pruning_line.push_back(y1);
	vec_pruning_line.push_back(x2);
	vec_pruning_line.push_back(y2);
	vec_pruning_line.push_back(a_out);
	vec_pruning_line.push_back(b_out);

	return vec_pruning_line;
}

vector<double> _pruning_line(double x_kstop, double y_kstop, int TOF1_gap){

	double x1, y1;
	double x2, y2;
	double b;
	double a_out = -1.;
	double b_out = -1.;

	//vector<double> vec_pruning_line;
	vec_pruning_line.clear();

	vec_pruning_line_coordinates.clear();

	//if(to_rotate) TOF1_gap = TOF1_rotated[TOF1_gap-1]+1;

	//Define min and max
	if(TOF1_gap==3 || TOF1_gap==9){

		x1 = x_kstop;
		y1 = -50.;
		
		x2 = x_kstop;
		y2 = 50.;

		a_out = 999.99;
		b_out = 999.99;
	}
	else{

		b = y_kstop - TOF1_line_slope[TOF1_gap-1]*x_kstop;

		x1 = -50.;
		y1 = TOF1_line_slope[TOF1_gap-1]*x1 + b;

		x2 = 50.;
		y2 = TOF1_line_slope[TOF1_gap-1]*x2 + b;

		a_out = TOF1_line_slope[TOF1_gap-1];
		b_out = b;
	}

	vec_pruning_line_coordinates.push_back(x1);
	vec_pruning_line_coordinates.push_back(y1);
	vec_pruning_line_coordinates.push_back(x2);
	vec_pruning_line_coordinates.push_back(y2);

	vec_pruning_line.push_back(x1);
	vec_pruning_line.push_back(y1);
	vec_pruning_line.push_back(x2);
	vec_pruning_line.push_back(y2);
	vec_pruning_line.push_back(a_out);
	vec_pruning_line.push_back(b_out);

	return vec_pruning_line;
}



vector<double> _pruning_line2(bool to_rotate, double distance, int gap_to_fit, double x_kstop, double y_kstop){
	int TOF1_gap = 99;
	double d_angle = 999.99;
	double ddx = 999.99;
	double ddy = 999.99;

	double a_out = 999.99;
	double b_out = 999.99;

	vector<double> V_temp;		V_temp.clear(); 
	vector<double> vec_out; 	vec_out.clear();

	//if(to_rotate) TOF1_gap = TOF1_rotated[gap_to_fit-1]+1;
	//else TOF1_gap = gap_to_fit;
	TOF1_gap = gap_to_fit;
    
    if(TOF1_gap==1)  d_angle =  60. - 90.;
    if(TOF1_gap==2)  d_angle =  30. - 90.;
    if(TOF1_gap==3)  d_angle =   0. - 90.;
    if(TOF1_gap==4)  d_angle = 330. - 90.;
    if(TOF1_gap==5)  d_angle = 300. - 90.;
    if(TOF1_gap==6)  d_angle = 270. - 90.;
    if(TOF1_gap==7)  d_angle = 240. - 90.;
    if(TOF1_gap==8)  d_angle = 210. - 90.;
    if(TOF1_gap==9)  d_angle = 180. - 90.;
    if(TOF1_gap==10) d_angle = 150. - 90.;
    if(TOF1_gap==11) d_angle = 120. - 90.;
    if(TOF1_gap==12) d_angle =  90. - 90.;


    ddx = distance*sin(d_angle*PI/180.);
    ddy = distance*cos(d_angle*PI/180.);


    V_temp = _pruning_line(x_kstop, y_kstop, gap_to_fit);

    a_out = ((V_temp[1]-ddy)-(V_temp[3]-ddy))/((V_temp[0]+ddx)-(V_temp[2]+ddx));
    b_out = (V_temp[1]-ddy) - a_out*(V_temp[0]+ddx);


    vec_out.push_back(V_temp[0]+ddx);
    vec_out.push_back(V_temp[1]-ddy);
    vec_out.push_back(V_temp[2]+ddx);
    vec_out.push_back(V_temp[3]-ddy);
    vec_out.push_back(a_out);
    vec_out.push_back(b_out);


    return vec_out;
}

vector<double> V2_line;

//void _Linear_pruning(bool to_rotate, double x_kstop, double y_kstop, int gap_to_fit,
//	                 double distance, vector<int> V_lepton){
vector<int> _Linear_pruning(bool to_rotate, double x_kstop, double y_kstop, int gap_to_fit,
	                 double distance, vector<int> V_lepton){

	//Determine Pruning Line using _pruning_line function
	vector<double> V_pruning_line;  V_pruning_line.clear();
	V_pruning_line = _pruning_line(x_kstop, y_kstop, gap_to_fit);



	//Pruning the data
	int TOF1_gap = 99;
	vector<int> V_temp;		V_temp.clear();
	vector<int> V_TBR;		V_TBR.clear();
	vector<int> V_pruned;	V_pruned.clear();
	
/*	if(to_rotate){
		for(unsigned int i=0; i<V_lepton.size(); i++){
			V_temp.push_back(TARGET_Rotated_index[V_lepton[i]]);
		}
		TOF1_gap  = TOF1_rotated[gap_to_fit-1]+1;
	}
*/	//else{
		V_temp = V_lepton;
		TOF1_gap = gap_to_fit;
	//}


	if(TOF1_gap==1 || TOF1_gap==2 || TOF1_gap==10 || TOF1_gap==11 || TOF1_gap==12){
    	for(unsigned int i=0; i<V_temp.size(); i++){
        	if(Yloc[V_temp[i]] < Xloc[V_temp[i]]*V_pruning_line[4]+V_pruning_line[5]){
            	V_TBR.push_back(V_temp[i]);
          	}
        } 
	}
	else if(TOF1_gap==4 || TOF1_gap==5 || TOF1_gap==6|| TOF1_gap==7 || TOF1_gap==8){     
    	for(unsigned int i=0; i<V_temp.size(); i++){
    		if(Yloc[V_temp[i]] > Xloc[V_temp[i]]*V_pruning_line[4]+V_pruning_line[5]){
    			V_TBR.push_back(V_temp[i]);
    		}
    	}
    }
    else if(TOF1_gap==3){
    	for(unsigned int i=0; i<V_temp.size(); i++){
    		if(Xloc[V_temp[i]]<x_kstop){
    			V_TBR.push_back(V_temp[i]);
    		}
    	}
    }
    else if(TOF1_gap==9){
    	for(unsigned int i=0; i<V_temp.size(); i++){
    		if(Xloc[V_temp[i]]>x_kstop){
    			V_TBR.push_back(V_temp[i]);
    		}
    	}
    }

    vector<double> V_x_temp;	V_x_temp.clear();
    vector<double> V_y_temp;	V_y_temp.clear();


    for(unsigned int i=0; i<V_TBR.size(); i++){
      V_temp.erase(find(V_temp.begin(), V_temp.end(), V_TBR[i]));
    }

    for(unsigned j=0; j<V_temp.size(); j++){
    	V_x_temp.push_back(Xloc[V_temp[j]]); 
    	V_y_temp.push_back(Yloc[V_temp[j]]); 
    }

  	// TGraph *gr_temp;
  	// gr_temp = new TGraph(V_temp.size(), &V_x_temp[0], &V_y_temp[0]);
  	// gr_temp->SetMarkerStyle(21);
  	// gr_temp->SetMarkerSize(0.8);
  	// gr_temp->SetMarkerColor(3);
  	// gr_temp->Draw("sameP");


/////////////////////////////////////////////////////////////

    /// REPRUNING IF TOO FEW LEPTONS KEPT 
  	
  	//vector<double> V2_line;		V2_line.clear();
  	V2_line.clear();

  	if(V_temp.size() <= lepton_cluster_size){
  		
  		V_temp.clear();
  		V_TBR.clear();

  		V2_line = _pruning_line2(to_rotate, distance, gap_to_fit, x_kstop, y_kstop);

  // 		if(to_rotate){
		// 	for(unsigned int i=0; i<V_lepton.size(); i++){
		// 		V_temp.push_back(TARGET_Rotated_index[V_lepton[i]]);
		// 	}
		// 	TOF1_gap  = TOF1_rotated[gap_to_fit-1]+1;
		// }
		//else{
			V_temp = V_lepton;
			TOF1_gap = gap_to_fit;
		//}

		if(TOF1_gap==1 || TOF1_gap==2 || TOF1_gap==10 || TOF1_gap==11 || TOF1_gap==12){
    		for(unsigned int i=0; i<V_temp.size(); i++){
        		if(Yloc[V_temp[i]] < Xloc[V_temp[i]]*V2_line[4]+V2_line[5]){
            		V_TBR.push_back(V_temp[i]);
          		}
        	} 
		}
		else if(TOF1_gap==4 || TOF1_gap==5 || TOF1_gap==6|| TOF1_gap==7 || TOF1_gap==8){     
    		for(unsigned int i=0; i<V_temp.size(); i++){
    			if(Yloc[V_temp[i]] > Xloc[V_temp[i]]*V2_line[4]+V2_line[5]){
    				V_TBR.push_back(V_temp[i]);
    			}
    		}
    	}
    	else if(TOF1_gap==3){
    		for(unsigned int i=0; i<V_temp.size(); i++){
    			if(Xloc[V_temp[i]]<V2_line[0]){
    				V_TBR.push_back(V_temp[i]);
    			}
    		}
    	}
    	else if(TOF1_gap==9){
    		for(unsigned int i=0; i<V_temp.size(); i++){
    			if(Xloc[V_temp[i]]>V2_line[0]){
    				V_TBR.push_back(V_temp[i]);
    			}
    		}
    	}

    	for(unsigned int i=0; i<V_TBR.size(); i++){
      		V_temp.erase(find(V_temp.begin(), V_temp.end(), V_TBR[i]));
    	}

    
    	V_x_temp.clear();
    	V_y_temp.clear();
    	for(unsigned j=0; j<V_temp.size(); j++){
    		V_x_temp.push_back(Xloc[V_temp[j]]); 
    		V_y_temp.push_back(Yloc[V_temp[j]]); 
    	}
    }

//////////////////////////////////////////////////////////////////
    vector<double> V_x_pruned;	V_x_pruned.clear();
    vector<double> V_y_pruned;	V_y_pruned.clear();

	if(V_temp.size() <= lepton_cluster_size){
		V_pruned = V_lepton;
	}    
	else{
		V_pruned = V_temp;
	}


	for(unsigned int i=0; i<V_pruned.size(); i++){
		V_x_pruned.push_back(Xloc[V_pruned[i]]);
		V_y_pruned.push_back(Yloc[V_pruned[i]]);
	}

   //  //DEBUG
   //  TGraph *gr_pruned;
  	// gr_pruned = new TGraph(V_pruned.size(), &V_x_pruned[0], &V_y_pruned[0]);
  	// gr_pruned->SetMarkerStyle(21);
  	// gr_pruned->SetMarkerSize(0.8);
  	// gr_pruned->SetMarkerColor(3);
  	// gr_pruned->Draw("sameP");


	vector<int> vec_pruned_out;	vec_pruned_out.clear();

	// if(to_rotate){
	// 	for(unsigned int i=0; i<V_pruned.size(); i++){
	// 		vec_pruned_out.push_back(TARGET_Rotated_index_inverse[V_pruned[i]]);
	// 	}
	// }
	//else vec_pruned_out = V_pruned;
	vec_pruned_out = V_pruned;

	return vec_pruned_out;


}
