

void _Drawing_pruning(int Run_Number, int ievt,
					            vector<int> vec_kaons, vector<double> vec_kstop_0,
					            vector<int> vec_triangle_leptons, vector<int> vec_line_leptons,
					            TGraph *gr_kaons_final, TGraph *gr_kaons_unfitted, TLine *line_fit_kaons,
					            TGraph *gr_leptons_final, TLine *line_fit_leptons){
	

	TCanvas *c_pruning2;
  	c_pruning2 = new TCanvas("Pruning","PRUNING NEW2",50,500,1050,700);
  	c_pruning2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
  	c_pruning2->Divide(3,2);

	float R_TARGET = 29.0; 
  	TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0); 

    TH2F *h_Pruning2[2];   
    char Title_Pruning2[2][100];   
    char Name_Pruning2[2][100];

  	sprintf(Title_Pruning2[0],"Linear Pruning  |  Run %d, Evt %d", Run_Number, ievt); 
  	sprintf(Name_Pruning2[0],"Linear Pruning"); 

  	sprintf(Title_Pruning2[1],"Triangular Pruning  |  Run %d, Evt %d", Run_Number, ievt); 
  	sprintf(Name_Pruning2[1],"Triangular Pruning"); 

  	h_Pruning2[0] = new TH2F(Name_Pruning2[0],Title_Pruning2[0],100,-50,50,100,-50,50);
  	h_Pruning2[1] = new TH2F(Name_Pruning2[1],Title_Pruning2[1],100,-50,50,100,-50,50);



  	/////// KAONS
  	vector<double> vec_x_kaons;		vec_x_kaons.clear();
  	vector<double> vec_y_kaons;		vec_y_kaons.clear();
  	vector<double> vec_ex_kaons;	vec_ex_kaons.clear();
  	vector<double> vec_ey_kaons;	vec_ey_kaons.clear();

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


  	line_fit_kaons->SetLineWidth(2);
  	line_fit_kaons->SetLineColor(4);
  	
  	line_fit_leptons->SetLineWidth(2);
  	line_fit_leptons->SetLineColor(2);

  	//******* Kstop_0
  	vector<double> vec_x_kstop_0;		vec_x_kstop_0.clear();
  	vector<double> vec_y_kstop_0;		vec_y_kstop_0.clear();

  	vec_x_kstop_0.push_back(vec_kstop_0[0]);
  	vec_y_kstop_0.push_back(vec_kstop_0[1]);
  
  	TGraph *gr_Kstop_0;
  	gr_Kstop_0 = new TGraph(vec_x_kstop_0.size(),&vec_x_kstop_0[0],&vec_y_kstop_0[0]);
  	gr_Kstop_0->SetMarkerStyle(34);
  	gr_Kstop_0->SetMarkerColor(1);
  	gr_Kstop_0->SetMarkerSize(1.3);


  	/////// PRUNED LEPTONS
  	//****** Triangular pruning
  	vector<double> vec_x_triangle_leptons;		vec_x_triangle_leptons.clear();
  	vector<double> vec_y_triangle_leptons;		vec_y_triangle_leptons.clear();
  	vector<double> vec_ex_triangle_leptons;		vec_ex_triangle_leptons.clear();
  	vector<double> vec_ey_triangle_leptons;		vec_ey_triangle_leptons.clear();

	vector<double> vec_x_triangle_leptons_pr; 	vec_x_triangle_leptons_pr.clear();
	vector<double> vec_y_triangle_leptons_pr;	vec_y_triangle_leptons_pr.clear();	

	// if(to_rotate){
 //  		for(unsigned int i=0; i<vec_triangle_leptons.size(); i++){
 // 			vec_x_triangle_leptons.push_back(Xloc[vec_triangle_leptons[i]]);
 // 			vec_y_triangle_leptons.push_back(Yloc[vec_triangle_leptons[i]]);
 // 			vec_ex_triangle_leptons.push_back(TARGET_Errors_X);
 // 			vec_ey_triangle_leptons.push_back(TARGET_Errors_Y);
 // 			vec_x_triangle_leptons_pr.push_back(Xloc[TARGET_Rotated_index[vec_triangle_leptons[i]]]);
 // 			vec_y_triangle_leptons_pr.push_back(Yloc[TARGET_Rotated_index[vec_triangle_leptons[i]]]);			
 //  		}
 //  	}
  	//else{
  		for(unsigned int i=0; i<vec_triangle_leptons.size(); i++){
 			vec_x_triangle_leptons.push_back(Xloc[vec_triangle_leptons[i]]);
 			vec_y_triangle_leptons.push_back(Yloc[vec_triangle_leptons[i]]);
 			vec_ex_triangle_leptons.push_back(TARGET_Errors_X);
 			vec_ey_triangle_leptons.push_back(TARGET_Errors_Y);
 			vec_x_triangle_leptons_pr.push_back(Xloc[vec_triangle_leptons[i]]);
 			vec_y_triangle_leptons_pr.push_back(Yloc[vec_triangle_leptons[i]]);			
  		}
  	//}

  	TGraph *gr_triangle_leptons;
  	gr_triangle_leptons = new TGraphErrors(vec_triangle_leptons.size(),&vec_x_triangle_leptons[0],&vec_y_triangle_leptons[0],
  										   &vec_ex_triangle_leptons[0],&vec_ey_triangle_leptons[0]);
  	
  	gr_triangle_leptons->SetMarkerColor(2);
  	gr_triangle_leptons->SetMarkerStyle(21);
  	gr_triangle_leptons->SetMarkerSize(0.8);

  	TGraph *gr_triangle_leptons_pr;
  	gr_triangle_leptons_pr = new TGraph(vec_triangle_leptons.size(),&vec_x_triangle_leptons_pr[0],&vec_y_triangle_leptons_pr[0]);
  	
  	gr_triangle_leptons_pr->SetMarkerColor(3);
  	gr_triangle_leptons_pr->SetMarkerStyle(21);
  	gr_triangle_leptons_pr->SetMarkerSize(0.8);

  	//Pruning Area
  	TGraph *gr_pruning_area;
  	gr_pruning_area = new TGraph(V_X_pruning_area.size(), &V_X_pruning_area[0], &V_Y_pruning_area[0]);
  	gr_pruning_area->SetLineWidth(2);
  	gr_pruning_area->SetLineColor(3);


  	// ****** Linear pruning
  	vector<double> vec_x_line_leptons;		vec_x_line_leptons.clear();
  	vector<double> vec_y_line_leptons;		vec_y_line_leptons.clear();
  	vector<double> vec_ex_line_leptons;		vec_ex_line_leptons.clear();
  	vector<double> vec_ey_line_leptons;		vec_ey_line_leptons.clear();

	vector<double> vec_x_line_leptons_pr; 	vec_x_line_leptons_pr.clear();
	vector<double> vec_y_line_leptons_pr;	vec_y_line_leptons_pr.clear();	

	// if(to_rotate){
 //  		for(unsigned int i=0; i<vec_line_leptons.size(); i++){
 // 			vec_x_line_leptons.push_back(Xloc[vec_line_leptons[i]]);
 // 			vec_y_line_leptons.push_back(Yloc[vec_line_leptons[i]]);
 // 			vec_ex_line_leptons.push_back(TARGET_Errors_X);
 // 			vec_ey_line_leptons.push_back(TARGET_Errors_Y);			
 // 			vec_x_line_leptons_pr.push_back(Xloc[TARGET_Rotated_index[vec_line_leptons[i]]]);
 // 			vec_y_line_leptons_pr.push_back(Yloc[TARGET_Rotated_index[vec_line_leptons[i]]]);			
 //  		}
 //  	}
  	//else{
  		for(unsigned int i=0; i<vec_line_leptons.size(); i++){
 			vec_x_line_leptons.push_back(Xloc[vec_line_leptons[i]]);
 			vec_y_line_leptons.push_back(Yloc[vec_line_leptons[i]]);
 			vec_ex_line_leptons.push_back(TARGET_Errors_X);
 			vec_ey_line_leptons.push_back(TARGET_Errors_Y);			
 			vec_x_line_leptons_pr.push_back(Xloc[vec_line_leptons[i]]);
 			vec_y_line_leptons_pr.push_back(Yloc[vec_line_leptons[i]]);			
  		}
  	//}

  	TGraph *gr_line_leptons;
  	gr_line_leptons = new TGraphErrors(vec_line_leptons.size(),&vec_x_line_leptons[0],&vec_y_line_leptons[0],
  									   &vec_ex_line_leptons[0],&vec_ey_line_leptons[0]);
  	
  	gr_line_leptons->SetMarkerColor(2);
  	gr_line_leptons->SetMarkerStyle(21);
  	gr_line_leptons->SetMarkerSize(0.8);

  	TGraph *gr_line_leptons_pr;
  	gr_line_leptons_pr = new TGraph(vec_line_leptons.size(),&vec_x_line_leptons_pr[0],&vec_y_line_leptons_pr[0]);
  	
  	gr_line_leptons_pr->SetMarkerColor(3);
  	gr_line_leptons_pr->SetMarkerStyle(21);
  	gr_line_leptons_pr->SetMarkerSize(0.8);

  
  	gr_leptons_final->SetMarkerColor(2);
  	gr_leptons_final->SetMarkerStyle(21);
  	gr_leptons_final->SetMarkerSize(0.8);



  	// Pruning Line
  	TLine *pruning_line;
  	pruning_line = new TLine(vec_pruning_line_coordinates[0], vec_pruning_line_coordinates[1], 
                                 vec_pruning_line_coordinates[2], vec_pruning_line_coordinates[3]);
  	pruning_line->SetLineWidth(2);
  	pruning_line->SetLineColor(3);

  	// Second Pruning Line if needed
  	TLine *pruning_line2;
  	if(V2_line.size()>0){
  		pruning_line2 = new TLine(V2_line[0], V2_line[1], V2_line[2], V2_line[3]);
  		pruning_line2->SetLineWidth(2);
  		pruning_line2->SetLineColor(6);
  	}


  	/////// DRAWING 
  	// PANEL 1	
  	c_pruning2->cd(1);
  	h_Pruning2[0]->Draw();
  	ell_Target->Draw("same");
  	gr_kaons_final->Draw("sameP");
  	line_fit_kaons->Draw("same");
  	gr_leptons_final->Draw("sameP");
  	line_fit_leptons->Draw("same");
  	gr_Kstop_0->Draw("sameP");

  	// PANEL 2
  	c_pruning2->cd(2);
  	h_Pruning2[0]->Draw();
  	ell_Target->Draw("same");
  	gr_kaons_final->Draw("sameP");
  	line_fit_kaons->Draw("same");
  	gr_leptons_final->Draw("sameP");
  	line_fit_leptons->Draw("same");
  	gr_line_leptons_pr->Draw("sameP");
  	pruning_line->Draw("same");
  	if(V2_line.size()>0) pruning_line2->Draw("same");
  	gr_Kstop_0->Draw("sameP");

  	// PANEL 3
  	c_pruning2->cd(3);
  	h_Pruning2[0]->Draw();
  	ell_Target->Draw("same");
  	//gr_kaons->Draw("sameP");
  	gr_kaons_unfitted->Draw("sameP");
  	gr_line_leptons->Draw("sameP");

  	// PANEL 4
  	c_pruning2->cd(4);
  	h_Pruning2[0]->Draw();
  	ell_Target->Draw("same");
  	gr_kaons_final->Draw("sameP");
  	line_fit_kaons->Draw("same");
  	gr_leptons_final->Draw("sameP");
  	line_fit_leptons->Draw("same");
  	gr_Kstop_0->Draw("sameP");

  	// PANEL 5
  	c_pruning2->cd(5);
  	h_Pruning2[1]->Draw();
  	ell_Target->Draw("same");
  	gr_kaons_final->Draw("sameP");
  	line_fit_kaons->Draw("same");
  	gr_leptons_final->Draw("sameP");
  	line_fit_leptons->Draw("same");
  	gr_triangle_leptons_pr->Draw("sameP");
  	gr_pruning_area->Draw("same");
  	gr_Kstop_0->Draw("sameP");

  	// PANEL 6
  	c_pruning2->cd(6);
 	h_Pruning2[1]->Draw();
  	ell_Target->Draw("same");
  	//gr_kaons->Draw("sameP");
  	gr_kaons_unfitted->Draw("sameP");
  	gr_triangle_leptons->Draw("sameP");

}

