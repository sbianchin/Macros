#include "CommonParameters.h"

vector<int> _remove_singles(vector<int> v_input, bool is_kaon[256], bool is_kaon2[256], int kaon_substitute){

	bool has_neighbours = false;
	vector<int> vec_TBR;	vec_TBR.clear();

	for(unsigned int i=0; i<v_input.size(); i++){
		has_neighbours = false;
		for(int j=0; j<8; j++){
			if(is_kaon[TARGET_neighbours[v_input[i]][j]] || is_kaon2[TARGET_neighbours[v_input[i]][j]]){
				has_neighbours = true;
			}
		}

		if(!has_neighbours && v_input[i]!=kaon_substitute) vec_TBR.push_back(v_input[i]);
	}

	// Removal of singlets
    for(unsigned int i=0; i<vec_TBR.size(); i++){
    	v_input.erase(find(v_input.begin(), v_input.end(), vec_TBR[i]));
    }

	return v_input;
}

vector<double> _get_kaon_centroid(vector<int> V_input){

	double X_weight = 0.;
	double Y_weight = 0.;
	double total_energy = 0.;

	vector<double> K_centroid;	K_centroid.clear();



	if(V_input.size()>0){
		for(unsigned int i=0; i<V_input.size(); i++){
			X_weight += ADC_Low_TARGET[V_input[i]]*Xloc[V_input[i]];
			Y_weight += ADC_Low_TARGET[V_input[i]]*Yloc[V_input[i]];
        	total_energy += ADC_Low_TARGET[V_input[i]];
		}	
	}

    if(V_input.size()>0 && total_energy!=0){
		K_centroid.push_back(X_weight/total_energy);
		K_centroid.push_back(Y_weight/total_energy);		
    }
    else{
		K_centroid.push_back(0.);
		K_centroid.push_back(0.);		
    }

    return K_centroid;
}

vector<double> _get_kstop_0(vector<int> V_lepton, vector<int> V_kaon, vector<double> V_centroid,
				  double a1, double b1, double a2, double b2){

    vector<double> V_kstop_0;		V_kstop_0.clear();
	vector<double> V_intersect_pr;	V_intersect_pr.clear();

	double distance_to_kaon;
	double distance_min = 999.;
	int K_closest = 999;

	int N_kaon_adc_LG_max = -999;
	int N_kaon_LG_max = -1; 

	if(V_kaon.size()>0){
		if(V_lepton.size() == 0){
			V_kstop_0.push_back(V_centroid[0]);
			V_kstop_0.push_back(V_centroid[1]);		
		}
		else if(V_kaon.size()<=kaon_cluster_size){
			for(unsigned int i=0; i<V_kaon.size(); i++){
				if(ADC_Low_TARGET[V_kaon[i]]+TARGET_ADC_Thr_LG_Offset>N_kaon_adc_LG_max){
					N_kaon_adc_LG_max = ADC_Low_TARGET[V_kaon[i]]+TARGET_ADC_Thr_LG_Offset;
           			N_kaon_LG_max = V_kaon[i];
				}
			}
			V_kstop_0.push_back(Xloc[N_kaon_LG_max]);
			V_kstop_0.push_back(Yloc[N_kaon_LG_max]);
		}
		else{

			V_intersect_pr = _2lines_intersect(a1,b1,a2,b2);         		

      		for(unsigned int j=0; j<V_kaon.size(); j++){
       			distance_to_kaon = sqrt(pow((Xloc[V_kaon[j]]-V_intersect_pr[0]),2) + pow((Yloc[V_kaon[j]]-V_intersect_pr[1]),2));
       			if(distance_to_kaon < distance_min){
           			distance_min = distance_to_kaon;
           			K_closest = V_kaon[j];
       			}
       		}
      		if(distance_min < DISTANCE_MAX_TO_KAON_PRUNED){
       			V_kstop_0.push_back(V_intersect_pr[0]);
       			V_kstop_0.push_back(V_intersect_pr[1]);
       		}
       		else{
       			V_kstop_0.push_back(V_centroid[0]);
       			V_kstop_0.push_back(V_centroid[1]);
       		}        		
		}
	}
	else{
  		V_kstop_0.push_back(999.);
    	V_kstop_0.push_back(999.);
	}

	return V_kstop_0;
}

vector<int> _select_leptons_to_fit(bool to_rotate, int gap_to_fit, vector<int> vec_leptons, 
							vector<double> vec_line_parameters, int EB_add_weight){

	vector<int> vec_temp;	vec_temp.clear();
	vector<int> vec_TBR;	vec_TBR.clear();
	vector<int> vec_EB_pr;	vec_EB_pr.clear();
	int TOF1_gap = -1;


	vector<int> vec_out; 	vec_out.clear();

	if(to_rotate){
		TOF1_gap = TOF1_rotated[gap_to_fit-1]+1;
		for(unsigned int i=0; i<vec_leptons.size(); i++){
			vec_temp.push_back(TARGET_Rotated_index[vec_leptons[i]]);
		}
	}
	else{
		TOF1_gap = gap_to_fit;
		vec_temp = vec_leptons;
	}


	/// Remove leptons to far from the fit linevec_temp
	if(vec_line_parameters.size()>0){
 		for(unsigned int i=0; i<vec_temp.size(); i++){

 			if(distance_to_line(Xloc[vec_temp[i]],Yloc[vec_temp[i]],vec_line_parameters[0],vec_line_parameters[1]) >= max_dist){
 				vec_TBR.push_back(vec_temp[i]);
			}
 		}

 		for(unsigned int i=0; i<vec_TBR.size(); i++){
      		vec_temp.erase(find(vec_temp.begin(), vec_temp.end(), vec_TBR[i]));
      	}
    }


    /// Add weight on Edge Bars
    for(unsigned int i=0; i<vec_temp.size(); i++){
    	if(IsIn(vec_temp[i], channel[TOF1_gap-1][0], channel[TOF1_gap-1][1], channel[TOF1_gap-1][2],
    		channel[TOF1_gap-1][3], channel[TOF1_gap-1][4], channel[TOF1_gap-1][5], 
    		channel[TOF1_gap-1][6], channel[TOF1_gap-1][7], channel[TOF1_gap-1][8], 
    		channel[TOF1_gap-1][9], channel[TOF1_gap-1][10], channel[TOF1_gap-1][11])){

    		vec_EB_pr.push_back(vec_temp[i]);
    	}
	}

	for(unsigned int i=0; i<vec_EB_pr.size(); i++){
		for(int j=0; j<EB_add_weight; j++){
			vec_temp.push_back(vec_EB_pr[i]);
		}
	}

	return vec_temp;
}

vector<int> _rotate_vector(vector<int> vec_input){

	vector<int> vec_output;		vec_output.clear();
	vector<int>	vec_temp;		vec_temp.clear();

	for(unsigned int i=0; i<vec_input.size(); i++){
		vec_temp.push_back(TARGET_Rotated_index[vec_input[i]]);
	}

	vec_input.clear();
	vec_input = vec_temp;

	return vec_input;
}

int _which_TOF1_section(int TOF1_gap, vector<double> vec_intersect){

	double d0, d1, d2, d3, d4;
	double dist_min = 999.99;
	double dist = 999.99;
	int closest_section = -1;

	for(int i=0; i<5; i++){
		dist = distance(vec_intersect[0], vec_intersect[1], TOF1_Xloc[TOF1_gap-1][i], TOF1_Yloc[TOF1_gap-1][i]);

		if(dist < dist_min){
			dist_min = dist;
			closest_section = i;
		}
	}

	return closest_section;
}

vector<double> _restore_fit_line(double a, double b){

	vector<double> vec_out;
	vec_out.clear();

	vec_out.push_back(50.*a + b);
	vec_out.push_back(-50.);
	vec_out.push_back(-50*a + b);
	vec_out.push_back(50.);

	return vec_out;

}

vector<double> _get_kstop_final(bool to_restore,
	vector<int> vec_leptons, vector<int> vec_kaons, vector<double> vec_kaon_centroid,
					  double a_kaons, double b_kaons, double a_fit3, double b_fit3,
					  bool has_kaon_sub, int kaon_substitute){

	vector<double> vec_kstop;		vec_kstop.clear();
	vector<double> vec_intersect;	vec_intersect.clear();

	double A, B;

	int kaon_LG_max = -1;
	int kaon_ADC_LG_max = -1;

	if(to_restore){
		A = -1/a_fit3;
		B = b_fit3/a_fit3;
	}


	if(vec_leptons.size()==0){
		vec_intersect.push_back(999.99);
		vec_intersect.push_back(999.99);

		if(vec_kaons.size()==0){
			vec_kstop.push_back(999.99);
			vec_kstop.push_back(999.99);
		}
		else if(vec_kaons.size() > kaon_cluster_size){
        	vec_kstop.push_back(vec_kaon_centroid[0]);
        	vec_kstop.push_back(vec_kaon_centroid[1]);
		}
		else{
			if(has_kaon_sub){
				kaon_LG_max = kaon_substitute;
			}
			else{
        		for(unsigned int i=0; i<vec_kaons.size(); i++){
          			if(ADC_Low_TARGET[vec_kaons[i]]+TARGET_ADC_Thr_LG_Offset>kaon_ADC_LG_max){
            			kaon_ADC_LG_max = ADC_Low_TARGET[vec_kaons[i]]+TARGET_ADC_Thr_LG_Offset;
            			kaon_LG_max = vec_kaons[i];
          			}
   				}
   			}

   			vec_kstop.push_back(Xloc[kaon_LG_max]);
   			vec_kstop.push_back(Yloc[kaon_LG_max]);
   		}
   	}
   	else{
   		if(vec_kaons.size()==0){
			vec_kstop.push_back(999.99);
			vec_kstop.push_back(999.99);
		}
		else if(vec_kaons.size() > kaon_cluster_size){
			if(to_restore) vec_intersect = _2lines_intersect(a_kaons, b_kaons, A, B);
			else vec_intersect = _2lines_intersect(a_kaons, b_kaons, a_fit3, b_fit3);

			vec_kstop = vec_intersect;
		}
		else{
			if(has_kaon_sub){
				kaon_LG_max = kaon_substitute;
			}
			else{
				for(unsigned int i=0; i<vec_kaons.size(); i++){
          			if(ADC_Low_TARGET[vec_kaons[i]]+TARGET_ADC_Thr_LG_Offset > kaon_ADC_LG_max){
            			kaon_ADC_LG_max = ADC_Low_TARGET[vec_kaons[i]]+TARGET_ADC_Thr_LG_Offset;
            			kaon_LG_max = vec_kaons[i];
          			}
				}
			}
   			
   			vec_kstop.push_back(Xloc[kaon_LG_max]);
   			vec_kstop.push_back(Yloc[kaon_LG_max]);
	
   		}
   	}


/*		if(vec_kaons.size() <= kaon_cluster_size){
			vec_kstop.push_back(vec_kaon_centroid[0]);
			vec_kstop.push_back(vec_kaon_centroid[1]);
		}
		else{
			if(to_restore) vec_intersect = _2lines_intersect(a_kaons, b_kaons, A, B);
			else vec_intersect = _2lines_intersect(a_kaons, b_kaons, a_fit3, b_fit3);

			vec_kstop = vec_intersect;
		}
	}
	else{  // <----- ONLY TEMPORARELY
		vec_kstop.push_back(999.99);
		vec_kstop.push_back(999.99);
	}
*/	

	return vec_kstop;
}

vector<int> _remove_doubles(vector<int> vec_in){
	vector<int> vec_out;	vec_out.clear();

	vec_out = vec_in;

   sort(vec_out.begin(), vec_out.end());
   vec_out.erase(unique(vec_out.begin(), vec_out.end()), vec_out.end());

   return vec_out; 
}