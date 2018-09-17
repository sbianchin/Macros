void _ED_printouts(char finput[200], bool NoFile, char walk_file[100],
	int nentries, char Version[],
	bool Good_Event, bool Switch, int enable_cout, int ievt,
	int kaon_TDC_min, int kaon_TDC_max, int TDC_min_Kstop, int TDC_max_Kstop,
	int TDC_average_NEW, char TDC_LE_TARGET_sel_string[256][16][20] , int gap_to_fit,
	double angle_final, double delta_phi_deg, double ChiS_final, int ndf_final,
	vector<double> vec_X_int_SFT, vector<double> vec_Y_int_SFT, 
	vector<double> vec_X_int_TOF1, vector<double> vec_Y_int_TOF1,
	vector<int> vec_TOF1_Gap, double Length_in_target){


  	int TDC_print[256];
  	for(int i=0; i<256; i++) TDC_print[i] = -1;


  	//int n_bar_selected = 999.;
  	//n_bar_selected = vec_TARGET_bar_selected.size();
  	//char TDC_LE_TARGET_sel_string[n_bar_selected][16][20];

    cout << "   " << endl;
    cout << "   " << endl;
    cout << "************************************************************************************************************" << endl;

    cout << "File opened:  " << finput << endl;
    if(!NoFile) cout << "Time Walk Parameter File:  " << walk_file << endl;
    if(NoFile) cout << "Time Walk Parameter File:  " << "NOT FOUND !" << endl;
    cout << "SFT Mapping File:  " << par_finput << endl;
    cout << "MWPC Mapping File:  " << par_finput2 << endl;
    cout << "Total Number of Events:  " << nentries <<endl;
    cout << "Event_Display.C -- " << Version << endl;
    cout << "************************************************************************************************************" << endl;
    cout << "  " << endl;
    cout << "  " << endl;

    if(Good_Event) cout << "Event: "<< ievt << "   --  GOOD EVENT!" << endl;
    cout << endl;
    cout << endl;

    cout << " **** RAW ADCs" << endl;
    cout << " ///////////   TOF1   ///////////   //////////////////////   TOF2   //////////////////////" << endl;
    for(Int_t j_TOF=0; j_TOF<12; j_TOF++){
      printf(" %2d   UP: %4d  |  DOWN: %4d    ||    %2d   AO: %4d  |  BO: %4d  |  AI: %4d  |  BI:%4d\n",
      j_TOF+1, ADC_tof1U[j_TOF], ADC_tof1D[j_TOF], j_TOF+1, ADC_tof2AO[j_TOF], ADC_tof2BO[j_TOF], ADC_tof2AI[j_TOF], ADC_tof2BI[j_TOF]);
    }

    cout << "NEW" << endl;
    select_kaon(true);

    cout << "TIME CUTS" << endl;
    cout << "# Letpon + Kaon : " << kaon_TDC_min << " <= " << "TDC" << " <= " << kaon_TDC_max << endl;
    cout << "# Kaon : " << TDC_min_Kstop << " <= " << "TDC" << " <= " << TDC_max_Kstop << endl;
    cout << endl;
    cout << endl;

    /// Print Values
    if(Good_Event && enable_cout==1){
    	if(!NoFile) cout << "//////  TARGET  (WALK CORRECTION APPLIED !)  //////" << endl;
        if(NoFile) cout << "//////  TARGET  //////" << endl;
       	cout << "Event: " << ievt << "     ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << "     ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
       	cout << "Fiber  HG-Ped  LG-Ped    TDC[0]   TDC[1]   TDC[2]   TDC[3]   TDC[4]   TDC[5]   Thr(HG)   Thr(LG)   HG/LG" << endl;  // SEBASTIEN
       	for(Int_t jj=0; jj<256; jj++){
         	if(TDC_LE_TARGET[jj]>-1 || tdc_le_target[jj][1]>-1 || tdc_le_target[jj][2]>-1 || tdc_le_target[jj][3]>-1){
           		printf("%3d     %4d    %4d     %4d     %4d     %4d     %4d     %4d     %4d      %4d      %4d     %2.2f\n",
           		jj, ADC_High_TARGET[jj] + TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset,
           		TDC_LE_TARGET[jj], tdc_le_target[jj][1], tdc_le_target[jj][2], tdc_le_target[jj][3], tdc_le_target[jj][4], tdc_le_target[jj][5],
           		HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Ped_LG[jj])));
         	}
       	}
    }

    if(Good_Event && (enable_cout==0 || enable_cout==9)){
    	if(!NoFile) cout << "//////  TARGET  (WALK CORRECTION APPLIED !)  //////" << endl;
    	if(NoFile) cout << "//////  TARGET  //////" << endl;
    	cout << "Event: " << ievt << "     ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << "     ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    	cout << "Fiber  HG-Ped   LG-Ped   T[0]-Av   T[1]-Av   T[2]-Av   T[3]-Av   ADC(HG)   ADC(LG)      TDC    Thr(HG)   Thr(LG)   HG/LG" << endl;
       	for(Int_t jj=0; jj<256; jj++){
        	if(ADC_High_TARGET[jj]>0 || (ADC_Low_TARGET[jj]>0 && Switch==1)){
            	for(int k=0; k<6; k++){
               		if(tdc_le_target[jj][k] > -1){
                 		TDC_print[jj] = tdc_le_target[jj][k];
             		}
           		}
          		printf("%3d     %4d    %4d      %s       %s      %s      %s      %4d      %4d      %4d     %4d      %4d     %2.2f\n",
           		jj, ADC_High_TARGET[jj]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset,
           		TDC_LE_TARGET_corr[jj][0], TDC_LE_TARGET_corr[jj][1], TDC_LE_TARGET_corr[jj][2], TDC_LE_TARGET_corr[jj][3], ADC_High_TARGET[jj], ADC_Low_TARGET[jj], TDC_print[jj],
           		HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], 
              float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Ped_LG[jj])));
              //float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))
              //);
          	}
       	}
    }

    if(Good_Event && enable_cout==2){
    	if(!NoFile) cout << "//////  TARGET  (WALK CORRECTION APPLIED !)  //////" << endl;
       	if(NoFile) cout << "//////  TARGET  //////" << endl;
       	cout << "Event: " << ievt << "     ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << "     ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
       	cout << "Fiber  HG-Ped   LG-Ped   T[0]-Av   T[1]-Av   T[2]-Av   T[3]-Av   T[4]-Av   T[5]-Av   Thr(HG)   Thr(LG)   HG/LG" << endl;   // SEBASTIEN
       	for(Int_t jj=0; jj<256; jj++){
        	printf("%3d     %4d     %4d      %4d      %4d      %4d      %4d      %4d      %4d     %4d      %4d     %2.2f\n",
         	jj, ADC_High_TARGET[jj]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset,
         	TDC_LE_TARGET_corrected[jj][0], TDC_LE_TARGET_corrected[jj][1], TDC_LE_TARGET_corrected[jj][2], TDC_LE_TARGET_corrected[jj][3], TDC_LE_TARGET_corrected[jj][4], TDC_LE_TARGET_corrected[jj][5],
         	HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], 
          float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Ped_LG[jj]))
          //float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))
          );
    	}
    }

    if(Good_Event && (enable_cout==0 || enable_cout==1 || enable_cout==9)){
    	cout << " " << endl;
       	cout << " " << endl;
       	cout << "//////  SFT  //////" << endl;
       	cout << "ADC Threshold Offset (HG): " << SFT_ADC_Thr_HG_Offset << endl;
       	cout << "ADC Threshold Offset (LG): " << SFT_ADC_Thr_LG_Offset << endl << endl;
       	cout << "Channel   Layer     Fiber     HG-Thr     LG-Thr   TDC[0]    TDC[1]    TDC[2]    TDC[3]   (HG) Thr   (LG) Thr" << endl;
       	for(Int_t jj=0; jj<128; jj++){
         	if((ADC_High_SFT_corr[jj] != 0) && (has_TDC_SFT_hit[jj] > 0)){
           		printf("%3d       %3d      %3d-%1c      %4d       %4d      %4d     %4d      %4d      %4d       %4d       %4d\n",
           		jj, SFT_channel_to_layer[jj], SFT_channel_to_fiber[jj], SFT_channel_to_US_DS[jj], ADC_High_SFT[jj], ADC_Low_SFT[jj],
           		TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3], HG_SFT_ADC_Thr[jj], LG_SFT_ADC_Thr[jj]);
         	}
       	}
   	}

    if(Good_Event && enable_cout==2){
    	cout << " " << endl;
       	cout << " " << endl;
       	cout << "//////  SFT  //////" << endl << endl;
       	cout << "ADC Threshold Offset (HG): " << SFT_ADC_Thr_HG_Offset << endl;
       	cout << "ADC Threshold Offset (LG): " << SFT_ADC_Thr_LG_Offset << endl << endl;
       	cout << "Channel   Layer     Fiber     HG-Thr     LG-Thr   TDC[0]    TDC[1]    TDC[2]    TDC[3]    Thr(HG)    Thr(LG)" << endl;
       	for(Int_t jj=0; jj<128; jj++){
        	printf("%3d       %3d      %3d-%1c      %4d       %4d      %4d     %4d      %4d      %4d       %4d       %4d\n",
         	jj, SFT_channel_to_layer[jj], SFT_channel_to_fiber[jj], SFT_channel_to_US_DS[jj], ADC_High_SFT[jj], ADC_Low_SFT[jj],
         	TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3], HG_SFT_ADC_Thr[jj], LG_SFT_ADC_Thr[jj]);
       	}
  	}

   	cout << "" << endl;
    cout << "" << endl;


    if(selected_TOF2 > 6){ //LEFT
    	for(int q=0; q<56; q++){ // C2 Counters
        	if(Good_Event && ADC_C2X_L[q]>0){
           		cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << int(ADC_C2X_L[q]) << endl;
         	}
       	}
	    for(int q=0; q<16; q++){
        	if(Good_Event && ADC_C2Y_L[q]>0){
        		cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << int(ADC_C2Y_L[q]) << endl;
         	}
       	}

       	cout << " " << endl;

       	for(int q=0; q<64; q++){ // C3 Counters
         	if(Good_Event && ADC_C3X_L[q]>0){
           		cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << int(ADC_C3X_L[q]) << endl;
         	}
       	}
       	for(int q=0; q<16; q++){
         	if(Good_Event && ADC_C3Y_L[q]>0){
           		cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << int(ADC_C3Y_L[q]) << endl;
         	}
       	}

       	cout << " " << endl;

       	for(int q=0; q<72; q++){ // C4 Counters
         	if(Good_Event && ADC_C4X_L[q]>0){
           		cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << int(ADC_C4X_L[q]) << endl;
         	}
       	}
       	for(int q=0; q<16; q++){
         	if(Good_Event && ADC_C4Y_L[q]>0){
           		cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << int(ADC_C4Y_L[q]) << endl;
         	}
       	}
    }
    else{
    	for(int q=0; q<56; q++){ // C2 Counters
        	if(Good_Event && ADC_C2X_R[q]>0){
           		cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << int(ADC_C2X_R[q]) << endl;
         	}
       	}
       	for(int q=0; q<16; q++){
         	if(Good_Event && ADC_C2Y_R[q]>0){
           		cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << int(ADC_C2Y_R[q]) << endl;
         	}
       	}

       	cout << " " << endl;

       	for(int q=0; q<64; q++){ // C3 Counters
         	if(Good_Event && ADC_C3X_R[q]>0){
           		cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << int(ADC_C3X_R[q]) << endl;
         	}
       	}
       	for(int q=0; q<16; q++){
         	if(Good_Event && ADC_C3Y_R[q]>0){
           		cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << int(ADC_C3Y_R[q]) << endl;
         	}
       	}

       	cout << " " << endl;

       	for(int q=0; q<72; q++){ // C4 Counters
         	if(Good_Event && ADC_C4X_R[q]>0){
           		cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << int(ADC_C4X_R[q]) << endl;
         	}
       	}
       	for(int q=0; q<16; q++){
         	if(Good_Event && ADC_C4Y_R[q]>0){
           		cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << int(ADC_C4Y_R[q]) << endl;
         	}
       	}
  	}

    cout << endl;
    if(selected_TOF2 > 6) cout << " LEFT" << endl;
    else cout << " RIGHT" << endl;

    cout << " C2X clusters = " << C2X_clusters << endl;
    cout << " C2Y clusters = " << C2Y_clusters << endl;
    cout << " C3X clusters = " << C3X_clusters << endl;
    cout << " C3Y clusters = " << C3Y_clusters << endl;
    cout << " C4X clusters = " << C4X_clusters << endl;
    cout << " C4Y clusters = " << C4Y_clusters << endl;

    cout << "   " << endl;
    cout << "   " << endl;

    cout << "//////  C_k / C_pi  //////" << endl;
    cout << "Counter       TDC[0]           TDC[1]           TDC[2]           TDC[3]           TDC[4]           TDC[5]" << endl;
    for(int ic = 0; ic < 14; ic++){
    	printf("%3d        %4d / %4d      %4d / %4d      %4d / %4d      %4d / %4d      %4d / %4d      %4d / %4d\n",ic,
       	tdc_ck[ic][0], tdc_cpi[ic][0], tdc_ck[ic][1], tdc_cpi[ic][1], tdc_ck[ic][2], tdc_cpi[ic][2],
       	tdc_ck[ic][3], tdc_cpi[ic][3], tdc_ck[ic][4], tdc_cpi[ic][4], tdc_ck[ic][5], tdc_cpi[ic][5]);
    }

    cerr << endl;
    cerr << "N(Ck) = " << vec_Ck.size() << "   " << "TDC Ck Avg. = " << TDC_ck_avg << "   " << "Var(Ck) = " << TDC_ck_sigma << "   " << "Sigma(Ck) = " << TDC_ck_sigma/sqrt(vec_Ck.size()) << endl;
    cerr << "N(Cpi) =  " << vec_Cpi.size() << "   " << "TDC Cpi Avg. = " << TDC_cpi_avg << "   " << "Var(Cpi) = " << TDC_cpi_sigma << "   " << "Sigma(Cpi) = " << TDC_cpi_sigma/sqrt(vec_Cpi.size()) << endl;
    cerr << endl;
    cerr << "N'(Ck) = " << TDC_ck_counter << "   " << "TDC'(Ck) = " << TDC_ck_avg2 << endl;
    cerr << "N'(Cpi) = " << TDC_cpi_counter << "   " << "TDC'(Cpi) = " << TDC_cpi_avg2 << endl;
    cerr << endl;
    cerr << "TDC Mean Ck - TDC Average = " << TDC_ck_avg2 - TDC_average_NEW << endl;
    cerr << endl;

    cout << endl;
    cout << "////// TDC TARGET - TDC Ck (corrected)  //////" << endl;
    cout << "Fiber  HG-Ped   LG-Ped     T[0]-Ck    T[1]-Ck    T[2]-Ck    T[3]-Ck    T[4]-Ck    T[5]-Ck" << endl;
    for(unsigned int jj=0; jj<vec_TARGET_bar_selected.size(); jj++){
        printf("%3d     %4d    %4d        %s      %s      %s      %s      %s      %s\n",
        vec_TARGET_bar_selected[jj], ADC_High_TARGET[vec_TARGET_bar_selected[jj]]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[vec_TARGET_bar_selected[jj]] + TARGET_ADC_Thr_LG_Offset,
        TDC_LE_TARGET_sel_string[jj][0], TDC_LE_TARGET_sel_string[jj][1], TDC_LE_TARGET_sel_string[jj][2], TDC_LE_TARGET_sel_string[jj][3], TDC_LE_TARGET_sel_string[jj][4], TDC_LE_TARGET_sel_string[jj][5]);
    }
    cout << endl;
    cout << endl;
    cout << "//////  TOF1 SCORING  //////" << endl;
    for (int k=0; k<12; k++) cout << "TOF1 score " << k+1 << ": " << vec_TOF1_Gap[k+1] << endl;

    cout << "" << endl;
    cout << "TOF1 SELECTED FOR FITTING:  TOF1 " << gap_to_fit << endl;

    cout << "" << endl;
    cout << "" << endl;
    for (Int_t i = 0; i<12; i++){
    	if(Good_Event){
        	if (ADC_TOF1[i] > 0 || ADC_TOF1[i+12] > 0){
           		if((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])){
            		printf("ADC TOF1 Up-%2d: %5d -- Down-%2d: %5d  |  TDC TOF1 Up-%2d: %5d -- Down-%2d: %5d\n", i+1, ADC_TOF1[i], i+1, ADC_TOF1[i+12], i+1, TDC_TOF1U[i], i+1, TDC_TOF1D[i]);
           		}
           		else {
             		printf("ADC TOF1 Up-%2d: %5d -- Down-%2d: %5d  |\n", i+1, ADC_TOF1[i], i+1, ADC_TOF1[i+12]);
           		}
        	}
        	else {
        		if((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {
            		printf("                                         |  TDC TOF1 Up-%2d: %5d -- Down-%2d: %5d\n", i+1, TDC_TOF1U[i], i+1, TDC_TOF1D[i]);
           		}
        	}
     	}
    }

   	cout << "" << endl;

   	for (Int_t i = 0; i<12; i++) {
      	if (ADC_TOF2AO[i] > 0 || ADC_TOF2AI[i] > 0){
        	if ((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) || (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])) {
           		printf("ADC TOF2 OutA-%2d: %5d -- InA-%2d: %5d  |  TDC TOF2 OutA-%2d: %5d -- InA-%2d: %5d\n", i+1, ADC_TOF2AO[i], i+1, ADC_TOF2AI[i], i+1, TDC_TOF2AO[i], i+1, TDC_TOF2AI[i]);
         	}
         	else{
           		printf("ADC TOF2 OutA-%2d: %5d -- InA-%2d: %5d  |\n", i+1, ADC_TOF2AO[i], i+1, ADC_TOF2AI[i]);
         	}
       	}
       	else {
         	if((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) || (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])) {
           		printf("                                        |  TDC TOF2 OutA-%2d: %5d -- InA-%2d: %5d\n", i+1, TDC_TOF2AO[i], i+1, TDC_TOF2AI[i]);
         	}
       	}
   	}


   	for (Int_t i = 0; i<12; i++){
    	if (i == 6) {
        	if (ADC_TOF2BO[i] > 0 || ADC_TOF2BO[i] > 0) {
           		if ((TDC_TOF2BO[6]>TOF2BO_TDC_min[6] && TDC_TOF2BO[6]<TOF2BO_TDC_max[6]) || (TDC_TOF2BI[6]>TOF2BI_TDC_min[6] && TDC_TOF2BI[6]<TOF2BI_TDC_max[6])) {
             		printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i], i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
           		}
           		else{
             		cout << "ADC TOF2 OutB-" << i+1 << ": " << ADC_TOF2BO[i] << " -- InB-" << i+1 << ": " << ADC_TOF2BI[i] << endl;
             		printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i]);
           		}
         	}
         	else {
           		if((TDC_TOF2BO[6]>TOF2BO_TDC_min[6] && TDC_TOF2BO[6]>TOF2BO_TDC_max[6]) || (TDC_TOF2BI[6]>TOF2BI_TDC_min[6] && TDC_TOF2BI[6]>TOF2BI_TDC_max[6])) {
             		printf("                                        |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
           		}
         	}
       	}
       	else {
        	if(ADC_TOF2BO[i] > 0 || ADC_TOF2BI[i] > 0) {
           		if((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) || (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])) {
             		printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i], i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
           		}
           		else {
             		printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i]);
           		}
         	}
         	else {
           		if((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) || (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])) {
           			printf("                                        |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
           		}
         	}
       	}	
   	}

    cerr << "" << endl;
    cerr << "" << endl;
    cerr << "Final Angle (deg.) = " << angle_final << " +- " << delta_phi_deg << "  ( GUIDE )" << endl;
    cerr << "Chi Square = " << ChiS_final << endl;
    cerr << "NdF = " << ndf_final << endl;
    cerr << "Reduced Chi Square = " << ChiS_final/ndf_final << endl;
    cerr << "" << endl;
    cerr << "" << endl;

   	cerr << "Intersect with SFT Layer 1:  " << "x = " << vec_X_int_SFT[0] << "   " << "y = " << vec_Y_int_SFT[0] << endl;
   	cerr << "Intersect with TOF1:  " << "x = " << vec_X_int_TOF1[0] << "   " << "y = " << vec_Y_int_TOF1[0] << endl;
    cerr << endl;



  	printf("Length of track in target (xy plane) = %5.2f\n", Length_in_target);
    cout << endl;
    cout << endl;
    cout << endl;
}