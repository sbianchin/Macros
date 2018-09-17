#include "TArrow.h"

void _Draw_Target_Frame(int Run_Number, int ievt,
	bool has_ADC_TOF1_hit[12], bool has_TDC_TOF1_hit[12],
	bool has_ADC_TOF2_hit[12], bool has_TDC_TOF2_hit[12]){

	  TLine *hline1 = new TLine(0.38,0.14,0.62,0.14);    hline1->SetLineWidth(2);
  	TLine *hline2 = new TLine(0.30,0.18,0.70,0.18);    hline2->SetLineWidth(2);
  	TLine *hline3 = new TLine(0.26,0.22,0.74,0.22);    hline3->SetLineWidth(2);
  	TLine *hline4 = new TLine(0.22,0.26,0.78,0.26);    hline4->SetLineWidth(2);
  	TLine *hline5 = new TLine(0.18,0.30,0.82,0.30);    hline5->SetLineWidth(2);
  	TLine *hline6 = new TLine(0.18,0.34,0.82,0.34);    hline6->SetLineWidth(2);
  	TLine *hline7 = new TLine(0.14,0.38,0.86,0.38);    hline7->SetLineWidth(2);
  	TLine *hline8 = new TLine(0.14,0.42,0.86,0.42);    hline8->SetLineWidth(2);
 	  TLine *hline9 = new TLine(0.14,0.46,0.86,0.46);    hline9->SetLineWidth(2);
  	TLine *hline10 = new TLine(0.14,0.50,0.86,0.50);   hline10->SetLineWidth(2);
  	TLine *hline11 = new TLine(0.14,0.54,0.86,0.54);   hline11->SetLineWidth(2);
  	TLine *hline12 = new TLine(0.14,0.58,0.86,0.58);   hline12->SetLineWidth(2);
  	TLine *hline13 = new TLine(0.14,0.62,0.86,0.62);   hline13->SetLineWidth(2);
  	TLine *hline14 = new TLine(0.18,0.66,0.82,0.66);   hline14->SetLineWidth(2);
  	TLine *hline15 = new TLine(0.18,0.70,0.82,0.70);   hline15->SetLineWidth(2);
  	TLine *hline16 = new TLine(0.22,0.74,0.78,0.74);   hline16->SetLineWidth(2);
  	TLine *hline17 = new TLine(0.26,0.78,0.74,0.78);   hline17->SetLineWidth(2);
  	TLine *hline18 = new TLine(0.30,0.82,0.70,0.82);   hline18->SetLineWidth(2);
  	TLine *hline19 = new TLine(0.38,0.86,0.62,0.86);   hline19->SetLineWidth(2);

  	TLine *vline1 = new TLine(0.14,0.38,0.14,0.62);    vline1->SetLineWidth(2);
  	TLine *vline2 = new TLine(0.18,0.30,0.18,0.70);    vline2->SetLineWidth(2);
  	TLine *vline3 = new TLine(0.22,0.26,0.22,0.74);    vline3->SetLineWidth(2);
  	TLine *vline4 = new TLine(0.26,0.22,0.26,0.78);    vline4->SetLineWidth(2);
  	TLine *vline5 = new TLine(0.30,0.18,0.30,0.82);    vline5->SetLineWidth(2);
  	TLine *vline6 = new TLine(0.34,0.18,0.34,0.82);    vline6->SetLineWidth(2);
  	TLine *vline7 = new TLine(0.38,0.14,0.38,0.86);    vline7->SetLineWidth(2);
  	TLine *vline8 = new TLine(0.42,0.14,0.42,0.86);    vline8->SetLineWidth(2);
  	TLine *vline9 = new TLine(0.46,0.14,0.46,0.86);    vline9->SetLineWidth(2);
  	TLine *vline10 = new TLine(0.50,0.14,0.50,0.86);   vline10->SetLineWidth(2);
  	TLine *vline11 = new TLine(0.54,0.14,0.54,0.86);   vline11->SetLineWidth(2);
  	TLine *vline12 = new TLine(0.58,0.14,0.58,0.86);   vline12->SetLineWidth(2);
  	TLine *vline13 = new TLine(0.62,0.14,0.62,0.86);   vline13->SetLineWidth(2);
  	TLine *vline14 = new TLine(0.66,0.18,0.66,0.82);   vline14->SetLineWidth(2);
  	TLine *vline15 = new TLine(0.70,0.18,0.70,0.82);   vline15->SetLineWidth(2);
  	TLine *vline16 = new TLine(0.74,0.22,0.74,0.78);   vline16->SetLineWidth(2);
  	TLine *vline17 = new TLine(0.78,0.26,0.78,0.74);   vline17->SetLineWidth(2);
  	TLine *vline18 = new TLine(0.82,0.30,0.82,0.70);   vline18->SetLineWidth(2);
  	TLine *vline19 = new TLine(0.86,0.38,0.86,0.62);   vline19->SetLineWidth(2);

  	TLine *vblue1 = new TLine(0.38,0.70,0.38,0.74);      vblue1->SetLineWidth(5);   vblue1->SetLineColor(kBlue-9);
  	TLine *vblue2 = new TLine(0.42,0.62,0.42,0.66);      vblue2->SetLineWidth(5);   vblue2->SetLineColor(kBlue-9);
  	TLine *vblue3 = new TLine(0.30,0.54,0.30,0.58);      vblue3->SetLineWidth(5);   vblue3->SetLineColor(kBlue-9);
  	TLine *vblue4 = new TLine(0.86,0.50,0.86,0.54);      vblue4->SetLineWidth(5);   vblue4->SetLineColor(kBlue-9);
  	TLine *vblue5 = new TLine(0.70,0.42,0.70,0.46);      vblue5->SetLineWidth(5);   vblue5->SetLineColor(kBlue-9);
  	TLine *vblue6 = new TLine(0.58,0.34,0.58,0.38);      vblue6->SetLineWidth(5);   vblue6->SetLineColor(kBlue-9);
  	TLine *vblue7 = new TLine(0.62,0.26,0.62,0.30);      vblue7->SetLineWidth(5);   vblue7->SetLineColor(kBlue-9);
  	TLine *vblue8 = new TLine(0.62,0.14,0.62,0.18);      vblue8->SetLineWidth(5);   vblue8->SetLineColor(kBlue-9);

  	hline1->Draw();     vline1->Draw();
  	hline2->Draw();     vline2->Draw();
  	hline3->Draw();     vline3->Draw();
  	hline4->Draw();     vline4->Draw();
  	hline5->Draw();     vline5->Draw();
  	hline6->Draw();     vline6->Draw();
  	hline7->Draw();     vline7->Draw();
  	hline8->Draw();     vline8->Draw();
  	hline9->Draw();     vline9->Draw();
  	hline10->Draw();    vline10->Draw();
  	hline11->Draw();    vline11->Draw();
  	hline12->Draw();    vline12->Draw();
  	hline13->Draw();    vline13->Draw();
  	hline14->Draw();    vline14->Draw();
  	hline15->Draw();    vline15->Draw();
  	hline16->Draw();    vline16->Draw();
  	hline17->Draw();    vline17->Draw();
  	hline18->Draw();    vline18->Draw();
  	hline19->Draw();    vline19->Draw();

  	vblue2->Draw();
  	vblue4->Draw();
  	vblue6->Draw();
  	vblue8->Draw();

    TLine *TOF_line1 = new TLine(0.8278461,0.8278461,0.62,0.9478461);	TOF_line1->SetLineWidth(10);	TOF_line1->SetLineColor(17);
    TLine *TOF_line2 = new TLine(0.9478461,0.62,0.8278461,0.8278461);	TOF_line2->SetLineWidth(10);	TOF_line2->SetLineColor(17);
    TLine *TOF_line3 = new TLine(0.9478461,0.38,0.9478461,0.62);		TOF_line3->SetLineWidth(10);	TOF_line3->SetLineColor(17);
    TLine *TOF_line4 = new TLine(0.8278461,0.172154,0.9478461,0.38);	TOF_line4->SetLineWidth(10);	TOF_line4->SetLineColor(17);
    TLine *TOF_line5 = new TLine(0.62,0.0521539,0.8278461,0.172154);	TOF_line5->SetLineWidth(10);	TOF_line5->SetLineColor(17);
    TLine *TOF_line6 = new TLine(0.38,0.0521539,0.62,0.0521539);		TOF_line6->SetLineWidth(10);	TOF_line6->SetLineColor(17);
    TLine *TOF_line7 = new TLine(0.172539,0.172154,0.38,0.0521539);		TOF_line7->SetLineWidth(10);	TOF_line7->SetLineColor(17);
    TLine *TOF_line8 = new TLine(0.172539,0.172154,0.052539,0.38);		TOF_line8->SetLineWidth(10);	TOF_line8->SetLineColor(17);
    TLine *TOF_line9 = new TLine(0.052539,0.38,0.052539,0.62);			TOF_line9->SetLineWidth(10);	TOF_line9->SetLineColor(17);
	  TLine *TOF_line10 = new TLine(0.052539,0.62,0.172539,0.8278461);	TOF_line10->SetLineWidth(10);	TOF_line10->SetLineColor(17);
    TLine *TOF_line11 = new TLine(0.38,0.9478461,0.172539,0.8278461);	TOF_line11->SetLineWidth(10);	TOF_line11->SetLineColor(17);
    TLine *TOF_line12 = new TLine(0.62,0.9478461,0.38,0.9478461);		TOF_line12->SetLineWidth(10);	TOF_line12->SetLineColor(17);

	TLine *TOF_line13 = new TLine(0.8348214,0.8393338,0.625558,0.9620716);		TOF_line13->SetLineWidth(20);	TOF_line13->SetLineColor(15);
    TLine *TOF_line14 = new TLine(0.9603795,0.627551,0.8394717,0.8369272);		TOF_line14->SetLineWidth(20);	TOF_line14->SetLineColor(15);
    TLine *TOF_line15 = new TLine(0.961,0.38,0.961,0.62);                  		TOF_line15->SetLineWidth(20);	TOF_line15->SetLineColor(15);
    TLine *TOF_line16 = new TLine(0.8394717,0.1630728,0.9580543,0.372449);  	TOF_line16->SetLineWidth(20);	TOF_line16->SetLineColor(15);
    TLine *TOF_line17 = new TLine(0.6232329,0.040335,0.8324963,0.1606662);  	TOF_line17->SetLineWidth(20);	TOF_line17->SetLineColor(15);
    TLine *TOF_line18 = new TLine(0.38,0.039,0.62,0.039);                   	TOF_line18->SetLineWidth(20);	TOF_line18->SetLineColor(15);
    TLine *TOF_line19 = new TLine(0.1651786,0.1606662,0.3721168,0.040335);  	TOF_line19->SetLineWidth(20);	TOF_line19->SetLineColor(15);
    TLine *TOF_line20 = new TLine(0.1605283,0.1630728,0.03962054,0.372449); 	TOF_line20->SetLineWidth(20);	TOF_line20->SetLineColor(15);
    TLine *TOF_line21 = new TLine(0.04,0.38,0.04,0.62);                     	TOF_line21->SetLineWidth(20);	TOF_line21->SetLineColor(15);
    TLine *TOF_line22 = new TLine(0.03962054,0.6251444,0.1605283,0.8345206);	TOF_line22->SetLineWidth(20);	TOF_line22->SetLineColor(15);
    TLine *TOF_line23 = new TLine(0.3721168,0.9572584,0.1605283,0.8369272); 	TOF_line23->SetLineWidth(20);	TOF_line23->SetLineColor(15);
    TLine *TOF_line24 = new TLine(0.62,0.96,0.38,0.96);                     	TOF_line24->SetLineWidth(20);	TOF_line24->SetLineColor(15);

    /// TOF1 COLOURS
     if (has_ADC_TOF1_hit[0]) TOF_line1->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(3);

     if (has_ADC_TOF1_hit[1]) TOF_line2->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[1]) TOF_line2->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[1] && has_TDC_TOF1_hit[1]) TOF_line2->SetLineColor(3);

     if (has_ADC_TOF1_hit[2]) TOF_line3->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[2]) TOF_line3->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[2] && has_TDC_TOF1_hit[2]) TOF_line3->SetLineColor(3);

     if (has_ADC_TOF1_hit[3]) TOF_line4->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[3]) TOF_line4->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[3] && has_TDC_TOF1_hit[3]) TOF_line4->SetLineColor(3);

     if (has_ADC_TOF1_hit[4]) TOF_line5->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[4]) TOF_line5->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[4] && has_TDC_TOF1_hit[4]) TOF_line5->SetLineColor(3);

     if (has_ADC_TOF1_hit[5]) TOF_line6->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[5]) TOF_line6->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[5] && has_TDC_TOF1_hit[5]) TOF_line6->SetLineColor(3);

     if (has_ADC_TOF1_hit[6]) TOF_line7->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[6]) TOF_line7->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[6] && has_TDC_TOF1_hit[6]) TOF_line7->SetLineColor(3);

     if (has_ADC_TOF1_hit[7]) TOF_line8->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[7]) TOF_line8->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[7] && has_TDC_TOF1_hit[7]) TOF_line8->SetLineColor(3);

     if (has_ADC_TOF1_hit[8]) TOF_line9->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[8]) TOF_line9->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[8] && has_TDC_TOF1_hit[8]) TOF_line9->SetLineColor(3);

     if (has_ADC_TOF1_hit[9]) TOF_line10->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[9]) TOF_line10->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[9] && has_TDC_TOF1_hit[9]) TOF_line10->SetLineColor(3);

     if (has_ADC_TOF1_hit[10]) TOF_line11->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[10]) TOF_line11->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[10] && has_TDC_TOF1_hit[10]) TOF_line11->SetLineColor(3);

     if (has_ADC_TOF1_hit[11]) TOF_line12->SetLineColor(kOrange+10);
     if (has_TDC_TOF1_hit[11]) TOF_line12->SetLineColor(kYellow);
     if (has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) TOF_line12->SetLineColor(3);



	/// SET TOFs COLOURS
    if (has_TDC_TOF2_hit[0] || has_ADC_TOF2_hit[0]) TOF_line13->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[0] && has_ADC_TOF2_hit[0]) TOF_line13->SetLineColor(3);

    if (has_TDC_TOF2_hit[1] || has_ADC_TOF2_hit[1]) TOF_line14->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[1] && has_ADC_TOF2_hit[1]) TOF_line14->SetLineColor(3);

    if (has_TDC_TOF2_hit[2] || has_ADC_TOF2_hit[2]) TOF_line15->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[2] && has_ADC_TOF2_hit[2]) TOF_line15->SetLineColor(3);

    if (has_TDC_TOF2_hit[3] || has_ADC_TOF2_hit[3]) TOF_line16->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[3] && has_ADC_TOF2_hit[3]) TOF_line16->SetLineColor(3);

    if (has_TDC_TOF2_hit[4] || has_ADC_TOF2_hit[4]) TOF_line17->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[4] && has_ADC_TOF2_hit[4]) TOF_line17->SetLineColor(3);

    if (has_TDC_TOF2_hit[5] || has_ADC_TOF2_hit[5]) TOF_line18->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[5] && has_ADC_TOF2_hit[5]) TOF_line18->SetLineColor(3);

    if (has_TDC_TOF2_hit[6] || has_ADC_TOF2_hit[6]) TOF_line19->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[6] && has_ADC_TOF2_hit[6]) TOF_line19->SetLineColor(3);

    if (has_TDC_TOF2_hit[7] || has_ADC_TOF2_hit[7]) TOF_line20->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[7] && has_ADC_TOF2_hit[7]) TOF_line20->SetLineColor(3);

    if (has_TDC_TOF2_hit[8] || has_ADC_TOF2_hit[8]) TOF_line21->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[8] && has_ADC_TOF2_hit[8]) TOF_line21->SetLineColor(3);

    if (has_TDC_TOF2_hit[9] || has_ADC_TOF2_hit[9]) TOF_line22->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[9] && has_ADC_TOF2_hit[9]) TOF_line22->SetLineColor(3);

    if (has_TDC_TOF2_hit[10] || has_ADC_TOF2_hit[10]) TOF_line23->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[10] && has_ADC_TOF2_hit[10]) TOF_line23->SetLineColor(3);

    if (has_TDC_TOF2_hit[11] || has_ADC_TOF2_hit[11]) TOF_line24->SetLineColor(kOrange+10);
    if (has_TDC_TOF2_hit[11] && has_ADC_TOF2_hit[11]) TOF_line24->SetLineColor(3);


   	TOF_line13->Draw();		TOF_line1->Draw();
   	TOF_line14->Draw();		TOF_line2->Draw();
   	TOF_line15->Draw();		TOF_line3->Draw();
   	TOF_line16->Draw();		TOF_line4->Draw();
   	TOF_line17->Draw();		TOF_line5->Draw();
   	TOF_line18->Draw();		TOF_line6->Draw();
   	TOF_line19->Draw();		TOF_line7->Draw();
   	TOF_line20->Draw();		TOF_line8->Draw();
   	TOF_line21->Draw();		TOF_line9->Draw();
   	TOF_line22->Draw();		TOF_line10->Draw();
   	TOF_line23->Draw();		TOF_line11->Draw();
   	TOF_line24->Draw();		TOF_line12->Draw();




  	TLatex *tex_Legend_TARGET[36];
  	tex_Legend_TARGET[0] = new TLatex(0.36,0.83,"0");     	tex_Legend_TARGET[9] = new TLatex(0.09,0.47,"128");
    tex_Legend_TARGET[1] = new TLatex(0.28,0.79,"6");     	tex_Legend_TARGET[10] = new TLatex(0.09,0.43,"146");
    tex_Legend_TARGET[2] = new TLatex(0.225,0.75,"16");     tex_Legend_TARGET[11] = new TLatex(0.09,0.39,"164");
    tex_Legend_TARGET[3] = new TLatex(0.185,0.71,"28");     tex_Legend_TARGET[12] = new TLatex(0.13,0.35,"182");
    tex_Legend_TARGET[4] = new TLatex(0.145,0.67,"42");     tex_Legend_TARGET[13] = new TLatex(0.13,0.31,"198");
    tex_Legend_TARGET[5] = new TLatex(0.145,0.63,"58");     tex_Legend_TARGET[14] = new TLatex(0.17,0.27,"214");
    tex_Legend_TARGET[6] = new TLatex(0.1025,0.59,"74");    tex_Legend_TARGET[15] = new TLatex(0.21,0.23,"228");
    tex_Legend_TARGET[7] = new TLatex(0.1025,0.55,"92");    tex_Legend_TARGET[16] = new TLatex(0.25,0.19,"240");
    tex_Legend_TARGET[8] = new TLatex(0.09,0.51,"110");     tex_Legend_TARGET[17] = new TLatex(0.33,0.15,"250");

    tex_Legend_TARGET[18] = new TLatex(0.635,0.83,"5");     tex_Legend_TARGET[27] = new TLatex(0.87,0.47,"145");
    tex_Legend_TARGET[19] = new TLatex(0.71,0.79,"15");     tex_Legend_TARGET[28] = new TLatex(0.87,0.43,"163");
    tex_Legend_TARGET[20] = new TLatex(0.75,0.75,"27");     tex_Legend_TARGET[29] = new TLatex(0.87,0.39,"181");
    tex_Legend_TARGET[21] = new TLatex(0.79,0.71,"41");     tex_Legend_TARGET[30] = new TLatex(0.83,0.35,"197");
    tex_Legend_TARGET[22] = new TLatex(0.83,0.67,"57");     tex_Legend_TARGET[31] = new TLatex(0.83,0.31,"213");
    tex_Legend_TARGET[23] = new TLatex(0.83,0.63,"73");     tex_Legend_TARGET[32] = new TLatex(0.79,0.27,"227");
    tex_Legend_TARGET[24] = new TLatex(0.885,0.59,"91");    tex_Legend_TARGET[33] = new TLatex(0.75,0.23,"239");
    tex_Legend_TARGET[25] = new TLatex(0.87,0.55,"109");    tex_Legend_TARGET[34] = new TLatex(0.71,0.19,"249");
    tex_Legend_TARGET[26] = new TLatex(0.87,0.51,"127");    tex_Legend_TARGET[35] = new TLatex(0.64,0.15,"255");

   	for(int i=0; i<36; i++){
       	tex_Legend_TARGET[i]->SetTextSize(0.03);
       	tex_Legend_TARGET[i]->SetLineWidth(2);
   		tex_Legend_TARGET[i]->Draw();
   	}

    TMarker *palette_TARGET[10];
    double flag_size_palette=1.6;

    palette_TARGET[0] = new TMarker(0.54,0.075,21);     palette_TARGET[0]->SetMarkerColor(kOrange+10);    palette_TARGET[0]->SetMarkerSize(flag_size_palette);
    palette_TARGET[1] = new TMarker(0.58,0.075,21);     palette_TARGET[1]->SetMarkerColor(kOrange+7);     palette_TARGET[1]->SetMarkerSize(flag_size_palette);
    palette_TARGET[2] = new TMarker(0.62,0.075,21);     palette_TARGET[2]->SetMarkerColor(kOrange+1);     palette_TARGET[2]->SetMarkerSize(flag_size_palette);
    palette_TARGET[3] = new TMarker(0.66,0.075,21);     palette_TARGET[3]->SetMarkerColor(kOrange-4);     palette_TARGET[3]->SetMarkerSize(flag_size_palette);
    palette_TARGET[4] = new TMarker(0.70,0.075,21);     palette_TARGET[4]->SetMarkerColor(kYellow-9);     palette_TARGET[4]->SetMarkerSize(flag_size_palette);
    palette_TARGET[5] = new TMarker(0.74,0.075,21);     palette_TARGET[5]->SetMarkerColor(kYellow-7);     palette_TARGET[5]->SetMarkerSize(flag_size_palette);
    palette_TARGET[6] = new TMarker(0.78,0.075,21);     palette_TARGET[6]->SetMarkerColor(kYellow-0);     palette_TARGET[6]->SetMarkerSize(flag_size_palette);
    palette_TARGET[7] = new TMarker(0.82,0.075,21);     palette_TARGET[7]->SetMarkerColor(kSpring-4);     palette_TARGET[7]->SetMarkerSize(flag_size_palette);
    palette_TARGET[8] = new TMarker(0.86,0.075,21);     palette_TARGET[8]->SetMarkerColor(kSpring-2);     palette_TARGET[8]->SetMarkerSize(flag_size_palette);
    palette_TARGET[9] = new TMarker(0.90,0.075,21);     palette_TARGET[9]->SetMarkerColor(kGreen-0);      palette_TARGET[9]->SetMarkerSize(flag_size_palette);
  
  	for(int i=0; i<10; i++) palette_TARGET[i]->Draw();

    TLatex *tex_palette_TARGET[10];
 	int ADC_palette_TARGET[10]={0,3,6,9,12,15,18,21,24,27};
	char ADC_palette_string_TARGET[10][100];
	for(int j=0; j<10; j++) sprintf(ADC_palette_string_TARGET[j],"%d",ADC_palette_TARGET[j]);
	sprintf(ADC_palette_string_TARGET[9],"%d+",ADC_palette_TARGET[9]);

    tex_palette_TARGET[0] = new TLatex(0.510,0.04,ADC_palette_string_TARGET[0]);      tex_palette_TARGET[0]->SetTextSize(0.02);
    tex_palette_TARGET[1] = new TLatex(0.545,0.04,ADC_palette_string_TARGET[1]);      tex_palette_TARGET[1]->SetTextSize(0.02);
    tex_palette_TARGET[2] = new TLatex(0.585,0.04,ADC_palette_string_TARGET[2]);      tex_palette_TARGET[2]->SetTextSize(0.02);
    tex_palette_TARGET[3] = new TLatex(0.625,0.04,ADC_palette_string_TARGET[3]);      tex_palette_TARGET[3]->SetTextSize(0.02);
    tex_palette_TARGET[4] = new TLatex(0.665,0.04,ADC_palette_string_TARGET[4]);      tex_palette_TARGET[4]->SetTextSize(0.02);
    tex_palette_TARGET[5] = new TLatex(0.705,0.04,ADC_palette_string_TARGET[5]);      tex_palette_TARGET[5]->SetTextSize(0.02);
    tex_palette_TARGET[6] = new TLatex(0.745,0.04,ADC_palette_string_TARGET[6]);      tex_palette_TARGET[6]->SetTextSize(0.02);
    tex_palette_TARGET[7] = new TLatex(0.785,0.04,ADC_palette_string_TARGET[7]);      tex_palette_TARGET[7]->SetTextSize(0.02);
    tex_palette_TARGET[8] = new TLatex(0.825,0.04,ADC_palette_string_TARGET[8]);      tex_palette_TARGET[8]->SetTextSize(0.02);
    tex_palette_TARGET[9] = new TLatex(0.865,0.04,ADC_palette_string_TARGET[9]);      tex_palette_TARGET[9]->SetTextSize(0.02);
    TLatex *tex_palette_TARGET_scale = new TLatex(0.905,0.04,"x 100");        		  tex_palette_TARGET_scale->SetTextSize(0.02);
  	
  	for(int j=0; j<10; j++) tex_palette_TARGET[j]->Draw();
	tex_palette_TARGET_scale->Draw();

	char event_string[100];
    sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt);
    TLatex *tex_event_TARGET;
    tex_event_TARGET = new TLatex(0.036,0.0,event_string);
    tex_event_TARGET->SetTextSize(0.038);
    tex_event_TARGET->SetLineWidth(2);
    tex_event_TARGET->Draw();
}


void _Drawing_Event_Display(int Run_Number, int ievt, bool to_restore,
	  int TOF1_gap, int selected_section,
	  TGraph* gr_leptons_final, TGraph *gr_kaon_pr, TLine *line_fit_kaons,
	  TLine *restored_line,
	  int ADC_Offset, int TDC_min, int TDC_max,
	  int HG_KAON, int LG_KAON, int TDC_min_Kstop, int TDC_max_Kstop, int kaon_bar, 
	  bool has_kaon_sub, bool has_TDC_hit_Kstop[256],
	  bool has_ADC_TOF1_hit[12], bool has_TDC_TOF1_hit[12],
	  bool has_ADC_TOF2_hit[12], bool has_TDC_TOF2_hit[12],
	  bool has_TDC_hit[256], int Switch,
	  int ADC_High_TARGET[256], int ADC_Low_TARGET[256],
	  double angle_final, double delta_phi_deg,
	  double ChiS_final, int ndf_final,
	  vector<double> vec_X_int_SFT, vector<double> vec_Y_int_SFT,
	  vector<double> vec_X_int_TARGET, vector<double> vec_Y_int_TARGET,
	  vector<double> vec_X_int_TOF1, vector<double>	vec_Y_int_TOF1,
	  vector<double> vec_kstop_FINAL,
	  vector<double> vec_intersect_TARGET_final,
	  bool b_bbar, vector<double> vec_bbar){

	TOF1_Detector_Attributes();

	char ch_Title1[100];
	char ch_Title2[100];
	char ch_Title3[100];
	sprintf(ch_Title1,"ADC HG Cut");
	sprintf(ch_Title2,"ADC HG & TDC Cut");
	sprintf(ch_Title3,"K Stop Cuts");

	char ch_Subtitle1[100];
	char ch_Subtitle2[100];
	char ch_Subtitle3[100];
	sprintf(ch_Subtitle1,"(ADC offset = %d)", ADC_Offset);
	sprintf(ch_Subtitle2,"(ADC offset = %d | %d #leq TDC #leq %d)",ADC_Offset,TDC_min,TDC_max);
	sprintf(ch_Subtitle3,"(ADC: HG #geq %d, LG #geq %d | %d #leq TDC K Stop #leq %d)",HG_KAON,LG_KAON,TDC_min_Kstop,TDC_max_Kstop);

  	TLatex *tex_Title_panel1;
  	tex_Title_panel1 = new TLatex(0.01759134,0.9295171,ch_Title1);
  	tex_Title_panel1->SetTextSize(0.07645875);
  	tex_Title_panel1->SetLineWidth(2);

  	TLatex *tex_Title_panel2;
  	tex_Title_panel2 = new TLatex(0.01759134,0.9295171,ch_Title2);
  	tex_Title_panel2->SetTextSize(0.07645875);
  	tex_Title_panel2->SetLineWidth(2);

  	TLatex *tex_Title_panel3 = new TLatex(0.01759134,0.9295171,ch_Title3);
  	tex_Title_panel3->SetTextSize(0.07);
  	tex_Title_panel3->SetLineWidth(2);

  	char ch_Kbar[100];
	 sprintf(ch_Kbar, "K_{stop} Bar = %d",kaon_bar);


  	TLatex *tex_Subtitle_panel1;
  	tex_Subtitle_panel1 = new TLatex(0.01759134,0.88,ch_Subtitle1);
    tex_Subtitle_panel1->SetTextSize(0.04);
  	tex_Subtitle_panel1->SetLineWidth(2);

  	TLatex *tex_Subtitle_panel2;
  	tex_Subtitle_panel2 = new TLatex(0.01759134,0.88,ch_Subtitle2);
  	tex_Subtitle_panel2->SetTextSize(0.04);
  	tex_Subtitle_panel2->SetLineWidth(2);

	TLatex *tex_Subtitle_panel3;
	tex_Subtitle_panel3 = new TLatex(0.01759134,0.88,ch_Subtitle3);
  	tex_Subtitle_panel3->SetTextSize(0.04);
  	tex_Subtitle_panel3->SetLineWidth(2);

	TLatex *tex_Kbar;
	tex_Kbar = new TLatex(0.64,0.9295171,ch_Kbar);
    tex_Kbar->SetTextSize(0.05);
    tex_Kbar->SetLineWidth(2);
    if(has_kaon_sub) tex_Kbar->SetTextColor(6);

  	char Version[100] = "Version X";
	TLatex *tex_version;
	tex_version = new TLatex(0.75,0.0,Version);
    tex_version->SetTextSize(0.04);
	tex_version->SetLineWidth(2);

	TCanvas *c_ED;
  	c_ED = new TCanvas("Event Display","EVENT DISPLAY",0,200,1050,700);
  	c_ED->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
  	c_ED->Divide(3,2);

  	/// MAX INDEXES
  	int max_index = 0;
  	int max_index2 = 0;
  	int max_index3 = 0;
  	int max_index4 = 0;

  	int max_ADC = -100000000;
  	int max_ADC2 = -100000000;
  	int max_ADC3 = -100000000;
  	int max_ADC4 = -100000000;

  	for(Int_t q=0; q<256; q++){
    	if(ADC_High_TARGET[q]>max_ADC){
      	max_index = q;
      	max_ADC = ADC_High_TARGET[q];
    	}
  	}

  	for(Int_t q=0; q<256; q++){
    	if (q == max_index) continue;
    	else{
      		if(ADC_High_TARGET[q]>max_ADC2){
        		max_index2 = q;
        		max_ADC2 = ADC_High_TARGET[q];
      		}
    	}
  	}

  	for(Int_t q=0; q<256; q++){
    	if ((q == max_index) || (q == max_index2)) continue;
   		else{
      		if(ADC_High_TARGET[q]>max_ADC3) {
       			max_index3 = q;
        		max_ADC3 = ADC_High_TARGET[q];
      		}
    	}
  	}

  	for(Int_t q=0; q<256; q++){
    	if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;
   		else{
      		if(ADC_High_TARGET[q]>max_ADC4) {
        		max_index4 = q;
        		max_ADC4 = ADC_High_TARGET[q];
      		}
    	}
  	}



  	//TARGET MARKERS
	TMarker *marker_ADC_TARGET[256];   
	TMarker *marker_ADCL_TARGET[256];   
    double flag_size_TARGET=1.35;

   	for(Int_t i1=0; i1<6; i1++) marker_ADC_TARGET[i1] = new TMarker(0.40+0.04*i1,0.84,21);
    for(Int_t i2=0; i2<10; i2++) marker_ADC_TARGET[i2+6] = new TMarker(0.32+0.04*i2,0.80,21);
    for(Int_t i3=0; i3<12; i3++) marker_ADC_TARGET[i3+16] = new TMarker(0.28+0.04*i3,0.76,21);
    for(Int_t i4=0; i4<14; i4++) marker_ADC_TARGET[i4+28] = new TMarker(0.24+0.04*i4,0.72,21);
    for(Int_t i5=0; i5<16; i5++) marker_ADC_TARGET[i5+42] = new TMarker(0.20+0.04*i5,0.68,21);
    for(Int_t i6=0; i6<16; i6++) marker_ADC_TARGET[i6+58] = new TMarker(0.20+0.04*i6,0.64,21);
    for(Int_t i7=0; i7<18; i7++) marker_ADC_TARGET[i7+74] = new TMarker(0.16+0.04*i7,0.60,21);
    for(Int_t i8=0; i8<18; i8++) marker_ADC_TARGET[i8+92] = new TMarker(0.16+0.04*i8,0.56,21);
    for(Int_t i9=0; i9<18; i9++) marker_ADC_TARGET[i9+110] = new TMarker(0.16+0.04*i9,0.52,21);
    for(Int_t i10=0; i10<18; i10++) marker_ADC_TARGET[i10+128] = new TMarker(0.16+0.04*i10,0.48,21);
    for(Int_t i11=0; i11<18; i11++) marker_ADC_TARGET[i11+146] = new TMarker(0.16+0.04*i11,0.44,21);
    for(Int_t i12=0; i12<18; i12++) marker_ADC_TARGET[i12+164] = new TMarker(0.16+0.04*i12,0.40,21);
    for(Int_t i13=0; i13<16; i13++) marker_ADC_TARGET[i13+182] = new TMarker(0.20+0.04*i13,0.36,21);
    for(Int_t i14=0; i14<16; i14++) marker_ADC_TARGET[i14+198] = new TMarker(0.20+0.04*i14,0.32,21);
    for(Int_t i15=0; i15<14; i15++) marker_ADC_TARGET[i15+214] = new TMarker(0.24+0.04*i15,0.28,21);
    for(Int_t i16=0; i16<12; i16++) marker_ADC_TARGET[i16+228] = new TMarker(0.28+0.04*i16,0.24,21);
    for(Int_t i17=0; i17<10; i17++) marker_ADC_TARGET[i17+240] = new TMarker(0.32+0.04*i17,0.20,21);
    for(Int_t i18=0; i18<6; i18++) marker_ADC_TARGET[i18+250] = new TMarker(0.40+0.04*i18,0.16,21);
    for(Int_t iSize=0; iSize<256; iSize++)  marker_ADC_TARGET[iSize]->SetMarkerSize(flag_size_TARGET);


	for(Int_t i1=0; i1<6; i1++) marker_ADCL_TARGET[i1] = new TMarker(0.40+0.04*i1,0.84,21);     // TO MOVE
    for(Int_t i2=0; i2<10; i2++) marker_ADCL_TARGET[i2+6] = new TMarker(0.32+0.04*i2,0.80,21);
    for(Int_t i3=0; i3<12; i3++) marker_ADCL_TARGET[i3+16] = new TMarker(0.28+0.04*i3,0.76,21);
    for(Int_t i4=0; i4<14; i4++) marker_ADCL_TARGET[i4+28] = new TMarker(0.24+0.04*i4,0.72,21);
    for(Int_t i5=0; i5<16; i5++) marker_ADCL_TARGET[i5+42] = new TMarker(0.20+0.04*i5,0.68,21);
    for(Int_t i6=0; i6<16; i6++) marker_ADCL_TARGET[i6+58] = new TMarker(0.20+0.04*i6,0.64,21);
    for(Int_t i7=0; i7<18; i7++) marker_ADCL_TARGET[i7+74] = new TMarker(0.16+0.04*i7,0.60,21);
    for(Int_t i8=0; i8<18; i8++) marker_ADCL_TARGET[i8+92] = new TMarker(0.16+0.04*i8,0.56,21);
    for(Int_t i9=0; i9<18; i9++) marker_ADCL_TARGET[i9+110] = new TMarker(0.16+0.04*i9,0.52,21);
    for(Int_t i10=0; i10<18; i10++) marker_ADCL_TARGET[i10+128] = new TMarker(0.16+0.04*i10,0.48,21);
    for(Int_t i11=0; i11<18; i11++) marker_ADCL_TARGET[i11+146] = new TMarker(0.16+0.04*i11,0.44,21);
    for(Int_t i12=0; i12<18; i12++) marker_ADCL_TARGET[i12+164] = new TMarker(0.16+0.04*i12,0.40,21);
    for(Int_t i13=0; i13<16; i13++) marker_ADCL_TARGET[i13+182] = new TMarker(0.20+0.04*i13,0.36,21);
    for(Int_t i14=0; i14<16; i14++) marker_ADCL_TARGET[i14+198] = new TMarker(0.20+0.04*i14,0.32,21);
    for(Int_t i15=0; i15<14; i15++) marker_ADCL_TARGET[i15+214] = new TMarker(0.24+0.04*i15,0.28,21);
    for(Int_t i16=0; i16<12; i16++) marker_ADCL_TARGET[i16+228] = new TMarker(0.28+0.04*i16,0.24,21);
    for(Int_t i17=0; i17<10; i17++) marker_ADCL_TARGET[i17+240] = new TMarker(0.32+0.04*i17,0.20,21);
    for(Int_t i18=0; i18<6; i18++) marker_ADCL_TARGET[i18+250] = new TMarker(0.40+0.04*i18,0.16,21);
    for(Int_t iSize=0; iSize<256; iSize++)  marker_ADCL_TARGET[iSize]->SetMarkerSize(flag_size_TARGET);


  	char Angle_string[30];
    char ChiS_string[30];
    //char ChiS_COS[30];

    sprintf(Angle_string,"#phi = %3.2f#circ#pm%3.2f", angle_final, delta_phi_deg);
    sprintf(ChiS_string,"#chi^{2}/ndf = %3.2f", ChiS_final/ndf_final);

    TLatex *tex_Angle;
    TLatex *tex_ChiS;
    //TLatex *tex_ChiS_COS;

    tex_Angle = new TLatex(-45.,43.,Angle_string);
    tex_Angle->SetTextSize(0.05);
    tex_Angle->SetLineWidth(2);

    tex_ChiS = new TLatex(-45.,37.,ChiS_string);
    tex_ChiS->SetTextSize(0.05);
    tex_ChiS->SetLineWidth(2);


    TLatex *x_sft;
    TLatex *y_sft;
    TLatex *x_target;
    TLatex *y_target;

    char X_SFT_String[30];
    char Y_SFT_String[30];
    char X_TARGET_String[30];
    char Y_TARGET_String[30];


    sprintf(X_SFT_String,"X(41.6) = %4.2f", vec_X_int_SFT[0]);    //  X' = Y  ;  Y' = -X
    sprintf(Y_SFT_String,"Y(41.6) = %4.2f", vec_Y_int_SFT[0]);     //  X' = Y  ;  Y' = -X
    sprintf(X_TARGET_String,"X(29) = %4.2f", vec_X_int_TARGET[0]);               //  X' = Y  ;  Y' = -X
    sprintf(Y_TARGET_String,"Y(29) = %4.2f", vec_Y_int_TARGET[0]);               //  X' = Y  ;  Y' = -X


    x_sft = new TLatex(10.,42.,X_SFT_String);
    x_sft->SetTextSize(0.05);
    x_sft->SetLineWidth(2);

    y_sft = new TLatex(10.,37.,Y_SFT_String);
    y_sft->SetTextSize(0.05);
    y_sft->SetLineWidth(2);

    x_target = new TLatex(-45.,-40.,X_TARGET_String);
    x_target->SetTextSize(0.05);
    x_target->SetLineWidth(2);

    y_target = new TLatex(-45.,-45.,Y_TARGET_String);
    y_target->SetTextSize(0.05);
    y_target->SetLineWidth(2);


    char str_Kstop_X[200];
    char str_Kstop_Y[200];

    sprintf(str_Kstop_X,"X_{Ks} = %5.2f",vec_kstop_FINAL[0]);
    sprintf(str_Kstop_Y,"Y_{Ks} = %5.2f",vec_kstop_FINAL[1]);

    TLatex *tex_Kstop_X = new TLatex(18.,-40.,str_Kstop_X);
    tex_Kstop_X->SetTextSize(0.05);
    tex_Kstop_X->SetLineWidth(2);

    TLatex *tex_Kstop_Y = new TLatex(18.,-45.,str_Kstop_Y);
    tex_Kstop_Y->SetTextSize(0.05);
    tex_Kstop_Y->SetLineWidth(2);


    /////
    TGraph *gr_int_SFT = new TGraph(vec_Y_int_SFT.size(),&vec_X_int_SFT[0],&vec_Y_int_SFT[0]);
    gr_int_SFT->SetMarkerStyle(20);
    gr_int_SFT->SetMarkerColor(4);
    gr_int_SFT->SetMarkerSize(0.8);

    TGraph *gr_int_TARGET = new TGraph(vec_X_int_TARGET.size(),&vec_X_int_TARGET[0],&vec_Y_int_TARGET[0]);
    gr_int_TARGET->SetMarkerStyle(20);
    gr_int_TARGET->SetMarkerColor(4);
    gr_int_TARGET->SetMarkerSize(0.8);

    TGraph *gr_int_TOF1 = new TGraph(vec_X_int_TOF1.size(),&vec_X_int_TOF1[0],&vec_Y_int_TOF1[0]);
    gr_int_TOF1->SetMarkerStyle(20);
    gr_int_TOF1->SetMarkerColor(4);
    gr_int_TOF1->SetMarkerSize(0.8);


    /// BLUE ARROWS (GUIDE)
    TArrow *X_guide;
    TArrow *Y_guide;

    double X_arrow = -1.;
    double Y_arrow = -1.;
  	int Axis_Vector_Length = 10;

    X_arrow = vec_intersect_TARGET_final[2];
    Y_arrow = vec_intersect_TARGET_final[3];

    X_guide = new TArrow(X_arrow, Y_arrow, X_arrow + Axis_Vector_Length, Y_arrow, 0.005, "|>");
    Y_guide = new TArrow(X_arrow, Y_arrow, X_arrow, Y_arrow + Axis_Vector_Length, 0.005, "|>");

    X_guide->SetLineWidth(2);     Y_guide->SetLineWidth(2);
    X_guide->SetLineColor(4);     Y_guide->SetLineColor(4);


    /// IF NO LEPTON LEFT
    cout << "////// BBABR : " << b_bbar << endl;
    TLine *line_bbar;
    if(b_bbar){
    	//TLine *line_bbar;
		line_bbar = new TLine(vec_bbar[0], vec_bbar[1], vec_bbar[2], vec_bbar[3]);  
		line_bbar->SetLineWidth(2);
      	line_bbar->SetLineColor(3);
  	}



  	TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  	TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");


  	/// TARGET CIRCLE
  	TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0); 
 	ell_Target->SetFillStyle(0);
  	ell_Target->SetLineColor(1);
  	ell_Target->SetLineWidth(1);

  	/// TOF1 CIRCLE
  	TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  	ell->SetFillStyle(0);
  	ell->SetLineColor(6);
  	ell->SetLineWidth(1);

  	/// SFT CIRCLE
  	TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  	ell_L1->SetFillStyle(0);
  	ell_L1->SetLineColor(4);
  	ell_L1->SetLineWidth(1);


  	/// TARGET CENTER
  	vector<double> vec_xx_Target_Center;  vec_xx_Target_Center.clear();
  	vector<double> vec_yy_Target_Center;  vec_yy_Target_Center.clear();
  	vec_xx_Target_Center.push_back(0.);
  	vec_yy_Target_Center.push_back(0.);

  	TGraph *gr_Target_Center = new TGraph(vec_xx_Target_Center.size(),&vec_xx_Target_Center[0],&vec_yy_Target_Center[0]);
  	gr_Target_Center->SetMarkerStyle(5);
  	gr_Target_Center->SetMarkerColor(1);
  	gr_Target_Center->SetMarkerSize(1);

  	///TOF1
  	vector<double> vec_x_TOF1;			vec_x_TOF1.clear();
  	vector<double> vec_y_TOF1;			vec_y_TOF1.clear();
  	vector<double> vec_x_TOF1_selected;	vec_x_TOF1_selected.clear();
  	vector<double> vec_y_TOF1_selected;	vec_y_TOF1_selected.clear();

  	if(selected_section != -1){
  		for(int i=0; i<5; i++){
  			if(i != selected_section){
  				vec_x_TOF1.push_back(TOF1_Xloc[TOF1_gap-1][i]);
  				vec_y_TOF1.push_back(TOF1_Yloc[TOF1_gap-1][i]);
  			}
  			else{
  				vec_x_TOF1_selected.push_back(TOF1_Xloc[TOF1_gap-1][i]);
  				vec_y_TOF1_selected.push_back(TOF1_Yloc[TOF1_gap-1][i]);
  			}
  		}
  	}

  	TGraph *gr_TOF1;
  	gr_TOF1 = new TGraph(vec_x_TOF1.size(), &vec_x_TOF1[0], &vec_y_TOF1[0]);
    gr_TOF1->SetMarkerStyle(20);
    gr_TOF1->SetMarkerColor(2);
    gr_TOF1->SetMarkerSize(1.5);

  	TGraph *gr_TOF1_selected;
  	gr_TOF1_selected = new TGraph(vec_x_TOF1_selected.size(), &vec_x_TOF1_selected[0], &vec_y_TOF1_selected[0]);
    gr_TOF1_selected->SetMarkerStyle(20);
    gr_TOF1_selected->SetMarkerColor(3);
    gr_TOF1_selected->SetMarkerSize(1.5);

  	/// LEPTON FINAL FIT
    char str_lepton_fit_final[100]; 	
    sprintf(str_lepton_fit_final,"Lepton Fit  |  Run %d  --  Event %d",Run_Number,ievt);
  	TH2F *h_lepton_fit_final;
    h_lepton_fit_final = new TH2F("Lepton Fit", str_lepton_fit_final, 500, -50, 50, 500, -50, 50);

    gr_leptons_final->SetMarkerStyle(21);
  	gr_leptons_final->SetMarkerColor(2);
  	gr_leptons_final->SetMarkerSize(0.8);
  	gr_leptons_final->GetXaxis()->SetLimits(-50.,50.);
  	gr_leptons_final->GetYaxis()->SetRangeUser(-50.,50.);


    /// KAON FINAL FIT
    char str_kaon_fit_final[100];
    sprintf(str_kaon_fit_final,"Kaon Fit  |  Run %d  --  Event %d",Run_Number,ievt);
  	TH2F *h_kaon_fit_final;
    h_kaon_fit_final = new TH2F("Kaon Fit", str_kaon_fit_final, 500, -50, 50, 500, -50, 50);

    /// KSTOP FINAL
    vector<double>	vec_x_kstop_final;	vec_x_kstop_final.clear();
    vector<double>	vec_y_kstop_final;	vec_y_kstop_final.clear();
    vec_x_kstop_final.push_back(vec_kstop_FINAL[0]);
    vec_y_kstop_final.push_back(vec_kstop_FINAL[1]);

    TGraph *gr_kstop_final;
    gr_kstop_final = new TGraph(vec_x_kstop_final.size(), &vec_x_kstop_final[0], &vec_y_kstop_final[0]);
  	gr_kstop_final->SetMarkerStyle(34);
  	gr_kstop_final->SetMarkerColor(1);
  	gr_kstop_final->SetMarkerSize(1.3);

    /// FINAL FIT !
    char str_final[100];
    sprintf(str_final,"Final Fit  |  Run %d  --  Event %d",Run_Number,ievt);
	TH2F *h_final;  
    h_final = new TH2F("Final Fit", str_final, 500, -50, 50, 500, -50, 50);

  	//////////////////////////////////////////////////////////////////////////////

    /// PANEL 1
    c_ED->cd(1);
    _Draw_Target_Frame(Run_Number, ievt,
    	has_ADC_TOF1_hit, has_TDC_TOF1_hit,
    	has_ADC_TOF2_hit, has_TDC_TOF2_hit
    	);
	tex_Title_panel1->Draw();
	tex_Subtitle_panel1->Draw();

  	for(Int_t icol=0; icol<256; icol++){
  		if(ADC_High_TARGET[icol]>=0 && ADC_High_TARGET[icol]<300){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=300 && ADC_High_TARGET[icol]<600){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=600 && ADC_High_TARGET[icol]<900){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=900 && ADC_High_TARGET[icol]<1200){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=1200 && ADC_High_TARGET[icol]<1500){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=1500 && ADC_High_TARGET[icol]<1800){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=1800 && ADC_High_TARGET[icol]<2100){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=2100 && ADC_High_TARGET[icol]<2400){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=2400 && ADC_High_TARGET[icol]<2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADC_TARGET[icol]->Draw();}
    	if(ADC_High_TARGET[icol]>=2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kGreen-0);	marker_ADC_TARGET[icol]->Draw();}

    	if(ADC_High_TARGET[icol]<0 && ADC_Low_TARGET[icol]>0){
      		if(has_TDC_hit[icol]) marker_ADC_TARGET[icol]->SetMarkerColor(11);
     	 	else marker_ADC_TARGET[icol]->SetMarkerColor(kBlue);
      		marker_ADC_TARGET[icol]->Draw();
    	}
  	}

  	if(ADC_High_TARGET[max_index] > 0){
   		marker_ADC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    	marker_ADC_TARGET[max_index]->Draw();
  	}

  	if(ADC_High_TARGET[max_index2] > 0){
    	marker_ADC_TARGET[max_index2]->SetMarkerColor(1);
    	marker_ADC_TARGET[max_index2]->Draw();
  	}

  	if(ADC_High_TARGET[max_index3] > 0){
    	marker_ADC_TARGET[max_index3]->SetMarkerColor(1);
    	marker_ADC_TARGET[max_index3]->Draw();
  	}

  	if(ADC_High_TARGET[max_index4] > 0){
    	marker_ADC_TARGET[max_index4]->SetMarkerColor(1);
    	marker_ADC_TARGET[max_index4]->Draw();
  	}


    /// PANEL 2
    c_ED->cd(2);
    _Draw_Target_Frame(Run_Number, ievt,
    	has_ADC_TOF1_hit, has_TDC_TOF1_hit,
    	has_ADC_TOF2_hit, has_TDC_TOF2_hit
    	);
  	tex_Title_panel2->Draw();
	  tex_Subtitle_panel2->Draw(); 	

  	for(Int_t icol=0; icol<256; icol++){           
    	if(has_TDC_hit[icol]){
      		if(ADC_High_TARGET[icol]>=0 && ADC_High_TARGET[icol]<300){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=300 && ADC_High_TARGET[icol]<600){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=600 && ADC_High_TARGET[icol]<900){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=900 && ADC_High_TARGET[icol]<1200){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=1200 && ADC_High_TARGET[icol]<1500){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=1500 && ADC_High_TARGET[icol]<1800){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=1800 && ADC_High_TARGET[icol]<2100){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=2100 && ADC_High_TARGET[icol]<2400){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=2400 && ADC_High_TARGET[icol]<2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADC_TARGET[icol]->Draw();}
      		if(ADC_High_TARGET[icol]>=2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kGreen-0); marker_ADC_TARGET[icol]->Draw();}
    	}

    	if(ADC_High_TARGET[icol]<0 && ADC_Low_TARGET[icol]>0 && TDC_min_TARGET && Switch==1){
      	if(has_TDC_hit[icol]) marker_ADC_TARGET[icol]->SetMarkerColor(11);
      		else marker_ADC_TARGET[icol]->SetMarkerColor(kBlue);
      		marker_ADC_TARGET[icol]->Draw();
    	}
  	}
  	
  	//if(ADC_High_TARGET[max_index] > 0){
    if(has_TDC_hit[max_index]){
   		marker_ADC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    	marker_ADC_TARGET[max_index]->Draw();
  	}

  	//if(ADC_High_TARGET[max_index2] > 0){
    if(has_TDC_hit[max_index2]){
    	marker_ADC_TARGET[max_index2]->SetMarkerColor(1);
    	marker_ADC_TARGET[max_index2]->Draw();
  	}

  	//if(ADC_High_TARGET[max_index3] > 0){
    if(has_TDC_hit[max_index3]){
    	marker_ADC_TARGET[max_index3]->SetMarkerColor(1);
    	marker_ADC_TARGET[max_index3]->Draw();
  	}

  	//if(ADC_High_TARGET[max_index4] > 0){
    if(has_TDC_hit[max_index4]){
    	marker_ADC_TARGET[max_index4]->SetMarkerColor(1);
    	marker_ADC_TARGET[max_index4]->Draw();
  	}


    
    /// PANEL 3
    c_ED->cd(3);
    _Draw_Target_Frame(Run_Number, ievt,
    	has_ADC_TOF1_hit, has_TDC_TOF1_hit,
    	has_ADC_TOF2_hit, has_TDC_TOF2_hit
    	);
  	
    tex_Title_panel3->Draw();
	  tex_Subtitle_panel3->Draw();
    tex_Kbar->Draw(); 
    tex_version->Draw();

    if(has_kaon_sub){
      marker_ADCL_TARGET[kaon_bar]->SetMarkerColor(6); 
      marker_ADCL_TARGET[kaon_bar]->Draw();
    }
    else{
      for(Int_t icol=0; icol<256; icol++){
        if(has_TDC_hit_Kstop[icol]){
          if(ADC_Low_TARGET[icol]>=LG_KAON && ADC_High_TARGET[icol]>=HG_KAON){
        		marker_ADCL_TARGET[icol]->SetMarkerColor(12); marker_ADCL_TARGET[icol]->Draw();
      	  }
    	  }
  	  }
    }


    /// PANEL 4
  	c_ED->cd(4);
  	h_lepton_fit_final->Draw();
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
  	
  	if(selected_section !=-1){
  		gr_TOF1->Draw("sameP");
  		gr_TOF1_selected->Draw("sameP");
  	}

  	gr_leptons_final->Draw("sameP");
  	if(to_restore) restored_line->Draw("same");

  	/// PANEL 5
  	c_ED->cd(5);
  	h_kaon_fit_final->Draw();
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

    gr_kaon_pr->Draw("sameP");
  	line_fit_kaons->Draw("same");


  	/// PANEL 6
    c_ED->cd(6);  
    h_final->Draw();
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

  	if(selected_section !=-1){
  		gr_TOF1->Draw("sameP");
  		gr_TOF1_selected->Draw("sameP");
  	}

  	gr_leptons_final->Draw("sameP");
  	if(to_restore) restored_line->Draw("same");
  	gr_kaon_pr->Draw("sameP");
  	gr_kstop_final->Draw("sameP");

  	gr_int_TARGET->Draw("sameP");
  	gr_int_SFT->Draw("sameP");
  	gr_int_TOF1->Draw("sameP");

  	tex_Angle->Draw("same");
  	tex_ChiS->Draw("same");
  	x_sft->Draw("same");
  	y_sft->Draw("same");
  	x_target->Draw("same");
  	y_target->Draw("same");
  	tex_Kstop_X->Draw("same");
  	tex_Kstop_Y->Draw("same");

  	X_guide->Draw();
  	Y_guide->Draw();

  	if(b_bbar) line_bbar->Draw("same");
}

