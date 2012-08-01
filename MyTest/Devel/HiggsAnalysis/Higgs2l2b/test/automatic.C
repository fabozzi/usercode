
void automatic(string var = "zllmass", string masspoint="500", string label="m_{ll}", bool set_logScale=false, float intLumi= 349.){
  //  gROOT->Reset();
  gROOT->SetStyle("Plain");

  enum selType{muChannel, elChannel};

  /*    START CUSTOMIZATION */

  //  string masspoint = "500";
  string output = masspoint +".txt";  
  //string label = "m_{ll} [GeV]";
  // string var = "zllmass";
  //  float intLumi = 349.;
  int rebin_fact = 2;
  // set log scale for plots
  // bool set_logScale = true;
  string lumi = "1 fb^{-1}";
  // set to true in order to set the range on the X axis
  bool custom = true;

  /*    END CUSTOMIZATION */

  if(custom){
    float rMin = 50.;
    float rMax = 150.;
}
  int ymax =20.;
  if(set_logScale) ymax = 1000.;
  string postfix  = "EdmNtp.root";
  ofstream f(output.c_str());
  float y_mu =200.;
  float y_el =200.;

  float x_mu = 50.;  
  float x_el = 50.;  
  int max = 100;
  int min = 800;
  float br_fact = 0.333;
  //float normFact = 1/1000;
  float normFact = intLumi/1000.;

  /* 250 */
  if(masspoint=="250"){
    int cut_value1=237;
    int cut_value2=260;
  }
  else if(masspoint=="300"){
    int cut_value1=280;
    int cut_value2=323;
  }
  else if(masspoint=="350"){
    int cut_value1=380;
    int cut_value2=330;
  }
  else if(masspoint=="4000"){
    int cut_value1=390;
    int cut_value2=460;
  }
  else if(masspoint=="450"){
    int cut_value1=420;
    int cut_value2=550;
  }
  else if(masspoint=="500"){
    int cut_value1=470;
    int cut_value2=1000;
  }
  
  vector<int> thisSel;
  thisSel.push_back(muChannel);
  thisSel.push_back(elChannel);

  // set selection type:
  //  selType this_selection = muChannel;
  // set rebin factor for histos


  // Set scale factors for lumi normalization
  float scaleZbb[]={0.00838281, 0.01604, 0.057394, 0.02507414};
  // right factors are those below, the previous ones are used just because a crash during the skim   
  //  float scaleZbb[]={0.00838281, 0.00812031, 0.057394, 0.02507414};
  float scaleZcc[]={0.0066633, 0.00883544, 0.05828496, 0.02760666};
  float scaleZjets_mu[]={1.6000847, 1.55144941, 1.060303407, 0.168596, 0.3089263, 0.11791534, 0.21228856, 0.15464988, 0.07783832, 0.15061466, 0.00185248, 0.00352303, 0.00567954, 0.01118749, 0.00511575, 0.00001615, 0.00007081, 0.00005568, 0.00004569, 0.00002192};
  float scaleZjets_el[]={0.0338, 0.098855, 0.1083411, 0.14836314,0.06671835, 0.00093808, 0.01061008, 0.02191539,0.020801, 0.0264994, 0.00001681, 0.00025568, 0.00072630,0.00135979, 0.00162018,0.00000036, 0.00000745, 0.00002377, 0.00003277, 0.00000393};
  float scaleZ0Jets =5.118;
  float scaleTT =0.0165;
  float scaleWZ =0.00863;
  float scaleZZ =0.002798;
  float scaleWW =0.0208;
  float scaleH400 =0.0005031;
  float scaleH450 =0.00032118;
  float scaleH350 = 0.00070501;
  float scaleH300 =0.00073510;
  float scaleH250 =0.00096653;
  float scaleH500 =0.00021778;

  float scaleFact250 =  scaleH250*normFact;
  float scaleFact300 =  scaleH300*normFact;
  float scaleFact350 =  scaleH350*normFact;
  float scaleFact400 =  scaleH400*normFact;
  float scaleFact450 =  scaleH450*normFact;
  float scaleFact500 =  scaleH500*normFact;
  
  float scaleFactZbb0 =  scaleZbb[1]*normFact;
  float scaleFactZbb1 =  scaleZbb[2]*normFact;
  float scaleFactZbb2 =  scaleZbb[3]*normFact;
  float scaleFactZbb3 =  scaleZbb[4]*normFact;

  float scaleFactZcc0 =  scaleZcc[0]*normFact;
  float scaleFactZcc1 =  scaleZcc[1]*normFact;
  float scaleFactZcc2 =  scaleZcc[2]*normFact;
  float scaleFactZcc3 =  scaleZcc[3]*normFact;

  float scaleFactZZ =  scaleZZ*normFact;
  float scaleFactWZ =  scaleWZ*normFact;
  float scaleFactTT =  scaleTT*normFact;
  float scaleFactWW =  scaleWW*normFact;
  
  float scaleFactZ0jet =  scaleZ0Jets*normFact*br_fact;

  float scaleFactZ1jet_mu =  scaleZjets_mu[0]*normFact*br_fact;
  float scaleFactZ2jet_mu =  scaleZjets_mu[1]*normFact*br_fact;
  float scaleFactZ3jet_mu =  scaleZjets_mu[2]*normFact*br_fact;
  float scaleFactZ4jet_mu =  scaleZjets_mu[3]*normFact*br_fact;
  float scaleFactZ5jet_mu =  scaleZjets_mu[4]*normFact*br_fact;

  float scaleFactZ1jet100_300_mu = scaleZjets_mu[5]*normFact*br_fact; 
  float scaleFactZ2jet100_300_mu = scaleZjets_mu[6]*normFact*br_fact; 
  float scaleFactZ3jet100_300_mu = scaleZjets_mu[7]*normFact*br_fact; 
  float scaleFactZ4jet100_300_mu = scaleZjets_mu[8]*normFact*br_fact; 
  float scaleFactZ5jet100_300_mu = scaleZjets_mu[9]*normFact*br_fact;

  float scaleFactZ1jet300_800_mu = scaleZjets_mu[10]*normFact*br_fact; 
  float scaleFactZ2jet300_800_mu = scaleZjets_mu[11]*normFact*br_fact; 
  float scaleFactZ3jet300_800_mu = scaleZjets_mu[12]*normFact*br_fact; 
  float scaleFactZ4jet300_800_mu = scaleZjets_mu[13]*normFact*br_fact; 
  float scaleFactZ5jet300_800_mu = scaleZjets_mu[14]*normFact*br_fact; 

  float scaleFactZ1jet800_1600_mu = scaleZjets_mu[15]*normFact*br_fact;
  float scaleFactZ2jet800_1600_mu = scaleZjets_mu[16]*normFact*br_fact;
  float scaleFactZ3jet800_1600_mu = scaleZjets_mu[17]*normFact*br_fact;
  float scaleFactZ4jet800_1600_mu = scaleZjets_mu[18]*normFact*br_fact;
  float scaleFactZ5jet800_1600_mu = scaleZjets_mu[19]*normFact*br_fact;


  float scaleFactZ1jet_el =  scaleZjets_el[0]*normFact*br_fact;
  float scaleFactZ2jet_el =  scaleZjets_el[1]*normFact*br_fact;
  float scaleFactZ3jet_el =  scaleZjets_el[2]*normFact*br_fact;
  float scaleFactZ4jet_el =  scaleZjets_el[3]*normFact*br_fact;
  float scaleFactZ5jet_el =  scaleZjets_el[4]*normFact*br_fact;

  float scaleFactZ1jet100_300_el = scaleZjets_el[5]*normFact*br_fact; 
  float scaleFactZ2jet100_300_el = scaleZjets_el[6]*normFact*br_fact; 
  float scaleFactZ3jet100_300_el = scaleZjets_el[7]*normFact*br_fact; 
  float scaleFactZ4jet100_300_el = scaleZjets_el[8]*normFact*br_fact; 
  float scaleFactZ5jet100_300_el = scaleZjets_el[9]*normFact*br_fact;

  float scaleFactZ1jet300_800_el = scaleZjets_el[10]*normFact*br_fact; 
  float scaleFactZ2jet300_800_el = scaleZjets_el[11]*normFact*br_fact; 
  float scaleFactZ3jet300_800_el = scaleZjets_el[12]*normFact*br_fact; 
  float scaleFactZ4jet300_800_el = scaleZjets_el[13]*normFact*br_fact; 
  float scaleFactZ5jet300_800_el = scaleZjets_el[14]*normFact*br_fact; 

  float scaleFactZ1jet800_1600_el = scaleZjets_el[15]*normFact*br_fact;
  float scaleFactZ2jet800_1600_el = scaleZjets_el[16]*normFact*br_fact;
  float scaleFactZ3jet800_1600_el = scaleZjets_el[17]*normFact*br_fact;
  float scaleFactZ4jet800_1600_el = scaleZjets_el[18]*normFact*br_fact;
  float scaleFactZ5jet800_1600_el = scaleZjets_el[19]*normFact*br_fact;


  for(int i=0; i<thisSel.size();++i){
    selType this_selection = thisSel[i];

    // INPUT DIR for the histo files
    string input_dir_base = "NoNorm_Mu/";
    if(this_selection == elChannel)
      input_dir_base = "NoNorm_El/";
    
    // OUTPUT DIR for the histo files
    string out_dir_base = "NoNorm_Mu/";
    if(this_selection == elChannel)
      out_dir_base = "NoNorm_El/";

    // OUTPUT .png file
    string name_eps = var+"_Lin_Mu.png";
    if(set_logScale) name_eps = var+"_Log_Mu.png";
    if(this_selection == elChannel){
      name_eps = var+"_Lin_El.png";
      if(set_logScale) name_eps = var+"_Log_El.png";    
    }


  // Print on the screen what selection you are applying
    if(this_selection == muChannel){
      cout << "####################### COMBINING MU HISTOS ##########################" << endl << endl;
      //     TLatex *t= new TLatex(x_mu,y_mu,"CMS Preliminary   ");
    }


    if(this_selection == elChannel){
      cout << "##################### COMBINING ELECTRON HISTOS ######################" << endl << endl;
      //TLatex *t= new TLatex(x_el,y_el,"CMS Preliminary");
    }
  
    //    t->SetTextFont(72);
    // t->SetTextSize(0.06);

    string create_out_dir = "mkdir "+out_dir_base;
    gSystem->Exec( create_out_dir.c_str() );
    
  /*h250*/

  /*h300*/
  string infile_name = input_dir_base+masspoint+postfix;
  TFile f_h500(infile_name.c_str());
  TH1F * var_h500_Raw = (TH1F*) f_h500.Get(var.c_str());
  TH1F * var_h500 = new TH1F(*var_h500_Raw);
  if(masspoint=="200") var_h500->Scale(scaleFact200);
  else if(masspoint=="250") var_h500->Scale(scaleFact250);
  else if(masspoint=="300") var_h500->Scale(scaleFact300);
  else if(masspoint=="350") var_h500->Scale(scaleFact350);
  else if(masspoint=="400") var_h500->Scale(scaleFact400);
  else if(masspoint=="450") var_h500->Scale(scaleFact450);
  else if(masspoint=="500") var_h500->Scale(scaleFact500);

  cout << "h500 raw entries" << var_h500_Raw->Integral() << endl;
  cout << "h500 1 fb entries" << var_h500->Integral() << endl;

  TH1F * var_All500 = new TH1F(*var_h500_Raw);
  var_All500->Reset();

  /*infile_name = input_dir_base+"h450VBF"+postfix;
  TFile f_h450vbf(infile_name.c_str());
  TH1F * var_h450VBF_Raw = (TH1F*) f_h450vbf.Get(var.c_str());

  TH1F * var_h450VBF = new TH1F(*var_h450VBF_Raw);
  var_h450VBF->Scale(scaleFactVBF450);
  cout << "h450VBF raw entries " << var_h450VBF_Raw->Integral() << endl;
  cout << "h450VBF 1 fb entries " << var_h450VBF->Integral() << endl;

  TH1F * var_h450 = new TH1F(*var_h450GF_Raw);
  var_h450->Add(var_h450GF,var_h450VBF);
  */
  cout << "---> h450ALL 1 fb entries " << var_h500->Integral() << endl;

  /*
  infile_name = input_dir_base+"data"+postfix;
  TFile f_data(infile_name.c_str());
  TH1F * var_data_Raw = (TH1F*) f_data.Get(var.c_str());

  TH1F * var_data = new TH1F(*var_data_Raw);
  cout << "data raw entries " << var_data_Raw->Integral() << endl;
  cout << "---> data "+lumi+" fb entries " << var_data->Integral() << endl;
  */


  infile_name = input_dir_base+"TT"+postfix;
  TFile f_tt(infile_name.c_str());
  TH1F * var_TT_Raw = (TH1F*) f_tt.Get(var.c_str());

  TH1F * var_TT = new TH1F(*var_TT_Raw);
  var_TT->Scale(scaleFactTT);
  cout << "TT raw entries " << var_TT_Raw->Integral() << endl;
  cout << "---> TT 1 fb entries " << var_TT->Integral() << endl;


  infile_name = input_dir_base+"ZZ"+postfix;
  TFile f_zz(infile_name.c_str());
  TH1F * var_ZZ_Raw = (TH1F*) f_zz.Get(var.c_str());
  TH1F * var_ZZ = new TH1F(*var_ZZ_Raw);
  var_ZZ->Scale(scaleFactZZ);
  cout << "ZZ raw entries " << var_ZZ_Raw->Integral() << endl;
  cout << "---> ZZ 1 fb entries " << var_ZZ->Integral() << endl;

  infile_name = input_dir_base+"WZ"+postfix;
  TFile f_wz(infile_name.c_str());
  TH1F * var_WZ_Raw = (TH1F*) f_wz.Get(var.c_str());
  TH1F * var_WZ = new TH1F(*var_WZ_Raw);
  var_WZ->Scale(scaleFactWZ);
  cout << "WZ raw entries " << var_WZ_Raw->Integral() << endl;
  cout << "---> WZ 1 fb entries " << var_WZ->Integral() << endl;

  infile_name = input_dir_base+"WW"+postfix;
  TFile f_ww(infile_name.c_str());
  TH1F * var_WW_Raw = (TH1F*) f_ww.Get(var.c_str());
  TH1F * var_WW = new TH1F(*var_WW_Raw);
  var_WW->Scale(scaleFactWW);
  cout << "WW raw entries " << var_WW_Raw->Integral() << endl;
  cout << "---> WW 1 fb entries " << var_WW->Integral() << endl;


  infile_name = input_dir_base+"ZBB0"+postfix;
  TFile f_zbb0(infile_name.c_str());
  TH1F * var_ZBB0_Raw = (TH1F*) f_zbb0.Get(var.c_str());
  TH1F * var_ZBB0 = new TH1F(*var_ZBB0_Raw);
  var_ZBB0->Scale(scaleFactZbb0);

  infile_name = input_dir_base+"ZBB1"+postfix;
  TFile f_zbb1(infile_name.c_str());
  TH1F * var_ZBB1_Raw = (TH1F*) f_zbb1.Get(var.c_str());
  TH1F * var_ZBB1 = new TH1F(*var_ZBB1_Raw);
  var_ZBB1->Scale(scaleFactZbb1);

  infile_name = input_dir_base+"ZBB2"+postfix;
  TFile f_zbb2(infile_name.c_str());
  TH1F * var_ZBB2_Raw = (TH1F*) f_zbb2.Get(var.c_str());
  TH1F * var_ZBB2 = new TH1F(*var_ZBB2_Raw);
  var_ZBB2->Scale(scaleFactZbb2);

  infile_name = input_dir_base+"ZBB3"+postfix;
  TFile f_zbb3(infile_name.c_str());
  TH1F * var_ZBB3_Raw = (TH1F*) f_zbb3.Get(var.c_str());
  TH1F * var_ZBB3 = new TH1F(*var_ZBB3_Raw);
  var_ZBB3->Scale(scaleFactZbb3);

  TH1F * var_ZBB = new TH1F(*var_ZBB0);
  var_ZBB->Add( var_ZBB1);
  var_ZBB->Add( var_ZBB2);
  var_ZBB->Add( var_ZBB3);

  cout << "---> ZBB 1 fb entries " << var_ZBB->Integral() << endl;

  infile_name = input_dir_base+"ZCC0"+postfix;
  TFile f_zcc0(infile_name.c_str());
  TH1F * var_ZCC0_Raw = (TH1F*) f_zcc0.Get(var.c_str());
  TH1F * var_ZCC0 = new TH1F(*var_ZCC0_Raw);
  var_ZCC0->Scale(scaleFactZcc0);

  infile_name = input_dir_base+"ZCC1"+postfix;
  TFile f_zcc1(infile_name.c_str());
  TH1F * var_ZCC1_Raw = (TH1F*) f_zcc1.Get(var.c_str());
  TH1F * var_ZCC1 = new TH1F(*var_ZCC1_Raw);
  var_ZCC1->Scale(scaleFactZcc1);

  infile_name = input_dir_base+"ZCC2"+postfix;
  TFile f_zcc2(infile_name.c_str());
  TH1F * var_ZCC2_Raw = (TH1F*) f_zcc2.Get(var.c_str());
  TH1F * var_ZCC2 = new TH1F(*var_ZCC2_Raw);
  var_ZCC2->Scale(scaleFactZcc2);

  infile_name = input_dir_base+"ZCC3"+postfix;
  TFile f_zcc3(infile_name.c_str());
  TH1F * var_ZCC3_Raw = (TH1F*) f_zcc3.Get(var.c_str());
  TH1F * var_ZCC3 = new TH1F(*var_ZCC3_Raw);
  var_ZCC3->Scale(scaleFactZcc3);

  TH1F * var_ZCC = new TH1F(*var_ZCC0);
  var_ZCC->Add(var_ZCC1);
  var_ZCC->Add(var_ZCC2);
  var_ZCC->Add(var_ZCC3);

  cout << "---> ZCC 1 fb entries " << var_ZCC->Integral() << endl;

  TH1F * var_ZBBZCC = new TH1F(*var_ZBB);
  var_ZBBZCC->Add(var_ZCC);
  cout << "---> ZBBZCC 1 fb entries " << var_ZBBZCC->Integral() << endl;


  infile_name = input_dir_base+"Z0Jets"+postfix;
  TFile f_z0jet(infile_name.c_str());
  TH1F * var_Z0JET_Raw = (TH1F*) f_z0jet.Get(var.c_str());
  TH1F * var_Z0JET = new TH1F(*var_Z0JET_Raw);
  cout<<"number events Z1Jet"<<endl;
  var_Z0JET->Scale(scaleFactZ0jet);
  

  infile_name = input_dir_base+"Z1Jets"+postfix;
  TFile f_z1jet(infile_name.c_str());
  TH1F * var_Z1JET_Raw = (TH1F*) f_z1jet.Get(var.c_str());
  TH1F * var_Z1JET = new TH1F(*var_Z1JET_Raw);

  if(this_selection == muChannel) var_Z1JET->Scale(scaleFactZ1jet_mu);
  if(this_selection == elChannel) var_Z1JET->Scale(scaleFactZ1jet_el);
  cout<<"number events Z1Jet "<<var_Z1JET->Integral()<<endl;


  infile_name = input_dir_base+"Z2Jets"+postfix;
  TFile f_z2jet(infile_name.c_str());
  TH1F * var_Z2JET_Raw = (TH1F*) f_z2jet.Get(var.c_str());
  TH1F * var_Z2JET = new TH1F(*var_Z2JET_Raw);
  if(this_selection == muChannel) var_Z2JET->Scale(scaleFactZ2jet_mu);
  if(this_selection == elChannel) var_Z2JET->Scale(scaleFactZ2jet_el);

  infile_name = input_dir_base+"Z3Jets"+postfix;
  TFile f_z3jet(infile_name.c_str());
  TH1F * var_Z3JET_Raw = (TH1F*) f_z3jet.Get(var.c_str());
  TH1F * var_Z3JET = new TH1F(*var_Z3JET_Raw);
  if(this_selection == muChannel) var_Z3JET->Scale(scaleFactZ3jet_mu);
  if(this_selection == elChannel) var_Z3JET->Scale(scaleFactZ3jet_el);

  infile_name = input_dir_base+"Z4Jets"+postfix;
  TFile f_z4jet(infile_name.c_str());
  TH1F * var_Z4JET_Raw = (TH1F*) f_z4jet.Get(var.c_str());
  TH1F * var_Z4JET = new TH1F(*var_Z4JET_Raw);
  if(this_selection == muChannel) var_Z4JET->Scale(scaleFactZ4jet_mu);
  if(this_selection == elChannel) var_Z4JET->Scale(scaleFactZ4jet_el);

  infile_name = input_dir_base+"Z5Jets"+postfix;
  TFile f_z5jet(infile_name.c_str());
  TH1F * var_Z5JET_Raw = (TH1F*) f_z5jet.Get(var.c_str());
  TH1F * var_Z5JET = new TH1F(*var_Z5JET_Raw);
  if(this_selection == muChannel) var_Z5JET->Scale(scaleFactZ5jet_mu);
  if(this_selection == elChannel) var_Z5JET->Scale(scaleFactZ5jet_el);


  infile_name = input_dir_base+"Z1Jets100_300"+postfix;
  TFile f_z1jet100_300(infile_name.c_str());
  TH1F * var_Z1JET100_300_Raw = (TH1F*) f_z1jet100_300.Get(var.c_str());
  TH1F * var_Z1JET100_300 = new TH1F(*var_Z1JET100_300_Raw);
  if(this_selection == muChannel) var_Z1JET100_300->Scale(scaleFactZ1jet100_300_mu);
  if(this_selection == elChannel) var_Z1JET100_300->Scale(scaleFactZ1jet100_300_el);
  cout<<"number events Z1Jet100_300 "<<var_Z1JET100_300->Integral()<<endl;

  infile_name = input_dir_base+"Z2Jets100_300"+postfix;
  TFile f_z2jet100_300(infile_name.c_str());
  TH1F * var_Z2JET100_300_Raw = (TH1F*) f_z2jet100_300.Get(var.c_str());
  TH1F * var_Z2JET100_300 = new TH1F(*var_Z2JET100_300_Raw);
  if(this_selection == muChannel) var_Z2JET100_300->Scale(scaleFactZ2jet100_300_mu);
  if(this_selection == elChannel) var_Z2JET100_300->Scale(scaleFactZ2jet100_300_el);

  infile_name = input_dir_base+"Z3Jets100_300"+postfix;
  TFile f_z3jet100_300(infile_name.c_str());
  TH1F * var_Z3JET100_300_Raw = (TH1F*) f_z3jet100_300.Get(var.c_str());
  TH1F * var_Z3JET100_300 = new TH1F(*var_Z3JET100_300_Raw);
  if(this_selection == muChannel) var_Z3JET100_300->Scale(scaleFactZ3jet100_300_mu);
  if(this_selection == elChannel) var_Z3JET100_300->Scale(scaleFactZ3jet100_300_el);

  infile_name = input_dir_base+"Z4Jets100_300"+postfix;
  TFile f_z4jet100_300(infile_name.c_str());
  TH1F * var_Z4JET100_300_Raw = (TH1F*) f_z4jet100_300.Get(var.c_str());
  TH1F * var_Z4JET100_300 = new TH1F(*var_Z4JET100_300_Raw);
  if(this_selection == muChannel) var_Z4JET100_300->Scale(scaleFactZ4jet100_300_mu);
  if(this_selection == elChannel) var_Z4JET100_300->Scale(scaleFactZ4jet100_300_el);

  infile_name = input_dir_base+"Z5Jets100_300"+postfix;
  TFile f_z5jet100_300(infile_name.c_str());
  TH1F * var_Z5JET100_300_Raw = (TH1F*) f_z5jet100_300.Get(var.c_str());
  TH1F * var_Z5JET100_300 = new TH1F(*var_Z5JET100_300_Raw);
  if(this_selection == muChannel) var_Z5JET100_300->Scale(scaleFactZ5jet100_300_mu);
  if(this_selection == elChannel) var_Z5JET100_300->Scale(scaleFactZ5jet100_300_el);

  infile_name = input_dir_base+"Z1Jets300_800"+postfix;
  TFile f_z1jet300_800(infile_name.c_str());
  TH1F * var_Z1JET300_800_Raw = (TH1F*) f_z1jet300_800.Get(var.c_str());
  TH1F * var_Z1JET300_800 = new TH1F(*var_Z1JET300_800_Raw);
  if(this_selection == muChannel) var_Z1JET300_800->Scale(scaleFactZ1jet300_800_mu);
  if(this_selection == elChannel) var_Z1JET300_800->Scale(scaleFactZ1jet300_800_el);
  cout<<"number events Z1Jet300_800 "<<var_Z1JET300_800->Integral()<<endl;

  infile_name = input_dir_base+"Z2Jets300_800"+postfix;
  TFile f_z2jet300_800(infile_name.c_str());
  TH1F * var_Z2JET300_800_Raw = (TH1F*) f_z2jet300_800.Get(var.c_str());
  TH1F * var_Z2JET300_800 = new TH1F(*var_Z2JET300_800_Raw);
  if(this_selection == muChannel) var_Z2JET300_800->Scale(scaleFactZ2jet300_800_mu);
  if(this_selection == elChannel) var_Z2JET300_800->Scale(scaleFactZ2jet300_800_el);

  infile_name = input_dir_base+"Z3Jets300_800"+postfix;
  TFile f_z3jet300_800(infile_name.c_str());
  TH1F * var_Z3JET300_800_Raw = (TH1F*) f_z3jet300_800.Get(var.c_str());
  TH1F * var_Z3JET300_800 = new TH1F(*var_Z3JET300_800_Raw);
  if(this_selection == muChannel) var_Z3JET300_800->Scale(scaleFactZ3jet300_800_mu);
  if(this_selection == elChannel) var_Z3JET300_800->Scale(scaleFactZ3jet300_800_el);

  infile_name = input_dir_base+"Z4Jets300_800"+postfix;
  TFile f_z4jet300_800(infile_name.c_str());
  TH1F * var_Z4JET300_800_Raw = (TH1F*) f_z4jet300_800.Get(var.c_str());
  TH1F * var_Z4JET300_800 = new TH1F(*var_Z4JET300_800_Raw);
  if(this_selection == muChannel) var_Z4JET300_800->Scale(scaleFactZ4jet300_800_mu);
  if(this_selection == elChannel) var_Z4JET300_800->Scale(scaleFactZ4jet300_800_el);

  infile_name = input_dir_base+"Z5Jets300_800"+postfix;
  TFile f_z5jet300_800(infile_name.c_str());
  TH1F * var_Z5JET300_800_Raw = (TH1F*) f_z5jet300_800.Get(var.c_str());
  TH1F * var_Z5JET300_800 = new TH1F(*var_Z5JET300_800_Raw);
  if(this_selection == muChannel) var_Z5JET300_800->Scale(scaleFactZ5jet300_800_mu);
  if(this_selection == elChannel) var_Z5JET300_800->Scale(scaleFactZ5jet300_800_el);


  infile_name = input_dir_base+"Z1Jets800_1600"+postfix;
  TFile f_z1jet800_1600(infile_name.c_str());
  TH1F * var_Z1JET800_1600_Raw = (TH1F*) f_z1jet800_1600.Get(var.c_str());
  TH1F * var_Z1JET800_1600 = new TH1F(*var_Z1JET800_1600_Raw);
  if(this_selection == muChannel)   var_Z1JET800_1600->Scale(scaleFactZ1jet800_1600_mu);
  if(this_selection == elChannel)   var_Z1JET800_1600->Scale(scaleFactZ1jet800_1600_el);
  cout<<"number events Z1Jet800_1600 "<<var_Z1JET800_1600->Integral()<<endl;

  infile_name = input_dir_base+"Z2Jets800_1600"+postfix;
  TFile f_z2jet800_1600(infile_name.c_str());
  TH1F * var_Z2JET800_1600_Raw = (TH1F*) f_z2jet800_1600.Get(var.c_str());
  TH1F * var_Z2JET800_1600 = new TH1F(*var_Z2JET800_1600_Raw);
  if(this_selection == muChannel)   var_Z2JET800_1600->Scale(scaleFactZ2jet800_1600_mu);
  if(this_selection == elChannel)   var_Z2JET800_1600->Scale(scaleFactZ2jet800_1600_el);

  infile_name = input_dir_base+"Z3Jets800_1600"+postfix;
  TFile f_z3jet800_1600(infile_name.c_str());
  TH1F * var_Z3JET800_1600_Raw = (TH1F*) f_z3jet800_1600.Get(var.c_str());
  TH1F * var_Z3JET800_1600 = new TH1F(*var_Z3JET800_1600_Raw);
  if(this_selection == muChannel)   var_Z3JET800_1600->Scale(scaleFactZ3jet800_1600_mu);
  if(this_selection == elChannel)   var_Z3JET800_1600->Scale(scaleFactZ3jet800_1600_el);

  infile_name = input_dir_base+"Z4Jets800_1600"+postfix;
  TFile f_z4jet800_1600(infile_name.c_str());
  TH1F * var_Z4JET800_1600_Raw = (TH1F*) f_z4jet800_1600.Get(var.c_str());
  TH1F * var_Z4JET800_1600 = new TH1F(*var_Z4JET800_1600_Raw);
  if(this_selection == muChannel)   var_Z4JET800_1600->Scale(scaleFactZ4jet800_1600_mu);
  if(this_selection == elChannel)   var_Z4JET800_1600->Scale(scaleFactZ4jet800_1600_el);

  infile_name = input_dir_base+"Z5Jets800_1600"+postfix;
  TFile f_z5jet800_1600(infile_name.c_str());
  TH1F * var_Z5JET800_1600_Raw = (TH1F*) f_z5jet800_1600.Get(var.c_str());
  TH1F * var_Z5JET800_1600 = new TH1F(*var_Z5JET800_1600_Raw);
  if(this_selection == muChannel)   var_Z5JET800_1600->Scale(scaleFactZ5jet800_1600_mu);
  if(this_selection == elChannel)   var_Z5JET800_1600->Scale(scaleFactZ5jet800_1600_el);

  TH1F * var_ZJET = new TH1F(*var_Z0JET);
    var_ZJET->Add( var_Z1JET);
  var_ZJET->Add( var_Z2JET);
  var_ZJET->Add( var_Z3JET);
  var_ZJET->Add( var_Z4JET);
  var_ZJET->Add( var_Z5JET);

  var_ZJET->Add( var_Z1JET100_300);
  var_ZJET->Add( var_Z2JET100_300);
  var_ZJET->Add( var_Z3JET100_300);
  var_ZJET->Add( var_Z4JET100_300);
  var_ZJET->Add( var_Z5JET100_300);

  var_ZJET->Add( var_Z1JET300_800);
  var_ZJET->Add( var_Z2JET300_800);
  var_ZJET->Add( var_Z3JET300_800);
  var_ZJET->Add( var_Z4JET300_800);
  var_ZJET->Add( var_Z5JET300_800);

  var_ZJET->Add( var_Z1JET800_1600);
  var_ZJET->Add( var_Z2JET800_1600);
  var_ZJET->Add( var_Z3JET800_1600);
  var_ZJET->Add( var_Z4JET800_1600);
  var_ZJET->Add( var_Z5JET800_1600);
  
  /*  var_ZJET->Add( var_Z1JET);
  var_ZJET->Add( var_Z1JET100_300);
  var_ZJET->Add( var_Z1JET300_800);
  var_ZJET->Add( var_Z1JET800_1600);
  */
  cout << "---> ZJet 1 fb entries after add" << var_ZJET->Integral() << endl;

  infile_name = input_dir_base+"MuRun2011A"+postfix;
  TFile f_dataMu(infile_name.c_str());
  TH1F * data_Mu = (TH1F*) f_dataMu.Get(var.c_str());
  infile_name = input_dir_base+"ElRun2011A"+postfix;
  TFile f_dataEl(infile_name.c_str());
  TH1F * data_El = (TH1F*) f_dataEl.Get(var.c_str());
  data_Mu->Add(data_El);
  cout<<"data"<<endl;

  // Eventually REBIN:
  if(rebin_fact > 1)
    {    
      //var_All250->Rebin(rebin_fact);
      //var_All300->Rebin(rebin_fact);
    var_All500->Rebin(rebin_fact);
    var_WW->Rebin(rebin_fact); 
    var_WZ->Rebin(rebin_fact); 
    var_TT->Rebin(rebin_fact); 
    var_ZCC->Rebin(rebin_fact); 
    var_ZBB->Rebin(rebin_fact); 
    var_ZJET->Rebin(rebin_fact); 
    var_ZZ->Rebin(rebin_fact); 
    //var_h250->Rebin(rebin_fact); 
    //var_h300->Rebin(rebin_fact); 
    var_h500->Rebin(rebin_fact); 
    data_Mu->Rebin(rebin_fact);
  }


  int n = var_h500->GetNbinsX();

  min = var_h500->GetXaxis()->GetXmin();
  max = var_h500->GetXaxis()->GetXmax();
  cout<<"min and max: "<<min<<" "<<max<<endl;
  cout<<"nBin: "<<n<<endl;
  float binsize = (max - min)/n;

  int bincut1 = (cut_value1 -min )/binsize +1;
  int bincut2 = (cut_value2 -min) /binsize +1;
  cout<<"bin cut "<<bincut1<<" --> "<<var_h500->GetBinCenter(bincut1)<<endl;
  cout<<"bin cut "<<bincut2<<" --> "<<var_h500->GetBinCenter(bincut2)<<endl;
 
  if(this_selection == elChannel)  f<<"Electron Channel "<<endl;
  if(this_selection == muChannel)  f<<"Muon Channel "<<endl;
  cout << "---> Higgs "<<masspoint<<" 1 fb entries in Signal region " << var_h500->Integral(bincut1,bincut2) << endl;
  cout << "---> ZJet 1 fb entries in Signal region " << var_ZJET->Integral(bincut1,bincut2) << endl;
  cout << "---> ZZ 1 fb entries in Signal region " << var_ZZ->Integral(bincut1,bincut2) << endl;
  cout << "---> WZ 1 fb entries in Signal region " << var_WZ->Integral(bincut1,bincut2) << endl;
  cout << "---> WW 1 fb entries in Signal region " << var_WW->Integral(bincut1,bincut2) << endl;
  cout << "---> Zbb 1 fb entries in Signal region " << var_ZBB->Integral(bincut1,bincut2) << endl;
  cout << "---> Zcc 1 fb entries in Signal region " << var_ZCC->Integral(bincut1,bincut2) << endl;
  cout << "---> TT 1 fb entries in Signal region " << var_TT->Integral(bincut1,bincut2) << endl;
 
  //  TH1F * S0 = new TH1F("S0","H mass",nbins,min,max);
  //  S0->SetMinimum(minimum);
  //  S0->SetMaximum(maximum);
  //  S0->SetTitle(var.c_str());

  // sum histos for stack plot
  TH1F * S1 = new TH1F(* var_WW);
  TH1F * S12 = new TH1F(* S1);
  //  TH1F * S12 = new TH1F(* var_TT);
  S12->Add(var_WZ);
  TH1F * S123 = new TH1F(* S12);
  //  TH1F * S12 = new TH1F(* var_TT);
  S123->Add(var_ZZ);
  TH1F * S1234 = new TH1F(* S123);
  //  TH1F * S123 = new TH1F(* var_ZBB);
  S1234->Add(var_TT);
  TH1F * S12345 = new TH1F(* S1234);
  //  TH1F * S123 = new TH1F(* var_ZBB);
  S12345->Add(var_ZBB);
  TH1F * S123456 = new TH1F(* S12345);
  //  TH1F * S1234 = new TH1F(* var_ZJET);
  S123456->Add(var_ZCC);
  TH1F * S1234567 = new TH1F(* S123456);
  //  TH1F * S12345 = new TH1F(* var_ZZ);
  S1234567->Add(var_ZJET);
  //  TH1F * S1234567 = new TH1F(* S123456);
  //  S1234567->Add(var_h250);

  //TH1F * S12345678 = new TH1F(* S123456);
  //  S12345678->Add(var_h500);

  TH1F * S12345678 = new TH1F(* S1234567);
  //  TH1F * S123456 = new TH1F(* var_h500);
  S12345678->Add(var_h500);


  // write data histo, stack histo and single histos into the output file

  TH1F *data = new TH1F("data","data", var_h500->GetNbinsX(), 100, 1000);
  data->SetName("data_obs");
  data->SetTitle("data_obs");

  
  if(this_selection == muChannel)  string outfile_name = out_dir_base+"Histo_"+masspoint+"_mu.root";
  if(this_selection == elChannel)  string outfile_name = out_dir_base+"Histo_"+masspoint+"_el.root";
  TFile outfile(outfile_name.c_str(), "RECREATE");
  //var_data->Write();

  data_Mu->SetMarkerStyle(22);
  data_Mu->Write();
  S12345678->SetName("data_obs");
  S12345678->Write();
  //  data->SetName("data_obs");
  var_WZ->SetName("WZtoany");
  var_WZ->Write();
  var_WW->SetName("WWtoany");
  var_WW->Write();
  var_TT->SetName("ttbar2L2Nu2B");
  var_TT->Write();
  var_ZBB->SetName("zbbbar");
  var_ZBB->Write();
  var_ZCC->SetName("zccbar");
  var_ZCC->Write();
  var_ZJET->SetName("zjets");
  var_ZJET->Write();
  var_ZZ->SetName("ZZtoany");
  var_ZZ->Write();
  //var_h250->Write();
  //var_h500->Write();
  var_h500->SetName("higgs");
  var_h500->Write();
  outfile.Close();



  // Draw stack plot into an .png file

  S12345678->SetStats(kFALSE);
  S1234567->SetStats(kFALSE);
  S123456->SetStats(kFALSE);
  S12345->SetStats(kFALSE);
  S1234->SetStats(kFALSE);
  S123->SetStats(kFALSE);
  S12->SetStats(kFALSE);
  S1->SetStats(kFALSE);

  S12345678->SetLineWidth(2);
  S1234567->SetLineWidth(2);
  S123456->SetLineWidth(2);
  S12345->SetLineWidth(2);
  S1234->SetLineWidth(2);
  S123->SetLineWidth(2);
  S12->SetLineWidth(2);
  S1->SetLineWidth(2);

  S1->SetFillColor(kYellow-7);
 //    S1->SetLineColor(kMagenta+4);
  S12->SetFillColor(kMagenta -3);
 //    S12->SetLineColor(kTeal+3);
  S123->SetFillColor(kTeal-5);
 //    S123->SetLineColor(kViolet+3);
  S1234->SetFillColor(kMagenta +3);
 //    S1234->SetLineColor(kPink-7);
  S12345->SetFillColor(kViolet+4);
 //    S12345->SetLineColor(kOrange+8);
  S123456->SetFillColor(kBlue-7);
  S1234567->SetFillColor(kOrange-2);
  S12345678->SetFillColor(kOrange+7);
  //  S123456->SetLineColor(kPink+3);


  //  TH1F * S1234567 = new TH1F(* S123456);
  //  S1234567->Add(h_svbf500);
  //  TH1F * S12345678 = new TH1F(* S1234567);
  //  S12345678->Add(h_svbf450);
 

  leg = new TLegend(.75,.75,1.0,1.0);
  leg->SetFillColor(0);
  leg->AddEntry(S1,"WW","f");
  leg->AddEntry(S12,"WZ","f");
  leg->AddEntry(S123,"ZZ","f");
  leg->AddEntry(S1234,"t #bar{t}","f");
  leg->AddEntry(S12345,"Zbb","f");
  leg->AddEntry(S123456,"Zcc","f");
  leg->AddEntry(S1234567,"ZJets","f");
  leg->AddEntry(S12345678,("H"+masspoint).c_str(),"f");
  //leg->AddEntry(S12345678,"H250 ","f");
  //leg->AddEntry(S12345678,"H500 ","f");
  //leg->AddEntry( var_data_Sumch,"data ","P");
  leg->AddEntry(data_Mu,"data ","P");

  float ymax_histosum = S12345678->GetMaximum();
  float ymax_data = data_Mu->GetMaximum();


  var_All500->SetMinimum(0.00001);
  if(ymax_histosum > ymax_data) var_All500->SetMaximum(ymax_histosum+ymax);
  else var_All500->SetMaximum(ymax_data+ymax);

  var_All500->SetTitle("");

  //  if(this_selection == muChannel)  var_All500->SetTitle("MuonChannel");
  //  if(this_selection == elChannel)  var_All500->SetTitle("ElectronChannel");


  if(this_selection == muChannel)  var_All500->GetXaxis()->SetTitle((label+" (#mu^{+} #mu^{-} channel)").c_str());
  if(this_selection == elChannel)  var_All500->GetXaxis()->SetTitle((label+" (e^{+} e^{-} channel)").c_str());

  float Xmin = var_h500->GetXaxis()->GetXmin();
  float Xmax = var_h500->GetXaxis()->GetXmax();

  float nEvts = (Xmax-Xmin)/n;
  cout<<"events"<<nEvts<<endl;
  //string Ylabel = "number of events/"+nEvts.str()+" GeV"
  var_All500->GetYaxis()->SetTitle(Form("number of events/ %.0f GeV",nEvts));
  if(custom) var_All500->GetXaxis()->SetRangeUser(rMin,rMax);
  //var_data->SetMarkerStyle(21);

  var_All500->Draw();
  S12345678->Draw("HISTsame");
  S1234567->Draw("HISTsame");
  S123456->Draw("HISTsame");
  S12345->Draw("HISTsame");
  S1234->Draw("HISTsame");
  S123->Draw("HISTsame");
  S12->Draw("HISTsame");
  S1->Draw("HISTsame");
  zllmass->Draw("HISTsame");
  data_Mu->Draw("*same");

  if(set_logScale){
    c1->SetLogy();
    var_All500->SetMinimum(1.);
  }
  leg->Draw("SAME");
  // t->Draw();


    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);

    latex.SetTextAlign(31); // align right
    //    latex.DrawLatex(0.65,0.90,"#sqrt{s} = 7 TeV");
    if (intLumi > 0.) {
      latex.SetTextAlign(31); // align right
      //      latex.DrawLatex(0.24,0.80,Form("#int #font[12]{L} dt = %.0f nb^{-1}",intLumi));
    }
    latex.SetTextAlign(11); // align left
    latex.DrawLatex(0.12, 0.92, "CMS preliminary 2011");
    latex.DrawLatex(0.45,0.92,Form("%.0f pb^{-1} at #sqrt{s} = 7 TeV",intLumi));

  c1->SaveAs(name_eps.c_str());


  }

}
