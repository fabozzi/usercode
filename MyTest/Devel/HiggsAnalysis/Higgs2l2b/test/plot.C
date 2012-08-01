double integralInRange( TH1F * var_h, float lower, float upper ){
  TAxis *axis = var_h->GetXaxis();
  int bmin = axis->FindBin(lower);
  int bmax = axis->FindBin(upper);
  double integr = var_h->Integral(bmin,bmax);
  integr -= ( (var_h->GetBinContent(bmin)) * (lower - (axis->GetBinLowEdge(bmin))) )/
    axis->GetBinWidth(bmin);
  integr -= ( (var_h->GetBinContent(bmax)) * ( (axis->GetBinUpEdge(bmax)) - upper) )/
    axis->GetBinWidth(bmax);
  return integr;
}

void plot(string runPeriod = "A", string var = "zllmass",  string label="m_{ll}", int rebin_fact = 1, float ymax = 50., bool set_logScale=false, bool custom = false, bool drawSig = false, string masspoint="500", bool drawSigStack = false, float lowcut = 70., float upcut = 110.){

  //  gROOT->Reset();

  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gROOT->SetStyle("tdrStyle");

  //  gROOT->SetStyle("Plain");

  enum selType{muChannel, elChannel};

  string output = masspoint;  
  //  int rebin_fact = 1;
  // set log scale for plots
  string lumi = "1 fb^{-1}";
  // set to true in order to set the range on the X axis

  // lumi for datat taking periods
  const float lumiA = 2132.6;
  const float lumiB = 2508.0;
    
  if(custom && var == "zllmass"){
    float rMin = 60.;
    float rMax = 120.;
  }
  if(custom && var == "npv"){
    float rMin = -0.5;
    float rMax = 24.5;
  }
  if(set_logScale) ymax = 1000.;
  string postfix  = "EdmNtp.root";
  ofstream f(("reports/"+output+var+".txt").c_str());
  float y_mu =200.;
  float y_el =200.;

  float x_mu = 50.;  
  float x_el = 50.;  
  int max = 100;
  int min = 800;
  float br_fact = 0.333;
  //float normFact = 1/1000;

  float intLumi = 1000.0;
  if(runPeriod == "A")
    intLumi = lumiA;
  else if(runPeriod == "B")
    intLumi = lumiB;
  else if(runPeriod == "All")
    intLumi = lumiA+lumiB;
  else
    cout << "Period not existent!!! Normalization set to 1 fb-1" << endl;
    
  float normFact = intLumi/1000.;

  //  if(this_selection == muChannel) 
  //    normFact *= 47531./46514.;
  //  if(this_selection == elChannel) 
  //    normFact *= 43508./40853.;

  float cut_value1f = 0.94 * atof(masspoint.c_str());
  float cut_value2f = 1.10 * atof(masspoint.c_str());
  if( (var!="lljjmass") && (var!="lljjmass_1btag") && (var!="lljjmass_0btag") ){
    cut_value1f = lowcut;
    cut_value2f = upcut;
  }

  f << "s.r. = [" << cut_value1f <<", " << cut_value2f << "]" << endl;
  f << endl;
  f.precision(4);

  // counters to store e+mu final yields  
  double integralSignalCombined(0), allSignalCombined(0);
  double integralDYCombined(0), allDYCombined(0);
  double integralTTCombined(0), allTTCombined(0);
  double integralZZCombined(0), allZZCombined(0);
  double integralWZCombined(0), allWZCombined(0);
  double integralWWCombined(0), allWWCombined(0);
  double integralDataCombined(0), allDataCombined(0);
  
  vector<int> thisSel;
  thisSel.push_back(muChannel);
  thisSel.push_back(elChannel);

  float scaleDYJets = 0.0840179;
  float scaleTT =0.0165;
  float scaleWZ =0.00863;
  float scaleZZ =0.002798;
  float scaleWW =0.0208;
  float scaleH200 =0.00074280;
  float scaleH250 =0.00052212;
  float scaleH300 =0.00039763;
  float scaleH350 =0.00037046;
  float scaleH400 =0.00027725;
  float scaleH450 =0.00017129;
  float scaleH500 =0.00011800;

  float scaleFactDYJets = scaleDYJets*normFact;
  float scaleFactTT =  scaleTT*normFact;
  float scaleFactWZ =  scaleWZ*normFact;
  float scaleFactZZ =  scaleZZ*normFact;
  float scaleFactWW =  scaleWW*normFact;
  float scaleFact200 =  scaleH200*normFact;
  float scaleFact250 =  scaleH250*normFact;
  float scaleFact300 =  scaleH300*normFact;
  float scaleFact350 =  scaleH350*normFact;
  float scaleFact400 =  scaleH400*normFact;
  float scaleFact450 =  scaleH450*normFact;
  float scaleFact500 =  scaleH500*normFact;
  
  for(int i=0; i<thisSel.size();++i){
    selType this_selection = thisSel[i];

    // INPUT DIR for the histo files
    string input_dir_base = "NoNorm_Mu/";
    if(this_selection == elChannel)
      input_dir_base = "NoNorm_El/";
    
    // OUTPUT DIR for the histo files
    string out_dir_base = "histos_Mu/";
    if(this_selection == elChannel)
      out_dir_base = "histos_El/";

    // OUTPUT .png file
    string name = var+"_Lin_Mu.";
    if(set_logScale) name = var+"_Log_Mu.";
    if(this_selection == elChannel){
      name = var+"_Lin_El.";
      if(set_logScale) name= var+"_Log_El.";    
    }

    if(drawSigStack)
      name = masspoint+"_"+name;
      

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

    //    cout<<"get histos from file"<<endl;
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

  TH1F * var_All500 = new TH1F(*var_h500_Raw);
  var_All500->Reset();

  infile_name = input_dir_base+"TT"+postfix;
  TFile f_tt(infile_name.c_str());
  TH1F * var_TT_Raw = (TH1F*) f_tt.Get(var.c_str());

  TH1F * var_TT = new TH1F(*var_TT_Raw);
  var_TT->Scale(scaleFactTT);
  //cout << "TT raw entries " << var_TT_Raw->Integral() << endl;
  //cout << "---> TT 1 fb entries " << var_TT->Integral() << endl;


  infile_name = input_dir_base+"ZZ"+postfix;
  TFile f_zz(infile_name.c_str());
  TH1F * var_ZZ_Raw = (TH1F*) f_zz.Get(var.c_str());
  TH1F * var_ZZ = new TH1F(*var_ZZ_Raw);
  var_ZZ->Scale(scaleFactZZ);
  //cout << "ZZ raw entries " << var_ZZ_Raw->Integral() << endl;
  //cout << "---> ZZ 1 fb entries " << var_ZZ->Integral() << endl;

  infile_name = input_dir_base+"WZ"+postfix;
  TFile f_wz(infile_name.c_str());
  TH1F * var_WZ_Raw = (TH1F*) f_wz.Get(var.c_str());
  TH1F * var_WZ = new TH1F(*var_WZ_Raw);
  var_WZ->Scale(scaleFactWZ);
  //cout << "WZ raw entries " << var_WZ_Raw->Integral() << endl;
  //cout << "---> WZ 1 fb entries " << var_WZ->Integral() << endl;

  infile_name = input_dir_base+"WW"+postfix;
  TFile f_ww(infile_name.c_str());
  TH1F * var_WW_Raw = (TH1F*) f_ww.Get(var.c_str());
  TH1F * var_WW = new TH1F(*var_WW_Raw);
  var_WW->Scale(scaleFactWW);
  //cout << "WW raw entries " << var_WW_Raw->Integral() << endl;
  //cout << "---> WW 1 fb entries " << var_WW->Integral() << endl;


  infile_name = input_dir_base+"DYJets_Summ11"+postfix;
  TFile f_dyjets(infile_name.c_str());
  TH1F * var_DYJets_Raw = (TH1F*) f_dyjets.Get(var.c_str());
  TH1F * var_DYJets = new TH1F(*var_DYJets_Raw);
  var_DYJets->Scale(scaleFactDYJets);

  if(this_selection == muChannel){
    infile_name = input_dir_base+"MuRun2011"+runPeriod+postfix;
  }
  if(this_selection == elChannel){
    infile_name = input_dir_base+"ElRun2011"+runPeriod+postfix;
  }
  TFile f_dataMu(infile_name.c_str());
  TH1F * data_Mu = (TH1F*) f_dataMu.Get(var.c_str());
  //  TFile f_dataEl(infile_name.c_str());
  //  TH1F * data_El = (TH1F*) f_dataEl.Get(var.c_str());
  //  data_Mu->Add(data_El);
  //  cout<<"data"<<endl;


  // getting integrals in signal region
  double integralSignal = integralInRange(var_h500, cut_value1f, cut_value2f);
  double integralDY = integralInRange(var_DYJets, cut_value1f, cut_value2f);
  double integralTT = integralInRange(var_TT, cut_value1f, cut_value2f);
  double integralZZ = integralInRange(var_ZZ, cut_value1f, cut_value2f);
  double integralWZ = integralInRange(var_WZ, cut_value1f, cut_value2f);
  double integralWW = integralInRange(var_WW, cut_value1f, cut_value2f);
  double integralBkg = integralDY + integralTT + integralZZ + integralWZ + integralWW;
  double integralData = integralInRange(data_Mu, cut_value1f, cut_value2f);

  double allSignal = var_h500->Integral();
  double allDY = var_DYJets->Integral();
  double allTT = var_TT->Integral();
  double allZZ = var_ZZ->Integral();
  double allWZ = var_WZ->Integral();
  double allWW = var_WW->Integral();
  double allBkg = allDY + allTT + allZZ + allWZ + allWW;
  double allData = data_Mu->Integral();

  integralSignalCombined += integralSignal;
  integralDYCombined += integralDY;
  integralTTCombined += integralTT;
  integralZZCombined += integralZZ;
  integralWZCombined += integralWZ;
  integralWWCombined += integralWW;
  integralDataCombined += integralData;
  allSignalCombined += allSignal;
  allDYCombined += allDY;
  allTTCombined += allTT;
  allZZCombined += allZZ;
  allWZCombined += allWZ;
  allWWCombined += allWW;
  allDataCombined += allData;

  cout.precision(5);
  cout << "Signal region = [" << cut_value1f <<", " << cut_value2f << "]" << endl;
  cout << "sample" << "\t" << "s.r." << "\t" << "all" << endl;
  cout << "DY" << "\t" << integralDY << "\t" << allDY << endl;
  cout << "TT" << "\t" << integralTT << "\t" << allTT << endl;
  cout << "ZZ" << "\t" << integralZZ << "\t" << allZZ << endl;
  cout << "WZ" << "\t" << integralWZ << "\t" << allWZ << endl;
  cout << "WW" << "\t" << integralWW << "\t" << allWW << endl;
  cout << "AllBkg" << "\t" << integralBkg << "\t" << allBkg << endl;
  cout << masspoint << "\t" << integralSignal << "\t" << allSignal << endl;
  cout << "data" << "\t" << integralData << "\t" << allData << endl;

  if(this_selection == elChannel)  f<<"Electron Channel, "<< intLumi << "pb-1" << endl;
  if(this_selection == muChannel)  f<<"Muon Channel, "<< intLumi << "pb-1" <<endl;        
  f << "\t" << "s.r." << "\t" << "all" << endl;
  f.precision(5);
  f << "DY" << "\t" << integralDY << "\t" << allDY << endl;
  f << "TT" << "\t" << integralTT << "\t" << allTT << endl;
  f << "ZZ" << "\t" << integralZZ << "\t" << allZZ << endl;
  f << "WZ" << "\t" << integralWZ << "\t" << allWZ << endl;
  f << "WW" << "\t" << integralWW << "\t" << allWW << endl;
  f << "AllBkg" << "\t" << integralBkg << "\t" << allBkg << endl;
  f << masspoint << "\t" << integralSignal << "\t" << allSignal << endl;
  f << "data" << "\t" << integralData << "\t" << allData << endl;
  f << endl;
  
  // Eventually REBIN:
  if(rebin_fact > 1){    
    //var_All250->Rebin(rebin_fact);
    //var_All300->Rebin(rebin_fact);
    var_All500->Rebin(rebin_fact);
    var_WW->Rebin(rebin_fact); 
    var_WZ->Rebin(rebin_fact); 
    var_TT->Rebin(rebin_fact); 
    //    var_ZCC->Rebin(rebin_fact); 
    //    var_ZBB->Rebin(rebin_fact); 
    //    var_ZJET->Rebin(rebin_fact); 
    var_ZZ->Rebin(rebin_fact); 
    var_DYJets->Rebin(rebin_fact);
    //var_h250->Rebin(rebin_fact); 
    //var_h300->Rebin(rebin_fact); 
    var_h500->Rebin(rebin_fact); 
    data_Mu->Rebin(rebin_fact);
  }


  int n = var_h500->GetNbinsX();
  min = var_h500->GetXaxis()->GetXmin();
  max = var_h500->GetXaxis()->GetXmax();
  //  cout<<"min and max: "<<min<<" "<<max<<endl;
  //    cout<<"nBin: "<<n<<endl;
  float binsize = (max - min)/n;
    
  //  int bincut1 = (cut_value1 -min )/binsize +1;
  //  int bincut2 = (cut_value2 -min) /binsize +1;
  //  cout<<"bin cut "<<bincut1<<" --> "<<var_h500->GetBinCenter(bincut1)<<endl;
  //  cout<<"bin cut "<<bincut2<<" --> "<<var_h500->GetBinCenter(bincut2)<<endl;
    
  //  TH1F * S0 = new TH1F("S0","H mass",nbins,min,max);
  //  S0->SetMinimum(minimum);
  //  S0->SetMaximum(maximum);
  //  S0->SetTitle(var.c_str());

  // sum histos for stack plot
  TH1F * S1 = new TH1F(* var_WW);
  TH1F * S12 = new TH1F(* S1);
  S12->Add(var_WZ);
  TH1F * S123 = new TH1F(* S12);
  S123->Add(var_ZZ);
  TH1F * S1234 = new TH1F(* S123);
  S1234->Add(var_TT);
  /*
    TH1F * S12345 = new TH1F(* S1234);
    S12345->Add(var_ZBB);
    TH1F * S123456 = new TH1F(* S12345);
    S123456->Add(var_ZJET);
  */
  TH1F * S12345 = new TH1F(* S1234);
  S12345->Add(var_DYJets);
  TH1F * S123456 = new TH1F(* S12345);
  //  var_h500->Scale(10);
  S123456->Add(var_h500);


  // write data histo, stack histo and single histos into the output file

  TH1F *data = new TH1F("data","data", var_h500->GetNbinsX(), 100, 1000);
  data->SetName("data_obs");
  data->SetTitle("data_obs");

  
  if(this_selection == muChannel)  string outfile_name = out_dir_base+"Histo_"+masspoint+"_mu.root";
  if(this_selection == elChannel)  string outfile_name = out_dir_base+"Histo_"+masspoint+"_el.root";
  TFile outfile(outfile_name.c_str(), "RECREATE");
  //var_data->Write();

  data_Mu->SetMarkerStyle(22);
  /*
  data_Mu->SetMarkerStyle(22);
  data_Mu->Write();
  S123456->SetName("data_obs");
  S123456->Write();
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

  */

  // Draw stack plot into an .png file


  //  S1234567->SetStats(kFALSE);
  S123456->SetStats(kFALSE);
  S12345->SetStats(kFALSE);
  S1234->SetStats(kFALSE);
  S123->SetStats(kFALSE);
  S12->SetStats(kFALSE);
  S1->SetStats(kFALSE);


  //  S1234567->SetLineWidth(2);
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
  //  S12345->SetFillColor(kViolet+4);
 //    S12345->SetLineColor(kOrange+8);
  S12345->SetFillColor(kOrange-2);
  S123456->SetFillColor(kOrange+7);
  //  S123456->SetLineColor(kPink+3);


  //  TH1F * S1234567 = new TH1F(* S123456);
  //  S1234567->Add(h_svbf500);
  //  TH1F * S12345678 = new TH1F(* S1234567);
  //  S12345678->Add(h_svbf450);
 
  //  pad1 = new TPad("pad","This is pad1",0.02,0.02,0.48,0.83,33);
  //  pad1->Draw();

  leg = new TLegend(.68,.67,0.94,0.94);
  leg->SetFillColor(0);
  leg->SetTextSize(0.040);
  leg->AddEntry(S1,"WW","f");
  leg->AddEntry(S12,"WZ","f");
  leg->AddEntry(S123,"ZZ","f");
  leg->AddEntry(S1234,"t #bar{t}","f");
  //  leg->AddEntry(S12345,"Zbb","f");
  leg->AddEntry(S12345,"DY Jets","f");
  //  leg->AddEntry(S123456,("H"+masspoint).c_str(),"f");
  //leg->AddEntry(S12345678,"H250 ","f");
  //leg->AddEntry(S12345678,"H500 ","f");
  //leg->AddEntry( var_data_Sumch,"data ","P");
  if(drawSig==true)  leg->AddEntry(var_h500,("H"+masspoint+" x 20").c_str(),"l");
  if(drawSigStack==true) leg->AddEntry(S123456,("H"+masspoint+"").c_str(),"f");
  leg->AddEntry(data_Mu,"data ","P");

  float ymax_histosum = S123456->GetMaximum();
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

  //string Ylabel = "number of events/"+nEvts.str()+" GeV"

  var_All500->GetYaxis()->SetTitleOffset(1.5);
  var_All500->GetYaxis()->SetTitleSize(0.05);
  var_All500->GetYaxis()->SetLabelSize(0.045);
  var_All500->GetXaxis()->SetTitleSize(0.05);
  var_All500->GetXaxis()->SetLabelSize(0.040);

  var_All500->GetYaxis()->SetTitle(Form("number of events/ %.2f GeV",nEvts));
  //var_All500->GetYaxis()->SetTitle(Form("number of events/ %.2f ",nEvts));
  if(custom) var_All500->GetXaxis()->SetRangeUser(rMin,rMax);
  //var_data->SetMarkerStyle(21);

  var_h500->SetLineColor(kOrange+7);
  var_h500->SetLineWidth(3);
  //  var_h500->Scale(200);
  var_h500->Scale(20.);

  var_All500->Draw();
  //  S1234567->Draw("HISTsame");

  if(drawSigStack)  S123456->Draw("HISTsame");
  S12345->Draw("HISTsame");
  S1234->Draw("HISTsame");
  S123->Draw("HISTsame");
  S12->Draw("HISTsame");
  S1->Draw("HISTsame");
  if(drawSig==true)  var_h500->Draw("HISTsame");
  data_Mu->Draw("*same");


  if(set_logScale){
    c1->SetLogy();
    var_All500->SetMinimum(1.);
  }
  leg->Draw("SAME");
  // t->Draw();


  gPad->RedrawAxis();

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
  latex.DrawLatex(0.17, 0.96, "CMS preliminary 2011");
  latex.DrawLatex(0.57,0.96,Form("%.0f pb^{-1} at #sqrt{s} = 7 TeV",intLumi));

  c1->SaveAs(("plots/"+name+"eps").c_str());
  c1->SaveAs(("plots/"+name+"png").c_str());
  
  }

  f<<"Total yields, "<< intLumi << "pb-1" << endl;
  f << "\t" << "s.r." << "\t" << "all" << endl;
  f << "DY" << "\t" << integralDYCombined << "\t" << allDYCombined << endl;
  f << "TT" << "\t" << integralTTCombined << "\t" << allTTCombined << endl;
  f << "ZZ" << "\t" << integralZZCombined << "\t" << allZZCombined << endl;
  f << "WZ" << "\t" << integralWZCombined << "\t" << allWZCombined << endl;
  f << "WW" << "\t" << integralWWCombined << "\t" << allWWCombined << endl;
  f << "AllBkg" << "\t" << integralDYCombined+integralTTCombined+integralZZCombined+integralWZCombined+integralWWCombined; 
  f  << "\t" << allDYCombined+allTTCombined+allZZCombined+allWZCombined+allWWCombined << endl;
  f << masspoint << "\t" << integralSignalCombined << "\t" << allSignalCombined << endl;
  f << "data" << "\t" << integralDataCombined << "\t" << allDataCombined << endl;

  f.close();
  gApplication->Terminate();
  
}
