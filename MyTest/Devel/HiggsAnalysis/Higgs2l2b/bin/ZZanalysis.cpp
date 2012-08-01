
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <Math/VectorUtil.h>
#include "DataFormats/Common/interface/Wrapper.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "HiggsAnalysis/Higgs2l2b/interface/BTagSFUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


using namespace std;

typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;

#define BRANCH(name) \
  edm::Wrapper<vfloat> * name = new edm::Wrapper<vfloat>(); \
  string alias_##name = events->GetAlias(#name); \
  alias_##name.erase(alias_##name.size() - 3);	\
  TBranch * b_##name = events->GetBranch(alias_##name.c_str());		\
  if (b_##name == 0) { cerr << "missing branch: " << alias_##name << endl; } \
  b_##name->SetAddress(& name) 

#define BRANCHVUINT(name) \
  edm::Wrapper< vector <unsigned int> > * name = new edm::Wrapper< vector <unsigned int> >(); \
  string alias_##name = events->GetAlias(#name); \
  alias_##name.erase(alias_##name.size() - 3);	\
  TBranch * b_##name = events->GetBranch(alias_##name.c_str());		\
  if (b_##name == 0) { cerr << "missing branch: " << alias_##name << endl; } \
  b_##name->SetAddress(& name) 

#define BRANCHINT(name) \
  edm::Wrapper<int> * name = new edm::Wrapper<int>(); \
  string alias_##name = events->GetAlias(#name); \
  alias_##name.erase(alias_##name.size() - 3);	\
  TBranch * b_##name = events->GetBranch(alias_##name.c_str());		\
  if (b_##name == 0) { cerr << "missing branch: " << alias_##name << endl; } \
  b_##name->SetAddress(& name) 

#define BRANCHFLOAT(name) \
  edm::Wrapper<float> * name = new edm::Wrapper<float>(); \
  string alias_##name = events->GetAlias(#name); \
  alias_##name.erase(alias_##name.size() - 3);	\
  TBranch * b_##name = events->GetBranch(alias_##name.c_str());		\
  if (b_##name == 0) { cerr << "missing branch: " << alias_##name << endl; } \
  b_##name->SetAddress(& name) 

#define BRANCHBOOL(name) \
  edm::Wrapper<bool> * name = new edm::Wrapper<bool>(); \
  string alias_##name = events->GetAlias(#name); \
  alias_##name.erase(alias_##name.size() - 3);	\
  TBranch * b_##name = events->GetBranch(alias_##name.c_str());		\
  if (b_##name == 0) { cerr << "missing branch: " << alias_##name << endl; } \
  b_##name->SetAddress(& name) 

template<typename T>
const T & get(const edm::Wrapper<std::vector<T> > * name, unsigned int idx) { return (* name->product())[idx]; }

template<typename T>
const T & getInt(const edm::Wrapper<T> * name) { return (* name->product()); }

#define GETENTRY(name, idx) \
  b_##name->GetEntry(idx); \
  assert(name != 0); \
  assert(name->isPresent())

void progress(float percent) {
  const int size = 50;
  std::cout<<"\r[";
  for(int i = 0; i < size; ++i) {
    if(i <= percent * size) cout << '=';
    else cout << ".";
  }
  std::cout <<"] ";
  std::cout.width(4);
  std::cout<<int(percent*100)+1<<"%"<<std::flush;
}

bool EleIDTightCuts( bool isEleBarrel, float DPhiAtVtx, float HOverE) {
  if( isEleBarrel ) {
    // cuts for barrel electron
    if( fabs(DPhiAtVtx) > 0.15 )
      return false;
  } else {
    // cuts for endcap electrons
    if( fabs(DPhiAtVtx) > 0.1 ) 
      return false;
    if( HOverE > 0.1 ) 
      return false;
  }
  return true;
}
 
// macro usage:
// 2l2b [input_ntuple_file] [output_histo_file] [DATA/MC] [SF/noSF] [Mu/El] [tree/notree] [txt/notxt]
// if MC --> apply PUreweight
// if SF --> apply b-tag scale factor
// Mu/El --> selection channel
// tree  --> write tree output
// txt  --> write txt output

int main(int argc, char **argv) {
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  const char * path = argv[1];
  const char * outFile = argv[2];
  string isData(argv[3]);
  string applybTagSF(argv[4]);
  string selChan(argv[5]);
  string writeTree(argv[6]);
  string writeTxt(argv[7]);
  string PUPeriod(argv[8]);
  string fixMuRunB(argv[9]);
  string outFileName(outFile);
  bool data = false;
  bool btagScale = false;
  bool fixMu = false;
  cout<<"filePath: "<<path<<endl;
  if (isData=="DATA") {
    data = true; 
    cout<<"Data is set to: "<<data<<endl;
    cout<<"No PUreweight applied!!!"<<endl;
  }

  if(applybTagSF == "SF") {
    btagScale = true;
    cout<<"SF will be applied!!!"<<endl;
  }

  if(fixMuRunB == "fixMu"){
    fixMu=true;
    cout << "Mu FIX for runB applied" << endl;
  }

  bool muChannel(false);
  if(selChan == "Mu")
    muChannel = true;
  
  bool writeT(false);
  if(writeTree == "tree")
    writeT = true;

  bool writeText(false);
  if(writeTxt == "txt")
    writeText = true;

  //  bool signalSel(false);
  //    if(selType == "sig")
  //      signalSel = true;

  vfloat bestHiggsMass, bestHiggsMass_1btag, bestHiggsMass_0btag;
  
  TFile *f = TFile::Open(path, "READ");
  TTree * events = (TTree*) f->Get("Events");
  unsigned int nEvents = events->GetEntries();
  cout << "events: " << nEvents << endl;

  // GET BRANCHES
  #include "HiggsAnalysis/Higgs2l2b/bin/2l2b_branches.h"

  // txt output file
  string reportName = "SelectedEvents"+selChan+".txt";
  ofstream fileout(reportName.c_str(),ios::app);

  // define histograms
  // histo for event variables

  TH1F cutFlow_0btag("cutFlow_0btag", "cutFlow_0btag", 12, 0., 12.);
  TH1F cutFlow_1btag("cutFlow_1btag", "cutFlow_1btag", 12, 0., 12.);
  TH1F cutFlow_2btag("cutFlow_2btag", "cutFlow_2btag", 12, 0., 12.);



  TH1F npv_woReweight("npv_woReweight", "npvwoReweight", 31, -0.5, 30.5);
  TH1F npv("npv", "npv", 31, -0.5, 30.5);
  TH1F npv1("npv1", "npv", 31, -0.5, 30.5);
  TH1F npv2("npv2", "npv", 31, -0.5, 30.5);
  TH1F npvfin("npvfin", "npv", 31, -0.5, 30.5);
  //  TH1F ngenint("ngenint", "ngenint", 36, -0.5, 35.5);
  // histo before lepton cuts
  TH1F db("db", "db", 100, 0, 0.2);
  // histo before Z mass cuts
  TH1F leptpt1("leptpt1", "P_{T} ", 300, 0, 300);
  TH1F leptpt2("leptpt2", "P_{T} ", 300, 0, 300);
  TH1F jetpt1("jetpt1", "P_{T} ", 300, 0, 300);
  TH1F jetpt2("jetpt2", "P_{T} ", 300, 0, 300);
  TH1F zjjmass("zjjmass", "m_{jj}", 500, 0, 500);
  TH1F zllmass("zllmass", "m_{ll}", 500, 0, 500);

  TH1F jetpt1_0btag("jetpt1_0btag", "P_{T} ", 300, 0, 300);
  TH1F jetpt1_2btag("jetpt2_2btag", "P_{T} ", 300, 0, 300);

  TH1F bb_0btag("bb_0btag", "bb shape m_{jj}", 500, 0, 500);
  TH1F notbb_0btag("notbb_0btag", "not bs shape m_{jj}", 500, 0, 500);
  TH1F bb_2btag("bb_2btag", "bb shape m_{jj}", 500, 0, 500);
  TH1F notbb_2btag("notbb_2btag", "not bs shape m_{jj}", 500, 0, 500);

  TH1F gg_2btag("gg_2btag", "gg shape m_{jj}", 500, 0, 500);

  TH1F stdCandle_0btag("stdCandle_0btag", "m_{ll}", 500, 0, 500);
  TH1F stdCandle_1btag("stdCandle_1btag", "m_{ll}", 500, 0, 500);
  TH1F stdCandle_2btag("stdCandle_2btag", "m_{ll}", 500, 0, 500);

  TH1F stdCandle_final_0("stdCandle_final_0", "m_{ll}", 500, 0, 500);
  TH1F stdCandle_final_1("stdCandle_final_1", "m_{ll}", 500, 0, 500);
  TH1F stdCandle_final_2("stdCandle_final_2", "m_{ll}", 500, 0, 500);

  TH1F zjjmass_0btag("zjjmass_0btag", "m_{jj}", 500, 0, 500);
  TH1F zjjmass_1btag("zjjmass_1btag", "m_{jj}", 500, 0, 500);
  TH1F zjjmass_2btag("zjjmass_2btag", "m_{jj}", 500, 0, 500);

  TH1F zjjmass_sbTT_0btag("zjjmass_sbTT_0btag", "m_{jj}", 500, 0, 500);
  TH1F zjjmass_sbTT_1btag("zjjmass_sbTT_1btag", "m_{jj}", 500, 0, 500);
  TH1F zjjmass_sbTT_2btag("zjjmass_sbTT_2btag", "m_{jj}", 500, 0, 500);

  TH1F jjm_0btag("jjm_0btag", "m_{jj}", 100, 50, 150);
  TH1F jjm_1btag("jjm_1btag", "m_{jj}", 100, 50, 150);
  TH1F jjm_2btag("jjm_2btag", "m_{jj}", 100, 50, 150);

  TH1F lljjmass_0btag("lljjmass_0btag", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_1btag("lljjmass_1btag", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_2btag("lljjmass_2btag", "m(lljj)", 1000, 0, 1000);

  TH1F llpt_2btag("llpt_2btag", "P_{T} Z_{ll}", 300, 0, 300);
  TH1F jjpt_2btag("jjpt_2btag", "P_{T} Z_{jj}", 300, 0, 300);
  TH1F jjdphi_2btag("jjdphi_2btag", "#Delta #Phi_{jj}", 100, -3., 3.);

  TH1F met_0btag("met_0btag", "MET", 100, 0, 140);
  TH1F met_1btag("met_1btag", "MET", 100, 0, 140);
  TH1F met_2btag("met_2btag", "MET", 100, 0, 140);

  TH1F t1corrMET_0btag("t1corrMET_0btag", "t1corrMET", 100, 0, 140);
  TH1F t1corrMET_1btag("t1corrMET_1btag", "t1corrMET", 100, 0, 140);
  TH1F t1corrMET_2btag("t1corrMET_2btag", "t1corrMET", 100, 0, 140);

  TH1F trkMET_0btag("trkMET_0btag", "trkMET", 100, 0, 140);
  TH1F trkMET_1btag("trkMET_1btag", "trkMET", 100, 0, 140);
  TH1F trkMET_2btag("trkMET_2btag", "trkMET", 100, 0, 140);

  TH2F trkMET_MET("trkMET_MET", "trkMET vs MET", 100, 0, 140, 100, 0, 140);

  TH1F minMET_2btag("minMET_2btag", "minMET", 100, 0, 140);
  TH1F minCorrMET_2btag("minCorrMET_2btag", "minCorrMET", 100, 0, 140);


  TH1F jzb_0btag("jzb_0btag", "JZB",800, -800., 800.);
  TH1F jzb_1btag("jzb_1btag", "JZB",800, -800., 800.);
  TH1F jzb_2btag("jzb_2btag", "JZB",800, -800., 800.);

  TH2F jzbmjj("jzbmjj", "JZB vs m_{jj}",500, 0., 500., 800, -800., 800. );
  TH2F zzdptmjj("zzdptmjj", "#Delta p_{T}_{ZZ} vs m_{jj}",500, 0., 500., 300, -300., 800. );
  TH2F jzbmjjcorr("jzbmjjcorr", "(JZB-m_{jj}) vs m_{jj}",500, 0., 500., 800, -800., 800. );

  TH2F jptmjj("jptmjj", "Jet p_{T} vs m_{jj}",500, 0., 500., 300, 0., 300. );

  TH1F drJ1L_0btag("drJ1L_0btag", "DR Leading Jet - Lepton", 100, 0, 6);
  TH1F drJ1L_1btag("drJ1L_1btag", "DR Leading Jet - Lepton", 100, 0, 6);
  TH1F drJ1L_2btag("drJ1L_2btag", "DR Leading Jet - Lepton", 100, 0, 6);

  TH1F drJ2L_0btag("drJ2L_0btag", "DR Leading Jet - Lepton", 100, 0, 6);
  TH1F drJ2L_1btag("drJ2L_1btag", "DR Leading Jet - Lepton", 100, 0, 6);
  TH1F drJ2L_2btag("drJ2L_2btag", "DR Leading Jet - Lepton", 100, 0, 6);


  TH2F drJ1L_mjj("drJ1Lmjj", "#Delta R{jl} (leading jet) vs m_{jj}", 500, 0., 500., 100, 0., 6. );
  TH2F drJ2L_mjj("drJ2Lmjj", "#Delta R{jl} (non leading jet) vs m_{jj}",  500, 0., 500., 100, 0., 6. );

  TH2F mjjmll("mjjmll", "m_{jj} vs m_{ll}", 500, 0, 500, 500, 0, 500);
  TH2F ptjjmjj("ptjjmjj", "P_{T} vs m_{jj}",  500., 0., 500., 300., 0., 300.);
  TH2F drllmjj("drllmjj", "#Delta R{ll} vs m_{jj}",500., 0., 500., 100, 0, 4);
  TH2F csvmjj("csvmjj", "CSV vs m_{jj}", 500, 0, 500, 100, 0, 1);
  TH2F drzzmjj("drzzmjj", "#Delta R_{ZZ} vs m_{jj}", 500, 0, 500, 100, 0, 4 );
  TH2F ptllptjj("ptllptjj", "P_{T} Z_{ll} vs P_{T} Z_{jj}", 300, 0, 300, 300, 0, 300 ); 
  TH2F dphijjmjj("dphijjmjj", "#Delta #Phi_{jj} vs m_{jj}", 300, 0, 300, 100, 0, 3. );
  TH2F detajjmjj("detajjmjj", "#Delta #Eta_{jj} vs m_{jj}", 300, 0, 300, 100, -3, 3. );



  TH1F zzdr("zzdr", "#Delta R_{ZZ}", 100, 0, 4);
  TH1F zzdphi("zzdphi", "#Delta #Phi_{ZZ}", 100, 0, 6);
  TH1F zzdeta("zzdeta", "#Delta #Eta_{ZZ}", 100, -3, 3);
  TH1F zzdr_btag("zzdr_btag", "#Delta R_{ZZ}", 100, 0, 4);
  TH1F lldr("lldr", "#Delta R_{ll}", 100, 0, 4);

  TH1F lljjpt("lljjpt", "P_{T} ", 300, 0, 300);
  TH1F lljjY("lljjY", "Y", 100, -4, 4);




  // histo before b-tagging cut
  TH1F btagJet1("TCHEJet1", "TCHE", 100, 0, 20);
  TH1F btagJet2("TCHEJet2", "TCHE", 100, 0, 20);
  TH1F btagJet("TCHEJet", "TCHE", 100, 0, 20);
  TH1F btagJet1_CSV("CSVJet1", "CSV", 100, 0, 1);
  TH1F btagJet2_CSV("CSVJet2", "CSV", 100, 0, 1);
  TH1F btagJet_CSV("CSVJet", "CSV", 100, 0, 1);
  //TH2F btagJet_TCHECSV("TCHE_VS_CSV", "TCHECSV", 100, 0, 1 ,100, 0, 20);
  //TH2F btagJet1_TCHECSV("TCHE_VS_CSV_Jet1", "TCHECSV", 100, 0, 1 ,100, 0, 20);
  //TH2F btagJet2_TCHECSV("TCHE_VS_CSV_Jet2", "TCHECSV", 100, 0, 1 ,100, 0, 20);



  TH1F njets_0btag("njets_0btag", "njets", 11, -0.5, 10.5);
  TH1F njets_2btag("njets_2btag", "njets", 11, -0.5, 10.5);


  TH1F h_met("met", "MET", 100, 0, 140);
  TH1F metSig("metSig", "METSig", 100, 0, 15);
  // histo for 2b candidates after b-tagging
  TH1F zllpt("zllpt", "P_{T}(Z_{ll})", 300, 0, 300);
  TH1F zjjpt("zjjpt", "P_{T}(Z_{jj})", 300, 0, 300);
  TH1F zzdpt("zzdpt", "P_{T}(Z_{ll})- P_{T}(Z_{jj})", 300,-300, 300);



  TH1F drjj("DRjj", "#Delta R_{jj}", 100, 0, 4);
  TH1F jjdphi("jjdphi", "#Delta #Phi_{jj}", 100, -3., 3.);
  TH1F jjdeta("jjdeta", "#Delta #Eta_{jj}", 100, -3., 3.);
  TH1F drjj_btag("DRjj_btag", "#Delta R_{jj}", 100, 0, 4);
  TH1F metSignif("metSignif", "METSignif", 100, 0, 15);
  TH1F lljjmass2tags("lljjmass2tags", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_btag("lljjmass_btag", "m(lljj)", 1000, 0, 1000);
  TH1F helyLD_RF2tags("HelyLDRefit2tags", "HelyLDRefit", 100, 0, 1);
  // histo after met significance cut
  TH1F zllptBeforeLD("zllptBefLD", "P_{T}(Z_{ll})", 300, 0, 300);
  TH1F drjjBefLD("DRjjBefLD", "#Delta R_{jj}", 100, 0, 4);
  TH1F lljjmassBeforeLD("lljjmassBefLD", "m(lljj)", 1000, 0, 1000);

  TH1F helyLD("HelyLD", "HelyLD", 100, 0, 1);
  TH1F cosT1("cosTheta1", "cosTheta1", 100, -1, 1);
  TH1F cosT2("cosTheta2", "cosTheta2", 100, 0, 1);
  TH1F cosT1Star("cosTheta1Star", "cosTheta1Star", 100, -1, 1);
  TH1F phi("phi", "phi", 100, -3, 3);
  TH1F phiStar("phiStar", "phiStar", 100, -3, 3);
  
  TH1F helyLD_RF("HelyLDRefit", "HelyLDRefit", 100, 0, 1);
  TH1F cosT1_RF("cosTheta1Refit", "cosTheta1Refit", 100, -1, 1);
  TH1F cosT2_RF("cosTheta2Refit", "cosTheta2Refit", 100, 0, 1);
  TH1F cosT1Star_RF("cosTheta1StarRefit", "cosTheta1StarRefit", 100, -1, 1);
  TH1F phi_RF("phiRefit", "phiRefit", 100, -3, 3);
  TH1F phiStar_RF("phiStarRefit", "phiStarRefit", 100, -3, 3);
  
  // final inv mass histo
  TH1F lljjmass("lljjmass", "m(lljj)", 1000, 0, 1000);
  //  TH1F lljjmass_0btag("lljjmass_0btag", "m(lljj)", 1000, 0, 1000);
  //  TH1F lljjmass_1btag("lljjmass_1btag", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_noRefit("lljjmassNoRefit", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_noRefit_1btag("lljjmassNoRefit_1btag", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_noRefit_0btag("lljjmassNoRefit_0btag", "m(lljj)", 1000, 0, 1000);
  //  TH1F lljjmassCands("lljjmassCands", "m(lljj)", 1000, 0, 1000);

  // inv mass histo for sb selection
  TH1F hmass_sb_all("lljjmass_sb_all", "m(lljj)", 1000, 0, 1000);
  TH1F hmassnorefit_sb_all("lljjmassNoRefit_sb_all", "m(lljj)", 1000, 0, 1000);
  TH1F hmass_sb_2b("lljjmass_sb_2b", "m(lljj)", 1000, 0, 1000);
  TH1F hmassnorefit_sb_2b("lljjmassNoRefit_sb_2b", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_sb("lljjmass_sb", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_noRefit_sb("lljjmassNoRefit_sb", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmassFinal_sb("lljjmassFinal_sb", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmassFinal_noRefit_sb("lljjmassFinalNoRefit_sb", "m(lljj)", 1000, 0, 1000);


  const unsigned int counts = 12;
  vint count(counts, 0);
  vint count_1(counts, 0);
  vint count_2(counts, 0);
  vint count_test(counts, 0);

  vint cand(counts, 0);


  //b-tag code initialization 
  int CandBTag[3];
  for(int kkk = 0; kkk<3 ; kkk++)
    CandBTag[kkk] = 0;
  // BTag SF utility
  BTagSFUtil* btsfutil = new BTagSFUtil(13);
  //TCHE SF
  TFile* file_loose = TFile::Open("BTagPayloads_TCHEL.root");
  TFile* file_medium = TFile::Open("BTagPayloads_TCHEM.root");

   //CSV SF
  //TFile* file_loose = TFile::Open("BTagPayloads_CSVL.root");
  //TFile* file_medium = TFile::Open("BTagPayloads_CSVM.root");
  //
  btsfutil->set_fileLoose(file_loose);
  btsfutil->set_fileMedium(file_medium);

  cand[1]=0;

  // PU reweighting util
  //  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting LumiWeights_;
  if(!data){
    //  --> Run2011A <--
    if(PUPeriod == "A"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011A" << endl;
      LumiWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011A_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    //  --> Run2011B <--
    else if(PUPeriod == "B"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011B" << endl;
      LumiWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011B_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    //  --> Run2011All <--
    else if(PUPeriod == "All"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011" << endl;
      LumiWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011All_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    else
      cout << "ERROR!!!!! Run period not existent!!!!!!" << endl;

    TString stringPUWei("Weight3D.root");
    const char * filePUwei = gSystem->FindFile("./", stringPUWei);
    if( filePUwei != NULL )
      LumiWeights_.weight3D_init("Weight3D.root");
    else
      LumiWeights_.weight3D_init(1.0);
  }

  int nEvts_0 = 0;
  int nEvts_1 = 0;
  int nEvts_2 = 0;

  int nCands_0 = 0;
  int nCands_1 = 0;
  int nCands_2 = 0;
  // Loop on EVENTS


  count[0] = nEvents;  
  count_1[0] = nEvents;  
  count_2[0] = nEvents;
  count_test[0] =nEvents;
  


  for(unsigned int i = 0; i <nEvents; ++i) {
    vbool pass(counts, false);
    vbool pass_1b(counts, false);
    vbool pass_2b(counts, false);
    //    if(i%100 == 0) progress(double(i)/double(nEvents));
    //    ++count[0];


    // variables for HLT selection
    bool passTrigPath(true);
    // NUMBER OF CANDIDATES IN THE EVENT
    unsigned int cands(0);
    unsigned int njets(1000);

    // GET ENTRIES FOR EVENT VARIABLES
    GETENTRY(numPV,i);
    GETENTRY(met,i);
    GETENTRY(t1corrMet,i);
    GETENTRY(metPhi,i);
    GETENTRY(metSignificance,i);
    GETENTRY(passSingleMuTrig,i);
    GETENTRY(passDoubleMuTrig,i);
    GETENTRY(passSingleElTrig,i);
    GETENTRY(passDoubleElTrig,i);
    // assign event weight for PU reweighting
    double evtWeight = 1.0;
    if(!data){
      //      BRANCHFLOAT(nGenInt);
      //      GETENTRY(nGenInt,i);
      BRANCHINT(nGenIntBXm1);
      BRANCHINT(nGenIntBX0);
      BRANCHINT(nGenIntBXp1);
      GETENTRY(nGenIntBXm1,i);
      GETENTRY(nGenIntBX0,i);
      GETENTRY(nGenIntBXp1,i);

      evtWeight = LumiWeights_.weight3D( getInt(nGenIntBXm1), getInt(nGenIntBX0), getInt(nGenIntBXp1) );
    }      



    
    if(muChannel) {
      if(data)
	passTrigPath = getInt(passDoubleMuTrig);
      //     if(fixMu && data)
      //	passTrigPath = true;
      //      passSingleTrig = getInt(passSingleMuTrig);
      //      passDoubleTrig = getInt(passDoubleMuTrig) && !(getInt(passDoubleElTrig));
      GETENTRY(muHiggsLeptDau1Pt,i);

      cands = muHiggsLeptDau1Pt->product()->size();

    } else {
      if(data)
      	passTrigPath = !(getInt(passDoubleMuTrig)) && getInt(passDoubleElTrig);
      //      passSingleTrig = getInt(passSingleElTrig);
      //      passDoubleTrig = getInt(passDoubleElTrig);
      GETENTRY(elHiggsLeptDau1Pt,i);
      cands = elHiggsLeptDau1Pt->product()->size();
    }
    
    // HLT selection: exclusive double trig selection
    if(!passTrigPath)
      continue;

    // fill histo of numPV (reweighted and not)
    npv.Fill(getInt(numPV),evtWeight);
    npv_woReweight.Fill(getInt(numPV));
    
    fill(pass.begin(), pass.end(), false);
    GETENTRY(CleanJetPt,i);
    njets = CleanJetPt->product()->size();
    if(cands > 0) {
      pass[1] = true;
      pass_1b[1] = true;
      pass_2b[1] = true;


      //GET ENTRIES FOR MU/EL VARIABLES
      if(muChannel){
#include "HiggsAnalysis/Higgs2l2b/bin/2l2bMu_getent.h"
      } else {
#include "HiggsAnalysis/Higgs2l2b/bin/2l2bEl_getent.h"
      }

      //#include "HiggsAnalysis/Higgs2l2b/bin/2l2b_getent.h"
      
      // variables to store the best H candidate in the event
      //      float met_, metSig_;
      float lljjmass_(0), lljjmass_noRefit_(0), lljjmass_1btag_(0), lljjmass_noRefit_1btag_(0);
      float lljjpt_(0), lljjY_(-10);
      float lljjmass_0btag_(0), lljjmass_noRefit_0btag_(0);
      float SumPtBest_(0), SumPtBest_Final(0), SumPtBest_Final_1btag(0);
      bool exist2tag(false), exist1tag(false), exist0tag(false);
      float ResidZllMassBest_(100.), ResidZllMassBest_Final(100.);
      float ResidZllMassBest1tag_(100.);
      float ResidZllMassBest0tag_(100.);
      float hLD_Best(-10), cos1_Best(-10), cos2_Best(-10), cosStar_Best(-10);
      float phi_Best(-10), phiStar_Best(-10);
      float hLD_rfBest(-10), cos1_rfBest(-10), cos2_rfBest(-10), cosStar_rfBest(-10);
      float phi_rfBest(-10), phiStar_rfBest(-10);
      float zzdphi_(-10), zzdeta_(-10), jjdphi_(-10), jjdeta_(-10);
      int EvtNum_(0), RunNum_(0);
      float metPhi_(-1);
      float minMET_(-1);
      float minCorrMET_(-1);
      float trkMet_(-1);
      float notbb_0btagBest(-1), bb_0btagBest(-1), notbb_2btagBest(-1), bb_2btagBest(-1);
      float SumPtBest_SB(0), SumPtBest_FinalSB(0), SumPtBest_FinalSB_1btag(0);
      float zjjmassBest_0(10000),zjjmassBest_1(10000), zjjmassBest_2(10000);
      float zjjmassBestSB_0(10000),zjjmassBestSB_1(10000), zjjmassBestSB_2(10000);

      float jetpt1Best_0(-1), jetpt1Best_1(-1), jetpt1Best_2(-1),trkMETBest_2(-1);  
      float lljjmassBest_2(-1), llptBest_2(-1000), jjptBest_2(-1000), jjdphiBest_2(-1000); 
      float trkMetBest_2(-1);

      float lljjmass_SB(0), lljjmass_noRefit_SB(0);
      float lljjmassFinal_SB(0), lljjmassFinal_noRefit_SB(0), lljjmassFinal_SB_1btag_(0), lljjmassFinal_noRefit_SB_1btag_(0);
      float zllcharge_(-10);

      // auxiliary variables for channel-independent selection
      float Jet1pt_(0), Jet2pt_(0), Jet1pt_rf(0), Jet2pt_rf(0);
      float Jet1eta_(0), Jet2eta_(0), Jet1phi_(0), Jet2phi_(0);
      float Lept1eta_(0), Lept2eta_(0), Lept1phi_(0), Lept2phi_(0);
      float lept1pt_(0), lept2pt_(0);
      float zllmass_(0), zjjmass_(0);
      float jet1tkhe_(0), jet2tkhe_(0);
      float jet1csv_(0), jet2csv_(0);
      int jet1flav_(0), jet2flav_(0);
      float zllpt_(-1.), zjjpt_(-1.), drjj_(-1.), zzdr_(-1), lldr_(-1) ;
      float hLD_(-10), cos1_(-10), cos2_(-10), cosStar_(-10), phi_(-10), phiStar_(-10);
      float hLD_rf(-10), cos1_rf(-10), cos2_rf(-10), cosStar_rf(-10), phi_rf(-10), phiStar_rf(-10);
      float higgsrefitMass(0), higgsMass(0);
      int higgsEvtNum(0), higgsRunNum(0);
      float zllphi_(-10);

      int flavourJet1(-1), flavourJet2(-1);


      // muon kine / isolation cut / cosmic rejection
      bool kineLepCut(false), dbLepCut(false), isoIDLepCut(false);

      // flag to dump event variables
      bool isEventSelected(false);

      bool hasCand_0(false);
      bool hasCand_1(false);
      bool hasCand_2(false);
      // LOOP on CANDIDATES IN THE EVENT


      int njets_0(1000), njets_2(1000);
     
      //  if (njets<4){ 
	pass[2] = true;
	pass_1b[2] = true;
	pass_2b[2] = true;

	for(unsigned int j = 0; j < cands; ++j) {
	++cand[1];

	// fill histo of db of the muon candidates
	if(muChannel) {
	  db.Fill(get(muHiggsLeptDau1dB,j),evtWeight);
	  db.Fill(get(muHiggsLeptDau2dB,j),evtWeight);
	}

	// assign variables for selection
	if(muChannel) {
	  zllcharge_= get(muHiggszllCharge,j);
	  Jet1pt_=get(muHiggsJetDau1Pt,j);
	  Jet2pt_=get(muHiggsJetDau2Pt,j);
	  Jet1pt_rf=get(muHiggsJetDau1RefitPt,j);
	  Jet2pt_rf=get(muHiggsJetDau2RefitPt,j);
	  Jet1eta_=get(muHiggsJetDau1Eta,j);
	  Jet2eta_=get(muHiggsJetDau2Eta,j);
	  Jet1phi_=get(muHiggsJetDau1Phi,j);
	  Jet2phi_=get(muHiggsJetDau2Phi,j);
	  Lept1eta_=get(muHiggsLeptDau1Eta,j);
	  Lept2eta_=get(muHiggsLeptDau2Eta,j);
	  Lept1phi_=get(muHiggsLeptDau1Phi,j);
	  Lept2phi_=get(muHiggsLeptDau2Phi,j);
	  lept1pt_=get(muHiggsLeptDau1Pt,j);
	  lept2pt_=get(muHiggsLeptDau2Pt,j);
	  dbLepCut = (fabs(get(muHiggsLeptDau1dB,j)) <0.02) && (fabs(get(muHiggsLeptDau2dB,j)) <0.02);
	  isoIDLepCut = (get(muHiggsLeptDau1CombRelIso,j) < 0.15) && (get(muHiggsLeptDau2CombRelIso,j) < 0.15);
	  zllmass_ = get(muHiggszllMass,j);
	  zjjmass_ = get(muHiggszjjMass,j); 
	  jet1tkhe_ = get(muHiggsJet1TKHE,j);
	  jet2tkhe_ = get(muHiggsJet2TKHE,j);
	  jet1csv_ = get(muHiggsJet1CSV,j);
	  jet2csv_ = get(muHiggsJet2CSV,j);
	  jet1flav_ = get(muHiggsJetDau1PartonFlavour,j);
	  jet2flav_ = get(muHiggsJetDau2PartonFlavour,j);
	  zllpt_ = get(muHiggszllPt,j);
	  zllphi_ = get(muHiggszllPhi,j);
	  zjjpt_ = get(muHiggszjjPt,j);
	  drjj_ = get(muHiggsjjdr,j);
	  jjdphi_ = get(muHiggsjjdPhi,j);
	  jjdeta_ = get(muHiggsjjdEta,j);
	  lldr_ = get(muHiggslldr,j);
	  zzdr_ = get(muHiggszzdr,j);
	  zzdeta_ = get(muHiggszzdEta,j);
	  zzdphi_ = get(muHiggszzdPhi,j);
	  cos1_=get(muHiggscosthetaNT1,j);
	  cos2_=get(muHiggscosthetaNT2,j);
	  cosStar_= get(muHiggscosthetastarNT,j);
	  phi_= get(muHiggsphiNT,j);
	  phiStar_= get(muHiggsphiNT1,j);
	  hLD_= get(muHiggsHelyLD,j);
	  cos1_rf=get(muHiggscosthetaNT1Refit,j);
	  cos2_rf=get(muHiggscosthetaNT2Refit,j);
	  cosStar_rf= get(muHiggscosthetastarNTRefit,j);
	  phi_rf= get(muHiggsphiNTRefit,j);
	  phiStar_rf= get(muHiggsphiNT1Refit,j);
	  hLD_rf= get(muHiggsHelyLDRefit,j);
	  higgsrefitMass = get(muHiggsRefitMass,j);
	  higgsMass = get(muHiggsMass,j);
	  lljjpt_ = get(muHiggsPt,j);
	  lljjY_ = get(muHiggsY,j);
	  higgsEvtNum = getInt(muHiggsEventNumber);
	  higgsRunNum = getInt(muHiggsRunNumber);
	  trkMet_ = get(muHiggstrkMet, j);
	  flavourJet1 = get(muHiggsJetDau1PartonFlavour, j);
	  flavourJet2 = get(muHiggsJetDau2PartonFlavour, j);
	} else {
	  zllcharge_= get(elHiggszllCharge,j);
	  Jet1pt_=get(elHiggsJetDau1Pt,j);
	  Jet2pt_=get(elHiggsJetDau2Pt,j);
	  Jet1pt_rf=get(elHiggsJetDau1RefitPt,j);
	  Jet2pt_rf=get(elHiggsJetDau2RefitPt,j);
	  Jet1eta_=get(elHiggsJetDau1Eta,j);
	  Jet2eta_=get(elHiggsJetDau2Eta,j);
	  Lept1eta_=get(elHiggsLeptDau1Eta,j);
	  Lept2eta_=get(elHiggsLeptDau2Eta,j);
	  Lept1phi_=get(elHiggsLeptDau1Phi,j);
	  Lept2phi_=get(elHiggsLeptDau2Phi,j);
	  lept1pt_=get(elHiggsLeptDau1Pt,j);
	  lept2pt_=get(elHiggsLeptDau2Pt,j);
	  dbLepCut = true;

	  // new electronID
	  bool EleDau1VBTFID = get(elHiggsEleDau1VBTF95CombID,j) == 7;
	  bool isEle1Barrel = get(elHiggsLeptDau1isEB,j);
	  bool EleDau2VBTFID = get(elHiggsEleDau2VBTF95CombID,j) == 7;
	  bool isEle2Barrel = get(elHiggsLeptDau2isEB,j);

	  bool EleDau1TightCuts = EleIDTightCuts( isEle1Barrel, get(elHiggsLeptDau1DeltaPhiAtVtx,j), 
						  get(elHiggsLeptDau1HOverE,j) );
	  bool EleDau2TightCuts = EleIDTightCuts( isEle2Barrel, get(elHiggsLeptDau2DeltaPhiAtVtx,j), 
						  get(elHiggsLeptDau2HOverE,j) );

	  isoIDLepCut = EleDau1VBTFID && EleDau1TightCuts && EleDau2VBTFID && EleDau2TightCuts;

	  zllmass_ = get(elHiggszllMass,j);
	  zjjmass_ = get(elHiggszjjMass,j); 
	  jet1tkhe_ = get(elHiggsJet1TKHE,j);
	  jet2tkhe_ = get(elHiggsJet2TKHE,j);
	  jet1csv_ = get(elHiggsJet1CSV,j);
	  jet2csv_ = get(elHiggsJet2CSV,j);
	  jet1flav_ = get(elHiggsJetDau1PartonFlavour,j);
	  jet2flav_ = get(elHiggsJetDau2PartonFlavour,j);
	  zllpt_ = get(elHiggszllPt,j);
	  zllphi_ = get(elHiggszllPhi,j);
	  zjjpt_ = get(elHiggszjjPt,j);
	  drjj_ = get(elHiggsjjdr,j);
	  jjdphi_ = get(elHiggsjjdPhi,j);
	  jjdeta_ = get(elHiggsjjdEta,j);
	  lldr_ = get(elHiggslldr,j);
	  zzdr_ = get(elHiggszzdr,j);
	  zzdeta_ = get(elHiggszzdEta,j);
	  zzdphi_ = get(elHiggszzdPhi,j);
	  cos1_=get(elHiggscosthetaNT1,j);
	  cos2_=get(elHiggscosthetaNT2,j);
	  cosStar_= get(elHiggscosthetastarNT,j);
	  phi_= get(elHiggsphiNT,j);
	  phiStar_= get(elHiggsphiNT1,j);
	  hLD_= get(elHiggsHelyLD,j);
	  cos1_rf=get(elHiggscosthetaNT1Refit,j);
	  cos2_rf=get(elHiggscosthetaNT2Refit,j);
	  cosStar_rf= get(elHiggscosthetastarNTRefit,j);
	  phi_rf= get(elHiggsphiNTRefit,j);
	  phiStar_rf= get(elHiggsphiNT1Refit,j);
	  hLD_rf= get(elHiggsHelyLDRefit,j);
	  higgsrefitMass = get(elHiggsRefitMass,j);
	  higgsMass = get(elHiggsMass,j);
	  lljjpt_ = get(elHiggsPt,j);
	  lljjY_ = get(elHiggsY,j);
	  higgsEvtNum = getInt(elHiggsEventNumber);
	  higgsRunNum = getInt(elHiggsRunNumber);
	  trkMet_ = get(elHiggstrkMet, j);
	  flavourJet1 = get(elHiggsJetDau1PartonFlavour, j);
	  flavourJet2 = get(elHiggsJetDau2PartonFlavour, j);
	}

	//kineLepCut = (lept1pt_ > 20 && lept2pt_ > 20  && zllcharge_ == 0) ;

	//kineLepCut = (lept1pt_ > 20 && lept2pt_ > 20 && Jet1pt_ > 40 && Jet2pt_ > 40 && zllcharge_ == 0) ;
	//	kineLepCut = (lept1pt_ > 40 && lept2pt_ > 20 && Jet1pt_ > 40 && zllcharge_ == 0) ;
	//	kineLepCut = (lept1pt_ > 40 && lept2pt_ > 20 && Jet1pt_ > 40 && zllcharge_ == 0) ;
	
	//	kineLepCut = (lept1pt_ > 40 && lept2pt_ > 20 && Jet1pt_ > 40 && Jet2pt_ > 30 && zllcharge_ == 0) ;
	
	//	kineLepCut = (lept1pt_ > 20 && lept2pt_ > 20 && Jet1pt_ > 30 && Jet2pt_ > 30 && zllcharge_ == 0) ;
	kineLepCut = (lept1pt_ > 40 && lept2pt_ > 20 && Jet1pt_ > 30 && Jet2pt_ > 30 ) ;



	float met_ = getInt(met);
	float t1corrMET_ = getInt(t1corrMet);
	float metPhi_ = getInt(metPhi);
	minMET_ = TMath::Min(met_, trkMet_);
	minCorrMET_ = TMath::Min(trkMet_, t1corrMET_);

	float dPhi = (TMath::Cos(metPhi_)*TMath::Cos(zllphi_)) + (TMath::Sin(metPhi_)*TMath::Sin(zllphi_)); 
	float JZB_0 = TMath::Sqrt( TMath::Power(met_,2) + TMath::Power(zllpt_,2) + (2*met_*zllpt_*dPhi) ) ;
	float JZB_1 = TMath::Abs(JZB_0);
	float JZB_2 = TMath::Abs(zllpt_);
	float JZB = JZB_1 - JZB_2;
	float DRlep1_Jet1 = deltaR(Jet1eta_, Jet1phi_, Lept1eta_, Lept1phi_);
	float DRlep2_Jet1 = deltaR(Jet1eta_, Jet1phi_, Lept2eta_, Lept2phi_);
	float DRlep1_Jet2 = deltaR(Jet2eta_, Jet2phi_, Lept1eta_, Lept1phi_);
	float DRlep2_Jet2 = deltaR(Jet2eta_, Jet2phi_, Lept2eta_, Lept2phi_);
	
	float minDRJet1 = TMath::Min(DRlep1_Jet1, DRlep2_Jet1);
	float minDRJet2 = TMath::Min(DRlep1_Jet2, DRlep2_Jet2);

	bool JZB_SB = JZB>50;

	jptmjj.Fill(zjjmass_, Jet1pt_, evtWeight );
	
	if(zllcharge_ == 0){
	  pass[3] = true;
	  pass_1b[3] = true;
	  pass_2b[3] = true;

	  if(dbLepCut && isoIDLepCut){
	    leptpt1.Fill(lept1pt_,evtWeight);
	    leptpt2.Fill(lept2pt_,evtWeight);
	    jetpt1.Fill(Jet1pt_,evtWeight);
	    jetpt2.Fill(Jet2pt_,evtWeight);
	  }

	  if(kineLepCut){
	    pass[4] = true;
	    pass_1b[4] = true;
	    pass_2b[4] = true;


	    if(dbLepCut){
	      pass[5] = true;
	      pass_1b[5] = true;
	      pass_2b[5] = true;
	      
	      if(isoIDLepCut){
		pass[6] = true;
		pass_1b[6] = true;
		pass_2b[6] = true;
		  

		zllmass.Fill(zllmass_,evtWeight);
		zllpt.Fill(zllpt_,evtWeight);	  
		zjjpt.Fill(zjjpt_,evtWeight);	  
		zjjmass.Fill(zjjmass_,evtWeight);
		mjjmll.Fill(zllmass_,  zjjmass_,  evtWeight);
		
		bool zllcut = (zllmass_ > 70) && (zllmass_ < 110);
		bool zjjcut_sigSel = (zjjmass_ > 75) && (zjjmass_ < 105);
		bool jet1_tagged_tight(false), jet1_tagged_medium(false), jet1_tagged_loose(false);
		bool jet2_tagged_tight(false), jet2_tagged_medium(false), jet2_tagged_loose(false);
		
		// jet1_tagged_medium = jet1tkhe_ > 3.3;
		// jet1_tagged_loose  = jet1tkhe_ > 1.7;	    
		// jet2_tagged_medium = jet2tkhe_ > 3.3;
		// jet2_tagged_loose  = jet2tkhe_ > 1.7;
		
		jet1_tagged_tight  = jet1csv_ > 0.898;
		jet1_tagged_medium = jet1csv_ > 0.679;
		jet1_tagged_loose  = jet1csv_ > 0.244;	    
		jet2_tagged_tight  = jet2csv_ > 0.898;
		jet2_tagged_medium = jet2csv_ > 0.679;
		jet2_tagged_loose  = jet2csv_ > 0.244;

		
		// eventually apply SF for MC
		if( btagScale ) {
		  btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, Jet1pt_, Jet1eta_, jet1flav_ );
		  btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, Jet2pt_, Jet2eta_, jet2flav_ );
		}
		
		bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_medium  ) ;
		bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_medium || jet2_tagged_medium );
		bool zeroTagCategory = ( !twoTagCategory ) && ( !oneTagCategory );

		// bool twoTagCategory  = ( jet1_tagged_tight && jet2_tagged_tight  ) ;
		// bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_tight || jet2_tagged_tight );
		// bool zeroTagCategory = ( !twoTagCategory ) && ( !oneTagCategory );
		
		csvmjj.Fill(zjjmass_, jet1csv_, evtWeight);
		csvmjj.Fill(zjjmass_, jet2csv_, evtWeight);
		
		if(zllcut){
		  pass[7] = true;
		  pass_1b[7] = true;
		  pass_2b[7] = true;

		  if(twoTagCategory ) {
		  
		  dphijjmjj.Fill(zjjmass_, jjdphi_, evtWeight );
		  detajjmjj.Fill(zjjmass_, jjdeta_, evtWeight );
		  drzzmjj.Fill(zjjmass_, zzdr_,   evtWeight);
		  ptjjmjj.Fill(zjjmass_, zjjpt_,evtWeight);
		  drllmjj.Fill(zjjmass_, lldr_,evtWeight);	   
		  }
		  
		  if(zzdphi_>1.){
		    pass[8] = true;
		    pass_1b[8] = true;
		    pass_2b[8] = true;

		    /* 0 btag category*/ 
		    if ( zeroTagCategory){
		      pass[9] = true;
		      lljjmass_0btag.Fill(higgsMass,evtWeight);
		      met_0btag.Fill(getInt(met),evtWeight);
		      jzb_0btag.Fill(JZB,evtWeight);
		      drJ1L_0btag.Fill(minDRJet1, evtWeight);
		      drJ2L_0btag.Fill(minDRJet2, evtWeight);
		      stdCandle_0btag.Fill(zllmass_,evtWeight);
		      t1corrMET_0btag.Fill(t1corrMET_,evtWeight);
		      
		      if( (zjjpt_ - zllpt_)>-50 && (zjjpt_ - zllpt_)<40  ){
			pass[10] = true;
			//if( (zjjpt_ - zllpt_)>-50 && (zjjpt_ - zllpt_)<40 &&  (trkMet_<60) ){
			//	    if( JZB>-40 && JZB<25 && zzdphi_>1. ){
			if(minMET_<55){
			  pass[11] = true;
			  //if(TMath::Abs(zjjmass_-91)< TMath::Abs(zjjmassBest_0 - 91)) {
			  if(Jet1pt_ > jetpt1Best_0){
			    if(flavourJet2== 5 && flavourJet2 == 5)bb_0btagBest = zjjmass_;
			    if(flavourJet2 != 5 && flavourJet2 != 5) notbb_0btagBest =zjjmass_;
			    njets_0 = njets;
			    zjjmassBest_0 = zjjmass_;
			    jetpt1Best_0 =Jet1pt_;
			  }
			  ++nCands_0;
			  if (!hasCand_0){
			    hasCand_0 = true;
			    ++nEvts_0;
			  }
			  
			  if(zjjcut_sigSel){
			    stdCandle_final_0.Fill(zllmass_,evtWeight);
			    jjm_0btag.Fill(zjjmass_ ,evtWeight);
			  }
			}
		      }
		      if( JZB_SB ){
			if(TMath::Abs(zjjmass_-91)< TMath::Abs(zjjmassBestSB_0 - 91)) zjjmassBestSB_0 = zjjmass_;
		      }
		    }
		    /* 1 btag category*/ 
		    if ( oneTagCategory){
		      pass_1b[9] = true;
		      lljjmass_1btag.Fill(higgsMass,evtWeight);
		      met_1btag.Fill(getInt(met),evtWeight);
		      jzb_1btag.Fill(JZB,evtWeight);
		      drJ1L_1btag.Fill(minDRJet1, evtWeight);
		      drJ2L_1btag.Fill(minDRJet2, evtWeight);
		      stdCandle_1btag.Fill(zllmass_,evtWeight);
		      t1corrMET_1btag.Fill(t1corrMET_,evtWeight);
		      trkMET_1btag.Fill(trkMet_,evtWeight);
		      
		      if( (zjjpt_ - zllpt_)>-50 && (zjjpt_ - zllpt_)<40 ){
			pass_1b[10] = true;
			//if( (zjjpt_ - zllpt_)>-50 && (zjjpt_ - zllpt_)<40 &&  (trkMet_<60)){
			//	    if((JZB>-40) && (JZB<25) && (zzdphi_>1.)){
			
			if(minMET_<55){
			  pass_1b[11] = true;
			  //if(TMath::Abs(zjjmass_-91)< TMath::Abs(zjjmassBest_1 - 91)) {
			  if(Jet1pt_ > jetpt1Best_1){
			    zjjmassBest_1 = zjjmass_;
			    jetpt1Best_1 =Jet1pt_;
			  }
			  ++nCands_1;
			  if (!hasCand_1){
			  hasCand_1 = true;
			  ++nEvts_1;
			  }
			  
			  if(zjjcut_sigSel){
			    stdCandle_final_1.Fill(zllmass_,evtWeight);
			    jjm_1btag.Fill(zjjmass_ ,evtWeight);
			  }
			}
		      }
		      if( JZB_SB ){
			if(TMath::Abs(zjjmass_-91)< TMath::Abs(zjjmassBestSB_1 - 91)) zjjmassBestSB_1 = zjjmass_;
		      }
		    }
		    /* 2 btag category*/ 
		    if (twoTagCategory){
		      pass_2b[9] = true;
		      lljjmass_2btag.Fill(higgsMass,evtWeight);
		      met_2btag.Fill(getInt(met),evtWeight);
		      jzb_2btag.Fill(JZB,evtWeight);
		      zzdpt.Fill(zjjpt_ - zllpt_, evtWeight);
		      drJ1L_2btag.Fill(minDRJet1, evtWeight);
		      drJ2L_2btag.Fill(minDRJet2, evtWeight);
		      drJ1L_mjj.Fill(zjjmass_, minDRJet1,evtWeight);
		      drJ2L_mjj.Fill(zjjmass_, minDRJet2,evtWeight);
		      jzbmjj.Fill(zjjmass_, JZB ,evtWeight);
		      jzbmjjcorr.Fill(zjjmass_, (JZB - zjjmass_) ,evtWeight);
		      zzdptmjj.Fill(zjjmass_, zllpt_ - zjjpt_ ,evtWeight);
		      stdCandle_2btag.Fill(zllmass_,evtWeight);
		      t1corrMET_2btag.Fill(t1corrMET_,evtWeight);
		      //	    trkMET_2btag.Fill(trkMet_,evtWeight);
		      if( (zjjpt_ - zllpt_)>-50 && (zjjpt_ - zllpt_)<40 ){
			pass_2b[10] = true;
		      //if( (zjjpt_ - zllpt_)>-50 && (zjjpt_ - zllpt_)<40 &&  (trkMet_<60)){
		      //	    if( (JZB>-40) && (JZB<25) && (zzdphi_>1.) ){
			minMET_2btag.Fill(minMET_,evtWeight);
			minCorrMET_2btag.Fill(minCorrMET_,evtWeight);
			
			if(minMET_<55){
			  pass_2b[11] = true;
			  //if(trkMet_<60){
			//if(TMath::Abs(zjjmass_ - 91)< TMath::Abs(zjjmassBest_2 - 91)) {
			  if(flavourJet2 == 21 && flavourJet2 == 21) gg_2btag.Fill(zjjmass_, evtWeight);
			  if(Jet1pt_ > jetpt1Best_2){
			    njets_2 = njets;
			    if(flavourJet2 == 5 && flavourJet2 == 5) bb_2btagBest = zjjmass_;
			    if(flavourJet2 != 5 && flavourJet2 != 5) notbb_2btagBest =zjjmass_;
			    jetpt1Best_2 =Jet1pt_;
			    zjjmassBest_2 = zjjmass_;	     
			    trkMetBest_2 =trkMet_;
			    lljjmassBest_2 = higgsMass;
			    llptBest_2 =  zllpt_;
			    jjptBest_2 =  zjjpt_;
			    jjdphiBest_2 = jjdphi_ ;
			  }
			  ++nCands_2;
			  if (!hasCand_2){
			    hasCand_2 = true;
			    ++nEvts_2;
			  }
			  
			  if( zjjcut_sigSel  ){
			    stdCandle_final_2.Fill(zllmass_,evtWeight);
			    jjm_2btag.Fill(zjjmass_ ,evtWeight);
			  }
			}
		      }	  
		      if( JZB_SB  ){
			if(TMath::Abs(zjjmass_-91)< TMath::Abs(zjjmassBestSB_2 - 91)) zjjmassBestSB_2 = zjjmass_;
		      }
		      
		  }
		  } // end cut on delta Phi zz
		} // end of zllcut 
		
		if(zjjcut_sigSel && zllcut){
		  //		  zjjpt.Fill(zjjpt_,evtWeight);	  
		  
		  //jetpt1.Fill(Jet1pt_,evtWeight);
		  //jetpt2.Fill(Jet2pt_,evtWeight);
		  drjj.Fill(drjj_,evtWeight);
		  btagJet1.Fill(jet1tkhe_,evtWeight);
		  btagJet2.Fill(jet2tkhe_,evtWeight);
		  btagJet.Fill(jet1tkhe_,evtWeight);
		  btagJet.Fill(jet2tkhe_,evtWeight);
		  btagJet1_CSV.Fill(jet1csv_,evtWeight);
		  btagJet2_CSV.Fill(jet2csv_,evtWeight);
		  btagJet_CSV.Fill(jet1csv_,evtWeight);
		  btagJet_CSV.Fill(jet2csv_,evtWeight);
		  h_met.Fill(getInt(met),evtWeight);
		  metSig.Fill(getInt(metSignificance),evtWeight);
		  lljjmass.Fill(higgsMass,evtWeight);
		  zzdr.Fill(zzdr_,evtWeight);
		  helyLD.Fill(hLD_,evtWeight);
		  cosT1.Fill(cos1_,evtWeight);
		  cosT2.Fill(cos2_,evtWeight);
		  cosT1Star.Fill(cosStar_,evtWeight);
		  phi.Fill(phi_,evtWeight);
		  phiStar.Fill(phiStar_,evtWeight);
		  metSignif.Fill(getInt(metSignificance),evtWeight);
		  trkMET_MET.Fill(met_, trkMet_, evtWeight);
		  // select only 2-tag candidates
		  if(twoTagCategory){
		    
	      
		    ptllptjj.Fill(zjjpt_, zllpt_, evtWeight);		
		    drjj_btag.Fill(drjj_,evtWeight);
		    lljjmass_btag.Fill(higgsMass,evtWeight);
		    zzdr_btag.Fill(zzdr_,evtWeight);
		    //	      zzdpt.Fill(zllpt_ - zjjpt_, evtWeight);
		    zzdphi.Fill(zzdphi_,evtWeight);
		    zzdeta.Fill(zzdeta_,evtWeight);
		    jjdphi.Fill(jjdphi_,evtWeight);
		    jjdeta.Fill(jjdeta_,evtWeight);
		    lldr.Fill(lldr_,evtWeight);
		    lljjpt.Fill(lljjpt_,evtWeight);
		    lljjY.Fill(lljjY_,evtWeight);		
		    
		  }
		}
		
	      }
	    } // isolation requirement
	  } // db cut
	}// end kinematic selection
	} // Leptonic Z charge == 0
	if(zjjmassBest_0!=10000) zjjmass_0btag.Fill(zjjmassBest_0,evtWeight);
	if(zjjmassBest_1!=10000) zjjmass_1btag.Fill(zjjmassBest_1,evtWeight);
	if(zjjmassBest_2!=10000) zjjmass_2btag.Fill(zjjmassBest_2,evtWeight);
	//if(jetpt1Best_0!=-1) jetpt1_0btag.Fill(jetpt1Best_0,evtWeight);
	//if(jetpt1Best_2!=-1) jetpt1_2btag.Fill(jetpt1Best_2,evtWeight);
	if(zjjmassBestSB_0!=10000)zjjmass_sbTT_0btag.Fill(zjjmassBestSB_0,evtWeight);
	if(zjjmassBestSB_1!=10000)zjjmass_sbTT_1btag.Fill(zjjmassBestSB_1,evtWeight);
	if(zjjmassBestSB_2!=10000)zjjmass_sbTT_2btag.Fill(zjjmassBestSB_2,evtWeight);
	if(trkMetBest_2!=-1)trkMET_2btag.Fill(trkMetBest_2,evtWeight);   
	if(lljjmassBest_2!=0)lljjmass_2btag.Fill(lljjmassBest_2,evtWeight);
	if(llptBest_2!=0)llpt_2btag.Fill(llptBest_2,evtWeight); 
	if(jjptBest_2!=0)jjpt_2btag.Fill(jjptBest_2,evtWeight);
	if(jjdphiBest_2!=0)jjdphi_2btag.Fill(jjdphiBest_2,evtWeight); 
	if(bb_2btagBest!= -1) bb_2btag.Fill(bb_2btagBest,evtWeight);
	if(notbb_2btagBest!= -1) notbb_2btag.Fill(notbb_2btagBest,evtWeight);
	if(bb_0btagBest!= -1) bb_0btag.Fill(bb_0btagBest,evtWeight);
	if(notbb_0btagBest!= -1) notbb_0btag.Fill(notbb_0btagBest,evtWeight);
	if(njets_2!=1000) njets_2btag.Fill(njets_2,evtWeight);
	if(njets_0!=1000) njets_0btag.Fill(njets_0,evtWeight);
      }//end loop on candidates
    // }//end if cands>0
    
    for(unsigned int k = 0; k < counts; ++k){
      if(pass[k]) ++count[k];
      if(pass_1b[k]) ++count_1[k];
      if(pass_2b[k]) ++count_2[k];
    }
    
  } // event loop
      
  // for(unsigned int k = 0; k < counts; ++k){
  //    cout<< "TEST---0btag--1btag--2btag-- "<<count[k]<<endl;  
  //    cout<< "TEST---"<<count[k]<<" "<<count_1[k]<<" "<<count_2[k]<<endl;  
  // }


   
   //--- write tree--
    // TTree output file
    if(writeT) {
      TFile *outTree = TFile::Open(("tree_"+outFileName).c_str(), "RECREATE");
      TTree lljjmassTree("lljjmassTree","lljjmass tree");
      // mzz branch for the output tree
      double mzz;
      double mzz_sb;
      lljjmassTree.Branch("mzz", &mzz, "mzz/D");
      lljjmassTree.Branch("mzz_sb", &mzz_sb, "mzz_sb/D");
      for(vfloat::iterator it=bestHiggsMass.begin(); it!= bestHiggsMass.end(); ++it ){
	mzz = *it;
      lljjmassTree.Fill();
      }
      for(vfloat::iterator it=bestHiggsMass_1btag.begin(); it!= bestHiggsMass_1btag.end(); ++it ){
	mzz_sb = *it;
	lljjmassTree.Fill();
      }
      lljjmassTree.Write();
      outTree->Close();
    }
    //----------------
  
  
  cout << endl << "--> Number of candidates/category after common selection <--" << endl;
  for (int kkk = 0; kkk<3; ++kkk)
    cout << kkk << "-Tags candidates = " << CandBTag[kkk] << endl;

  TH1F h_cuts("h_cuts", "ProgressiveCuts", 11, 0, 11.);
  
  for (unsigned int k = 1; k < counts+1; ++k){
    cutFlow_0btag.SetBinContent(k, count[k-1]);
    cutFlow_1btag.SetBinContent(k, count_1[k-1]);
    cutFlow_2btag.SetBinContent(k, count_2[k-1]);
   }

  // cutFlow_0btag.SetBinContent(counts, zjjmass_0btag.GetEntries());
  // cutFlow_1btag.SetBinContent(counts, zjjmass_1btag.GetEntries());
  // cutFlow_2btag.SetBinContent(counts, zjjmass_2btag.GetEntries());

  //  cout<<"test on cutFlow - number of candidate after best cand selection for 0btag: "<< cutFlow_0btag[counts-1] <<endl;
  // cout<<"test on cutFlow "<< endl;
  // for (unsigned int k = 0; k < counts; ++k)
  //   cout<<cutFlow_0btag.GetBinContent(k+1)<<endl;
  // cout<<"test on cutFlow 2btag: "<< endl;
  // for (unsigned int k = 0; k < counts; ++k)
  //   cout<<cutFlow_1btag.GetBinContent(k+1)<<endl; 
  // for (unsigned int k = 0; k < counts; ++k)
  //   cout<<cutFlow_2btag.GetBinContent(k+1)<<endl; 

  // output histogram file
  //  TFile histos(outFile, "RECREATE");
  TFile histos(outFile, "RECREATE");
  histos.cd();

  const char* labels[] = {"events", "events with at least 1 cand", "veto on #jets", "Zll charge = 0", "pt cut", "db cut", "iso cut", "zll mass cut", "zzdphi cut", "btag", "zzdpt cut", "minMET cut", "bestCandChoice"};
  for (int k = 1; k < int(counts +1); ++k){
    cutFlow_0btag.GetXaxis()->SetBinLabel(k,labels[k-1]);	 
    cutFlow_1btag.GetXaxis()->SetBinLabel(k,labels[k-1]);	 
    cutFlow_2btag.GetXaxis()->SetBinLabel(k,labels[k-1]);	 
  } 

  cutFlow_0btag.Write();
  cutFlow_1btag.Write();
  cutFlow_2btag.Write();

  mjjmll.Write();
  lljjmass.Write();  
  lljjmass_0btag.Write();  
  lljjmass_1btag.Write();  
  lljjmass_2btag.Write();  


  lljjpt.Write();  
  lljjY.Write();  
  lljjmass_btag.Write();  

  /*  lljjmass_noRefit.Write();  
  lljjmass_1btag.Write();  
  lljjmass_noRefit_1btag.Write();
  lljjmass_0btag.Write();  
  lljjmass_noRefit_0btag.Write();  
  */
  //  lljjmassCands.Write();  
  
  jetpt1_2btag.Write();
  jetpt1_0btag.Write();

  leptpt1.Write();  
  leptpt2.Write();  
  jetpt1.Write();  
  jetpt2.Write();  
  h_met.Write();  
  metSig.Write();  
  metSignif.Write();  
  trkMET_MET.Write();

  bb_0btag.Write();
  notbb_0btag.Write();
  bb_2btag.Write();
  notbb_2btag.Write();
  gg_2btag.Write();

  minMET_2btag.Write();
  minCorrMET_2btag.Write();

  /*
  lljjmass2tags.Write();
  helyLD_RF2tags.Write();
  
  zllptBeforeLD.Write();
  drjjBefLD.Write();
  lljjmassBeforeLD.Write();
  */

  cosT1.Write();
  cosT2.Write();
  cosT1Star.Write();
  phi.Write();
  phiStar.Write();
  helyLD.Write();

  njets_0btag.Write();
  njets_2btag.Write();

  
  /*
    cosT1_RF.Write();
    cosT2_RF.Write();
    cosT1Star_RF.Write();
    phi_RF.Write();
    phiStar_RF.Write();
    helyLD_RF.Write();
  */
  btagJet1.Write();  
  btagJet2.Write();  
  btagJet.Write();
  btagJet1_CSV.Write();  
  btagJet2_CSV.Write();  
  btagJet_CSV.Write();

  drjj.Write();  
  drjj_btag.Write();  
  lldr.Write();  
  zzdr.Write(); 
  zzdeta.Write();  
  zzdphi.Write();  
  jjdeta.Write();  
  jjdphi.Write();  
  zzdr_btag.Write();  
  zllpt.Write();
  zjjpt.Write();
  zzdpt.Write();
  zjjmass.Write();
  zllmass.Write();
  zjjmass_0btag.Write();
  zjjmass_1btag.Write();
  zjjmass_2btag.Write();
  zjjmass_sbTT_0btag.Write();
  zjjmass_sbTT_1btag.Write();
  zjjmass_sbTT_2btag.Write();
  met_0btag.Write();
  met_1btag.Write();
  met_2btag.Write();
  t1corrMET_0btag.Write();
  t1corrMET_1btag.Write();
  t1corrMET_2btag.Write();

  trkMET_0btag.Write();
  trkMET_1btag.Write();
  trkMET_2btag.Write();

  jzb_0btag.Write();
  jzb_1btag.Write();
  jzb_2btag.Write();

  llpt_2btag.Write();
  jjpt_2btag.Write();
  jjdphi_2btag.Write();

  jzbmjj.Write();
  jzbmjjcorr.Write();
  zzdptmjj.Write();
  jptmjj.Write();

  jjm_0btag.Write();
  jjm_1btag.Write();
  jjm_2btag.Write();

  drJ1L_0btag.Write();
  drJ1L_1btag.Write();
  drJ1L_2btag.Write();
  drJ2L_0btag.Write();
  drJ2L_1btag.Write();
  drJ2L_2btag.Write();

  drJ1L_mjj.Write();
  drJ2L_mjj.Write();

  db.Write();

  csvmjj.Write();

  drzzmjj.Write();
  ptjjmjj.Write();
  drllmjj.Write();
  ptllptjj.Write();
  dphijjmjj.Write();
  detajjmjj.Write();


  stdCandle_final_0.Write();
  stdCandle_final_1.Write();
  stdCandle_final_2.Write();
  stdCandle_0btag.Write();
  stdCandle_1btag.Write();
  stdCandle_2btag.Write();

  //  ngenint.Write();
  npv_woReweight.Write();
  npv.Write();
  //  npv1.Write();
  //  npv2.Write();
  //  npvfin.Write();
  
   /*
	hmass_sb_all.Write();
	hmassnorefit_sb_all.Write();
	hmass_sb_2b.Write();
	hmassnorefit_sb_2b.Write();
	lljjmass_sb.Write();
	lljjmass_noRefit_sb.Write();
	lljjmassFinal_sb.Write();
	lljjmassFinal_noRefit_sb.Write();
  */

  cout << endl;

  cout<<"------------------------------------------------------------------------------------"<<endl;
  cout<<"# of evts with at least 1 cand in 0btag category / # of Candidates in 0btag category"<<endl;
  cout<<nEvts_0<<" / "<<nCands_0<<endl;
  cout<<"Events in the histo: "<<zjjmass_0btag.Integral()<<endl;
  cout<<"Events in the histo(not scaled): "<<zjjmass_0btag.GetEntries()<<endl;
  cout<<"------------------------------------------------------------------------------------"<<endl;
  cout<<"# of evts with at least 1 cand in 1btag category / # of Candidates in 1btag category"<<endl;
  cout<<nEvts_1<<" / "<<nCands_1<<endl;
  cout<<"Events in the histo: "<<zjjmass_1btag.Integral()<<endl;
  cout<<"Events in the histo(not scaled): "<<zjjmass_1btag.GetEntries()<<endl;
  cout<<"------------------------------------------------------------------------------------"<<endl;
  cout<<"# of evts with at least 1 cand in 2btag category / # of Candidates in 2btag category"<<endl;
  cout<<nEvts_2<<" / "<<nCands_2<<endl;
  cout<<"Events in the histo: "<<zjjmass_2btag.Integral()<<endl;
  cout<<"Events in the histo(not scaled): "<<zjjmass_2btag.GetEntries()<<endl;
  cout<<"------------------------------------------------------------------------------------"<<endl;
  /*  cout<<"#TotEvts EvtsAtLeast1Cand Pt+Iso+db zll/zjjMmass   btag  HelyLD  metSig"<<endl;
  string names[]={"#TotEvents","basicSel", "pt/Eta/Iso", "dB", "zll/zjjMass", "btag",  "HelyLD", "met"};

  cout << "selected events: "<<endl;
  for(unsigned int k = 0; k < counts; ++k) {
  cout << count[k] << " ";
    h_cuts.SetBinContent(k,count[k]);
    h_cuts.GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    }
  h_cuts.Write();
  cout<<endl;
  cout << "selected candidates: "<<endl;
  for(unsigned int k = 0; k < counts; ++k) {
  cout << cand[k] << " ";
  }
  cout << endl;
  cout<<"number of candidates passing entire selection: "<<lljjmass.Integral();
  cout<<endl;
  cout<<"number of events passing entire selection: "<<count[counts-1]<< endl;
  */


  //Efficiencies
  cout << "efficiencies: ";
  for(unsigned int k = 0; k < counts; ++k) {
    cout << double(count[k])/double(count[0]) << " ";
  }
  cout << endl;
  //

  //Writing efficiencies to a txt file
  /*
    string EffName = "Efficiencies"+selChan+".txt";
    ofstream fileeffout(EffName.c_str(),ios::app);
    fileeffout <<"Ntupla: "<< path << endl;
    fileeffout <<"Progressive Efficiencies"<< endl;
    fileeffout <<"EfficiencyPgr - 1 Candidate:                  "<< double(count[1])/double(count[0])<< endl;
    fileeffout <<"EfficiencyPgr - Kine, db and Iso on Leptons:  "<< double(count[2])/double(count[0])<< endl;
    fileeffout <<"EfficiencyPgr - zll and zjj:                  "<< double(count[3])/double(count[0])<< endl;
    fileeffout <<"EfficiencyPgr - 2-Tag:                        "<< double(count[4])/double(count[0])<< endl;
    fileeffout <<"EfficiencyPgr - Met Significance:             "<< double(count[5])/double(count[0])<< endl;
    fileeffout <<"EfficiencyPgr - LD:                           "<< double(count[6])/double(count[0])<< endl;
    fileeffout <<"Factorized Efficiencies"<< endl;
    fileeffout <<"EfficiencyFact - 1 Candidate:                 "<< double(count[1])/double(count[0])<< endl;
    fileeffout <<"EfficiencyFact - Kine, db and Iso on Leptons: "<< double(count[2])/double(count[1])<< endl;
    fileeffout <<"EfficiencyFact - zll and zjj:                 "<< double(count[3])/double(count[2])<< endl;
    fileeffout <<"EfficiencyFact - 2-Tag:                       "<< double(count[4])/double(count[3])<< endl;
    fileeffout <<"EfficiencyFact - Met Significance:            "<< double(count[5])/double(count[4])<< endl;
    fileeffout <<"EfficiencyFact - LD:                          "<< double(count[6])/double(count[5])<< endl;
    fileeffout <<"Counting"<< endl;
  */


  //fileeffout <<"Events passing selection:                     "<<count[counts-1]<< endl;
  //fileeffout <<"Candidates passing selection:                 "<<lljjmass.Integral()<< endl;
  //fileeffout fileeffout.close();
  //end Efficiencies to txt write

  return 0;
}

