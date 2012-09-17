#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include "DataFormats/Common/interface/Wrapper.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "HiggsAnalysis/Higgs2l2b/interface/BTagSFUtil.h"
#include "HiggsAnalysis/Higgs2l2b/interface/QGLikelihoodCalculator.h"
#include "HiggsAnalysis/Higgs2l2b/interface/MuonEffectiveArea.h"
#include "HiggsAnalysis/Higgs2l2b/interface/ElectronEffectiveArea.h"
#include "HiggsAnalysis/Higgs2l2b/interface/LineshapeWeight.h"
#include "../interface/LeptonID.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HiggsAnalysis/Higgs2l2b/interface/table.h"
#include "HiggsAnalysis/Higgs2l2b/interface/Weights.h"


using namespace std;

typedef vector<double> vdouble;
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

#define BRANCHDOUBLE(name) \
  edm::Wrapper<double> * name = new edm::Wrapper<double>(); \
  string alias_##name = events->GetAlias(#name); \
  alias_##name.erase(alias_##name.size() - 3);  \
  TBranch * b_##name = events->GetBranch(alias_##name.c_str());         \
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

// LR --> apply Lineshape Reweighting (for powheg MC only)

int main(int argc, char **argv) {
  std::cout<<"Let's start"<<endl;
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
  //  string fixMuRunB(argv[9]);
  string scaleFactor(argv[9]);
  string applyLR_(argv[10]);
  string hmasshyp(argv[11]);
  string outFileName(outFile);

  bool data = false;
  bool btagScale = false;
  //  bool fixMu = false;
  std::stringstream ss;
  ss<<scaleFactor;
  float scaleFact;
  ss>>scaleFact;
  cout<< "scaleFactor "<<scaleFact << endl;

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

  bool applyLR(false);
  if( applyLR_ == "LR"){ 
    applyLR = true;
    cout << "Lineshape Reweighting will be applied!!! (Powheg signal MC only)" << endl;
 }

  //  if(fixMuRunB == "fixMu"){
  //    fixMu=true;
  //    cout << "Mu FIX for runB applied" << endl;
  //  }

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

  vdouble bestHiggsMass, bestmjj, bestHiggsMass_1btag, bestHiggsMass_0btag, weight, btagFlag, channel, puWeight, lumiWeight;


  TFile *f = TFile::Open(path, "READ");
  TTree * events = (TTree*) f->Get("Events");
  unsigned int nEvents = events->GetEntries();
  cout << "events: " << nEvents << endl;

  // GET BRANCHES
  #include "HiggsAnalysis/Higgs2l2b/bin/2l2b_branches_2012.h"
  
  cout << " DONE BRANCHES" << endl;

  // txt output file
  string reportName = "SelectedEvents"+selChan+".txt";
  ofstream fileout;
  fileout.open(reportName.c_str(),ios::in | ios::out |ios::app);
  fileout<<"RunNumber "<<"  "
	 <<"EventNumber"<< "  " 
	 <<"LumiSect"<< "  " 
	 <<"leppt_1"<< "  " 
	 <<"leppt_2"<< "  " 
	 <<"zllmass"<< "  " 
	 <<"jetpt_1_ref"<< "  " 
	 <<"jetpt_2_ref"<< "  " 
	 <<"zjjmass"<< "  " 
	 << "lljj_m "<< "  " 
	 << "qgprod " <<  "  " 
	 << "metsig" << "  " 
	 << "LD" << "  " 
	 << "LD_cut" << "  " 
	 << "b-cat" << endl;

  // string reportNameCheck = "SelectedEventsCheck"+selChan+".txt";
  // ofstream fileoutCheck(reportName.c_str(),ios::trunc);
  // fileoutCheck<< "Event "<<"  "
  // 	 <<"leppt_1"<< "  " 
  // 	 <<"leppt_2"<< "  " 
  // 	 <<"zllmass"<< "  " 
  // 	 <<"jetpt_1_ref"<< "  " 
  // 	 <<"jetpt_2_ref"<< "  " 
  // 	 <<"zjjmass"<< "  " 
  // 	 <<"b-cat"<< endl;



  // define histograms
  // control histograms for lineshape reweighting
  TH1F genhmass("genh_mass", "genh_mass", 500, 200, 1200);
  TH1F genhmass_rew("genh_mass_rew", "genh_mass_rew", 500, 200, 1200);
  TH1F LR_weights("LR_wei", "LR_wei", 100, 0, 5);
  TH2F LR_weights_mass("wei_mass", "wei_mass", 250, 200, 1200, 50, 0, 5);
  // histo for event variables
  TH1F npv_woReweight("npv_woReweight", "npvwoReweight", 31, -0.5, 30.5);
  TH1F npv("npv", "npv", 31, -0.5, 30.5);
  TH1F npv1("npv1", "npv", 31, -0.5, 30.5);
  TH1F npv2("npv2", "npv", 31, -0.5, 30.5);
  TH1F npvfin("npvfin", "npv", 31, -0.5, 30.5);
  //  TH1F ngenint("ngenint", "ngenint", 36, -0.5, 35.5);
  // histo before lepton cuts
  TH1F db("db", "db", 100, 0, 0.2);
  // histo before Z mass cuts
  TH1F zllmass("zllmass", "m_{ll}", 500, 0, 500);
  TH1F zllmass_BB("zllmass_BB", "m_{ll}", 500, 0, 500);
  // histo after Z mass cuts
  TH1F stdCandle("stdCandle", "m_{ll}", 50, 65, 115);
  TH1F leptpt1("leptpt1", "P_{T} ", 300, 0, 300);
  TH1F leptpt2("leptpt2", "P_{T} ", 300, 0, 300);
  TH1F jetpt1("jetpt1", "P_{T} ", 300, 0, 300);
  TH1F jetpt2("jetpt2", "P_{T} ", 300, 0, 300);
  TH1F lepteta1("lepteta1", "#eta ", 100, -6., 6.);
  TH1F lepteta2("lepteta2", "#eta ", 100,-6., 6.);
  TH1F leptphi1("leptphi1", "#phi ", 100, -3., 3.);
  TH1F leptphi2("leptphi2", "#phi ", 100,-3., 3.);
  TH1F zjjmass("zjjmass", "m_{jj}", 500, 0, 500);

  TH1F elMVAID("ElMVAID", "MVA ID", 1000, -1.1, 1.1);
  TH1F beta("beta", "beta", 1000, -0.1, 1.1);

  TH2F eta_pt("pt_eta", "pt vs eta", 100, -6., 6., 300, 0., 300. );
  TH2F eta_phi("phi_eta", "phi vs eta", 100, -6., 6., 100, -3., 3. );
  TH2F pt_phi("pt_phi", "pt vs phi", 100, -3., 3., 300, 0., 300. );
  TH2F eta_zllmass("zllmass_eta", "Z_{#mu#mu} mass vs #eta",100, -6., 6., 500, 0, 500);
  // histo before b-tagging cut
  TH1F qgLD("qgd", "Quark-Gluon Discriminant", 100, 0, 1);
  TH1F njets("njets", "njets", 11, -0.5, 10.5);
  TH1F btagJet1("TCHEJet1", "TCHE", 100, 0, 20);
  TH1F btagJet2("TCHEJet2", "TCHE", 100, 0, 20);
  TH1F btagJet("TCHEJet", "TCHE", 100, 0, 20);
  TH1F btagJet1_CSV("CSVJet1", "CSV", 100, 0, 1);
  TH1F btagJet2_CSV("CSVJet2", "CSV", 100, 0, 1);
  TH1F btagJet_CSV("CSVJet", "CSV", 100, 0, 1);
  TH1F btagJet1_CSVMVA("CSVMVAJet1", "CSVMVA", 100, 0, 1);
  TH1F btagJet2_CSVMVA("CSVMVAJet2", "CSVMVA", 100, 0, 1);
  TH1F btagJet_CSVMVA("CSVMVAJet", "CSVMVA", 100, 0, 1);
  TH1F btagJet1_JP("JPJet1", "JP", 100, 0, 1);
  TH1F btagJet2_JP("JPJet2", "JP", 100, 0, 1);
  TH1F btagJet_JP("JPJet", "JP", 100, 0, 1);

  //TH2F btagJet_TCHECSV("TCHE_VS_CSV", "TCHECSV", 100, 0, 1 ,100, 0, 20);
  //TH2F btagJet1_TCHECSV("TCHE_VS_CSV_Jet1", "TCHECSV", 100, 0, 1 ,100, 0, 20);
  //TH2F btagJet2_TCHECSV("TCHE_VS_CSV_Jet2", "TCHECSV", 100, 0, 1 ,100, 0, 20);

  TH1F h_met("met", "MET", 100, 0, 140);
  TH1F metSig("metSig", "METSig", 100, 0, 15);
  // histo for 2b candidates after b-tagging
  TH1F zllpt("zllpt", "P_{T}(Z_{ll})", 300, 0, 300);
  TH1F drjj("DRjj", "#Delta R_{jj}", 100, 0, 4);
  TH1F metSignif("metSignif", "METSignif", 100, 0, 15);
  TH1F lljjmass2tags("lljjmass2tags", "m(lljj)", 1000, 0, 1000);
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
  
  TH1F helyLD_RF_zjj("HelyLDRefit_zjj", "HelyLDRefit", 100, 0, 1);
  TH1F cosT1_RF_zjj("cosTheta1Refit_zjj", "cosTheta1Refit", 100, -1, 1);
  TH1F cosT2_RF_zjj("cosTheta2Refit_zjj", "cosTheta2Refit", 100, 0, 1);
  TH1F cosT1Star_RF_zjj("cosTheta1StarRefit_zjj", "cosTheta1StarRefit", 100, -1, 1);
  TH1F phi_RF_zjj("phiRefit_zjj", "phiRefit", 100, -3, 3);
  TH1F phiStar_RF_zjj("phiStarRefit_zjj", "phiStarRefit", 100, -3, 3);
  

  // final inv mass histo
  TH1F lljjmass("lljjmass", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_0btag("lljjmass_0btag", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_1btag("lljjmass_1btag", "m(lljj)", 1000, 0, 1000);

  TH1F lljjmass_2btagSB("lljjmass_2btagSB", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_0btagSB("lljjmass_0btagSB", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_1btagSB("lljjmass_1btagSB", "m(lljj)", 1000, 0, 1000);

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
 TH1F lljjmassFinal_0btagsb("lljjmassFinal_0btagsb", "m(lljj)", 1000, 0, 1000);
 TH1F lljjmassFinal_1btagsb("lljjmassFinal_1btagsb", "m(lljj)", 1000, 0, 1000);


  const unsigned int counts = 7;
  vint count(counts, 0);
  vbool pass(counts, false);
  vint cand(counts, 0);
  //b-tag code initialization 
  int CandBTag[3];
  for(int kkk = 0; kkk<3 ; kkk++)
    CandBTag[kkk] = 0;
  // BTag SF utility

  BTagSFUtil* btsfutil = new BTagSFUtil(13);
  
  //TCHE SF
  // TFile* file_loose = TFile::Open("BTagPayloads_TCHEL.root");
  //  TFile* file_medium = TFile::Open("BTagPayloads_TCHEM.root");
  //CSV SF
  //TFile* file_loose = TFile::Open("BTagPayloads_CSVL.root");
  //TFile* file_medium = TFile::Open("BTagPayloads_CSVM.root");
  //
  //  btsfutil->set_fileLoose(file_loose);
  //  btsfutil->set_fileMedium(file_medium);

  cand[1]=0;

  // PU reweighting util
  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting Lumi3DWeights_;
  if(!data){
    //  --> Run2012A <--
    if( PUPeriod == "2012All"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2012" << endl;
      LumiWeights_ = edm::LumiReWeighting("PUDist_Summer12MC_Extended.root", "PUDist_Run2012AB_Truth_JSON_29Aug.root", "PUS7_Distr", "pileup");
    }
    //  --> Run2011A <--
    else if(PUPeriod == "2011A"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011A" << endl;
      Lumi3DWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011A_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    //  --> Run2011B <--
    else if(PUPeriod == "2011B"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011B" << endl;
      Lumi3DWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011B_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    //  --> Run2011All <--
    else if(PUPeriod == "2011All"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011" << endl;
      Lumi3DWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011All_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
   
      TString stringPUWei("Weight3D.root");
      
      const char * filePUwei = gSystem->FindFile("./", stringPUWei);
      if( filePUwei != NULL ) Lumi3DWeights_.weight3D_init("Weight3D.root");
      else Lumi3DWeights_.weight3D_init(1.0);
    }
    else cout << "ERROR!!!!! Run period not existent!!!!!!" << endl;

  
  
  }

  // Lineshape reweighting util
  LineshapeWeight * LRUtil;
  
  //QGLike discriminator
  string QGFilePDF = "QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root";
  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator( QGFilePDF );

  // LOOP on EVENTS

  vdouble runNumber, evtNumber,lumiBlock;
  vint runNumSkim;

  for(unsigned int i = 0; i <nEvents; ++i) {
   
    //    if(i%100 == 0) progress(double(i)/double(nEvents));
    ++count[0];
    // NUMBER OF CANDIDATES IN THE EVENT
    unsigned int cands(0);

    // GET ENTRIES FOR EVENT VARIABLES
    GETENTRY(numPV,i);
    GETENTRY(met,i);
    GETENTRY(metSignificance,i);
    GETENTRY(rhoRestrictedEta,i);
    GETENTRY(rho,i);

    // assign event weights 
    double evtWeight = 1.0;
    double PUWeight = 1.0;
    double IDweight_1 = 1.0;
    double IDweight_2 = 1.0;
    double IDweight = 1.0;
    double LRweight(1.0);
    double LRweightErrp(0.0);
    double LRweightErrm(0.0);
    
    // assign event weight for Lineshape Reweighting
    if(applyLR){
      LRUtil = new LineshapeWeight("../../../MMozer/powhegweight/data/mZZ_Higgs"+hmasshyp+"_8TeV_W.txt_I.txt");
      BRANCHFLOAT(genHiggsMass);    
      GETENTRY(genHiggsMass,i);    
    // generated higgs mass before Lineshape Reweighting
      genhmass.Fill(getInt(genHiggsMass));
      LRUtil->getWeight( getInt(genHiggsMass), LRweight, LRweightErrp, LRweightErrm);
      //      cout << "Hmass = " << getInt(genHiggsMass) << " LRweight = " << LRweight << endl;
      // generated higgs mass after Lineshape Reweighting
      genhmass_rew.Fill(getInt(genHiggsMass), LRweight);
      LR_weights.Fill(LRweight);
      LR_weights_mass.Fill(getInt(genHiggsMass), LRweight);
    }

    // get PU weight
    if(!data){
      //      BRANCHFLOAT(nGenInt);
      //      GETENTRY(nGenInt,i);
      BRANCHINT(nGenIntBXm1);
      BRANCHINT(nGenIntBX0);
      BRANCHINT(nGenIntBXp1);
      BRANCHFLOAT(trueNInt);
      //BRANCHFLOAT(rhoRestrictedEta);
      GETENTRY(nGenIntBXm1,i);
      GETENTRY(nGenIntBX0,i);
      GETENTRY(nGenIntBXp1,i);
      GETENTRY(trueNInt,i);
      //GETENTRY(rhoRestrictedEta, i);
      //      cout<<"numInt: "<< getInt(trueNInt) <<endl;
      if(PUPeriod == "2012All") {
	PUWeight = LumiWeights_.weight(getInt(trueNInt));
	//cout<<"WEIGHT: "<< evtWeight <<endl;
      }
      else PUWeight = Lumi3DWeights_.weight3D( getInt(nGenIntBXm1), getInt(nGenIntBX0), getInt(nGenIntBXp1) );
    }      

    evtWeight = PUWeight;

    // fill histo of numPV (reweighted and not)
    npv.Fill(getInt(numPV),evtWeight);
    npv_woReweight.Fill(getInt(numPV));
    //runNumber.push_back(get((muHiggsRunNumber));
    GETENTRY(CleanJetPt,i);
    unsigned int nj = CleanJetPt->product()->size();

    bool hasAlreadyCand = false;
    bool hasAlreadyCand_afterzjj = false;
    bool hasAlreadyCand_preZllcut = false;

    // including also LR to the event weight
    if(applyLR){
      evtWeight *= LRweight;
    }

    // get number of candidates in the event
    if(muChannel) {
      GETENTRY(muHiggsLeptDau1Pt,i);
      cands = muHiggsLeptDau1Pt->product()->size();
    } else {
      GETENTRY(elHiggsLeptDau1Pt,i);
      cands = elHiggsLeptDau1Pt->product()->size();
    }

    if(cands > 0) {
      pass[1] = true;
   
      //  std::cout<<"Entering cands >1 if statement"<<std::endl;

      //GET ENTRIES FOR MU/EL VARIABLES
      if(muChannel){
#include "HiggsAnalysis/Higgs2l2b/bin/2l2bMu_getent_2012.h"
      } else {
#include "HiggsAnalysis/Higgs2l2b/bin/2l2bEl_getent_2012.h"
      }
      //  std::cout<<"accessed to getent files"<<std::endl;

      //#include "HiggsAnalysis/Higgs2l2b/bin/2l2b_getent.h"
      
      // variables to store the best H candidate in the event
      //      float met_, metSig_;
      float lljjmass_(0), lljjmass_noRefit_(0), lljjmass_1btag_(0), lljjmass_noRefit_1btag_(0);
      float lljjmass_0btag_(0), lljjmass_noRefit_0btag_(0);
      float SumPtBest_(0), SumPtBest_Final(0), SumPtBest_Final_1btag(0);
      bool exist2tag(false), exist1tag(false), exist0tag(false);
      bool exist2tagSB(false), exist1tagSB(false), exist0tagSB(false);
      double ch2tag(-10), ch1tag(-10), ch0tag(-10);
      double ch2tagSB(-10), ch1tagSB(-10), ch0tagSB(-10);
      double puw2tag(-10), puw1tag(-10), puw0tag(-10);
      double puw2tagSB(-10), puw1tagSB(-10), puw0tagSB(-10);

      float ResidZllMassBest_(10000.), ResidZllMassBest_Final(10000.);
      float ResidZllMassBest1tag_(10000.);
      float ResidZllMassBest0tag_(10000.);
      float ResidZllMassBestSB0tag_(10000.);
      float ResidZllMassBestSB1tag_(10000.);
      float ResidZllMassBestSB2tag_(10000.);
      float hLD_Best(-10), cos1_Best(-10), cos2_Best(-10), cosStar_Best(-10);
      float phi_Best(-10), phiStar_Best(-10);
      float hLD_rfBest(-10), cos1_rfBest(-10), cos2_rfBest(-10), cosStar_rfBest(-10);
      float phi_rfBest(-10), phiStar_rfBest(-10);
      int EvtNum2Tag_(0), RunNum2Tag_(0), LumiBlock2Tag_(0);
      int EvtNum1Tag_(0), RunNum1Tag_(0), LumiBlock1Tag_(0);
      int EvtNum0Tag_(0), RunNum0Tag_(0), LumiBlock0Tag_(0);
      int EvtNum2TagSB_(0), RunNum2TagSB_(0), LumiBlock2TagSB_(0);
      int EvtNum1TagSB_(0), RunNum1TagSB_(0), LumiBlock1TagSB_(0);
      int EvtNum0TagSB_(0), RunNum0TagSB_(0), LumiBlock0TagSB_(0);
      float ld_2btag_(-10), ld_1btag_(-10), ld_0btag_(-10);
      float zllmass_2btag_(-10), zllmass_1btag_(-10), zllmass_0btag_(-10);
      float zjjmass_2btag_(-10), zjjmass_1btag_(-10), zjjmass_0btag_(-10);
      float zjjmass_SB2btag_(-10), zjjmass_SB1btag_(-10), zjjmass_SB0btag_(-10);
      float qgprod_2btag_(-10), qgprod_1btag_(-10), qgprod_0btag_(-10);
      float metsig_2btag_(-10), metsig_1btag_(-10), metsig_0btag_(-10);
      float jet1pt_2btag_(-10), jet1pt_1btag_(-10), jet1pt_0btag_(-10);
      float jet2pt_2btag_(-10), jet2pt_1btag_(-10), jet2pt_0btag_(-10);
      float lep1pt_2btag_(-10), lep1pt_1btag_(-10), lep1pt_0btag_(-10);
      float lep2pt_2btag_(-10), lep2pt_1btag_(-10), lep2pt_0btag_(-10);

      float SumPtBest_SB(0), SumPtBest_FinalSB(0), SumPtBest_FinalSB_1btag(0);
      float lljjmass_SB(0), lljjmass_noRefit_SB(0);
      float lljjmassFinal_SB(0), lljjmassFinal_noRefit_SB(0), lljjmassFinal_SB_1btag_(0), lljjmassFinal_noRefit_SB_1btag_(0), lljjmassFinal_SB_0btag_(0);


      // auxiliary variables for channel-independent selection
      float lept1ChHadIso(-1), lept2ChHadIso(-1), lept1NeuHadIso(-1), lept2NeuHadIso(-1), lept1PhotonIso(-1), lept2PhotonIso(-1), lept1PUIso(-1), lept2PUIso(-1);
      float iso1(-1), iso2(-1);
      float njets_(-1);
      float beta_(-1), beta2_(-1);
      float Jet1pt_(0), Jet2pt_(0), Jet1pt_rf(0), Jet2pt_rf(0);;
      float Jet1eta_(0), Jet2eta_(0);
      float Jet1phi_(0), Jet2phi_(0);
      float lept1pt_(0), lept2pt_(0);
      float lept1eta_(0), lept2eta_(0);
      float lept1phi_(0), lept2phi_(0);

      float zllcharge_(0), zllmass_(0), zjjmass_(0);
      float jet1tkhe_(0), jet2tkhe_(0);
      float jet1jp_(0), jet2jp_(0);
      float jet1csv_(0), jet2csv_(0);
      float jet1csvmva_(0), jet2csvmva_(0);
      int jet1flav_(0), jet2flav_(0);

      int jet1Nch_(-1), jet2Nch_(-1);
      int jet1Nneu_(-1), jet2Nneu_(-1);
      float jet1ptD_(-1), jet2ptD_(-1);

      float zllpt_(-1.), drjj_(-1.);
      float qgLD_(-10);
      float hLD_(-10), cos1_(-10), cos2_(-10), cosStar_(-10), phi_(-10), phiStar_(-10);
      float hLD_rf(-10), cos1_rf(-10), cos2_rf(-10), cosStar_rf(-10), phi_rf(-10), phiStar_rf(-10);
      float higgsrefitMass(0), higgsMass(0);
      int higgsEvtNum(0), higgsRunNum(0), higgsLumiBlock(0);
      
      bool kineLepCut(false), fiducialCut(false), kineJetCut(false);
      bool dbLepCut(false);
      bool IDLep1Cut(false), IDLep2Cut(false), isoIDLepCut(false);

      // flag to dump event variables
      bool isEventSelected(false);

      table ElTable("sFGsfIdLoose.txt"); 
      //table MuTable("PfToId.txt"); 

      Weights elEff("Weights.root", "ElIDEff");
      // LOOP on CANDIDATES IN THE EVENT
      for(unsigned int j = 0; j < cands; ++j) {
	++cand[1];
	double ch;
	
	// std::cout<<"candidates loop"<<std::endl;
	  
	// fill histo of db of the muon candidates
	if(muChannel) {

	  //std::cout<<"muons db"<<std::endl;

	  db.Fill(get(muHiggsLeptDau1dB,j),evtWeight);
	  db.Fill(get(muHiggsLeptDau2dB,j),evtWeight);
	}

	// assign variables for selection
	if(muChannel) {
	  ch = 1.;

	  //std::cout<<"muons loop"<<std::endl;
	  
	  Jet1pt_=get(muHiggsJetDau1Pt,j);
	  Jet2pt_=get(muHiggsJetDau2Pt,j);
	  Jet1pt_rf=get(muHiggsJetDau1RefitPt,j);
	  Jet2pt_rf=get(muHiggsJetDau2RefitPt,j);
	  Jet1eta_=get(muHiggsJetDau1Eta,j);
	  Jet2eta_=get(muHiggsJetDau2Eta,j);
	  Jet1phi_=get(muHiggsJetDau1Phi,j);
	  Jet2phi_=get(muHiggsJetDau2Phi,j);
	  lept1pt_=get(muHiggsLeptDau1Pt,j);
	  lept2pt_=get(muHiggsLeptDau2Pt,j);
	  lept1eta_=get(muHiggsLeptDau1Eta,j);
	  lept2eta_=get(muHiggsLeptDau2Eta,j);
	  lept1phi_=get(muHiggsLeptDau1Phi,j);
	  lept2phi_=get(muHiggsLeptDau2Phi,j);

	  fiducialCut = true;
	  // db cut
	  dbLepCut = (fabs(get(muHiggsLeptDau1dB,j)) <0.02) && (fabs(get(muHiggsLeptDau2dB,j)) <0.02);


	  /* muonID cut */


	  bool MuID1(false);
	  bool MuID2(false);
	  LeptonID muID(PUPeriod);

	  /* 2012 Muon ID */
	  if((PUPeriod == "2012A") || (PUPeriod == "2012B") ||  (PUPeriod == "2012All") ){
	    MuID1 = muID.muonID2012(get(muHiggsLeptDau1GlobalMuonBit,j), get(muHiggsLeptDau1TrackerMuonBit,j), get(muHiggsMuDau1PFMuonBit,j), get(muHiggsLeptDau1NormChi2,j), get(muHiggsLeptDau1NofMuonHits,j), get(muHiggsLeptDau1NofMatchedStations,j), get(muHiggsMuDau1DzVtx,j), get(muHiggsLeptDau1NofPixelHits,j), get(muHiggsLeptDau1NofTrackerLayers,j), get(muHiggsLeptDau1dB,j) );
								      
	    MuID2 = muID.muonID2012(get(muHiggsLeptDau2GlobalMuonBit,j), get(muHiggsLeptDau2TrackerMuonBit,j), get(muHiggsMuDau2PFMuonBit,j), get(muHiggsLeptDau2NormChi2,j), get(muHiggsLeptDau2NofMuonHits,j), get(muHiggsLeptDau2NofMatchedStations,j), get(muHiggsMuDau2DzVtx,j),  get(muHiggsLeptDau2NofPixelHits,j), get(muHiggsLeptDau2NofTrackerLayers,j), get(muHiggsLeptDau2dB,j) );
	  }
	  
	  /* 2011 Muon ID */

	  else {
	    MuID1 = muID.muonID2011(get(muHiggsLeptDau1GlobalMuonBit,j),get(muHiggsLeptDau1TrackerMuonBit,j), get(muHiggsLeptDau1NormChi2,j), get(muHiggsLeptDau1NofTrackerHits,j) ,   get(muHiggsLeptDau1NofPixelHits,j), get(muHiggsLeptDau1NofMuonHits,j), get(muHiggsLeptDau1NofMatches,j), get(muHiggsLeptDau1dB,j)   );	   
	    MuID2 = muID.muonID2011(get(muHiggsLeptDau2GlobalMuonBit,j),get(muHiggsLeptDau2TrackerMuonBit,j), get(muHiggsLeptDau2NormChi2,j), get(muHiggsLeptDau2NofTrackerHits,j) ,   get(muHiggsLeptDau2NofPixelHits,j), get(muHiggsLeptDau2NofMuonHits,j), get(muHiggsLeptDau2NofMatches,j), get(muHiggsLeptDau2dB,j)   );

	  }


	  /* Muon PF Isolation rho*EA correction */
	  
	  //	  MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = MuonEffectiveArea::kMuEAFall11MC;
	  MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = MuonEffectiveArea::kMuEAData2012;
	  MuonEffectiveArea::MuonEffectiveAreaType effAreaType_   = MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
	  

	
	  
	  double EffectiveArea1 = 0.;
	  double EffectiveArea2 = 0.;
	  float abseta1=fabs(lept1eta_);
	  float abseta2=fabs(lept2eta_);
	  
	  EffectiveArea1 = MuonEffectiveArea::GetMuonEffectiveArea(effAreaType_, abseta1, effAreaTarget_);
	  EffectiveArea2 = MuonEffectiveArea::GetMuonEffectiveArea(effAreaType_, abseta2, effAreaTarget_);
	  

	  //std::cout<<"effective Area"<<std::endl;

	  lept1ChHadIso =  (get(muHiggsLeptDau1ChHadIso,j));
	  lept2ChHadIso =  (get(muHiggsLeptDau2ChHadIso,j));
	  lept1NeuHadIso =  (get(muHiggsLeptDau1NeuHadIso,j));
	  lept2NeuHadIso =  (get(muHiggsLeptDau2NeuHadIso,j));
	  lept1PhotonIso =  (get(muHiggsLeptDau1PhotonIso,j));
	  lept2PhotonIso =  (get(muHiggsLeptDau2PhotonIso,j));
	  lept1PUIso     =  (get(muHiggsLeptDau1PUChHadIso,j));
	  lept2PUIso     =  (get(muHiggsLeptDau2PUChHadIso,j));
	  //std::cout<<"PU iso muons: "<<lept1PUIso<<std::endl;

	  
	  iso1 = lept1ChHadIso;
	  iso1 += max<float>(0.,(lept1NeuHadIso + lept1PhotonIso) - EffectiveArea1*max<float>(0.0,getInt(rho)) );

	  iso2 = lept2ChHadIso;
	  iso2 += max<float>(0.,(lept2NeuHadIso + lept2PhotonIso) - EffectiveArea2*max<float>(0.0,getInt(rho)) );
	  

	  /******double beta correction********/
	  /*  iso1 = lept1ChHadIso;
	      iso1 += max<float>(0.,(lept1NeuHadIso + lept1PhotonIso - 0.5* lept1PUIso ) );
	      
	      iso2 = lept2ChHadIso;
	      iso2 += max<float>(0.,(lept2NeuHadIso + lept2PhotonIso - 0.5* lept2PUIso ) );
	  */
	  /***********************************/
	  iso1 = iso1/lept1pt_;
	  iso2 = iso2/lept2pt_;
	  
	  if(PUPeriod == "2011A"|| PUPeriod == "2011B" || PUPeriod == "2011All") isoIDLepCut = (get(muHiggsLeptDau1CombRelIso,j) < 0.15) && (get(muHiggsLeptDau2CombRelIso,j) < 0.15) && MuID1 && MuID2 ;

	  else isoIDLepCut = (iso1< 0.12) && (iso2 < 0.12) && MuID1 && MuID2 ;

	  //std::cout<<"muon ID performed: "<<lept1PUIso<<std::endl;

	  zllcharge_ = get(muHiggszllCharge,j);
	  zllmass_ = get(muHiggszllMass,j);
	  zjjmass_ = get(muHiggszjjMass,j); 
	  jet1tkhe_ = get(muHiggsJet1TKHE,j);
	  jet2tkhe_ = get(muHiggsJet2TKHE,j);
	  jet1jp_ = get(muHiggsJet1JProb,j);
	  jet2jp_ = get(muHiggsJet2JProb,j);
	  jet1csv_ = get(muHiggsJet1CSV,j);
	  jet2csv_ = get(muHiggsJet2CSV,j);
	  jet1csvmva_ = get(muHiggsJet1CSVMVA,j);
	  jet2csvmva_ = get(muHiggsJet2CSVMVA,j);
	  jet1flav_ = get(muHiggsJetDau1PartonFlavour,j);
	  jet2flav_ = get(muHiggsJetDau2PartonFlavour,j);

	  beta_ = get(muHiggsJetDau1puBeta, j);
	  beta2_ = get(muHiggsJetDau2puBeta, j);

	  jet1Nch_ = int( get(muHiggsJetDau1ChHadMult,j) );
	  jet2Nch_ = int( get(muHiggsJetDau2ChHadMult,j) );
	  jet1Nneu_ = int( get(muHiggsJetDau1NeuHadMult,j) ) + int( get(muHiggsJetDau1PhotMult,j) );
	  jet2Nneu_ = int( get(muHiggsJetDau2NeuHadMult,j) ) + int( get(muHiggsJetDau2PhotMult,j) );
	  jet1ptD_ = get(muHiggsJetDau1PtDJet,j); 
	  jet2ptD_ = get(muHiggsJetDau2PtDJet,j);

	  zllpt_ = get(muHiggszllPt,j);
	  drjj_ = get(muHiggsjjdr,j);
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
	  higgsEvtNum = getInt(muHiggsEventNumber);
	  higgsRunNum = getInt(muHiggsRunNumber);
	  higgsLumiBlock = getInt(muHiggsLumiblock);

	  //std::cout<<"variables gotten: "<<lept1PUIso<<std::endl;

	} else {

	  // std::cout<<"electron channel"<<std::endl;
      
	  ch = 0.;
	  Jet1pt_=get(elHiggsJetDau1Pt,j);
	  Jet2pt_=get(elHiggsJetDau2Pt,j);
	  Jet1pt_rf=get(elHiggsJetDau1RefitPt,j);
	  Jet2pt_rf=get(elHiggsJetDau2RefitPt,j);
	  Jet1eta_=get(elHiggsJetDau1Eta,j);
	  Jet2eta_=get(elHiggsJetDau2Eta,j);
	  Jet1phi_=get(elHiggsJetDau1Phi,j);
	  Jet2phi_=get(elHiggsJetDau2Phi,j);
	  lept1pt_=get(elHiggsLeptDau1Pt,j);
	  lept2pt_=get(elHiggsLeptDau2Pt,j);
	  lept1eta_=get(elHiggsLeptDau1Eta,j);
	  lept2eta_=get(elHiggsLeptDau2Eta,j);
	  lept1phi_=get(elHiggsLeptDau1Phi,j);
	  lept2phi_=get(elHiggsLeptDau2Phi,j);

	  fiducialCut = (fabs(get(elHiggsLeptDau1Eta,j)) < 2.5 ) &&
	    (fabs(get(elHiggsLeptDau2Eta,j)) < 2.5 ) &&
	    !(  (fabs(get(elHiggsLeptDau1EtaSC,j)) < 1.566) &&
		(fabs(get(elHiggsLeptDau1EtaSC,j)) > 1.4442) ) &&  
	    !(  (fabs(get(elHiggsLeptDau2EtaSC,j)) < 1.566) &&
		(fabs(get(elHiggsLeptDau2EtaSC,j)) > 1.4442) );

	  dbLepCut = true;
	  //	  isoIDLepCut =(get(elHiggsEleDau1VBTF80CombID,j) ==7) || (get(elHiggsEleDau2VBTF80CombID,j) ==7);
	  // new electronID


	  bool EleDau1VBTFID = get(elHiggsEleDau1VBTF95CombID,j) == 7;
	  bool isEle1Barrel = fabs(get(elHiggsLeptDau1EtaSC,j)) <= 1.4442;
	  bool EleDau2VBTFID = get(elHiggsEleDau2VBTF95CombID,j) == 7;
	  bool isEle2Barrel = fabs(get(elHiggsLeptDau2EtaSC,j)) <= 1.4442;

	  // bool EleDau1TightCuts = EleIDTightCuts( isEle1Barrel, get(elHiggsLeptDau1DeltaPhiAtVtx,j), 
	  // 					  get(elHiggsLeptDau1HOverE,j) );
	  // bool EleDau2TightCuts = EleIDTightCuts( isEle2Barrel, get(elHiggsLeptDau2DeltaPhiAtVtx,j), 
	  // 					  get(elHiggsLeptDau2HOverE,j) );



	  //	  if(getInt(elHiggsEventNumber) == 266333 || getInt(elHiggsEventNumber) == 140798 || 
	  //	     getInt(elHiggsEventNumber) == 201809 || getInt(elHiggsEventNumber) == 202476 ) {
	  /*
	  if(getInt(elHiggsEventNumber) == 270594) {
	  cout << "-------------------------------------------------------------" << endl;
	  cout << "EVENT = " << getInt(elHiggsEventNumber) << endl;
	  cout << "CANDIDATE " << j << " of " << cands-1 << endl;
	  cout << "Lept1 eta " << fabs(get(elHiggsLeptDau1Eta,j)) << endl; 
	  cout << "Lept2 eta " << fabs(get(elHiggsLeptDau2Eta,j)) << endl; 
	  cout << "Lept1 etaSC " << fabs(get(elHiggsLeptDau1EtaSC,j)) << endl; 
	  cout << "Lept2 etaSC " << fabs(get(elHiggsLeptDau2EtaSC,j)) << endl; 
	  cout << "EleDau1VBTFID95 " << get(elHiggsEleDau1VBTF95CombID,j) << endl;
	  cout << "EleDau2VBTFID95 " << get(elHiggsEleDau2VBTF95CombID,j) << endl;
	  cout << "EleDau1DPhiSC " << get(elHiggsLeptDau1DeltaPhiAtVtx,j) << endl;
	  cout << "EleDau2DPhiSC " << get(elHiggsLeptDau2DeltaPhiAtVtx,j) << endl;
	  cout << "EleDau1H/E " << get(elHiggsLeptDau1HOverE,j) << endl;
	  cout << "EleDau2H/E " << get(elHiggsLeptDau2HOverE,j) << endl;
	  cout << "isEle1Barrel " << isEle1Barrel << endl;
	  cout << "isEle2Barrel " << isEle2Barrel << endl;
	  cout << "isEle1Tight " << EleIDTightCuts( isEle1Barrel, get(elHiggsLeptDau1DeltaPhiAtVtx,j),get(elHiggsLeptDau1HOverE,j) ) << endl;
	  cout << "isEle2Tight " << EleIDTightCuts( isEle2Barrel, get(elHiggsLeptDau2DeltaPhiAtVtx,j),get(elHiggsLeptDau2HOverE,j) ) << endl;
	}
	  */

	  //ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAFall11MC;
	  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
	  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaTypeGammaNeuHad_   = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;
  
	  float EffectiveAreaGammaNeuHad1 = 0.;
	  float EffectiveAreaGammaNeuHad2 = 0.;
	  float abseta1=fabs(get(elHiggsLeptDau1EtaSC,j));
	  float abseta2=fabs(get(elHiggsLeptDau2EtaSC,j));

	
	  EffectiveAreaGammaNeuHad1 = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaTypeGammaNeuHad_, abseta1, effAreaTarget_);
	  EffectiveAreaGammaNeuHad2 = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaTypeGammaNeuHad_, abseta2, effAreaTarget_);
	  
	  lept1ChHadIso =  (get(elHiggsLeptDau1ChHadIso,j));
	  lept2ChHadIso =  (get(elHiggsLeptDau2ChHadIso,j));
	  lept1NeuHadIso =  (get(elHiggsLeptDau1NeuHadIso,j));
	  lept2NeuHadIso =  (get(elHiggsLeptDau2NeuHadIso,j));
	  lept1PhotonIso =  (get(elHiggsLeptDau1PhotonIso,j));
	  lept2PhotonIso =  (get(elHiggsLeptDau2PhotonIso,j));
	  lept1PUIso =  (get(elHiggsLeptDau1PUChHadIso,j));
	  lept2PUIso =  (get(elHiggsLeptDau2PUChHadIso,j));

	  //  //std::cout<<"PU iso electrons: "<<lept1PUIso<<std::endl;

	  
	  iso1 = lept1ChHadIso;
	  iso1 += max<float>(0.,(lept1NeuHadIso + lept1PhotonIso) - EffectiveAreaGammaNeuHad1*getInt(rhoRestrictedEta) );

	  iso2 = lept2ChHadIso;
	  iso2 += max<float>(0.,(lept2NeuHadIso + lept2PhotonIso) - EffectiveAreaGammaNeuHad2*getInt(rhoRestrictedEta) );
	  

	  //std::cout<<"isolation"<<std::endl;

	  float combRelIso1 = iso1/ lept1pt_;
	  float combRelIso2 = iso2/ lept2pt_;

	  bool EleID1 = 0;
	  bool EleID2 = 0;
	  //float combRelIso1 = ( ( get(elHiggsEleDau1TrkIso03 ,j) +  get(elHiggsEleDau1EcalIso03 ,j) +  get(elHiggsEleDau1HcalIso03 ,j)) /  lept1pt_);
	  //float combRelIso2 = ( get(elHiggsEleDau2TrkIso03 ,j) +  get(elHiggsEleDau2EcalIso03 ,j) +  get(elHiggsEleDau2HcalIso03 ,j)) /  lept2pt_;
	  // float a = fabs( (1/get( elHiggsLeptDau1EcalEn, j)) - (get(elHiggsLeptDau1EPin,j)/get( elHiggsLeptDau1EcalEn, j))   )
	  
	  bool isEB  = get(elHiggsLeptDau1isEB, j) ? true : false;
	  
	  LeptonID eleID(PUPeriod);

	  if((PUPeriod == "2012A") || (PUPeriod == "2012B") ||  (PUPeriod == "2012All") ){	      
	    EleID1 = eleID.electronID2012(get(elHiggsEleDau1dxy, j), get(elHiggsEleDau1dz, j), get(elHiggsEleDau1TrkIso03, j), get(elHiggsEleDau1EcalIso03, j), get(elHiggsEleDau1HcalIso03, j), lept1pt_, combRelIso1, get(elHiggsLeptDau1DeltaEtaAtVtx, j), get(elHiggsLeptDau1DeltaPhiAtVtx, j), get(elHiggsLeptDau1Sigmaee, j) , get(elHiggsLeptDau1HOverE, j), get(elHiggsLeptDau1EcalEn, j) , get(elHiggsLeptDau1EPin, j) , get(elHiggsEleDau1mHits, j), get(elHiggsEleDau1hasMatchConv, j), get(elHiggsLeptDau1isEB, j)  );
	     
	    EleID2 = eleID.electronID2012(get(elHiggsEleDau2dxy, j), get(elHiggsEleDau2dz, j), get(elHiggsEleDau2TrkIso03, j), get(elHiggsEleDau2EcalIso03, j), get(elHiggsEleDau2HcalIso03, j), lept2pt_, combRelIso2, get(elHiggsLeptDau2DeltaEtaAtVtx, j), get(elHiggsLeptDau2DeltaPhiAtVtx, j), get(elHiggsLeptDau2Sigmaee, j) , get(elHiggsLeptDau2HOverE, j),  get(elHiggsLeptDau2EcalEn, j) , get(elHiggsLeptDau2EPin, j) , get(elHiggsEleDau2mHits, j), get(elHiggsEleDau2hasMatchConv, j), get(elHiggsLeptDau2isEB, j)  );
	      
	     }
	  else{
	  
	  EleID1 = eleID.electronID2011(get(elHiggsLeptDau1DeltaPhiAtVtx, j), get(elHiggsLeptDau1HOverE, j), get(elHiggsEleDau1VBTF95CombID,j), get(elHiggsLeptDau1EtaSC,j) , lept1pt_, combRelIso1);
	  EleID2 = eleID.electronID2011(get(elHiggsLeptDau2DeltaPhiAtVtx, j), get(elHiggsLeptDau2HOverE, j), get(elHiggsEleDau2VBTF95CombID,j), get(elHiggsLeptDau2EtaSC,j) , lept2pt_, combRelIso2);
	   }
	  isoIDLepCut =  EleID1 && EleID2;

	  //std::cout<<"id boolean: "<<isoIDLepCut<<std::endl;

	  //IDweight_1 = ElTable.Val(lept1pt_, lept1eta_);
	  IDweight_1 = elEff.getEff( lept1eta_, lept1pt_);
	  //std::cout<<"El eff1: "<<IDweight_1<<std::endl;
	  //IDweight_2 = ElTable.Val(lept2pt_, lept2eta_);
	  IDweight_2 = elEff.getEff( lept2eta_, lept2pt_);

	  //std::cout<<"El eff2: "<<IDweight_2<<std::endl;



	  zllcharge_ = get(elHiggszllCharge,j);
	  zllmass_ = get(elHiggszllMass,j);
	  zjjmass_ = get(elHiggszjjMass,j); 
	  jet1tkhe_ = get(elHiggsJet1TKHE,j);
	  jet2tkhe_ = get(elHiggsJet2TKHE,j);
	  jet1jp_ = get(elHiggsJet1JProb,j);
	  jet2jp_ = get(elHiggsJet2JProb,j);
	  jet1csv_ = get(elHiggsJet1CSV,j);
	  jet2csv_ = get(elHiggsJet2CSV,j);
	  jet1csvmva_ = get(elHiggsJet1CSVMVA,j);
	  jet2csvmva_ = get(elHiggsJet2CSVMVA,j);
	  jet1flav_ = get(elHiggsJetDau1PartonFlavour,j);
	  jet2flav_ = get(elHiggsJetDau2PartonFlavour,j);

	  beta_ = get(elHiggsJetDau1puBeta, j);
	  beta2_ = get(elHiggsJetDau2puBeta, j);

	  jet1Nch_ = int( get(elHiggsJetDau1ChHadMult,j) );
	  jet2Nch_ = int( get(elHiggsJetDau2ChHadMult,j) );
	  jet1Nneu_ = int( get(elHiggsJetDau1NeuHadMult,j) ) + int( get(elHiggsJetDau1PhotMult,j) );
	  jet2Nneu_ = int( get(elHiggsJetDau2NeuHadMult,j) ) + int( get(elHiggsJetDau2PhotMult,j) );
	  jet1ptD_ = get(elHiggsJetDau1PtDJet,j); 
	  jet2ptD_ = get(elHiggsJetDau2PtDJet,j);

	  zllpt_ = get(elHiggszllPt,j);
	  drjj_ = get(elHiggsjjdr,j);
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
	  higgsEvtNum = getInt(elHiggsEventNumber);
	  higgsRunNum = getInt(elHiggsRunNumber);
	  higgsLumiBlock = getInt(elHiggsLumiblock);
	  
	}



	kineLepCut = (lept1pt_ > 40 && lept2pt_ > 20) || ( lept1pt_ > 20 &&  lept2pt_ > 40 );
	kineJetCut = (Jet1pt_ > 30) && (Jet2pt_ > 30);

	//std::cout<<"kinematic definition: "<<std::endl;

	bool zchargecut = (zllcharge_ == 0) ;

	//elMVAID.Fill(get(elHiggsEleDau1mvaTrigV0,j),evtWeight);


	float DRlep1_Jet1 = deltaR(Jet1eta_, Jet1phi_, lept1eta_, lept1phi_);
	float DRlep2_Jet1 = deltaR(Jet1eta_, Jet1phi_, lept2eta_, lept2phi_);
	float DRlep1_Jet2 = deltaR(Jet2eta_, Jet2phi_, lept1eta_, lept1phi_);
	float DRlep2_Jet2 = deltaR(Jet2eta_, Jet2phi_, lept2eta_, lept2phi_);

	bool unclean = ( (DRlep1_Jet1 < 0.5)  || (DRlep2_Jet1 < 0.5) || (DRlep1_Jet2 < 0.5) ||  (DRlep2_Jet2 < 0.5));


	//BTagSFUtil* btsfutil = new BTagSFUtil(TMath::Sin(Jet1phi_*1000000));
	//BTagSFUtil* btsfutil2 = new BTagSFUtil(TMath::Sin(Jet2phi_*1000000));
	if(!data) IDweight = IDweight_1 * IDweight_2;
	//	evtWeight = PUWeight*IDweight;
	evtWeight *= IDweight;


	if(kineLepCut && fiducialCut && dbLepCut && isoIDLepCut && kineJetCut &&  zchargecut && !unclean) {
	  //if(kineLepCut  && isoIDLepCut   ) {


	  //std::cout<<"kinematic selection: "<<std::endl;


	  //	if(1==1) {
	  pass[2] = true;
	  ++cand[2];

	  // fill plots after lepton selection
	  if (!hasAlreadyCand_preZllcut){
	    zllmass.Fill(zllmass_,evtWeight);
	    if(abs(lept1eta_)< 0.8 && abs(lept2eta_)< 0.8)zllmass_BB.Fill(zllmass_,evtWeight);
	    hasAlreadyCand_preZllcut = true;
	    }
	  npv1.Fill(getInt(numPV),evtWeight);
	  njets_ = nj;
	  // Z inv. mass selection	  
	  bool zllcut = (zllmass_ > 70) && (zllmass_ < 110);
	  // bool zllcut = (zllmass_ > 60) && (zllmass_ < 120);
	  // Z charge selection	  

	  // zjj cut different for signal vs. sideband
	  bool zjjcut_sigSel = (zjjmass_ > 75) && (zjjmass_ < 105);
	  //	  bool zjjcut_sbSel = zjjmass_ > 120;
	  bool zjjcut_sbSel =(  ((zjjmass_ > 60) &&  (zjjmass_ < 75))  || (  (zjjmass_ > 105) &&  (zjjmass_ < 130))  );


	  float QGLikeJet1 = qglikeli->computeQGLikelihoodPU( Jet1pt_, getInt(rhoRestrictedEta), jet1Nch_, jet1Nneu_, jet1ptD_ );
	  float QGLikeJet2 = qglikeli->computeQGLikelihoodPU( Jet2pt_, getInt(rhoRestrictedEta), jet2Nch_, jet2Nneu_, jet2ptD_ );
	  qgLD_ = QGLikeJet1 * QGLikeJet1;


	  float ResidZZMass_ = (fabs(zllmass_-91.187) + fabs(zjjmass_-91.187));  


	  //SIGNAL SELECTION
	  if( zllcut) {
	    if (!hasAlreadyCand){
	      runNumSkim.push_back(higgsRunNum);
	      stdCandle.Fill(zllmass_,evtWeight);	   
	      zjjmass.Fill(zjjmass_,evtWeight);	   	  
	      jetpt1.Fill(Jet1pt_,evtWeight);
	      jetpt2.Fill(Jet2pt_,evtWeight);
	      leptpt1.Fill(lept1pt_,evtWeight);
	      leptpt2.Fill(lept2pt_,evtWeight);
	      lepteta1.Fill( lept1eta_, evtWeight);
	      lepteta2.Fill( lept2eta_, evtWeight);
	      leptphi1.Fill( lept1phi_, evtWeight);
	      leptphi2.Fill( lept2phi_, evtWeight);
	      eta_pt.Fill(  lept1eta_, lept1pt_, evtWeight);
	      eta_zllmass.Fill(lept1eta_, zllmass_, evtWeight);
	      eta_phi.Fill( lept1eta_, lept1phi_, evtWeight);
	      pt_phi.Fill(lept1phi_, lept1pt_, evtWeight);
  
	      beta.Fill(beta_, evtWeight);

	      qgLD.Fill(qgLD_, evtWeight);
	      helyLD_RF.Fill(hLD_rf,evtWeight);
	      cosT1_RF.Fill(cos1_rf,evtWeight);
	      cosT2_RF.Fill(cos2_rf,evtWeight);
	      cosT1Star_RF.Fill(cosStar_rf,evtWeight);
	      phi_RF.Fill(phi_rf,evtWeight);
	      phiStar_RF.Fill(phiStar_rf,evtWeight);

	      btagJet1_JP.Fill(jet1jp_,evtWeight);
	      btagJet2_JP.Fill(jet2jp_,evtWeight);
	      btagJet_JP.Fill(jet1jp_,evtWeight);
	      btagJet_JP.Fill(jet2jp_,evtWeight);
	      btagJet1.Fill(jet1tkhe_,evtWeight);
	      btagJet2.Fill(jet2tkhe_,evtWeight);
	      btagJet.Fill(jet1tkhe_,evtWeight);
	      btagJet.Fill(jet2tkhe_,evtWeight);
	      btagJet1_CSV.Fill(jet1csv_,evtWeight);
	      btagJet2_CSV.Fill(jet2csv_,evtWeight);
	      btagJet_CSV.Fill(jet1csv_,evtWeight);
	      btagJet_CSV.Fill(jet2csv_,evtWeight);
	      btagJet1_CSVMVA.Fill(jet1csvmva_,evtWeight);
	      btagJet2_CSVMVA.Fill(jet2csvmva_,evtWeight);
	      btagJet_CSVMVA.Fill(jet1csvmva_,evtWeight);
	      btagJet_CSVMVA.Fill(jet2csvmva_,evtWeight);
	      
	      h_met.Fill(getInt(met),evtWeight);
	      metSig.Fill(getInt(metSignificance),evtWeight);
	      metSignif.Fill(getInt(metSignificance),evtWeight);
	      
	      npv2.Fill(getInt(numPV),evtWeight);
	      hasAlreadyCand = true;
	    }

	    if(zjjcut_sigSel){
	      if(!hasAlreadyCand_afterzjj){
		hasAlreadyCand_afterzjj = true;
		
		helyLD_RF_zjj.Fill(hLD_rf,evtWeight);
		cosT1_RF_zjj.Fill(cos1_rf,evtWeight);
		cosT2_RF_zjj.Fill(cos2_rf,evtWeight);
		cosT1Star_RF_zjj.Fill(cosStar_rf,evtWeight);
		phi_RF_zjj.Fill(phi_rf,evtWeight);
		phiStar_RF_zjj.Fill(phiStar_rf,evtWeight);

	      }

	    }	    
	    if( zjjcut_sigSel && beta_>=0.2 && beta2_ >=0.2) {	      
	      pass[3] = true;
	      ++cand[3];
	      
	      // fill plots after Z mass selection
	   
	  
	      // end of common selection
	      // now jet classification according btag
	      
	      //Tag algorithm cut definition for signal
	      bool jet1_tagged_medium(false), jet1_tagged_loose(false);
	      bool jet2_tagged_medium(false), jet2_tagged_loose(false);

	      jet1_tagged_medium = jet1jp_ > 0.545;
	      jet1_tagged_loose  = jet1jp_ > 0.275;	    
	      jet2_tagged_medium = jet2jp_ > 0.545;
	      jet2_tagged_loose  = jet2jp_ > 0.275;
      
 	      // jet1_tagged_medium = jet1tkhe_ > 3.3;
 	      // jet1_tagged_loose  = jet1tkhe_ > 1.7;	    
 	      // jet2_tagged_medium = jet2tkhe_ > 3.3;
 	      // jet2_tagged_loose  = jet2tkhe_ > 1.7;

// 	      jet1_tagged_medium = jet1csv_ > 0.679;
// 	      jet1_tagged_loose  = jet1csv_ > 0.244;	    
// 	      jet2_tagged_medium = jet2csv_ > 0.679;
// 	      jet2_tagged_loose  = jet2csv_ > 0.244;


	      // eventually apply SF for MC



	      if( btagScale ) {
		btsfutil->modifyBTagsWithSF( "JP", jet1_tagged_loose, jet1_tagged_medium, Jet1pt_, Jet1eta_, jet1flav_ );
		btsfutil->modifyBTagsWithSF( "JP", jet2_tagged_loose, jet2_tagged_medium, Jet2pt_, Jet2eta_, jet2flav_ );
	      }
	      
// 	      bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_medium  );
	      bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_loose  )
	      	|| (jet1_tagged_loose  && jet2_tagged_medium );
// 	      bool twoTagCategory  = ( jet1_tagged_loose && jet2_tagged_loose  );
	      bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_loose || jet2_tagged_loose );
// 	      bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_medium || jet2_tagged_medium );
	      bool zeroTagCategory = ( !twoTagCategory ) && ( !oneTagCategory );
	      



	      // count the candidates per category	    
	      if(zeroTagCategory) 
		++CandBTag[0];
	      if(oneTagCategory) 
		++CandBTag[1];
	      if(twoTagCategory) 
		++CandBTag[2];
	      
	      //std::cout<<"sf applied: "<<std::endl;

	      // select only 2-tag candidates

	   

	      if(twoTagCategory){
		pass[4] = true;
		++cand[4];
		// fill plots for 2-tag candidates selection
		zllpt.Fill(zllpt_,evtWeight);
		drjj.Fill(drjj_,evtWeight);
	
		lljjmass2tags.Fill(higgsrefitMass,evtWeight);
		helyLD_RF2tags.Fill(hLD_rf,evtWeight);
		
		// apply selection for 2-tag candidates
		// met significance cut
		if( getInt(metSignificance) < 10) {
		  pass[5] = true;
		  ++cand[5];
		  zllptBeforeLD.Fill(zllpt_,evtWeight);
		  drjjBefLD.Fill(drjj_,evtWeight);
		  lljjmassBeforeLD.Fill(higgsrefitMass,evtWeight);
		  
		  // check if it is the best candidate
		  //		if ( ((Jet1pt_rf+ Jet2pt_rf)>SumPtBest_) && (higgsrefitMass>0) ){
		  if (  (ResidZZMass_ < ResidZllMassBest_) && (higgsrefitMass>0) ){
		    // if best candidate store its variables
		    //		  SumPtBest_= Jet1pt_rf+ Jet2pt_rf;
		    ResidZllMassBest_=  ResidZZMass_;
		    cos1_Best=cos1_;
		    cos2_Best=cos2_;
		    cosStar_Best= cosStar_;
		    phi_Best= phi_;
		    phiStar_Best= phiStar_;
		    hLD_Best= hLD_;
		    
		    cos1_rfBest=cos1_rf;
		    cos2_rfBest=cos2_rf;
		    cosStar_rfBest= cosStar_rf;
		    phi_rfBest= phi_rf;
		    phiStar_rfBest= phiStar_rf;
		    hLD_rfBest= hLD_rf;
		    
		  }
		  
		  // LD cut
		  if( hLD_rf > 0.5 ){
		    pass[6] = true;
		    ++cand[6];

		    // dump event variables
		    // check if it is the best candidate at this level
		    //		  if ( ( (Jet1pt_rf+ Jet2pt_rf)>SumPtBest_Final) && (higgsrefitMass>0) ){
		    if ( ( ResidZZMass_ < ResidZllMassBest_Final) && (higgsrefitMass>0) ){
		      // if best candidate store lljj mass
		      //		    SumPtBest_Final= Jet1pt_rf+ Jet2pt_rf;
		      exist2tag = true;
		      puw2tag = evtWeight;
		      ch2tag = ch;
		      //cout<<"lept channel "<<ch2tag<<endl;
		      ResidZllMassBest_Final=  ResidZZMass_;
		      lljjmass_= higgsrefitMass;
		      lljjmass_noRefit_= higgsMass;
		      EvtNum2Tag_=higgsEvtNum;
		      RunNum2Tag_=higgsRunNum;
		      LumiBlock1Tag_=higgsLumiBlock; 
		      ld_2btag_ = hLD_rf;
		      qgprod_2btag_ = -1;
		      metsig_2btag_ = getInt(metSignificance);
		      jet1pt_2btag_ = Jet1pt_rf;
		      jet2pt_2btag_ = Jet2pt_rf;
		      lep1pt_2btag_ = lept1pt_;
		      lep2pt_2btag_ = lept2pt_;
		      zllmass_2btag_ = zllmass_;
		      zjjmass_2btag_ = zjjmass_;
		    }
		    
		  }
		}
	      }//end 2-tag category

	      
	      // select only 1-tag candidates
	      if(oneTagCategory){
	
		// apply selection for 1-tag candidates
		// LD cut
		if( hLD_rf >(0.302+(0.000656*higgsrefitMass)) ){
		  // check if it is the best candidate at this level
		  // ( ( (Jet1pt_rf+ Jet2pt_rf)>SumPtBest_Final_1btag) && (higgsrefitMass>0) ){
		  if ( ( ResidZZMass_ < ResidZllMassBest1tag_) && (higgsrefitMass>0) ){
		    // if best candidate store lljj mass
		    exist1tag = true;
		    puw1tag = evtWeight;
		      
		    ch1tag = ch;
		    ResidZllMassBest1tag_ =  ResidZZMass_;
		    lljjmass_1btag_ = higgsrefitMass;
		    lljjmass_noRefit_1btag_ = higgsMass;
		    EvtNum1Tag_=higgsEvtNum;
		    RunNum1Tag_=higgsRunNum;
		    LumiBlock1Tag_=higgsLumiBlock; 
		    ld_1btag_ = hLD_rf;
		    qgprod_1btag_ = -1;
		    metsig_1btag_ = -1;
		    jet1pt_1btag_ = Jet1pt_rf;
		    jet2pt_1btag_ = Jet2pt_rf;
		    lep1pt_1btag_ = lept1pt_;
		    lep2pt_1btag_ = lept2pt_;
		    zllmass_1btag_ = zllmass_;
		    zjjmass_1btag_ = zjjmass_;
		  }
		}
	      } // end 1-tag category	      
	      
	      // select only 0-tag candidates
	      if(zeroTagCategory){
		// apply selection for 0-tag candidates

		// QG cut
		//		float QGLikeJet1 = qglikeli->computeQGLikelihoodPU( Jet1pt_, getInt(rhoRestrictedEta), jet1Nch_, jet1Nneu_, jet1ptD_ );
		//		float QGLikeJet2 = qglikeli->computeQGLikelihoodPU( Jet2pt_, getInt(rhoRestrictedEta), jet2Nch_, jet2Nneu_, jet2ptD_ );

		
		//		if( QGLikeJet1 * QGLikeJet2 > 0.10) {
		if(1) {
		  // LD cut
		  
	
		  if( hLD_rf >(0.55+(0.00025*higgsrefitMass)) ){
		    // check if it is the best candidate at this level
		    if ( ( ResidZZMass_< ResidZllMassBest0tag_) && (higgsrefitMass>0) ){
		      exist0tag = true;
		      puw0tag = evtWeight;
		      ch0tag = ch;
		      ResidZllMassBest0tag_ =  ResidZZMass_;
		      lljjmass_0btag_ = higgsrefitMass;
		      lljjmass_noRefit_0btag_ = higgsMass;
		      EvtNum0Tag_=higgsEvtNum;
		      RunNum0Tag_=higgsRunNum;
		      LumiBlock0Tag_=higgsLumiBlock; 
		      ld_0btag_ = hLD_rf;
		      qgprod_0btag_ = QGLikeJet1*QGLikeJet2;
		      metsig_0btag_ = -1;
		      jet1pt_0btag_ = Jet1pt_rf;
		      jet2pt_0btag_ = Jet2pt_rf;
		      lep1pt_0btag_ = lept1pt_;
		      lep2pt_0btag_ = lept2pt_;
		      zllmass_0btag_ = zllmass_;
		      zjjmass_0btag_ = zjjmass_;
		    }
		  } // end LD cut
		} // end QG cut
	      } // end 0-tag category
	    


	    }//end zjj cut

	  } // end of signal selection (zll cut)


	  // SIDEBAND SELECTION
	  if( zllcut && zjjcut_sbSel && beta_>0.2 && beta2_ >0.2 ) {
	    
	    // all candidates in the sideband
	    hmass_sb_all.Fill(higgsrefitMass);
	    hmassnorefit_sb_all.Fill(higgsMass);
	    //

	    //Tag algorithm cut definition for sideband
	    bool jet1_tagged_medium(false), jet1_tagged_loose(false);
	    bool jet2_tagged_medium(false), jet2_tagged_loose(false);
 	    jet1_tagged_medium = jet1jp_ > 0.545;
 	    jet1_tagged_loose  = jet1jp_ > 0.275;	    
 	    jet2_tagged_medium = jet2jp_ > 0.545;
 	    jet2_tagged_loose  = jet2jp_ > 0.275;

 	    // jet1_tagged_medium = jet1tkhe_ > 3.3;
 	    // jet1_tagged_loose  = jet1tkhe_ > 1.7;	    
 	    // jet2_tagged_medium = jet2tkhe_ > 3.3;
 	    // jet2_tagged_loose  = jet2tkhe_ > 1.7;

// 	    jet1_tagged_medium = jet1csv_ > 0.679;
// 	    jet1_tagged_loose  = jet1csv_ > 0.244;	    
// 	    jet2_tagged_medium = jet2csv_ > 0.679;
// 	    jet2_tagged_loose  = jet2csv_ > 0.244;


	    if( btagScale ) {
	      btsfutil->modifyBTagsWithSF( "JP", jet1_tagged_loose, jet1_tagged_medium, Jet1pt_, Jet1eta_, jet1flav_ );
	      btsfutil->modifyBTagsWithSF( "JP", jet2_tagged_loose, jet2_tagged_medium, Jet2pt_, Jet2eta_, jet2flav_ );
	    }
	    
//  	    bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_medium  );
	    bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_loose  )
	      || (jet1_tagged_loose  && jet2_tagged_medium );
	    // 	      bool twoTagCategory  = ( jet1_tagged_loose && jet2_tagged_loose  );
	    bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_loose || jet2_tagged_loose );
	    // 	      bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_medium || jet2_tagged_medium );
	    bool zeroTagCategory = ( !twoTagCategory ) && ( !oneTagCategory );
	      

	    if(twoTagCategory){
	      // candidates in the sideband after btagging
	      hmass_sb_2b.Fill(higgsrefitMass);
	      hmassnorefit_sb_2b.Fill(higgsMass);
	      //
	      if( getInt(metSignificance) < 10) {
		// best candidate in the sideband after metSig
		if ( ((Jet1pt_rf+ Jet2pt_rf)>SumPtBest_SB) && (higgsrefitMass>0) ){
		  SumPtBest_SB= Jet1pt_rf+ Jet2pt_rf;
		  lljjmass_SB= higgsrefitMass;
		  lljjmass_noRefit_SB= higgsMass;
		}		
		if( hLD_rf > 0.5 ){
		  // best candidate in the sideband after hLD
		  if ( ( ResidZZMass_ < ResidZllMassBestSB2tag_) && (higgsrefitMass>0) ){
		    //		  if ( ( (Jet1pt_rf+ Jet2pt_rf)>SumPtBest_FinalSB) && (higgsrefitMass>0) ){
		    // SumPtBest_FinalSB= Jet1pt_rf+ Jet2pt_rf;
		    exist2tagSB = true;
		    puw2tagSB = evtWeight;
		    ch2tagSB = ch;
		    ResidZllMassBestSB2tag_ =  ResidZZMass_; 
		    lljjmassFinal_SB = higgsrefitMass;
		    lljjmassFinal_noRefit_SB = higgsMass;
		    zjjmass_SB2btag_ = zjjmass_;
		    EvtNum2TagSB_=higgsEvtNum;
		    RunNum2TagSB_=higgsRunNum;
		    LumiBlock2TagSB_=higgsLumiBlock; 
		  }
		}
	      }
	      } // end of 2btag sideband selection

	    // select only 1-tag candidates
	      if(oneTagCategory){
			// apply selection for 1-tag candidates
		// LD cut
		if( hLD_rf >(0.302+(0.000656*higgsrefitMass)) ){
		  // check if it is the best candidate at this level
		  // ( ( (Jet1pt_rf+ Jet2pt_rf)>SumPtBest_Final_1btag) && (higgsrefitMass>0) ){
		  if ( ( ResidZZMass_< ResidZllMassBestSB1tag_) && (higgsrefitMass>0) ){
		    // if best candidate store lljj mass
		    exist1tagSB = true;
		    puw1tagSB = evtWeight;
		    ch1tagSB = ch;
		    ResidZllMassBestSB1tag_ =  ResidZZMass_;
		    lljjmassFinal_SB_1btag_= higgsrefitMass;
		    zjjmass_SB1btag_ = zjjmass_;
		    EvtNum1TagSB_=higgsEvtNum;
		    RunNum1TagSB_=higgsRunNum;
		    LumiBlock1TagSB_=higgsLumiBlock; 	
		  }
		}
	      } // end 1-tag category	      
	        
	    if(zeroTagCategory){
		// apply selection for 0-tag candidates
		
	      if(1){	
		//if( QGLikeJet1 * QGLikeJet2 > 0.10) {
		  // LD cut
	 
		  if( hLD_rf >(0.55+(0.00025*higgsrefitMass)) ){
		    // check if it is the best candidate at this level
		    if ( ( ResidZZMass_ < ResidZllMassBestSB0tag_) && (higgsrefitMass>0) ){
		      exist0tagSB = true;
		      puw0tagSB = evtWeight;
		      ch0tagSB = ch;
		      ResidZllMassBestSB0tag_ =  ResidZZMass_;
		      lljjmassFinal_SB_0btag_ = higgsrefitMass;
		      zjjmass_SB0btag_ = zjjmass_;
		      EvtNum0TagSB_=higgsEvtNum;
		      RunNum0TagSB_=higgsRunNum;
		      LumiBlock0TagSB_=higgsLumiBlock; 
		    }
		  } // end LD cut
		} // end QG cut
	      } // end 0-tag category in sideband


	  } // end of sideband selection


	}
      } // candidate loop
      
      //std::cout<<"writeText: "<<writeText<<std::endl;
      // fill best candidate plots
      if((lljjmass_!=0) && exist2tag) { 
	lljjmass.Fill(lljjmass_, evtWeight);
	// store the best cand mass
	bestHiggsMass.push_back(lljjmass_);
	bestmjj.push_back(zjjmass_2btag_);
	btagFlag.push_back(2.);
	if(muChannel) channel.push_back(1.);
	else channel.push_back(0.);
	weight.push_back(evtWeight*scaleFact);
	puWeight.push_back(PUWeight);
	lumiWeight.push_back(scaleFact);
	runNumber.push_back( RunNum2Tag_);
	evtNumber.push_back( EvtNum2Tag_);
	lumiBlock.push_back( LumiBlock2Tag_);
	if(writeText) 
	  fileout<< RunNum2Tag_<<"  "
		 << EvtNum2Tag_<<"  "
		 << LumiBlock2Tag_<<"  "
		 <<lep1pt_2btag_<< "  " 
		 <<lep2pt_2btag_<< "  " 
		 <<zllmass_2btag_<< "  " 
		 <<jet1pt_2btag_<< "  " 
		 <<jet2pt_2btag_<< "  " 
		 <<zjjmass_2btag_<< "  " 
		 <<lljjmass_<< "  " 
		 << qgprod_2btag_ <<  "  " 
		 << metsig_2btag_ << "  " 
		 << ld_2btag_ << "  " 
		 << "0.5" << "  " 
		 << "2b" << endl;
      }
      if((lljjmass_1btag_!=0) && !exist2tag && exist1tag) { 
	lljjmass_1btag.Fill(lljjmass_1btag_, evtWeight);
	// store the best cand mass
	bestHiggsMass.push_back(lljjmass_1btag_);
	bestmjj.push_back(zjjmass_1btag_);
	btagFlag.push_back(1.);
	if(muChannel) channel.push_back(1.);
	else channel.push_back(0.);

	weight.push_back(evtWeight*scaleFact);
	puWeight.push_back(PUWeight);
	lumiWeight.push_back(scaleFact);
	runNumber.push_back( RunNum1Tag_);
	evtNumber.push_back( EvtNum1Tag_);
	lumiBlock.push_back( LumiBlock1Tag_);
	if(writeText) 
	  fileout<< RunNum1Tag_<<"  "
		 << EvtNum1Tag_<<"  "
		 << LumiBlock1Tag_<<"  "
		 <<lep1pt_1btag_<< "  " 
		 <<lep2pt_1btag_<< "  " 
		 <<zllmass_1btag_<< "  " 
		 <<jet1pt_1btag_<< "  " 
		 <<jet2pt_1btag_<< "  " 
		 <<zjjmass_1btag_<< "  " 
		 <<lljjmass_1btag_<< "  " 
		 << qgprod_1btag_ <<  "  " 
		 << metsig_1btag_ << "  " 
		 << ld_1btag_ << "  " 
		 << 0.000656*lljjmass_1btag_+0.302 << "  " 
		 << "1b" << endl;
      }
      if((lljjmass_0btag_!=0) && !exist2tag && !exist1tag && exist0tag) { 
	lljjmass_0btag.Fill(lljjmass_0btag_, evtWeight);
	// store the best cand mass
	bestHiggsMass.push_back(lljjmass_0btag_);
	bestmjj.push_back(zjjmass_0btag_);
	btagFlag.push_back(0.);
	if(muChannel) channel.push_back(1.);
	else channel.push_back(0.);
	weight.push_back(evtWeight*scaleFact);
	puWeight.push_back(PUWeight);
	lumiWeight.push_back(scaleFact);
	runNumber.push_back( RunNum0Tag_);
	evtNumber.push_back( EvtNum0Tag_);
	lumiBlock.push_back( LumiBlock0Tag_);
	if(writeText) 
	  fileout<< RunNum0Tag_<<"  "
		 << EvtNum0Tag_<<"  "
		 << LumiBlock0Tag_<<"  "
		 <<lep1pt_0btag_<< "  " 
		 <<lep2pt_0btag_<< "  " 
		 <<zllmass_0btag_<< "  " 
		 <<jet1pt_0btag_<< "  " 
		 <<jet2pt_0btag_<< "  " 
		 <<zjjmass_0btag_<< "  " 
		 <<lljjmass_0btag_<< "  " 
		 << qgprod_0btag_ <<  "  " 
		 << metsig_0btag_ << "  " 
		 << ld_0btag_ << "  " 
		 << 0.00025 * lljjmass_0btag_+0.55 << "  " 
		 << "0b" << endl;
      }

      if((lljjmassFinal_SB!=0) && !exist2tag && !exist1tag && !exist0tag &&  exist2tagSB) {
      	lljjmass_2btagSB.Fill(lljjmassFinal_SB, evtWeight);
      	bestHiggsMass.push_back(lljjmassFinal_SB);
      	bestmjj.push_back(zjjmass_SB2btag_);
      	btagFlag.push_back(2.);
	if(muChannel) channel.push_back(1.);
	else channel.push_back(0.);
  
      	weight.push_back(evtWeight*scaleFact);
      	puWeight.push_back(PUWeight);
      	lumiWeight.push_back(scaleFact);
	runNumber.push_back( RunNum2TagSB_);
	evtNumber.push_back( EvtNum2TagSB_);
	lumiBlock.push_back( LumiBlock2TagSB_);
      }
      if((lljjmassFinal_SB_1btag_!=0) && !exist2tag && !exist1tag && !exist0tag &&  !exist2tagSB && exist1tagSB) {
      	lljjmass_1btagSB.Fill(lljjmassFinal_SB_1btag_, evtWeight);
      	bestHiggsMass.push_back(lljjmassFinal_SB_1btag_);
      	bestmjj.push_back(zjjmass_SB1btag_);
      	btagFlag.push_back(1.);
	if(muChannel) channel.push_back(1.);
	else channel.push_back(0.);
  
      	weight.push_back(evtWeight*scaleFact);
      	puWeight.push_back(PUWeight);
      	lumiWeight.push_back(scaleFact);
	runNumber.push_back( RunNum1TagSB_);
	evtNumber.push_back( EvtNum1TagSB_);
	lumiBlock.push_back( LumiBlock1TagSB_);
      }
      if((lljjmassFinal_SB_0btag_!=0) && !exist2tag && !exist1tag && !exist0tag &&  !exist2tagSB && !exist1tagSB && exist0tagSB) {
      	lljjmass_0btagSB.Fill(lljjmassFinal_SB_0btag_, evtWeight);
      	bestHiggsMass.push_back(lljjmassFinal_SB_0btag_);
      	bestmjj.push_back(zjjmass_SB0btag_);
      	btagFlag.push_back(0.);
	if(muChannel) channel.push_back(1.);
	else channel.push_back(0.);
  
      	weight.push_back(evtWeight*scaleFact);
      	puWeight.push_back(PUWeight);
      	lumiWeight.push_back(scaleFact);
	runNumber.push_back( RunNum0TagSB_);
	evtNumber.push_back( EvtNum0TagSB_);
	lumiBlock.push_back( LumiBlock0TagSB_);
      }


      if(njets_!=-1) njets.Fill(njets_,evtWeight);
      if((lljjmass_noRefit_!=0) && exist2tag) 
	lljjmass_noRefit.Fill(lljjmass_noRefit_,evtWeight);
      if(hLD_Best!=-10) helyLD.Fill(hLD_Best,evtWeight);
      if(cos1_Best!=-10) cosT1.Fill(cos1_Best,evtWeight);
      if(cos2_Best!=-10) cosT2.Fill(cos2_Best,evtWeight);
      if(cosStar_Best!=-10) cosT1Star.Fill(cosStar_Best,evtWeight);
      if(phi_Best!=-10) phi.Fill(phi_Best,evtWeight);
      if(phiStar_Best!=-10) phiStar.Fill(phiStar_Best,evtWeight);
 
      // if(hLD_rfBest!=-10) helyLD_RF.Fill(hLD_rfBest,evtWeight);
      // if(cos1_rfBest!=-10) cosT1_RF.Fill(cos1_rfBest,evtWeight);
      // if(cos2_rfBest!=-10) cosT2_RF.Fill(cos2_rfBest,evtWeight);
      // if(cosStar_rfBest!=-10) cosT1Star_RF.Fill(cosStar_rfBest,evtWeight);
      // if(phi_rfBest!=-10) phi_RF.Fill(phi_rfBest,evtWeight);
      // if(phiStar_rfBest!=-10) phiStar_RF.Fill(phiStar_rfBest,evtWeight);

      if(lljjmass_SB != 0) lljjmass_sb.Fill(lljjmass_SB,evtWeight);
      if(lljjmass_noRefit_SB != 0) lljjmass_noRefit_sb.Fill(lljjmass_noRefit_SB,evtWeight);
      if(lljjmassFinal_SB != 0) lljjmassFinal_sb.Fill(lljjmassFinal_SB,evtWeight);
      if(lljjmassFinal_noRefit_SB != 0) lljjmassFinal_noRefit_sb.Fill(lljjmassFinal_noRefit_SB,evtWeight);


    }  // if cands > 0 ...
    
    for(unsigned int k = 0; k < counts; ++k)
      if(pass[k]) ++count[k];
  } // event loop
  
  
  //--- write tree--
  // TTree output file
 

  if(writeT) {
    TFile *outTree = TFile::Open(("tree_"+outFileName).c_str(), "RECREATE");
    TTree lljjmassTree("lljjmassTree","lljjmass tree");
    // mzz branch for the output tree
    double mzz;
    double mjj;
    double w;
    double pu_w;
    double lumi_w;
    double btagCat;
    double leptChannel;
    double runNum;
    double evtNum;
    double lumiB;
    lljjmassTree.Branch("mZZ", &mzz, "mZZ/D");
    lljjmassTree.Branch("mJJ", &mjj, "mJJ/D");
    lljjmassTree.Branch("weight", &w, "weight/D");
    lljjmassTree.Branch("PUweight", &pu_w, "PUweight/D");
    lljjmassTree.Branch("lumiweight", &lumi_w, "lumiweight/D");
    lljjmassTree.Branch("nBTags", &btagCat, "nBTags/D");
    lljjmassTree.Branch("lep", &leptChannel, "lep/D");
    lljjmassTree.Branch("runNum", &runNum, "runNum/D");
    lljjmassTree.Branch("evtNum", &evtNum, "evtNum/D");
    lljjmassTree.Branch("lumiB", &lumiB, "lumiB/D");
 
    


    for(size_t i = 0; i<bestHiggsMass.size(); ++i){

      mzz = bestHiggsMass[i];
      mjj = bestmjj[i];
      btagCat = btagFlag[i];
      w = weight[i];
      pu_w = puWeight[i];
      lumi_w = lumiWeight[i];
      leptChannel = channel[i];
      runNum = runNumber[i];
      evtNum = evtNumber[i];
      lumiB = lumiBlock[i];

      lljjmassTree.Fill();
    }


    TTree runNumTree("runNumTree","runNum tree");
    double runNumSkim_;
    runNumTree.Branch("runNumSkim_", &runNumSkim_, "runNumSkim_/D");

    for(size_t i = 0; i<runNumSkim.size(); ++i){
      runNumSkim_=runNumSkim[i];
      runNumTree.Fill();
    }
  

    // }
    // for(vdouble::iterator it=bestHiggsMass.begin(); it!= bestHiggsMass.end(); ++it ){
    //   mzz = *it;
    //   lljjmassTree.Fill();
    //   // cout<<"filling tree"<<endl;
    // }
   //  for(vdouble::iterator it=bestmjj.begin(); it!= bestmjj.end(); ++it ){
  //     mjj = *it;
  //     lljjmassTree.Fill();
  //   }
  //   for(vdouble::iterator itbtag=btagFlag.begin(); itbtag!= btagFlag.end(); ++itbtag ){
  //     btagCat = *itbtag;
  //     lljjmassTree.Fill();
  //   }
  //   for(vdouble::iterator it=weight.begin(); it!= weight.end(); ++it ){
  //     w = *it;
  //     lljjmassTree.Fill();
  //   }
  // for(vdouble::iterator it=puWeight.begin(); it!= puWeight.end(); ++it ){
  //     pu_w = *it;
  //     lljjmassTree.Fill();
  //   }
  // for(vdouble::iterator it=lumiWeight.begin(); it!= lumiWeight.end(); ++it ){
  //     lumi_w = *it;
  //     lljjmassTree.Fill();
  //   }
  //  for(vdouble::iterator it=channel.begin(); it!= channel.end(); ++it ){
  //     leptChannel = *it;
  //     lljjmassTree.Fill();
  //   }

    lljjmassTree.Write();
    runNumTree.Write();
    outTree->Close();
  }
  //----------------
  
  bestHiggsMass.clear();
  cout << endl << "--> Number of candidates/category after common selection <--" << endl;
  for (int kkk = 0; kkk<3; ++kkk)
    cout << kkk << "-Tags candidates = " << CandBTag[kkk] << endl;

  TH1F h_cuts("h_cuts", "ProgressiveCuts", 11, 0, 11.);
  
  cout<<count[1]<<"   "<<cand[1]<<endl;  
  
  // output histogram file
  //  TFile histos(outFile, "RECREATE");
  TFile histos(outFile, "RECREATE");
  histos.cd();
  lljjmass.Write();  
  lljjmass_noRefit.Write();  
  lljjmass_1btag.Write();  
  lljjmass_noRefit_1btag.Write();
  lljjmass_0btag.Write();  
  lljjmass_noRefit_0btag.Write();  

  lljjmass_0btagSB.Write();
  lljjmass_1btagSB.Write();
  lljjmass_2btagSB.Write();

  //  lljjmassCands.Write();  
  leptpt1.Write();  
  leptpt2.Write(); 
  lepteta1.Write();
  lepteta2.Write(); 
  leptphi1.Write();
  leptphi2.Write();

  eta_phi.Write();
  pt_phi.Write();
  eta_zllmass.Write();
  eta_pt.Write();
  jetpt1.Write();  
  jetpt2.Write();  
  h_met.Write();  
  metSig.Write();  
  metSignif.Write();  
  lljjmass2tags.Write();
  helyLD_RF2tags.Write();

  zllptBeforeLD.Write();
  drjjBefLD.Write();
  lljjmassBeforeLD.Write();


  cosT1.Write();
  cosT2.Write();
  cosT1Star.Write();
  phi.Write();
  phiStar.Write();
  helyLD.Write();
  qgLD.Write();

  cosT1_RF.Write();
  cosT2_RF.Write();
  cosT1Star_RF.Write();
  phi_RF.Write();
  phiStar_RF.Write();
  helyLD_RF.Write();


  helyLD_RF_zjj.Write();
  cosT1_RF_zjj.Write();
  cosT2_RF_zjj.Write();
  cosT1Star_RF_zjj.Write();
  phi_RF_zjj.Write();
  phiStar_RF_zjj.Write();


  njets.Write();
  btagJet1.Write();  
  btagJet2.Write();  
  btagJet.Write();
  btagJet1_JP.Write();  
  btagJet2_JP.Write();  
  btagJet_JP.Write();
  btagJet1_CSV.Write();  
  btagJet2_CSV.Write();  
  btagJet_CSV.Write();
  btagJet1_CSVMVA.Write();  
  btagJet2_CSVMVA.Write();  
  btagJet_CSVMVA.Write();
  
  drjj.Write();  
  zllpt.Write();
  stdCandle.Write();
  zjjmass.Write();
  zllmass.Write();
  zllmass_BB.Write();
  db.Write();
  //  ngenint.Write();
  npv_woReweight.Write();
  npv.Write();
  npv1.Write();
  npv2.Write();
  npvfin.Write();
  genhmass.Write();
  genhmass_rew.Write();
  LR_weights.Write();
  LR_weights_mass.Write();
   
  hmass_sb_all.Write();
  hmassnorefit_sb_all.Write();
  hmass_sb_2b.Write();
  hmassnorefit_sb_2b.Write();
  lljjmass_sb.Write();
  lljjmass_noRefit_sb.Write();
  lljjmassFinal_sb.Write();
  lljjmassFinal_noRefit_sb.Write();
  elMVAID.Write();
  beta.Write();

  cout << endl;
  cout<<"#TotEvts EvtsAtLeast1Cand Pt+Iso+db zll/zjjMmass   btag  HelyLD  metSig"<<endl;
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

  //Efficiencies
  cout << "efficiencies: ";
  for(unsigned int k = 0; k < counts; ++k) {
    cout << double(count[k])/double(count[0]) << " ";
  }
  cout << endl;
  //
  //Writing efficiencies to a txt file
  string EffName = "Efficiencies"+selChan+".txt";
  ofstream fileeffout(EffName.c_str(),ios::trunc);
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
  //fileeffout <<"Events passing selection:                     "<<count[counts-1]<< endl;
  //fileeffout <<"Candidates passing selection:                 "<<lljjmass.Integral()<< endl;
  //fileeffout fileeffout.close();
  //end Efficiencies to txt write

  fileout.close();
  return 0;
}

