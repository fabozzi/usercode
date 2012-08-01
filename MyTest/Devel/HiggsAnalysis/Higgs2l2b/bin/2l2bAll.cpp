
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TSystem.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include "DataFormats/Common/interface/Wrapper.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "HiggsAnalysis/Higgs2l2b/interface/BTagSFUtil.h"
// QG discriminant to be added for next skim
//#include "HiggsAnalysis/Higgs2l2b/interface/QGLikelihoodCalculator.h"

using namespace std;

typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;

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
  ofstream fileout(reportName.c_str(),ios::trunc);
  fileout<< "Event "<<"  "
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

  // define histograms
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
  TH1F leptpt1("leptpt1", "P_{T} ", 300, 0, 300);
  TH1F leptpt2("leptpt2", "P_{T} ", 300, 0, 300);
  TH1F jetpt1("jetpt1", "P_{T} ", 300, 0, 300);
  TH1F jetpt2("jetpt2", "P_{T} ", 300, 0, 300);
  TH1F zjjmass("zjjmass", "m_{jj}", 500, 0, 500);
  TH1F zllmass("zllmass", "m_{ll}", 500, 0, 500);
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
  
  // final inv mass histo
  TH1F lljjmass("lljjmass", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_0btag("lljjmass_0btag", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_1btag("lljjmass_1btag", "m(lljj)", 1000, 0, 1000);
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

  // ----- AVAILABLE FOR FUTURE NTUPLES ----------
  //QGLike discriminator
  //  string QGFilePDF = "QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root";
  //  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator( QGFilePDF );
  // ---------------------------------------------

  // LOOP on EVENTS
  for(unsigned int i = 0; i <nEvents; ++i) {
    //    if(i%100 == 0) progress(double(i)/double(nEvents));
    ++count[0];

    // variables for HLT selection
    bool passTrigPath(true);
    // NUMBER OF CANDIDATES IN THE EVENT
    unsigned int cands(0);

    // GET ENTRIES FOR EVENT VARIABLES
    GETENTRY(numPV,i);
    GETENTRY(met,i);
    GETENTRY(metSignificance,i);
    // ------ AVAILABLE FOR FUTURE NTUPLES --------
    //    GETENTRY(rhoRestrictedEta,i);
    // --------------------------------------------

    // ------ TO DECIDE WHAT TO DO WITH TRIGGER BITS -------
    GETENTRY(passSingleMuTrig,i);
    GETENTRY(passDoubleMuTrig,i);
    GETENTRY(passSingleElTrig,i);
    GETENTRY(passDoubleElTrig,i);
    // -----------------------------------------------------

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
      // ------ TO DECIDE WHAT TO DO WITH TRIGGER BITS -------
      if(data)
      	passTrigPath = getInt(passDoubleMuTrig);
      // -----------------------------------------------------
      GETENTRY(muHiggsLeptDau1Pt,i);
      cands = muHiggsLeptDau1Pt->product()->size();
    } else {
      // ------ TO DECIDE WHAT TO DO WITH TRIGGER BITS -------
      if(data)
      	passTrigPath = !(getInt(passDoubleMuTrig)) && getInt(passDoubleElTrig);
      // -----------------------------------------------------
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
    if(cands > 0) {
      pass[1] = true;


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
      float lljjmass_0btag_(0), lljjmass_noRefit_0btag_(0);
      float SumPtBest_(0), SumPtBest_Final(0), SumPtBest_Final_1btag(0);
      bool exist2tag(false), exist1tag(false), exist0tag(false);
      float ResidZllMassBest_(10000.), ResidZllMassBest_Final(10000.);
      float ResidZllMassBest1tag_(10000.);
      float ResidZllMassBest0tag_(10000.);
      float hLD_Best(-10), cos1_Best(-10), cos2_Best(-10), cosStar_Best(-10);
      float phi_Best(-10), phiStar_Best(-10);
      float hLD_rfBest(-10), cos1_rfBest(-10), cos2_rfBest(-10), cosStar_rfBest(-10);
      float phi_rfBest(-10), phiStar_rfBest(-10);
      int EvtNum2Tag_(0), RunNum2Tag_(0);
      int EvtNum1Tag_(0), RunNum1Tag_(0);
      int EvtNum0Tag_(0), RunNum0Tag_(0);
      float ld_2btag_(-10), ld_1btag_(-10), ld_0btag_(-10);
      float zllmass_2btag_(-10), zllmass_1btag_(-10), zllmass_0btag_(-10);
      float zjjmass_2btag_(-10), zjjmass_1btag_(-10), zjjmass_0btag_(-10);
      float qgprod_2btag_(-10), qgprod_1btag_(-10), qgprod_0btag_(-10);
      float metsig_2btag_(-10), metsig_1btag_(-10), metsig_0btag_(-10);
      float jet1pt_2btag_(-10), jet1pt_1btag_(-10), jet1pt_0btag_(-10);
      float jet2pt_2btag_(-10), jet2pt_1btag_(-10), jet2pt_0btag_(-10);
      float lep1pt_2btag_(-10), lep1pt_1btag_(-10), lep1pt_0btag_(-10);
      float lep2pt_2btag_(-10), lep2pt_1btag_(-10), lep2pt_0btag_(-10);

      float SumPtBest_SB(0), SumPtBest_FinalSB(0), SumPtBest_FinalSB_1btag(0);
      float lljjmass_SB(0), lljjmass_noRefit_SB(0);
      float lljjmassFinal_SB(0), lljjmassFinal_noRefit_SB(0), lljjmassFinal_SB_1btag_(0), lljjmassFinal_noRefit_SB_1btag_(0);


      // auxiliary variables for channel-independent selection
      float Jet1pt_(0), Jet2pt_(0), Jet1pt_rf(0), Jet2pt_rf(0);;
      float Jet1eta_(0), Jet2eta_(0);
      float lept1pt_(0), lept2pt_(0);
      float zllcharge_(0), zllmass_(0), zjjmass_(0);
      float jet1tkhe_(0), jet2tkhe_(0);
      float jet1csv_(0), jet2csv_(0);
      int jet1flav_(0), jet2flav_(0);
      // ----------- variables for QG discriminant ----------
      int jet1Nch_(-1), jet2Nch_(-1);
      int jet1Nneu_(-1), jet2Nneu_(-1);
      float jet1ptD_(-1), jet2ptD_(-1);
      // ----------------------------------------------------
      float zllpt_(-1.), drjj_(-1.);
      float hLD_(-10), cos1_(-10), cos2_(-10), cosStar_(-10), phi_(-10), phiStar_(-10);
      float hLD_rf(-10), cos1_rf(-10), cos2_rf(-10), cosStar_rf(-10), phi_rf(-10), phiStar_rf(-10);
      float higgsrefitMass(0), higgsMass(0);
      int higgsEvtNum(0), higgsRunNum(0);
      
      bool kineLepCut(false), fiducialCut(false), kineJetCut(false);
      bool dbLepCut(false);
      bool IDLep1Cut(false), IDLep2Cut(false), isoIDLepCut(false);

      // flag to dump event variables
      bool isEventSelected(false);

      // LOOP on CANDIDATES IN THE EVENT
      for(unsigned int j = 0; j < cands; ++j) {
	++cand[1];

	// fill histo of db of the muon candidates
	if(muChannel) {
	  db.Fill(get(muHiggsLeptDau1dB,j),evtWeight);
	  db.Fill(get(muHiggsLeptDau2dB,j),evtWeight);
	}

	// assign variables for selection
	if(muChannel) {
	  Jet1pt_=get(muHiggsJetDau1Pt,j);
	  Jet2pt_=get(muHiggsJetDau2Pt,j);
	  Jet1pt_rf=get(muHiggsJetDau1RefitPt,j);
	  Jet2pt_rf=get(muHiggsJetDau2RefitPt,j);
	  Jet1eta_=get(muHiggsJetDau1Eta,j);
	  Jet2eta_=get(muHiggsJetDau2Eta,j);
	  lept1pt_=get(muHiggsLeptDau1Pt,j);
	  lept2pt_=get(muHiggsLeptDau2Pt,j);

	  fiducialCut = true;
	  // db cut
	  dbLepCut = (fabs(get(muHiggsLeptDau1dB,j)) <0.02) && (fabs(get(muHiggsLeptDau2dB,j)) <0.02);
	  // muonID cut
	  IDLep1Cut = (get(muHiggsLeptDau1GlobalMuonBit,j)) && 
	    (get(muHiggsLeptDau1TrackerMuonBit,j)) && 
	    (get(muHiggsLeptDau1NormChi2,j) < 10) && 
	    //            (get(muHiggsLeptDau1NofTrackerHits,j) > 10) &&
            (get(muHiggsLeptDau1NofPixelHits,j) > 0) &&   
            (get(muHiggsLeptDau1NofMuonHits,j) > 0); //&& 
	    //	    (get(muHiggsLeptDau1NofMatches,j) > 1);
	
	  IDLep2Cut = (get(muHiggsLeptDau2GlobalMuonBit,j)) && 
	    (get(muHiggsLeptDau2TrackerMuonBit,j)) && 
	    (get(muHiggsLeptDau2NormChi2,j) < 10) && 
	    //            (get(muHiggsLeptDau2NofTrackerHits,j) > 10) &&
            (get(muHiggsLeptDau2NofPixelHits,j) > 0) &&   
            (get(muHiggsLeptDau2NofMuonHits,j) > 0) ; //&& 
	    //	    (get(muHiggsLeptDau2NofMatches,j) > 1);
	
	  // isolation + ID cut
	  isoIDLepCut = (get(muHiggsLeptDau1CombRelIso,j) < 0.15) && (get(muHiggsLeptDau2CombRelIso,j) < 0.15) && IDLep1Cut && IDLep2Cut ;

	  zllcharge_ = get(muHiggszllCharge,j);
	  zllmass_ = get(muHiggszllMass,j);
	  zjjmass_ = get(muHiggszjjMass,j); 
	  jet1tkhe_ = get(muHiggsJet1TKHE,j);
	  jet2tkhe_ = get(muHiggsJet2TKHE,j);
	  jet1csv_ = get(muHiggsJet1CSV,j);
	  jet2csv_ = get(muHiggsJet2CSV,j);
	  jet1flav_ = get(muHiggsJetDau1PartonFlavour,j);
	  jet2flav_ = get(muHiggsJetDau2PartonFlavour,j);

	  // ------ AVAILABLE FOR FUTURE NTUPLES --------
	  //	  jet1Nch_ = int( get(muHiggsJetDau1ChHadMult,j) );
	  //	  jet2Nch_ = int( get(muHiggsJetDau2ChHadMult,j) );
	  //	  jet1Nneu_ = int( get(muHiggsJetDau1NeuHadMult,j) ) + int( get(muHiggsJetDau1PhotMult,j) );
	  //	  jet2Nneu_ = int( get(muHiggsJetDau2NeuHadMult,j) ) + int( get(muHiggsJetDau2PhotMult,j) );
	  //	  jet1ptD_ = get(muHiggsJetDau1PtDJet,j); 
	  //	  jet2ptD_ = get(muHiggsJetDau2PtDJet,j);
	  // --------------------------------------------

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
	} else {
	  Jet1pt_=get(elHiggsJetDau1Pt,j);
	  Jet2pt_=get(elHiggsJetDau2Pt,j);
	  Jet1pt_rf=get(elHiggsJetDau1RefitPt,j);
	  Jet2pt_rf=get(elHiggsJetDau2RefitPt,j);
	  Jet1eta_=get(elHiggsJetDau1Eta,j);
	  Jet2eta_=get(elHiggsJetDau2Eta,j);
	  lept1pt_=get(elHiggsLeptDau1Pt,j);
	  lept2pt_=get(elHiggsLeptDau2Pt,j);

	  // ----- FUTURE NTUPLES ------
	  //	  fiducialCut = (fabs(get(elHiggsLeptDau1Eta,j)) < 2.5 ) &&
	  //	    (fabs(get(elHiggsLeptDau2Eta,j)) < 2.5 ) &&
	  //	    !(  (fabs(get(elHiggsLeptDau1EtaSC,j)) < 1.566) &&
	  //		(fabs(get(elHiggsLeptDau1EtaSC,j)) > 1.4442) ) &&  
	  //	    !(  (fabs(get(elHiggsLeptDau2EtaSC,j)) < 1.566) &&
	  //		(fabs(get(elHiggsLeptDau2EtaSC,j)) > 1.4442) );
	  //-----------------------------
	  fiducialCut = true;
	  dbLepCut = true;
	  // new electronID
	  bool EleDau1VBTFID = get(elHiggsEleDau1VBTF95CombID,j) == 7;

	  // ----- FUTURE NTUPLES ------
	  //	  bool isEle1Barrel = fabs(get(elHiggsLeptDau1EtaSC,j)) <= 1.4442;
	  //----------------------------
	  bool isEle1Barrel = get(elHiggsLeptDau1isEB,j);

	  bool EleDau2VBTFID = get(elHiggsEleDau2VBTF95CombID,j) == 7;

	  // ----- FUTURE NTUPLES ------
	  //	  bool isEle2Barrel = fabs(get(elHiggsLeptDau2EtaSC,j)) <= 1.4442;
	  // ---------------------------
	  bool isEle2Barrel = get(elHiggsLeptDau2isEB,j);

	  bool EleDau1TightCuts = EleIDTightCuts( isEle1Barrel, get(elHiggsLeptDau1DeltaPhiAtVtx,j), 
						  get(elHiggsLeptDau1HOverE,j) );
	  bool EleDau2TightCuts = EleIDTightCuts( isEle2Barrel, get(elHiggsLeptDau2DeltaPhiAtVtx,j), 
						  get(elHiggsLeptDau2HOverE,j) );

	  isoIDLepCut = EleDau1VBTFID && EleDau1TightCuts && EleDau2VBTFID && EleDau2TightCuts;

	  zllcharge_ = get(elHiggszllCharge,j);
	  zllmass_ = get(elHiggszllMass,j);
	  zjjmass_ = get(elHiggszjjMass,j); 
	  jet1tkhe_ = get(elHiggsJet1TKHE,j);
	  jet2tkhe_ = get(elHiggsJet2TKHE,j);
	  jet1csv_ = get(elHiggsJet1CSV,j);
	  jet2csv_ = get(elHiggsJet2CSV,j);
	  jet1flav_ = get(elHiggsJetDau1PartonFlavour,j);
	  jet2flav_ = get(elHiggsJetDau2PartonFlavour,j);

	  // ------ AVAILABLE FOR FUTURE NTUPLES --------
	  //jet1Nch_ = int( get(elHiggsJetDau1ChHadMult,j) );
	  //jet2Nch_ = int( get(elHiggsJetDau2ChHadMult,j) );
	  //jet1Nneu_ = int( get(elHiggsJetDau1NeuHadMult,j) ) + int( get(elHiggsJetDau1PhotMult,j) );
	  //jet2Nneu_ = int( get(elHiggsJetDau2NeuHadMult,j) ) + int( get(elHiggsJetDau2PhotMult,j) );
	  //jet1ptD_ = get(elHiggsJetDau1PtDJet,j); 
	  //jet2ptD_ = get(elHiggsJetDau2PtDJet,j);
	  // --------------------------------------------

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
	  
	}

	kineLepCut = (lept1pt_ > 40 && lept2pt_ > 20) || ( lept1pt_ > 20 &&  lept2pt_ > 40 );
	kineJetCut = (Jet1pt_ > 30) && (Jet2pt_ > 30);
	
	if(kineLepCut && fiducialCut && dbLepCut && isoIDLepCut && kineJetCut) {
	  pass[2] = true;
	  ++cand[2];

	  // fill plots after lepton selection
	  leptpt1.Fill(lept1pt_,evtWeight);
	  leptpt2.Fill(lept2pt_,evtWeight);
	  zllmass.Fill(zllmass_,evtWeight);
	  npv1.Fill(getInt(numPV),evtWeight);
	  
	  // Z inv. mass selection	  
	  bool zllcut = (zllmass_ > 70) && (zllmass_ < 110);
	  // Z charge selection	  
	  bool zchargecut = (zllcharge_ == 0) ;
	  // zjj cut different for signal vs. sideband
	  bool zjjcut_sigSel = (zjjmass_ > 75) && (zjjmass_ < 105);
	  bool zjjcut_sbSel = zjjmass_ > 120;


	  //SIGNAL SELECTION
	  if( zllcut && zchargecut ) {
	    
	    zjjmass.Fill(zjjmass_,evtWeight);	   
	    jetpt1.Fill(Jet1pt_,evtWeight);
	    jetpt2.Fill(Jet2pt_,evtWeight);
	    
	    if( zjjcut_sigSel ) {	      
	      pass[3] = true;
	      ++cand[3];
	      
	      // fill plots after Z mass selection
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
	      npv2.Fill(getInt(numPV),evtWeight);
	  
	      // end of common selection
	      // now jet classification according btag
	      
	      //Tag algorithm cut definition for signal
	      bool jet1_tagged_medium(false), jet1_tagged_loose(false);
	      bool jet2_tagged_medium(false), jet2_tagged_loose(false);
	      
 	      jet1_tagged_medium = jet1tkhe_ > 3.3;
 	      jet1_tagged_loose  = jet1tkhe_ > 1.7;	    
 	      jet2_tagged_medium = jet2tkhe_ > 3.3;
 	      jet2_tagged_loose  = jet2tkhe_ > 1.7;

// 	      jet1_tagged_medium = jet1csv_ > 0.679;
// 	      jet1_tagged_loose  = jet1csv_ > 0.244;	    
// 	      jet2_tagged_medium = jet2csv_ > 0.679;
// 	      jet2_tagged_loose  = jet2csv_ > 0.244;

	      // eventually apply SF for MC
	      if( btagScale ) {
		btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, Jet1pt_, Jet1eta_, jet1flav_ );
		btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, Jet2pt_, Jet2eta_, jet2flav_ );
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
	      
	      // select only 2-tag candidates
	      if(twoTagCategory){
		pass[4] = true;
		++cand[4];
		// fill plots for 2-tag candidates selection
		zllpt.Fill(zllpt_,evtWeight);
		drjj.Fill(drjj_,evtWeight);
		metSignif.Fill(getInt(metSignificance),evtWeight);
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
		  if ( (fabs(zllmass_-91.187) < ResidZllMassBest_) && (higgsrefitMass>0) ){
		    // if best candidate store its variables
		    //		  SumPtBest_= Jet1pt_rf+ Jet2pt_rf;
		    ResidZllMassBest_= fabs(zllmass_-91.187);
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
		    if(!isEventSelected)
		      npvfin.Fill(getInt(numPV),evtWeight);

		    isEventSelected = true;
		    // check if it is the best candidate at this level
		    //		  if ( ( (Jet1pt_rf+ Jet2pt_rf)>SumPtBest_Final) && (higgsrefitMass>0) ){
		    if ( (fabs(zllmass_-91.187) < ResidZllMassBest_Final) && (higgsrefitMass>0) ){
		      // if best candidate store lljj mass
		      //		    SumPtBest_Final= Jet1pt_rf+ Jet2pt_rf;
		      exist2tag = true;
		      ResidZllMassBest_Final= fabs(zllmass_-91.187);
		      lljjmass_= higgsrefitMass;
		      lljjmass_noRefit_= higgsMass;
		      EvtNum2Tag_=higgsEvtNum;
		      RunNum2Tag_=higgsRunNum;
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
		  if ( (fabs(zllmass_-91.187) < ResidZllMassBest1tag_) && (higgsrefitMass>0) ){
		    // if best candidate store lljj mass
		    exist1tag = true;
		    ResidZllMassBest1tag_ = fabs(zllmass_-91.187);
		    lljjmass_1btag_ = higgsrefitMass;
		    lljjmass_noRefit_1btag_ = higgsMass;
		    EvtNum1Tag_=higgsEvtNum;
		    RunNum1Tag_=higgsRunNum;
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

	  // ------ AVAILABLE FOR FUTURE NTUPLES --------
		// QG cut
		//		float QGLikeJet1 = qglikeli->computeQGLikelihoodPU( Jet1pt_, getInt(rhoRestrictedEta), jet1Nch_, jet1Nneu_, jet1ptD_ );
		//		float QGLikeJet2 = qglikeli->computeQGLikelihoodPU( Jet2pt_, getInt(rhoRestrictedEta), jet2Nch_, jet2Nneu_, jet2ptD_ );

		//		if( QGLikeJet1 * QGLikeJet2 > 0.10) {
	  // --------------------------------------------

		// LD cut
		if( hLD_rf >(0.55+(0.00025*higgsrefitMass)) ){
		  // check if it is the best candidate at this level
		  if ( (fabs(zllmass_-91.187) < ResidZllMassBest0tag_) && (higgsrefitMass>0) ){
		    exist0tag = true;
		    ResidZllMassBest0tag_ = fabs(zllmass_-91.187);
		    lljjmass_0btag_ = higgsrefitMass;
		    lljjmass_noRefit_0btag_ = higgsMass;
		    EvtNum0Tag_=higgsEvtNum;
		    RunNum0Tag_=higgsRunNum;
		    ld_0btag_ = hLD_rf;
		    //		    qgprod_0btag_ = QGLikeJet1*QGLikeJet2;
		    metsig_0btag_ = -1;
		    jet1pt_0btag_ = Jet1pt_rf;
		    jet2pt_0btag_ = Jet2pt_rf;
		    lep1pt_0btag_ = lept1pt_;
		    lep2pt_0btag_ = lept2pt_;
		    zllmass_0btag_ = zllmass_;
		    zjjmass_0btag_ = zjjmass_;
		  }
		} // end LD cut
		//		} // end QG cut <------- FUTURE NTUPLES
	      } // end 0-tag category
	    


	    }//end zjj cut

	  } // end of signal selection (zll cut)


	  // SIDEBAND SELECTION
	  if( zllcut && zjjcut_sbSel ) {
	    
	    // all candidates in the sideband
	    hmass_sb_all.Fill(higgsrefitMass);
	    hmassnorefit_sb_all.Fill(higgsMass);
	    //

	    //Tag algorithm cut definition for sideband
	    bool jet1_tagged_medium(false), jet1_tagged_loose(false);
	    bool jet2_tagged_medium(false), jet2_tagged_loose(false);
 	    jet1_tagged_medium = jet1tkhe_ > 3.3;
 	    jet1_tagged_loose  = jet1tkhe_ > 1.7;	    
 	    jet2_tagged_medium = jet2tkhe_ > 3.3;
 	    jet2_tagged_loose  = jet2tkhe_ > 1.7;

// 	    jet1_tagged_medium = jet1csv_ > 0.679;
// 	    jet1_tagged_loose  = jet1csv_ > 0.244;	    
// 	    jet2_tagged_medium = jet2csv_ > 0.679;
// 	    jet2_tagged_loose  = jet2csv_ > 0.244;

	    if( btagScale ) {
	      btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, Jet1pt_, Jet1eta_, jet1flav_ );
	      btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, Jet2pt_, Jet2eta_, jet2flav_ );
	    }
	    
//  	    bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_medium  );
	    bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_loose  )
	      || (jet1_tagged_loose  && jet2_tagged_medium );
// 	    bool twoTagCategory  = ( jet1_tagged_loose && jet2_tagged_loose  );
// 	    bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_loose || jet2_tagged_loose );
// 	    bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_medium || jet2_tagged_medium );
// 	    bool zeroTagCategory = ( !twoTagCategory ) && ( !oneTagCategory );


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
		  if ( ( (Jet1pt_rf+ Jet2pt_rf)>SumPtBest_FinalSB) && (higgsrefitMass>0) ){
		    SumPtBest_FinalSB= Jet1pt_rf+ Jet2pt_rf;
		    lljjmassFinal_SB = higgsrefitMass;
		    lljjmassFinal_noRefit_SB = higgsMass;
		  }
		}
	      }
	    }

	  } // end of sideband selection


	}
      } // candidate loop
      

      // fill best candidate plots
      if((lljjmass_!=0) && exist2tag) { 
	lljjmass.Fill(lljjmass_, evtWeight);
	// store the best cand mass
	bestHiggsMass.push_back(lljjmass_);
	if(writeText) 
	  fileout<< EvtNum2Tag_<<"  "
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
	bestHiggsMass_1btag.push_back(lljjmass_1btag_);
	if(writeText) 
	  fileout<< EvtNum1Tag_<<"  "
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
	bestHiggsMass_0btag.push_back(lljjmass_0btag_);
	if(writeText) 
	  fileout<< EvtNum0Tag_<<"  "
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

      if((lljjmass_noRefit_!=0) && exist2tag) 
	lljjmass_noRefit.Fill(lljjmass_noRefit_,evtWeight);
      if(hLD_Best!=-10) helyLD.Fill(hLD_Best,evtWeight);
      if(cos1_Best!=-10) cosT1.Fill(cos1_Best,evtWeight);
      if(cos2_Best!=-10) cosT2.Fill(cos2_Best,evtWeight);
      if(cosStar_Best!=-10) cosT1Star.Fill(cosStar_Best,evtWeight);
      if(phi_Best!=-10) phi.Fill(phi_Best,evtWeight);
      if(phiStar_Best!=-10) phiStar.Fill(phiStar_Best,evtWeight);
 
      if(hLD_rfBest!=-10) helyLD_RF.Fill(hLD_rfBest,evtWeight);
      if(cos1_rfBest!=-10) cosT1_RF.Fill(cos1_rfBest,evtWeight);
      if(cos2_rfBest!=-10) cosT2_RF.Fill(cos2_rfBest,evtWeight);
      if(cosStar_rfBest!=-10) cosT1Star_RF.Fill(cosStar_rfBest,evtWeight);
      if(phi_rfBest!=-10) phi_RF.Fill(phi_rfBest,evtWeight);
      if(phiStar_rfBest!=-10) phiStar_RF.Fill(phiStar_rfBest,evtWeight);

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
  //  lljjmassCands.Write();  
  leptpt1.Write();  
  leptpt2.Write();  
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

  cosT1_RF.Write();
  cosT2_RF.Write();
  cosT1Star_RF.Write();
  phi_RF.Write();
  phiStar_RF.Write();
  helyLD_RF.Write();

  btagJet1.Write();  
  btagJet2.Write();  
  btagJet.Write();
  btagJet1_CSV.Write();  
  btagJet2_CSV.Write();  
  btagJet_CSV.Write();
  
  drjj.Write();  
  zllpt.Write();
  zjjmass.Write();
  zllmass.Write();
  db.Write();
  //  ngenint.Write();
  npv_woReweight.Write();
  npv.Write();
  npv1.Write();
  npv2.Write();
  npvfin.Write();

  hmass_sb_all.Write();
  hmassnorefit_sb_all.Write();
  hmass_sb_2b.Write();
  hmassnorefit_sb_2b.Write();
  lljjmass_sb.Write();
  lljjmass_noRefit_sb.Write();
  lljjmassFinal_sb.Write();
  lljjmassFinal_noRefit_sb.Write();


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

  return 0;
}

