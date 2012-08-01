
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
  TH1F npv_woReweight("npv_woReweight", "npvwoReweight", 31, -0.5, 30.5);
  TH1F npv("npv", "npv", 31, -0.5, 30.5);
  TH1F npv1("npv1", "npv", 31, -0.5, 30.5);
  TH1F npv2("npv2", "npv", 31, -0.5, 30.5);
  TH1F npvfin("npvfin", "npv", 31, -0.5, 30.5);
  TH1F db("db", "db", 100, -0.05, 0.05);
  // histo before Z mass cuts
  TH1F leptpt1("leptpt1", "P_{T} ", 300, 0, 300);
  TH1F leptpt2("leptpt2", "P_{T} ", 300, 0, 300);
  TH1F jetpt1("jetpt1", "P_{T} ", 600, 0, 600);
  TH1F jetpt2("jetpt2", "P_{T} ", 600, 0, 600);
  TH1F zjjmass("zjjmass", "m_{jj}", 1500, 0, 1500);
  TH1F zllmass("zllmass", "m_{ll}", 1500, 0, 1500);
  TH1F lljjmassallcands("lljjmassallcands", "m(lljj)", 2000, 0, 2000);

  // histo before b-tagging cut
  TH1F zllmass_zveto("zllmass_zveto", "m_{ll}", 1500, 0, 1500);
  TH1F btagJet1("TCHEJet1", "TCHE", 100, 0, 20);
  TH1F btagJet2("TCHEJet2", "TCHE", 100, 0, 20);
  TH1F btagJet("TCHEJet", "TCHE", 100, 0, 20);
  TH1F btagJet1_CSV("CSVJet1", "CSV", 100, 0, 1);
  TH1F btagJet2_CSV("CSVJet2", "CSV", 100, 0, 1);
  TH1F btagJet_CSV("CSVJet", "CSV", 100, 0, 1);
  TH1F jetpt1_zveto("jetpt1_zveto", "P_{T} ", 600, 0, 600);
  TH1F jetpt2_zveto("jetpt2_zveto", "P_{T} ", 600, 0, 600);
  TH1F zjjmass_zveto("zjjmass_zveto", "m_{jj}", 1500, 0, 1500);

  TH1F lljjmassalltags("lljjmassalltags", "m(lljj)", 2000, 0, 2000);

  // histo for 2b candidates
  TH1F zllpt("zllpt", "P_{T}(Z_{ll})", 600, 0, 600);
  TH1F drjj("DRjj", "#Delta R_{jj}", 100, 0, 4);
  TH1F metSignif("metSignif", "METSignif", 100, 0, 15);
  TH1F jetpt1_2tags("jetpt1_2tags", "P_{T} ", 600, 0, 600);
  TH1F jetpt2_2tags("jetpt2_2tags", "P_{T} ", 600, 0, 600);
  TH1F zjjmass_2tags("zjjmass_2tags", "m_{jj}", 1500, 0, 1500);
  TH1F lljjmass2tags("lljjmass2tags", "m(lljj)", 1000, 0, 2000);

  // histo for 1b candidates
  TH1F lljjmass1tag("lljjmass1tag", "m(lljj)", 2000, 0, 2000);
  TH1F jetpt1_1tag("jetpt1_1tag", "P_{T} ", 600, 0, 600);
  TH1F jetpt2_1tag("jetpt2_1tag", "P_{T} ", 600, 0, 600);
  TH1F zjjmass_1tag("zjjmass_1tag", "m_{jj}", 1500, 0, 1500);

  // histo for 1b+2b candidates
  TH1F lljjmass1and2tags("lljjmass1and2tags", "m(lljj)", 2000, 0, 2000);
  TH1F jetpt1_1and2tags("jetpt1_1and2tags", "P_{T} ", 600, 0, 600);
  TH1F jetpt2_1and2tags("jetpt2_1and2tags", "P_{T} ", 600, 0, 600);
  TH1F zjjmass_1and2tags("zjjmass_1and2tags", "m_{jj}", 1500, 0, 1500);


  // final inv mass histo
  TH1F lljjmass("lljjmass", "m(lljj)", 2000, 0, 2000);


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
      //      LumiWeights_ = edm::LumiReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011A_Truth_v2.root", "PUS4_Distr", "pileup");
      LumiWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011A_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    //  --> Run2011B <--
    else if(PUPeriod == "B"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011B" << endl;
      //      LumiWeights_ = edm::LumiReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011B_Truth_v2.root", "PUS4_Distr", "pileup");
      LumiWeights_ = edm::Lumi3DReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011B_Truth_v2_finebin.root", "PUS4_Distr", "pileup");
    }
    //  --> Run2011All <--
    else if(PUPeriod == "All"){
      cout << "PU reweighting ------> USING PU MODEL FOR Run2011" << endl;
      //      LumiWeights_ = edm::LumiReWeighting("PUDist_Summer11MC_Flat10.root", "PUDist_Run2011All_Truth_v2.root", "PUS4_Distr", "pileup");
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
    //      LumiWeights_.weight3D_init();
  }

  // LOOP on EVENTS
  for(unsigned int i = 0; i <nEvents; ++i) {
    //    if(i%100 == 0) progress(double(i)/double(nEvents));
    ++count[0];

    // variables for HLT selection
    //    bool passDoubleTrig(false), passSingleTrig(false);
    bool passTrigPath(false);
    // NUMBER OF CANDIDATES IN THE EVENT
    unsigned int cands(0);

    // GET ENTRIES FOR EVENT VARIABLES
    GETENTRY(numPV,i);
    GETENTRY(met,i);
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
      //      evtWeight = LumiWeights_.weight3BX(getInt(nGenInt));
      //      ngenint.Fill( getInt(nGenInt) );
    }      

    
    if(muChannel) {
      passTrigPath = getInt(passDoubleMuTrig);
      //      if(fixMu && data)
      //	passTrigPath = true;
      //      passSingleTrig = getInt(passSingleMuTrig);
      //      passDoubleTrig = getInt(passDoubleMuTrig) && !(getInt(passDoubleElTrig));
      GETENTRY(muHiggsLeptDau1Pt,i);
      cands = muHiggsLeptDau1Pt->product()->size();
    } else {
      passTrigPath = !(getInt(passDoubleMuTrig)) && getInt(passDoubleElTrig);
      //      passSingleTrig = getInt(passSingleElTrig);
      //      passDoubleTrig = getInt(passDoubleElTrig);
      GETENTRY(elHiggsLeptDau1Pt,i);
      cands = elHiggsLeptDau1Pt->product()->size();
    }
    
    // HLT selection: exclusive double trig selection
    //    if(!passTrigPath)
    //      continue;

    fill(pass.begin(), pass.end(), false);
    if(cands > 0) {
      pass[1] = true;


      // fill histo of numPV (reweighted and not)
      npv.Fill(getInt(numPV),evtWeight);
      npv_woReweight.Fill(getInt(numPV));
    
      //GET ENTRIES FOR MU/EL VARIABLES
      if(muChannel){
#include "HiggsAnalysis/Higgs2l2b/bin/2l2bMu_getent.h"
      } else {
#include "HiggsAnalysis/Higgs2l2b/bin/2l2bEl_getent.h"
      }

      // auxiliary variables for channel-independent selection
      float Jet1pt_(0), Jet2pt_(0), Jet1pt_rf(0), Jet2pt_rf(0);;
      float Jet1eta_(0), Jet2eta_(0);
      float lept1pt_(0), lept2pt_(0);
      float zllmass_(0), zjjmass_(0);
      float jet1tkhe_(0), jet2tkhe_(0);
      float jet1csv_(0), jet2csv_(0);
      int jet1flav_(0), jet2flav_(0);
      float zllpt_(-1.), drjj_(-1.);
      float higgsrefitMass(0), higgsMass(0);
      int higgsEvtNum(0), higgsRunNum(0);
      
      // muon kine / isolation cut / cosmic rejection
      bool kineLepCut(false), dbLepCut(false), isoIDLepCut(false);

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
	  drjj_ = get(muHiggsjjdr,j);
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
	  drjj_ = get(elHiggsjjdr,j);
	  higgsrefitMass = get(elHiggsRefitMass,j);
	  higgsMass = get(elHiggsMass,j);
	  higgsEvtNum = getInt(elHiggsEventNumber);
	  higgsRunNum = getInt(elHiggsRunNumber);
	  
	}

	kineLepCut = (lept1pt_ > 25 && lept2pt_ > 20) || ( lept1pt_ > 20 &&  lept2pt_ > 25 );

	leptpt1.Fill(lept1pt_,evtWeight);
	leptpt2.Fill(lept2pt_,evtWeight);
	
	if(kineLepCut && dbLepCut && isoIDLepCut) {

	  // jet pT selection
	  if( (Jet1pt_>120) || (Jet2pt_>120) ) {

	    pass[2] = true;
	    ++cand[2];
	    
	    // fill plots after lepton/jet selection
	    
	    zllmass.Fill(zllmass_,evtWeight);
	    npv1.Fill(getInt(numPV),evtWeight);
	    zjjmass.Fill(zjjmass_,evtWeight);	   
	    jetpt1.Fill(Jet1pt_,evtWeight);
	    jetpt2.Fill(Jet2pt_,evtWeight);
	    lljjmassallcands.Fill(higgsMass,evtWeight);
	    
	    // Z inv. mass selection -> MODIFIED FOR BUMP SEARCH
	    bool zllcut = (zllmass_ > 70) && (zllmass_ < 100);

	    //BUMP SEARCH -> INVERT Z->ll MASS CUT
	    if( !zllcut ) {
	      
	      // DO NOT APPLY ANY Z->jj mass cut
	      pass[3] = true;
	      ++cand[3];
	    
	      // fill plots after Z mass selection
	      zllmass_zveto.Fill(zllmass_,evtWeight);
	      btagJet1.Fill(jet1tkhe_,evtWeight);
	      btagJet2.Fill(jet2tkhe_,evtWeight);
	      btagJet.Fill(jet1tkhe_,evtWeight);
	      btagJet.Fill(jet2tkhe_,evtWeight);
	      btagJet1_CSV.Fill(jet1csv_,evtWeight);
	      btagJet2_CSV.Fill(jet2csv_,evtWeight);
	      btagJet_CSV.Fill(jet1csv_,evtWeight);
	      btagJet_CSV.Fill(jet2csv_,evtWeight);

	      zjjmass_zveto.Fill(zjjmass_,evtWeight);	   
	      jetpt1_zveto.Fill(Jet1pt_,evtWeight);
	      jetpt2_zveto.Fill(Jet2pt_,evtWeight);
	      npv2.Fill(getInt(numPV),evtWeight);
	  
	      lljjmassalltags.Fill(higgsMass,evtWeight);

	      // end of common selection
	      // now jet classification according btag
	      
	      //Tag algorithm cut definition for signal
	      bool jet1_tagged_medium(false), jet1_tagged_loose(false);
	      bool jet2_tagged_medium(false), jet2_tagged_loose(false);

	      // MODIFY B-TAG CUT DEFINITION FOR BUMP SEARCH
	      jet1_tagged_medium = jet1tkhe_ > 2.0;
	      jet1_tagged_loose  = jet1tkhe_ > 2.0;	    
	      jet2_tagged_medium = jet2tkhe_ > 2.0;
	      jet2_tagged_loose  = jet2tkhe_ > 2.0;
	    
	      // IN THE MODIFIED CUT CASE DO NOT APPLY ANY SF

	      bool twoTagCategory  = ( jet1_tagged_medium && jet2_tagged_loose  )
		|| (jet1_tagged_loose  && jet2_tagged_medium );
	      bool oneTagCategory  = ( !twoTagCategory ) && ( jet1_tagged_loose || jet2_tagged_loose );
	      bool zeroTagCategory = ( !twoTagCategory ) && ( !oneTagCategory );
	      
	      // count the candidates per category	    
	      if(zeroTagCategory) 
		++CandBTag[0];
	      if(oneTagCategory) 
		++CandBTag[1];
	      if(twoTagCategory) 
		++CandBTag[2];
	    
	      // if 1-tag or 2 tags candidates
	      if(oneTagCategory || twoTagCategory){
		lljjmass1and2tags.Fill(higgsMass,evtWeight);
		zjjmass_1and2tags.Fill(zjjmass_,evtWeight);	   
		jetpt1_1and2tags.Fill(Jet1pt_,evtWeight);
		jetpt2_1and2tags.Fill(Jet2pt_,evtWeight);
		// dump event variables
		if(!isEventSelected){
		  npvfin.Fill(getInt(numPV),evtWeight);
		  isEventSelected = true;
		}
	      }
	      
	      // select only 1-tag candidates
	      if(oneTagCategory){
		lljjmass1tag.Fill(higgsMass,evtWeight);
		zjjmass_1tag.Fill(zjjmass_,evtWeight);	   
		jetpt1_1tag.Fill(Jet1pt_,evtWeight);
		jetpt2_1tag.Fill(Jet2pt_,evtWeight);
	      }
	      
	      // select only 2-tag candidates
	      if(twoTagCategory){
		pass[4] = true;
		++cand[4];
		// fill plots for 2-tag candidates selection
		zllpt.Fill(zllpt_,evtWeight);
		drjj.Fill(drjj_,evtWeight);
		metSignif.Fill(getInt(metSignificance),evtWeight);
		lljjmass2tags.Fill(higgsMass,evtWeight);
		zjjmass_2tags.Fill(zjjmass_,evtWeight);	   
		jetpt1_2tags.Fill(Jet1pt_,evtWeight);
		jetpt2_2tags.Fill(Jet2pt_,evtWeight);
	      } // if two tag category


	    } // end of anti-Zll selection
	    
	  } // end of jet pT selection
	  
	} // end of lepton quality selection

      } // candidate loop
      
    }  // if cands > 0 ...
    
    for(unsigned int k = 0; k < counts; ++k)
      if(pass[k]) ++count[k];
  } // event loop
  
  
  cout << endl << "--> Number of candidates/category after common selection <--" << endl;
  for (int kkk = 0; kkk<3; ++kkk)
    cout << kkk << "-Tags candidates = " << CandBTag[kkk] << endl;

  TH1F h_cuts("h_cuts", "ProgressiveCuts", 11, 0, 11.);
  
  cout<<count[1]<<"   "<<cand[1]<<endl;  
  
  TFile histos(outFile, "RECREATE");
  histos.cd();
  lljjmass.Write();  
  leptpt1.Write();  
  leptpt2.Write();  
  jetpt1.Write();  
  jetpt2.Write();  
  metSignif.Write();  
  zjjmass_zveto.Write();
  jetpt1_zveto.Write();  
  jetpt2_zveto.Write();  

  lljjmassallcands.Write();
  lljjmassalltags.Write();
  lljjmass1tag.Write();
  zjjmass_1tag.Write();
  jetpt1_1tag.Write();  
  jetpt2_1tag.Write();  

  lljjmass2tags.Write();
  zjjmass_2tags.Write();
  jetpt1_2tags.Write();  
  jetpt2_2tags.Write();  

  lljjmass1and2tags.Write();
  zjjmass_1and2tags.Write();
  jetpt1_1and2tags.Write();  
  jetpt2_1and2tags.Write();  

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
  zllmass_zveto.Write();
  db.Write();
  npv_woReweight.Write();
  npv.Write();
  npv1.Write();
  npv2.Write();
  npvfin.Write();

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

  return 0;
}

