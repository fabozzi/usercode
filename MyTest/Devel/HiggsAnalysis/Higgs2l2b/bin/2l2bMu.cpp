
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
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
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

// macro usage:
// mu2l2b [input_ntuple_file] [output_histo_file] [DATA/MC] [SF/noSF]
// if MC --> apply PUreweight
// if SF --> apply b-tag scale factor

int main(int argc, char **argv) {
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  const char * path = argv[1];
  const char * outFile = argv[2];
  string isData(argv[3]);
  string applybTagSF(argv[4]);
  string outFileName(outFile);
  bool data = false;
  bool btagScale = false;
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
  
  TFile *f = TFile::Open(path, "READ");
  TTree * events = (TTree*) f->Get("Events");
  unsigned int nEvents = events->GetEntries();
  cout << "events: " << nEvents << endl;

  // GET BRANCHES
  BRANCH(muHiggsLeptDau1Pt);
  BRANCH(muHiggsLeptDau2Pt);
  BRANCH(muHiggsJetDau1Pt);
  BRANCH(muHiggsJetDau2Pt);
  BRANCH(muHiggsLeptDau1Eta);
  BRANCH(muHiggsLeptDau2Eta);
  BRANCH(muHiggsLeptDau1dB);
  BRANCH(muHiggsLeptDau2dB);
  BRANCH(muHiggsJetDau1Eta);
  BRANCH(muHiggsJetDau2Eta);
  BRANCH(muHiggsjjdr);
  BRANCH(muHiggszllPt);
  BRANCH(muHiggsMass);
  BRANCH(muHiggszllMass);
  BRANCH(muHiggszjjMass);
  BRANCH(muHiggsJet1TKHE);
  BRANCH(muHiggsJet2TKHE);
  BRANCH(muHiggsJetDau1PartonFlavour);
  BRANCH(muHiggsJetDau2PartonFlavour);
  BRANCH(muHiggsLeptDau1CombRelIso);
  BRANCH(muHiggsLeptDau2CombRelIso);
  //  BRANCH(muHiggsLeptDau1TrkIso);
  //  BRANCH(muHiggsLeptDau2TrkIso);
  //  BRANCH(muHiggsLeptDau1EcalIso);
  //  BRANCH(muHiggsLeptDau2EcalIso);
  //  BRANCH(muHiggsLeptDau1HcalIso);
  //  BRANCH(muHiggsLeptDau2HcalIso);
  BRANCHFLOAT(met);
  BRANCHFLOAT(metSignificance);
  BRANCH(muHiggsMet2);

  /* Angular variables  */
  BRANCH(muHiggscosthetaNT1);
  BRANCH(muHiggscosthetaNT2);
  BRANCH(muHiggscosthetastarNT);
  BRANCH(muHiggsphiNT);
  BRANCH(muHiggsphiNT1);
  BRANCH(muHiggsHelyLD);

  BRANCH(muHiggscosthetaNT1Refit);
  BRANCH(muHiggscosthetaNT2Refit);
  BRANCH(muHiggscosthetastarNTRefit);
  BRANCH(muHiggsphiNTRefit);
  BRANCH(muHiggsphiNT1Refit);
  BRANCH(muHiggsHelyLDRefit);

  BRANCH(muHiggsRefitMass);
  BRANCHINT(muHiggsEventNumber);
  BRANCHINT(muHiggsRunNumber);
  BRANCHINT(muHiggsLumiblock);
  BRANCHINT(numPV);

  // txt output file
  ofstream fileout("SelectedEventsMu.txt",ios::app);
  fileout<<"Muon channel "<<endl;
  fileout<<"EvtNumber  RunNumber  mlljj"<<endl<<endl;

  // TTree output file
  TFile *outTree = TFile::Open(("tree_"+outFileName).c_str(), "RECREATE");
  TTree lljjmassTree("lljjmassTree","lljjmass tree");
  // mzz branch for the output tree
  double mzz;
  lljjmassTree.Branch("mzz", &mzz, "mzz/D");

  // define histograms
  TH1F lljjmass("lljjmass", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmass_noReFit("lljjmassNoRefit", "m(lljj)", 1000, 0, 1000);
  TH1F lljjmassCands("lljjmassCands", "m(lljj)", 1000, 0, 1000);
  TH1F zllpt("zllpt", "P_{T}(Z_{ll})", 300, 0, 300);
  TH1F leptpt1("leptpt1", "P_{T} ", 300, 0, 300);
  TH1F leptpt2("leptpt2", "P_{T} ", 300, 0, 300);

  TH1F zjjmass("zjjmass", "m_{jj}", 500, 0, 500);
  TH1F zllmass("zllmass", "m_{ll}", 500, 0, 500);
  TH1F drjj("DRjj", "#Delta R_{jj}", 100, 0, 4);
  TH1F h_met("met", "MET", 100, 0, 140);
  TH1F met2("met2", "MET2", 100, 0, 140);
  TH1F metSig("metSig", "METSig", 100, 0, 15);
  TH1F btagJet1("TCHEJet1", "TCHE", 100, 0, 20);
  TH1F btagJet2("TCHEJet2", "TCHE", 100, 0, 20);

  TH1F helyLD("HelyLD", "HelyLD", 100, 0, 1);
  TH1F cosT1("cosTheta1", "cosTheta1", 100, -1, 1);
  TH1F cosT2("cosTheta2", "cosTheta2", 100, -1, 1);
  TH1F cosT1Star("cosTheta1Star", "cosTheta1Star", 100, -1, 1);
  TH1F phi("phi", "phi", 100, -3, 3);
  TH1F phiStar("phiStar", "phiStar", 100, -3, 3);

  TH1F helyLD_RF("HelyLDRefit", "HelyLDRefit", 100, 0, 1);
  TH1F cosT1_RF("cosTheta1Refit", "cosTheta1Refit", 100, -1, 1);
  TH1F cosT2_RF("cosTheta2Refit", "cosTheta2Refit", 100, -1, 1);
  TH1F cosT1Star_RF("cosTheta1StarRefit", "cosTheta1StarRefit", 100, -1, 1);
  TH1F phi_RF("phiRefit", "phiRefit", 100, -3, 3);
  TH1F phiStar_RF("phiStarRefit", "phiStarRefit", 100, -3, 3);

  TH1F db("db", "db", 100, 0, 0.2);
  TH1F npv_woReweight("npv_woReweight", "npvwoReweight", 100, 0, 25);
  TH1F npv("npv", "npv", 100, 0, 25);

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
  TFile* file_loose = TFile::Open("BTagPayloads_TCHEL.root");
  TFile* file_medium = TFile::Open("BTagPayloads_TCHEM.root");
  btsfutil->set_fileLoose(file_loose);
  btsfutil->set_fileMedium(file_medium);

  //  float zllMassNominal = 91.187;

  cand[1]=0;

  // PU reweighting util
  edm::LumiReWeighting LumiWeights_;
  if(!data){
    LumiWeights_ = edm::LumiReWeighting("PUDist_MCClone.root", "PUDist_Data.root", "MC_clone", "pileup");
  }

  // LOOP on EVENTS
  for(unsigned int i = 0; i <nEvents; ++i) {
    //    if(i%100 == 0) progress(double(i)/double(nEvents));
    ++count[0];
    // GET THE NUMBER OF CANDIDATES IN THE EVENT
    GETENTRY(muHiggsLeptDau1Pt,i);
    unsigned int cands = muHiggsLeptDau1Pt->product()->size();
    fill(pass.begin(), pass.end(), false);
    if(cands > 0) {
      pass[1] = true;

      // GET ENTRIES FOR VARIABLES IN THE EVENT
      GETENTRY(muHiggsLeptDau1CombRelIso,i);
      GETENTRY(muHiggsLeptDau2CombRelIso,i);
      //      GETENTRY(muHiggsLeptDau1TrkIso,i);
      //      GETENTRY(muHiggsLeptDau2TrkIso,i);
      //      GETENTRY(muHiggsLeptDau1EcalIso,i);
      //      GETENTRY(muHiggsLeptDau2EcalIso,i);
      //      GETENTRY(muHiggsLeptDau1HcalIso,i);
      //      GETENTRY(muHiggsLeptDau2HcalIso,i);
      GETENTRY(muHiggsLeptDau2Pt,i);	  
      GETENTRY(muHiggsLeptDau1dB,i);
      GETENTRY(muHiggsLeptDau2dB,i);
      GETENTRY(muHiggszllMass,i);
      GETENTRY(muHiggszjjMass,i);
      GETENTRY(muHiggszllPt,i);
      GETENTRY(muHiggsJet1TKHE,i);
      GETENTRY(muHiggsJet2TKHE,i);
      GETENTRY(muHiggsJetDau1PartonFlavour,i);
      GETENTRY(muHiggsJetDau2PartonFlavour,i);
      GETENTRY(muHiggsjjdr,i);
      GETENTRY(met,i);
      GETENTRY(muHiggsMet2,i);
      GETENTRY(metSignificance,i);
      
      GETENTRY(muHiggsHelyLD,i);	  
      GETENTRY(muHiggscosthetaNT1,i);
      GETENTRY(muHiggscosthetaNT2,i);
      GETENTRY(muHiggscosthetastarNT,i);
      GETENTRY(muHiggsphiNT,i);
      GETENTRY(muHiggsphiNT1,i);
      
      GETENTRY(muHiggsHelyLDRefit,i);	  
      GETENTRY(muHiggscosthetaNT1Refit,i);
      GETENTRY(muHiggscosthetaNT2Refit,i);
      GETENTRY(muHiggscosthetastarNTRefit,i);
      GETENTRY(muHiggsphiNTRefit,i);
      GETENTRY(muHiggsphiNT1Refit,i);
      
      GETENTRY(muHiggsJetDau1Pt,i);
      GETENTRY(muHiggsJetDau2Pt,i);
      GETENTRY(muHiggsJetDau1Eta,i);
      GETENTRY(muHiggsJetDau2Eta,i);
      GETENTRY(muHiggsEventNumber,i);
      GETENTRY(muHiggsLumiblock,i);
      GETENTRY(muHiggsRunNumber,i);
      GETENTRY(muHiggsRefitMass,i);
      GETENTRY(muHiggsMass,i);

      // GETENTRY(muHiggsLeptDau1Eta,i);
      // GETENTRY(muHiggsLeptDau2Eta,i);

      GETENTRY(numPV,i);
      double evtWeight = 1;

      // assign event weight for PU reweighting
      if(!data){
	BRANCHINT(nGenInt);
	GETENTRY(nGenInt,i);
	evtWeight = LumiWeights_.weight(getInt(nGenInt));
      }      

      // define variables of the best H candidate in the event
      float met_, met2_, metSig_;
      float lljjmass_(0), lljjmass_noReFit_(0);
      float zllpt_, drjj_;
      float btagJet1_, btagJet2_;
      float zllmass_, zllmassBest_;
      float zjjmass_; 
      float Jet1pt_(0), Jet2pt_(0), SumPtBest_(0), SumPtBest_Final(0);
      float lept1pt_(0), lept2pt_(0);
      float hLD_(-10), cos1_(-10), cos2_(-10), cosStar_(-10), phi_(-10), phiStar_(-10);
      float  hLD_rf(-10), cos1_rf(-10), cos2_rf(-10), cosStar_rf(-10), phi_rf(-10), phiStar_rf(-10);
      int EvtNum_(0), RunNum_(0);


      // LOOP on CANDIDATES IN THE EVENT
      for(unsigned int j = 0; j < cands; ++j) {
	++cand[1];

	// fill histo of db of the muon candidates
	db.Fill(get(muHiggsLeptDau1dB,j),evtWeight);
	db.Fill(get(muHiggsLeptDau2dB,j),evtWeight);

	// fill histo of numPV (reweighted and not)
	npv.Fill(getInt(numPV),evtWeight);
	npv_woReweight.Fill(getInt(numPV));

 	double isoDau1 = get(muHiggsLeptDau1CombRelIso,j);
	double isoDau2 = get(muHiggsLeptDau2CombRelIso,j);
	//	double isoDau1 = (get(muHiggsLeptDau1TrkIso,j) + get(muHiggsLeptDau1EcalIso,j) + get(muHiggsLeptDau1HcalIso,j))/get(muHiggsLeptDau1Pt,j); 
	//	double isoDau2 = (get(muHiggsLeptDau2TrkIso,j) + get(muHiggsLeptDau2EcalIso,j) + get(muHiggsLeptDau2HcalIso,j))/get(muHiggsLeptDau2Pt,j); 

	Jet1pt_=get(muHiggsJetDau1Pt,j);
	Jet2pt_=get(muHiggsJetDau2Pt,j);
	lept1pt_=get(muHiggsLeptDau1Pt,j);
	lept2pt_=get(muHiggsLeptDau2Pt,j);

	//	if( ((get(muHiggsLeptDau1Pt,j) > 40 && get(muHiggsLeptDau2Pt,j) > 20 ) || (get(muHiggsLeptDau1Pt,j) >20 && get(muHiggsLeptDau2Pt,j) > 40 )) && 
	//	    (fabs(get(muHiggsLeptDau1dB,j)) <0.02 && fabs(get(muHiggsLeptDau2dB,j)) <0.02)  ) {

	// muon kine && isolation cut && cosmic rejection
	bool kineLepCut = (get(muHiggsLeptDau1Pt,j) > 40 && get(muHiggsLeptDau2Pt,j) > 20) || (get(muHiggsLeptDau1Pt,j) >20 && get(muHiggsLeptDau2Pt,j) > 40 );
	bool dbLepCut = (fabs(get(muHiggsLeptDau1dB,j)) <0.02) && (fabs(get(muHiggsLeptDau2dB,j)) <0.02);
	bool isoLepCut = (get(muHiggsLeptDau1CombRelIso,j) < 0.15) && (get(muHiggsLeptDau2CombRelIso,j) < 0.15);
	
	if(kineLepCut && dbLepCut && isoLepCut) {
	  pass[2] = true;
	  ++cand[2];

	  // fill plots after muon selection
	  zllpt.Fill(get(muHiggszllPt,j),evtWeight);
	  zllmass.Fill(get(muHiggszllMass,j),evtWeight);
	  zjjmass.Fill(get(muHiggszjjMass,j),evtWeight);	   
	  drjj.Fill(get(muHiggsjjdr,j),evtWeight);
	  h_met.Fill(getInt(met),evtWeight);
	  met2.Fill(get(muHiggsMet2,j),evtWeight);
	  metSig.Fill(getInt(metSignificance),evtWeight);
	  btagJet1.Fill(get(muHiggsJet1TKHE,j),evtWeight);
	  btagJet2.Fill(get(muHiggsJet2TKHE,j),evtWeight);
	  leptpt1.Fill(lept1pt_,evtWeight);
	  leptpt2.Fill(lept2pt_,evtWeight);
	  
	  // Z inv. mass selection	  
	  if(get(muHiggszllMass,j)>70 &&  get(muHiggszllMass,j)<110  && get(muHiggszjjMass,j)>75 && get(muHiggszjjMass,j)<105 ) {
	    pass[3] = true;
	    ++cand[3];
	    
	    // end of common selection
	    // now jet classification according btag
	      
	    bool jet1_tagged_medium = get(muHiggsJet1TKHE,j)>3.3;
	    bool jet1_tagged_loose  = get(muHiggsJet1TKHE,j)>1.7;
	    
	    bool jet2_tagged_medium = get(muHiggsJet2TKHE,j)>3.3;
	    bool jet2_tagged_loose  = get(muHiggsJet2TKHE,j)>1.7;
	    
	    // eventually apply SF for MC
	    if( btagScale ) {
	      btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, get(muHiggsJetDau1Pt,j), get(muHiggsJetDau1Eta,j), (int) get(muHiggsJetDau1PartonFlavour,j) );
	      btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, get(muHiggsJetDau2Pt,j), get(muHiggsJetDau2Eta,j), (int) get(muHiggsJetDau2PartonFlavour,j) );
	    }
	    
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

	    // select only 2-tag candidates
	    if(twoTagCategory){
	      pass[4] = true;
	      ++cand[4];
	      
	      // apply selection for 2-tag candidates
	      // met significance cut
	    if(getInt(metSignificance)<10){
		pass[5] = true;
		++cand[5];
		
		// check if it is the best candidate
		if ( (Jet1pt_+ Jet2pt_)>SumPtBest_ && get(muHiggsRefitMass,j)>0 ){
		  // if best candidate store its variables
		  SumPtBest_= Jet1pt_+ Jet2pt_;
		  cos1_=get(muHiggscosthetaNT1,j);
		  cos2_=get(muHiggscosthetaNT2,j);
		  cosStar_= get(muHiggscosthetastarNT,j);
		  phi_= get(muHiggsphiNT,j);
		  phiStar_= get(muHiggsphiNT1,j);
		  hLD_= get(muHiggsHelyLD,j);
		  
		  cos1_rf=get(muHiggscosthetaNT1,j);
		  cos2_rf=get(muHiggscosthetaNT2,j);
		  cosStar_rf= get(muHiggscosthetastarNT,j);
		  phi_rf= get(muHiggsphiNT,j);
		  phiStar_rf= get(muHiggsphiNT1,j);
		  hLD_rf= get(muHiggsHelyLD,j);
		}
	     
		// LD cut
		if(get(muHiggsHelyLD,j)>0.5){

		  // check if it is the best candidate at this level
		  if ( (Jet1pt_+ Jet2pt_)>SumPtBest_Final && get(muHiggsRefitMass,j)>0 ){
		    // if best candidate store lljj mass
		    SumPtBest_Final= Jet1pt_+ Jet2pt_;
		    lljjmass_=get(muHiggsRefitMass,j);
		    lljjmass_noReFit_=get(muHiggsMass,j);
		    EvtNum_=getInt(muHiggsEventNumber);
		    RunNum_=getInt(muHiggsRunNumber);
		  }

		  pass[6] = true;
		  ++cand[6];

		}
	      }
	    }
	  }
	}
      } // candidate loop
      
      if(lljjmass_!=0) { 
	// fill best candidates plot
	lljjmass.Fill(lljjmass_, evtWeight);
	// dump the evt number for the best cand
	fileout<<EvtNum_<<"  "<<RunNum_<<"   "<<lljjmass_<<endl;
	// fill the tree with best cand mass
	mzz=lljjmass_;
	lljjmassTree.Fill();
      }
      
      if(lljjmass_noReFit_!=0)  lljjmass_noReFit.Fill(lljjmass_noReFit_,evtWeight);
      if(hLD_!=-10) helyLD.Fill(hLD_,evtWeight);
      if(cos1_!=-10) cosT1.Fill(cos1_,evtWeight);
      if(cos2_!=-10) cosT2.Fill(cos2_,evtWeight);
      if(cosStar_!=-10) cosT1Star.Fill(cosStar_,evtWeight);
      if(phi_!=-10) phi.Fill(phi_,evtWeight);
      if(phiStar_!=-10) phiStar.Fill(phiStar_,evtWeight);
 
      if(hLD_rf!=-10) helyLD_RF.Fill(hLD_rf,evtWeight);
      if(cos1_rf!=-10) cosT1_RF.Fill(cos1_rf,evtWeight);
      if(cos2_rf!=-10) cosT2_RF.Fill(cos2_rf,evtWeight);
      if(cosStar_rf!=-10) cosT1Star_RF.Fill(cosStar_rf,evtWeight);
      if(phi_rf!=-10) phi_RF.Fill(phi_rf,evtWeight);
      if(phiStar_rf!=-10) phiStar_RF.Fill(phiStar_rf,evtWeight);
   }  // if cands > 0 ...
    
    for(unsigned int k = 0; k < counts; ++k)
      if(pass[k]) ++count[k];
  } // event loop

 
  //--- write tree--
  //  outTree->cd();
  lljjmassTree.Write();
  outTree->Close();
  //----------------


  cout << endl << "--> Number of candidates/category after common selection <--" << endl;
  for (int kkk = 0; kkk<3; ++kkk)
    cout << kkk << "-Tags candidates = " << CandBTag[kkk] << endl;

  TH1F h_cuts("h_cuts", "ProgressiveCuts", 11, 0, 11.);

  cout<<count[1]<<"   "<<cand[1]<<endl;  

  // output histogram file
  TFile histos(outFile, "RECREATE");

  lljjmass.Write();  
  lljjmass_noReFit.Write();  
  //  lljjmassCands.Write();  
  leptpt1.Write();  
  leptpt2.Write();  
  h_met.Write();  
  met2.Write();  
  metSig.Write();  
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
  drjj.Write();  
  zllpt.Write();
  zjjmass.Write();
  zllmass.Write();
  db.Write();
  npv.Write();
  npv_woReweight.Write();

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
  
  //    cout << "efficiencies: ";
  //for(unsigned int k = 0; k < counts; ++k) {
  //cout << double(count[k])/double(count[0]) << " ";
  //}
  //cout << endl;
  return 0;
}

