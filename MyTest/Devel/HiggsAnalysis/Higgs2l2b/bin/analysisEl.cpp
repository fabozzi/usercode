
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

/*  first argument:  directory
    second argument: output fileName
*/

int main(int argc, char **argv) {
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  const char * path = argv[1];
  const char * outFile = argv[2];
  cout<<"filePath: "<<path<<endl;
  TFile *f = TFile::Open(path, "READ");
  TTree * events = (TTree*) f->Get("Events");
  unsigned int nEvents = events->GetEntries();
  cout << "events: " << nEvents << endl;
  BRANCH(elHiggsLeptDau1Pt);
  BRANCH(elHiggsLeptDau2Pt);
  BRANCH(elHiggsJetDau1Pt);
  BRANCH(elHiggsJetDau2Pt);
  BRANCH(elHiggsjjdr);
  BRANCH(elHiggszllPt);
  BRANCH(elHiggsMet);
  BRANCH(elHiggsMet2);
  BRANCH(elHiggsMass);
  BRANCH(elHiggszllMass);
  BRANCH(elHiggszjjMass);
  BRANCH(elHiggsJet1TKHE);
  BRANCH(elHiggsJet2TKHE);
  BRANCH(elHiggsMetSig);
  BRANCH(elHiggsEleDau1VBTF80CombID);
  BRANCH(elHiggsEleDau2VBTF80CombID);
  BRANCHINT(elHiggsEventNumber);
  BRANCHINT(elHiggsRunNumber);
  BRANCHINT(elHiggsLumiblock);

  TFile histos(outFile, "RECREATE");
  TH1F lljjmass("lljjmass", "m(lljj)", 100, 0, 1000);
  TH1F lljjmassCands("lljjmassCands", "m(lljj)", 100, 0, 1000);
  TH1F zllpt("zllpt", "p_{T}(Z_{ll})", 100, 0, 300);
  TH1F zllmass("zllmass", "m_{ll}", 100, 0, 200);
  TH1F zjjmass("zjjmass", "m_{jj}", 100, 0, 200);
  TH1F drjj("DRjj", "#Delta R_{jj}", 100, 0, 4);
  TH1F met("met", "MET", 100, 0, 140);
  TH1F met2("met2", "MET2", 100, 0, 140);
  TH1F metSig("metSig", "METSig", 100, 0, 15);
  TH1F btagJet1("TCHEJet1", "TCHE", 100, 0, 20);
  TH1F btagJet2("TCHEJet2", "TCHE", 100, 0, 20);
  const unsigned int counts = 8;
  vint count(counts, 0);
  vbool pass(counts, false);
  vint cand(counts, 0);
  float zllMassNominal = 91.187;
  ofstream fileout("selectedEvents.txt");
  cand[1]=0;
  for(unsigned int i = 0; i <nEvents; ++i) {
    if(i%100 == 0) progress(double(i)/double(nEvents));
    ++count[0];
    GETENTRY(elHiggsLeptDau1Pt,i);
    unsigned int cands = elHiggsLeptDau1Pt->product()->size();
    fill(pass.begin(), pass.end(), false);
    if(cands > 0) {
      pass[1] = true;
      GETENTRY(elHiggsLeptDau2Pt,i);
      GETENTRY(elHiggsEleDau1VBTF80CombID,i);
      GETENTRY(elHiggsEleDau2VBTF80CombID,i);

      float met_, met2_, metSig_, lljjmass_, zllpt_, drjj_, btagJet1_, btagJet2_, zllmass_, zllmassBest_, zjjmass_, Jet1pt_, Jet2pt_, SumPtBest_;
      Jet1pt_=0;
      Jet2pt_=0;
      SumPtBest_=0;
      lljjmass_=0;

      for(unsigned int j = 0; j < cands; ++j) {
	if( get(elHiggsLeptDau1Pt,j) > 20 && get(elHiggsLeptDau2Pt,j) > 20 && (get(elHiggsEleDau1VBTF80CombID,j)==7||get(elHiggsEleDau2VBTF80CombID,j)==7) ) {
	  pass[2] = true;
	  ++cand[2];
	  GETENTRY(elHiggszllMass,i);
	  GETENTRY(elHiggszjjMass,i);
	  GETENTRY(elHiggsJet1TKHE,i);
	  GETENTRY(elHiggsJet2TKHE,i);
	  GETENTRY(elHiggsMet,i);
	  GETENTRY(elHiggsMet2,i);
	  GETENTRY(elHiggsMetSig,i);
	  GETENTRY(elHiggsjjdr,i);
	  GETENTRY(elHiggszllPt,i);
	  GETENTRY(elHiggsjjdr,i);

	  drjj.Fill(get(elHiggsjjdr,j));
	  met.Fill(get(elHiggsMet,j));
	  met2.Fill(get(elHiggsMet2,j));
	  metSig.Fill(get(elHiggsMetSig,j));
	  btagJet1.Fill(get(elHiggsJet1TKHE,j));
	  btagJet2.Fill(get(elHiggsJet2TKHE,j));
	  zllpt.Fill(get(elHiggszllPt,j));
	  zllmass.Fill(get(elHiggszllMass,j));
	  zjjmass.Fill(get(elHiggszjjMass,j));


	  if(TMath::Abs(get(elHiggszllMass,j)-91.187) < 20 && TMath::Abs(get(elHiggszjjMass,j)-91.187)<15  ) {
	    pass[3] = true;
	    ++cand[3];	  

	    if((get(elHiggsJet1TKHE,j) >3.3 && get(elHiggsJet2TKHE,j) >1.7)||(get(elHiggsJet2TKHE,j) >3.3 && get(elHiggsJet1TKHE,j) >1.7)) {
	      pass[4] = true;
	      ++cand[4];
	      GETENTRY(elHiggsMass,i);
	      GETENTRY(elHiggsjjdr,i);
	      GETENTRY(elHiggszllPt,i);
	      zllpt.Fill(get(elHiggszllPt,j));
	      if(get(elHiggsjjdr,j) < (1.9 + (1./500.)*(250. - get(elHiggsMass,j)))) {
		pass[5] = true;
		++cand[5];
		GETENTRY(elHiggszllPt,i);
		if(get(elHiggszllPt,j) > (150. - (6./15.)*(500. + get(elHiggsMass,j)))) {
		  pass[6] = true;
		  ++cand[6];
		  GETENTRY(elHiggsMet,i);
		  GETENTRY(elHiggsMetSig,i);
		  GETENTRY(elHiggsEventNumber,i);
		  GETENTRY(elHiggsRunNumber,i);
		  GETENTRY(elHiggsLumiblock,i);
		  GETENTRY(elHiggsJetDau1Pt,i);
		  GETENTRY(elHiggsJetDau2Pt,i);
		  //cout<<endl;
		  //cout<<"mass: "<<get(elHiggsMass,j)<<endl;
		  //cout<<"met: "<<get(elHiggsMet,j)<<endl;
		  //cout<<"metSig: "<<get(elHiggsMetSig,j)<<endl;
		  //		  cout<<"eventNumber: "<<get(elHiggsEventNumber,j)<<endl;
		  //		  fileout<<get(elHiggsEventNumber,j)<<endl;


		  if(get(elHiggsMetSig,j) < 10) {

	
			Jet1pt_=get(elHiggsJetDau1Pt,j);
			Jet2pt_=get(elHiggsJetDau2Pt,j);
			if ( (Jet1pt_+ Jet2pt_)>SumPtBest_ ){
			  SumPtBest_= Jet1pt_+ Jet2pt_;
			  lljjmass_=get(elHiggsMass,j);
			}
		  
			//		    lljjmass.Fill(get(elHiggsMass,j));
		    /*cout<<"EventNumber "<<getInt(elHiggsEventNumber)<<endl;
		    cout<<"LumiNumber "<<getInt(elHiggsLumiblock)<<endl;
		    cout<<"RunNumber "<<getInt(elHiggsRunNumber)<<endl;
		    */pass[7] = true;
		    ++cand[7];


		  }
		}
	      }
	    }
	  }
	}
      } // candidate loop
      if(lljjmass_!=0)lljjmass.Fill(lljjmass_);
    }  // if cands > 0 ...
    for(unsigned int k = 0; k < counts; ++k)
    if(pass[k]) ++count[k];
  } // event loop

TH1F h_cuts("h_cuts", "ProgressiveCuts", 10, 0, 10);

  
 lljjmass.Write();
 lljjmassCands.Write();
 zllpt.Write();
 zjjmass.Write();
 zllmass.Write();
 drjj.Write();
 btagJet1.Write();
 btagJet2.Write();
 met.Write();
 met2.Write();
 metSig.Write();
 cout << endl;
 cout<<"basic   Pt/Eta  zll/zjj mass  btag   jjdr   zllpt   met"<<endl;
 string names[]={"basic", "lep/jet pt/Eta","Isolation", "zll/zjj mass", "btag", "jjdr", "zllpt", "met"};

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

  cout << "efficiencies: ";
  for(unsigned int k = 0; k < counts; ++k) {
    cout << double(count[k])/double(count[0]) << " ";
  }
  cout << endl;
  return 0;
}

