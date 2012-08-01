
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
  BRANCH(muHiggsLeptDau1Pt);
  BRANCH(muHiggsLeptDau2Pt);
  BRANCH(muHiggsJetDau1Pt);
  BRANCH(muHiggsJetDau2Pt);
  BRANCH(muHiggsLeptDau1Eta);
  BRANCH(muHiggsLeptDau2Eta);
  BRANCH(muHiggsLeptDau1dB);
  BRANCH(muHiggsLeptDau2dB);
  BRANCH(muHiggsjjdr);
  BRANCH(muHiggszllPt);
  BRANCH(muHiggsMet);
  BRANCH(muHiggsMass);
  BRANCH(muHiggszllMass);
  BRANCH(muHiggszjjMass);
  BRANCH(muHiggsJet1TKHE);
  BRANCH(muHiggsJet2TKHE);
  BRANCH(muHiggsLeptDau1TrkIso);
  BRANCH(muHiggsLeptDau2TrkIso);
  BRANCH(muHiggsLeptDau1EcalIso);
  BRANCH(muHiggsLeptDau2EcalIso);
  BRANCH(muHiggsLeptDau1HcalIso);
  BRANCH(muHiggsLeptDau2HcalIso);
  BRANCH(muHiggsMetSig);
  BRANCH(muHiggsMet2);
  BRANCH(muHiggsHelyLD);
  BRANCH(muHiggsRefitMass);
  BRANCHINT(muHiggsEventNumber);
  BRANCHINT(muHiggsRunNumber);
  BRANCHINT(muHiggsLumiblock);

  TFile histos(outFile, "RECREATE");
  TH1F lljjmass("lljjmass", "m(lljj)", 100, 0, 1000);
  TH1F lljjmass_noReFit("lljjmassNoRefit", "m(lljj)", 100, 0, 1000);
  TH1F lljjmassCands("lljjmassCands", "m(lljj)", 100, 0, 1000);
  TH1F zllpt("zllpt", "P_{T}(Z_{ll})", 100, 0, 300);
  TH1F zjjmass("zjjmass", "m_{jj}", 100, 0, 200);
  TH1F zllmass("zllmass", "m_{ll}", 100, 0, 200);
  TH1F drjj("DRjj", "#Delta R_{jj}", 100, 0, 4);
  TH1F met("met", "MET", 100, 0, 140);
  TH1F met2("met2", "MET2", 100, 0, 140);
  TH1F metSig("metSig", "METSig", 100, 0, 15);
  TH1F btagJet1("TCHEJet1", "TCHE", 100, 0, 20);
  TH1F btagJet2("TCHEJet2", "TCHE", 100, 0, 20);
  TH1F helyLD("HelyLD", "HelyLD", 100, 0, 1);
  const unsigned int counts = 11;
  vint count(counts, 0);
  vbool pass(counts, false);
  vint cand(counts, 0);
  float zllMassNominal = 91.187;
  ofstream fileout("SelectedEventsMu.txt");
  cand[1]=0;
  for(unsigned int i = 0; i <nEvents; ++i) {
    if(i%100 == 0) progress(double(i)/double(nEvents));
    ++count[0];
    GETENTRY(muHiggsLeptDau1Pt,i);
    unsigned int cands = muHiggsLeptDau1Pt->product()->size();
    fill(pass.begin(), pass.end(), false);
    if(cands > 0) {
      pass[1] = true;
      GETENTRY(muHiggsLeptDau2Pt,i);
      GETENTRY(muHiggsLeptDau1Eta,i);
      GETENTRY(muHiggsLeptDau2Eta,i);
      float met_, met2_, metSig_, lljjmass_, lljjmass_noReFit_, zllpt_, drjj_, btagJet1_, btagJet2_, zllmass_, zllmassBest_, zjjmass_, Jet1pt_,Jet2pt_, SumPtBest_;
      Jet1pt_=0;
      Jet2pt_=0;
      SumPtBest_=0;
      lljjmass_=0;
      lljjmass_noReFit_=0;

      for(unsigned int j = 0; j < cands; ++j) {
	++cand[1];
	GETENTRY(muHiggszllMass,i);
	GETENTRY(muHiggszjjMass,i);
	GETENTRY(muHiggszllPt,i);

	if((get(muHiggsLeptDau1Pt,j) > 40 && get(muHiggsLeptDau2Pt,j) > 20 ) || (get(muHiggsLeptDau1Pt,j) > 20 && get(muHiggsLeptDau2Pt,j) > 40 ) ) {
	  pass[2] = true;
	  ++cand[2];
	  GETENTRY(muHiggsLeptDau1dB,i);
	  GETENTRY(muHiggsLeptDau2dB,i);

	  if(fabs(get(muHiggsLeptDau1dB,j)) <0.02 && fabs(get(muHiggsLeptDau2dB,j)) <0.02) {
	    pass[3] = true;
	    ++cand[3];

	    GETENTRY(muHiggsLeptDau1TrkIso,i);
	    GETENTRY(muHiggsLeptDau2TrkIso,i);
	    GETENTRY(muHiggsLeptDau1EcalIso,i);
	    GETENTRY(muHiggsLeptDau2EcalIso,i);
	    GETENTRY(muHiggsLeptDau1HcalIso,i);
	    GETENTRY(muHiggsLeptDau2HcalIso,i);
	    if( ((get(muHiggsLeptDau1TrkIso,j) + get(muHiggsLeptDau1EcalIso,j) + get(muHiggsLeptDau1HcalIso,j)) < 0.15*get(muHiggsLeptDau1Pt,j) ) && ((get(muHiggsLeptDau2TrkIso,j) + get(muHiggsLeptDau2EcalIso,j) + get(muHiggsLeptDau2HcalIso,j)) < 0.15*get(muHiggsLeptDau2Pt,j) )){
              pass[4] = true;
	      ++cand[4];	  

	      GETENTRY(muHiggsJet1TKHE,i);
	      GETENTRY(muHiggsJet2TKHE,i);
	      GETENTRY(muHiggsjjdr,i);
	      GETENTRY(muHiggsMet,i);
	      GETENTRY(muHiggsMet2,i);
	      GETENTRY(muHiggsMetSig,i);

	      zllpt.Fill(get(muHiggszllPt,j));
	      zllmass.Fill(get(muHiggszllMass,j));
	      zjjmass.Fill(get(muHiggszjjMass,j));	   
	      drjj.Fill(get(muHiggsjjdr,j));
	      met.Fill(get(muHiggsMet,j));
	      met2.Fill(get(muHiggsMet2,j));
	      metSig.Fill(get(muHiggsMetSig,j));
	      btagJet1.Fill(get(muHiggsJet1TKHE,j));
	      btagJet2.Fill(get(muHiggsJet2TKHE,j));
	
	      //	      cout<<"zllmass"<<get(muHiggszllMass,j)<<endl;
	      if(TMath::Abs(get(muHiggszllMass,j)-91.187) < 20 && TMath::Abs(get(muHiggszjjMass,j)-91.187)<15  ) {
		pass[5] = true;
		++cand[5];
		GETENTRY(muHiggsHelyLD,i);	  
		helyLD.Fill(get(muHiggsHelyLD,j));

		if((get(muHiggsJet1TKHE,j) >3.3 && get(muHiggsJet2TKHE,j) >1.7)||(get(muHiggsJet2TKHE,j) >3.3 && get(muHiggsJet1TKHE,j) >1.7)) {
		  pass[6] = true;
		  ++cand[6];
		  GETENTRY(muHiggsMass,i);
		  GETENTRY(muHiggsjjdr,i);
		  GETENTRY(muHiggszllPt,i);

		  if(get(muHiggsjjdr,j) < (1.9 + (1./500.)*(250. - get(muHiggsMass,j)))) {
		    pass[7] = true;
		    ++cand[7];
		    //		    GETENTRY(muHiggszllPt,i);
		    if(get(muHiggszllPt,j) > (150. - (6./15.)*(500. + get(muHiggsMass,j)))) {
		      pass[8] = true;
		      ++cand[8];
		     
		      if(get(muHiggsHelyLD,j)>0.5){
			pass[9] = true;
			++cand[9];
			GETENTRY(muHiggsMet,i);
			GETENTRY(muHiggsMetSig,i);
			GETENTRY(muHiggsJetDau1Pt,i);
			GETENTRY(muHiggsJetDau2Pt,i);
			GETENTRY(muHiggsEventNumber,i);
			GETENTRY(muHiggsLumiblock,i);
			GETENTRY(muHiggsRunNumber,i);
			if(get(muHiggsMetSig,j)<10){
			  //cout<<endl;
			  //cout<<"mass: "<<get(muHiggsMass,j)<<endl;
			  //cout<<"met: "<<get(muHiggsMet,j)<<endl;
			  //cout<<"metSig: "<<get(muHiggsMetSig,j)<<endl;
			  //cout<<"eventNumber: "<<getInt(muHiggsEventNumber)<<endl;
			  fileout<<getInt(muHiggsEventNumber)<< endl;
			  //cout<<"runNumber: "<<get(muHiggsRunNumber,j)<<endl;
			  GETENTRY(muHiggsRefitMass,i);
			  GETENTRY(muHiggsMass,i);
			  Jet1pt_=get(muHiggsJetDau1Pt,j);
			  Jet2pt_=get(muHiggsJetDau2Pt,j);
			  if ( (Jet1pt_+ Jet2pt_)>SumPtBest_ && get(muHiggsRefitMass,j)>0 ){
			    SumPtBest_= Jet1pt_+ Jet2pt_;
			    lljjmass_=get(muHiggsRefitMass,j);
			    lljjmass_noReFit_=get(muHiggsMass,j);
			  }
			  
			  //lljjmass.Fill(get(muHiggsMass,j));
			  /*cout<<"met "<<get(muHiggsMet,j)<<endl;
			    cout<<"metSig "<<get(muHiggsMetSig,j)<<endl;
			    cout<<"mlljj "<<get(muHiggsMass,j)<<endl;
			    cout<<"EventNumber "<<getInt(muHiggsEventNumber)<<endl;
			    cout<<"LumiNumber "<<getInt(muHiggsLumiblock)<<endl;
			    cout<<"RunNumber "<<getInt(muHiggsRunNumber)<<endl;
			    cout<<"zllmass"<<get(muHiggszllMass,j)<<endl;
			    cout<<"zjjmass"<<get(muHiggszjjMass,j)<<endl;
			  */
			  pass[10] = true;
			  ++cand[10];
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      } // candidate loop

      if(lljjmass_!=0)  lljjmass.Fill(lljjmass_);
      if(lljjmass_noReFit_!=0)  lljjmass_noReFit.Fill(lljjmass_noReFit_);
    }  // if cands > 0 ...

    for(unsigned int k = 0; k < counts; ++k)
      if(pass[k]) ++count[k];
  } // event loop
  
  TH1F h_cuts("h_cuts", "ProgressiveCuts", 11, 0, 11.);
  
  
  lljjmass.Write();  
  lljjmass_noReFit.Write();  
  //  lljjmassCands.Write();  
  met.Write();  
  met2.Write();  
  metSig.Write();  
  btagJet1.Write();  
  btagJet2.Write();  
  drjj.Write();  
  zllpt.Write();
  zjjmass.Write();
  zllmass.Write();
  helyLD.Write();
  cout << endl;
  cout<<"#TotEvts basicSel   Pt/Eta   dB  Iso  zll/zjjMmass   btag   jjdr   zllpt HelyLD  metSig"<<endl;
  string names[]={"#TotEvents","basicSel", "lep/jet pt/Eta", "dB","Isolation", "zll/zjjMass", "btag", "jjdr", "zllpt", "HelyLD", "met"};

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

