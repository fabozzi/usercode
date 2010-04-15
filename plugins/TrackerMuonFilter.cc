#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include <string>
using namespace edm;
using namespace std;
using namespace reco;


class TrackerMuonFilter : public edm::EDFilter {
public:
  TrackerMuonFilter(const edm::ParameterSet& pset);
  virtual bool filter( edm::Event& event, const edm::EventSetup& setup );
private:
  InputTag muons_;
  string algoType_;
  muon::SelectionType seleMu_;
};

TrackerMuonFilter::TrackerMuonFilter(const edm::ParameterSet& pset) : 
  muons_(pset.getParameter<InputTag>("src") ),
  algoType_(pset.getParameter<string>("algo") )
{ 
}

bool TrackerMuonFilter::filter(edm::Event& event, const edm::EventSetup& setup) {
  bool isGood = false;
  if( algoType_ == "TMLastStationAngTight"){
    seleMu_ = muon::TMLastStationAngTight;
  }
  else{
    cout<< "UNKNOWN selection type!!!" << endl;
    return false;
  }
  Handle<CandidateView> muonCands;
  event.getByLabel(muons_, muonCands);
  for( unsigned int i = 0; i < muonCands->size(); i++ ) {
    const Candidate & goodMuon = (*muonCands)[i];
    const Muon * muCand = dynamic_cast<const Muon*>(&goodMuon);
    if(muon::isGoodMuon(*muCand, seleMu_)){
      isGood = true;
      return isGood;
    }
  }
  return isGood;
} 
   
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( TrackerMuonFilter );
  
