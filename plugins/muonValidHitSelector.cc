#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include <iterator>

using namespace edm;
using namespace std;
using namespace reco;


class MuonValidHitsSelector : public edm::EDFilter {
public:
  MuonValidHitsSelector(const edm::ParameterSet& pset);
  virtual bool filter( edm::Event& event, const edm::EventSetup& setup );
private:
  InputTag muons_;

};

MuonValidHitsSelector::MuonValidHitsSelector(const edm::ParameterSet& pset) : 
  muons_(pset.getParameter<InputTag>("src") )
{ 
}

bool MuonValidHitsSelector::filter(edm::Event& event, const edm::EventSetup& setup) {
  bool hasValidHits = false;
  Handle<CandidateView> muonCands;
  event.getByLabel(muons_, muonCands);
  
  ESHandle<MagneticField> theMGField;
  setup.get<IdealMagneticFieldRecord>().get(theMGField);
  
  for( unsigned int i = 0; i < muonCands->size(); i++ ) {
    const Candidate & goodMuon = (*muonCands)[i];
    const Muon * muCand = dynamic_cast<const Muon*>(&goodMuon);
    TrackRef staTrack = muCand->outerTrack();
    reco::TransientTrack muTrack(staTrack,&*theMGField); 
    for (trackingRecHit_iterator it = muTrack.recHitsBegin();  it != muTrack.recHitsEnd(); it++) {
      if ((*it)->isValid ()) { 
	hasValidHits = true;
	return hasValidHits;
      }
    }
  }
  return hasValidHits;
} 
   
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( MuonValidHitsSelector );
  
