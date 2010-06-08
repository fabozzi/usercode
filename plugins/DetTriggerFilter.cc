#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include <string>
using namespace edm;
using namespace std;

class DetTriggerFilter : public edm::EDFilter {
public:
  DetTriggerFilter(const edm::ParameterSet& pset);
  ~DetTriggerFilter(){};
  virtual bool filter( edm::Event& event, const edm::EventSetup& setup );
private:
  InputTag gtLbl;

  int bxLimit;

};

DetTriggerFilter::DetTriggerFilter(const edm::ParameterSet& pset):
  gtLbl(pset.getParameter<InputTag>("gtTag")),
  bxLimit(pset.getParameter<int>("absMaxBX"))
{
}

bool
DetTriggerFilter::filter(edm::Event& event, const edm::EventSetup& setup){

  Handle<L1MuGMTReadoutCollection> gmtCol;
  event.getByLabel(gtLbl,gmtCol);
  
  vector<L1MuGMTReadoutRecord> gmtRecs(gmtCol->getRecords());
  vector<bool> bxBools;bxBools.reserve(gmtRecs.size());
  for(vector<L1MuGMTReadoutRecord>::const_iterator rIt=gmtRecs.begin();rIt!=gmtRecs.end();rIt++){

    if(abs(rIt->getBxInEvent())>bxLimit)continue;
    
    bool trBx(true);
    vector<L1MuGMTExtendedCand> gmtCds(rIt->getGMTCands());
    for(vector<L1MuGMTExtendedCand>::const_iterator gmtIt=gmtCds.begin();gmtIt!=gmtCds.end();gmtIt++){
      if(!gmtIt->isRPC())continue;
      trBx=false;
    }
    bxBools.push_back(trBx);
  }

  bool aEv;
  for(unsigned int bls=0;bls<bxBools.size()-1;bls++)
    aEv=(bxBools[bls] || bxBools[bls+1]);

  return aEv;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DetTriggerFilter);
