#ifndef LeptonID_h
#define LeptonID_h

//#include "TF1.h"
//#include "TLegend.h"
//#include <TCanvas.h>
//#include "HiggsAnalysis/Higgs2l2b/bin/2l2b_branches_2012.h"
#include <string>
#include <cmath>

using namespace std;

class LeptonID{

 public:

  LeptonID(string period = "2012A");
  ~LeptonID();

  bool electronID2012(float dxy, float dz, float trkIso03, float ecalIso03, float hcalIso03, float pt, float combRelIso,  float dEtaVtx, float  dPhiVtx, float sigmaee, float HoverE, float EcalEn, float Epin, float mHits, bool hasMatchConv, bool isEB);

  bool electronID2011(float  dPhiVtx, float HoverE, float vbtf95, float etaSC, float pt, float combRelIso);

  bool electronIDMVA(float  mvaID,  bool isEB);
  
  bool muonID2012(bool isGlobal, bool isTrackMu, bool isPFMu, float normChi2, int NofMuHits, int NofMatchedStations, float dz, int NofPixelHits, int NofTrackerLayers, float db);
  bool muonID2011(bool isGlobal, bool isTrackMu, float normChi2, int NofTkrHits, int NofPixelHits, int NofMuHits, int NofMatches, float db);

 private: 


  
  string period_;
  
};
#endif
