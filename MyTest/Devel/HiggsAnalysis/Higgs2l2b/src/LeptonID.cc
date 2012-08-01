
//#include "HiggsAnalysis/Higgs2l2b/interface/SFlightFuncs.h"

#include "../interface/LeptonID.h"


LeptonID::LeptonID( string period ){
  period_ = period;
}

LeptonID::~LeptonID(){}

bool LeptonID::electronID2012(float dxy, float dz, float trkIso03, float ecalIso03, float hcalIso03, float pt, float combRelIso,  float dEtaVtx, float  dPhiVtx, float sigmaee, float HoverE, float EcalEn, float EPin, float mHits, bool hasMatchConv, bool isEB){
 
  bool EleID(false);
  
  
  /* Apply Cut-Based ID, Veto WP + Tight Trigger, Default*/
 /* Apply Cut-Based ID, Loose WP*/
  float ooemoop = (1.0/EcalEn - EPin/EcalEn);
  EleID =(fabs(dxy) < 0.02) && 
    (fabs(dz) < 0.2) && 
    (fabs(ooemoop) < 0.05) &&
    ((trkIso03 /pt )< 0.2) &&
    ((ecalIso03/pt )< 0.2) &&
    ((hcalIso03/pt) < 0.2)&&
    (!hasMatchConv) &&
    (mHits <= 1) &&
    (combRelIso < 0.15);    
  
  if(isEB){
    
    EleID = EleID && (fabs (dEtaVtx) < 0.007) &&
      (fabs(dPhiVtx)< 0.15) &&
      (sigmaee < 0.01) && 
      (HoverE < 0.12); 
  } 
    else{
      EleID = EleID && (fabs (dEtaVtx) < 0.009) &&
	(fabs(dPhiVtx)< 0.10) &&
	(sigmaee < 0.03) && 
	(HoverE < 0.10); 
    }
  
    
  
  return EleID;
}


bool LeptonID::electronID2011(float  dPhiVtx, float HoverE, float vbtf95, float etaSC, float pt, float combRelIso){
 
  bool EleID(false);
  bool VBTFID = (vbtf95 == 7 || vbtf95 == 5 ) ;
  bool isEleBarrel = fabs(etaSC) <= 1.4442;
  
  EleID =  VBTFID;
    
  if(isEleBarrel)EleID = EleID && (fabs(dPhiVtx )< 0.15) && (combRelIso < 0.15);
  else EleID = EleID && (fabs(dPhiVtx)< 0.1) && HoverE < 0.10 && (combRelIso < 0.15);
  return EleID;

}


bool LeptonID::electronIDMVA(float  mvaID,  bool isEB ){
 
  bool EleID(false);
  if(isEB)  EleID = (mvaID > 0.94);
  else  EleID = (mvaID > 0.91);
  return EleID;

}



bool LeptonID::muonID2012(bool isGlobal, bool isTrackMu, bool isPFMu, float normChi2, int NofMuHits, int NofMatchedStations, float dz, int NofPixelHits, int NofTrackerLayers, float db){ 

  bool MuID(false);

  MuID = (isGlobal) && 
    (isTrackMu) &&
    (isPFMu) && 
    (normChi2< 10) && 
    (NofMuHits > 0) && 
    (NofMatchedStations> 1) && 
    (fabs(dz) < 0.5) &&
    (NofPixelHits > 0) &&   
    (NofTrackerLayers > 5) &&
    (fabs(db) <0.2);

  return MuID;}

bool LeptonID::muonID2011(bool isGlobal, bool isTrackMu, float normChi2, int NofTrackerHits, int NofPixelHits, int NofMuHits, int NofMatches, float db){ 

  bool MuID(false);

  MuID = (isGlobal) && 
    (isTrackMu) &&
    (normChi2< 10) && 
    (NofTrackerHits > 10) &&
    (NofPixelHits > 0) &&   
    (NofMuHits > 0) && 
    (NofMatches > 1) &&
    (fabs(db) <0.2);
 
  return MuID;}

