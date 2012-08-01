  // EVENT BRANCHES
  BRANCHFLOAT(met);
  BRANCHFLOAT(metSignificance);
  BRANCHINT(numPV);
  BRANCHFLOAT(rhoRestrictedEta);
  BRANCHFLOAT(rho);
  std::cout<<"debug 0"<<std::endl;
  // MUON BRANCHES
  BRANCHBOOL(passSingleMuTrig);
  BRANCHBOOL(passDoubleMuTrig);
  BRANCH(CleanJetPt);
  std::cout<<"debug 0"<<std::endl;
  BRANCH(muHiggsLeptDau1Pt);
  BRANCH(muHiggsLeptDau2Pt);
  BRANCH(muHiggsJetDau1Pt);
  BRANCH(muHiggsJetDau1RefitPt);
  BRANCH(muHiggsJetDau2Pt);
  BRANCH(muHiggsJetDau2RefitPt);
  BRANCH(muHiggsLeptDau1Eta);
  BRANCH(muHiggsLeptDau2Eta);
  BRANCH(muHiggsLeptDau1Phi);
  BRANCH(muHiggsLeptDau2Phi);
  BRANCH(muHiggsLeptDau1dB);
  BRANCH(muHiggsLeptDau2dB);
  BRANCHVUINT(muHiggsLeptDau1GlobalMuonBit);
  BRANCHVUINT(muHiggsLeptDau2GlobalMuonBit);
  BRANCHVUINT(muHiggsLeptDau1TrackerMuonBit);
  BRANCHVUINT(muHiggsLeptDau2TrackerMuonBit);
  BRANCH(muHiggsLeptDau1NormChi2);
  BRANCH(muHiggsLeptDau2NormChi2);
  BRANCH(muHiggsLeptDau1NofTrackerHits);
  BRANCH(muHiggsLeptDau2NofTrackerHits);
  BRANCH(muHiggsLeptDau1NofPixelHits);
  BRANCH(muHiggsLeptDau2NofPixelHits);
  BRANCH(muHiggsLeptDau1NofMuonHits);
  BRANCH(muHiggsLeptDau2NofMuonHits);
  BRANCH(muHiggsLeptDau1NofMatches);
  BRANCH(muHiggsLeptDau2NofMatches);
 cout<<"lepton id branhces for 2012 run "<<endl;
  BRANCH(muHiggsMuDau1PFMuonBit);
  BRANCH(muHiggsMuDau2PFMuonBit);
  BRANCH(muHiggsLeptDau1NofMatchedStations);
  BRANCH(muHiggsLeptDau2NofMatchedStations);
  BRANCH(muHiggsMuDau1DzVtx);
  BRANCH(muHiggsMuDau2DzVtx);
  BRANCH(muHiggsLeptDau1NofTrackerLayers);
  BRANCH(muHiggsLeptDau2NofTrackerLayers);
 cout<<"lepton id branches for 2012 run "<<endl;
  BRANCH(muHiggsJetDau1Eta);
  BRANCH(muHiggsJetDau2Eta);
  BRANCH(muHiggsjjdr);
  BRANCH(muHiggslldr);
  BRANCH(muHiggszzdr);
  BRANCH(muHiggszllPt);
  BRANCH(muHiggszllCharge);
  BRANCH(muHiggsMass);
  BRANCH(muHiggsPt);
  BRANCH(muHiggsY);
  BRANCH(muHiggszllMass);
  BRANCH(muHiggszjjMass);
  BRANCH(muHiggsJet1TKHE);
  BRANCH(muHiggsJet2TKHE);
  BRANCH(muHiggsJet1JProb);
  BRANCH(muHiggsJet2JProb);

  BRANCH(muHiggsJet1CSV);
  BRANCH(muHiggsJet2CSV);
  BRANCH(muHiggsJet1CSVMVA);
  BRANCH(muHiggsJet2CSVMVA);
  BRANCH(muHiggsJetDau1PartonFlavour);
  BRANCH(muHiggsJetDau2PartonFlavour);

  BRANCH(muHiggsJetDau1ChHadMult);
  BRANCH(muHiggsJetDau2ChHadMult);
  BRANCH(muHiggsJetDau1NeuHadMult);
  BRANCH(muHiggsJetDau2NeuHadMult);
  BRANCH(muHiggsJetDau1PhotMult);
  BRANCH(muHiggsJetDau2PhotMult);
  BRANCH(muHiggsJetDau1PtDJet);
  BRANCH(muHiggsJetDau2PtDJet);

  BRANCH(muHiggsLeptDau1CombRelIso);
  BRANCH(muHiggsLeptDau2CombRelIso);
  BRANCH(muHiggsLeptDau1ChHadIso);
  BRANCH(muHiggsLeptDau2ChHadIso);
  BRANCH(muHiggsLeptDau1NeuHadIso);
  BRANCH(muHiggsLeptDau2NeuHadIso);
  BRANCH(muHiggsLeptDau1PhotonIso);
  BRANCH(muHiggsLeptDau2PhotonIso);
  BRANCH(muHiggsLeptDau1PUChHadIso);
  BRANCH(muHiggsLeptDau2PUChHadIso);
 cout<<"lepton id branches for 2012 run "<<endl;

BRANCH(muHiggsJetDau1puBeta);
BRANCH(muHiggsJetDau2puBeta);
 cout<<"lepton id branches for 2012 run "<<endl;
  //  BRANCH(muHiggsLeptDau1TrkIso);
  //  BRANCH(muHiggsLeptDau2TrkIso);
  //  BRANCH(muHiggsLeptDau1EcalIso);
  //  BRANCH(muHiggsLeptDau2EcalIso);
  //  BRANCH(muHiggsLeptDau1HcalIso);
  //  BRANCH(muHiggsLeptDau2HcalIso);

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

 cout<<"end muons "<<endl;


  // ELECTRON BRANCHES
 BRANCHBOOL(passSingleElTrig);
  BRANCHBOOL(passDoubleElTrig);
  BRANCH(elHiggsLeptDau1Pt);
  BRANCH(elHiggsLeptDau2Pt);
  BRANCH(elHiggsJetDau1Pt);
  BRANCH(elHiggsJetDau2Pt);
  BRANCH(elHiggsJetDau1RefitPt);
  BRANCH(elHiggsJetDau2RefitPt);
  BRANCH(elHiggsLeptDau1Eta);
  BRANCH(elHiggsLeptDau2Eta);
  BRANCH(elHiggsLeptDau1Phi);
  BRANCH(elHiggsLeptDau2Phi);
  BRANCH(elHiggsLeptDau1EtaSC);
  BRANCH(elHiggsLeptDau2EtaSC);
  BRANCH(elHiggsLeptDau1dB);
  BRANCH(elHiggsLeptDau2dB);
  BRANCH(elHiggsJetDau1Eta);
  BRANCH(elHiggsJetDau2Eta);
  BRANCH(elHiggsjjdr);
  BRANCH(elHiggslldr);
  BRANCH(elHiggszzdr);
  BRANCH(elHiggszllPt);
  BRANCH(elHiggszllCharge);
  BRANCH(elHiggsMass);
  BRANCH(elHiggsPt);
  BRANCH(elHiggsY);
  BRANCH(elHiggszllMass);
  BRANCH(elHiggszjjMass);
  BRANCH(elHiggsJet1TKHE);
  BRANCH(elHiggsJet2TKHE);
  BRANCH(elHiggsJet1CSV);
  BRANCH(elHiggsJet2CSV);
  BRANCH(elHiggsJet1JProb);
  BRANCH(elHiggsJet2JProb);
  BRANCH(elHiggsJet1CSVMVA);
  BRANCH(elHiggsJet2CSVMVA);
  BRANCH(elHiggsJetDau1PartonFlavour);
  BRANCH(elHiggsJetDau2PartonFlavour);

  BRANCH(elHiggsJetDau1ChHadMult);
  BRANCH(elHiggsJetDau2ChHadMult);
  BRANCH(elHiggsJetDau1NeuHadMult);
  BRANCH(elHiggsJetDau2NeuHadMult);
  BRANCH(elHiggsJetDau1PhotMult);
  BRANCH(elHiggsJetDau2PhotMult);
  BRANCH(elHiggsJetDau1PtDJet);
  BRANCH(elHiggsJetDau2PtDJet);
 cout<<"lepton id branches for 2012 run "<<endl;
  BRANCH(elHiggsEleDau1mvaTrigV0); 
  BRANCH(elHiggsEleDau2mvaTrigV0); 
  BRANCH(elHiggsJetDau1puBeta);
  BRANCH(elHiggsJetDau2puBeta);
 cout<<"lepton id branches: mva and beta "<<endl;
  //  BRANCH(elHiggsLeptDau1TrkIso);
  //  BRANCH(elHiggsLeptDau2TrkIso);
  //  BRANCH(elHiggsLeptDau1EcalIso);
  //  BRANCH(elHiggsLeptDau2EcalIso);
  //  BRANCH(elHiggsLeptDau1HcalIso);
  //  BRANCH(elHiggsLeptDau2HcalIso);

  /* Angular variables  */
  BRANCH(elHiggscosthetaNT1);
  BRANCH(elHiggscosthetaNT2);
  BRANCH(elHiggscosthetastarNT);
  BRANCH(elHiggsphiNT);
  BRANCH(elHiggsphiNT1);
  BRANCH(elHiggsHelyLD);


  BRANCH(elHiggscosthetaNT1Refit);
  BRANCH(elHiggscosthetaNT2Refit);
  BRANCH(elHiggscosthetastarNTRefit);
  BRANCH(elHiggsphiNTRefit);
  BRANCH(elHiggsphiNT1Refit);
  BRANCH(elHiggsHelyLDRefit);

  /*                    */
  BRANCH(elHiggsRefitMass);
  BRANCHINT(elHiggsEventNumber);
  BRANCHINT(elHiggsRunNumber);
  BRANCHINT(elHiggsLumiblock);

  /* electronID variables */ 
  BRANCH(elHiggsEleDau1VBTF80CombID);
  BRANCH(elHiggsEleDau2VBTF80CombID);
  BRANCH(elHiggsEleDau1VBTF95CombID);
  BRANCH(elHiggsEleDau2VBTF95CombID);

  BRANCH(elHiggsLeptDau1ChHadIso);
  BRANCH(elHiggsLeptDau2ChHadIso);
  BRANCH(elHiggsLeptDau1NeuHadIso);
  BRANCH(elHiggsLeptDau2NeuHadIso);
  BRANCH(elHiggsLeptDau1PhotonIso);
  BRANCH(elHiggsLeptDau2PhotonIso);
  BRANCH(elHiggsLeptDau1PUChHadIso);
  BRANCH(elHiggsLeptDau2PUChHadIso);

std::cout<<"isolation"<<std::endl;
  BRANCHVUINT(elHiggsLeptDau1isEB);
  BRANCHVUINT(elHiggsLeptDau2isEB);
  BRANCH(elHiggsLeptDau1DeltaPhiAtVtx);
  BRANCH(elHiggsLeptDau2DeltaPhiAtVtx);
  BRANCH(elHiggsLeptDau1HOverE);
  BRANCH(elHiggsLeptDau2HOverE);
std::cout<<"electron ID first part"<<std::endl;

BRANCH(elHiggsEleDau1dz);
BRANCH(elHiggsEleDau2dz);
BRANCH(elHiggsEleDau1dxy);
BRANCH(elHiggsEleDau2dxy);
std::cout<<"electron ID second part"<<std::endl;

BRANCH(elHiggsEleDau1TrkIso03);
BRANCH(elHiggsEleDau2TrkIso03);
BRANCH(elHiggsEleDau1EcalIso03);
BRANCH(elHiggsEleDau2EcalIso03);
BRANCH(elHiggsEleDau1HcalIso03);
BRANCH(elHiggsEleDau2HcalIso03);
std::cout<<"electron ID third part"<<std::endl;

BRANCH(elHiggsLeptDau1DeltaEtaAtVtx);
BRANCH(elHiggsLeptDau2DeltaEtaAtVtx);
BRANCH(elHiggsLeptDau1Sigmaee);
BRANCH(elHiggsLeptDau2Sigmaee);
BRANCH(elHiggsLeptDau1EPin);
BRANCH(elHiggsLeptDau2EPin);
BRANCH(elHiggsLeptDau1Fbrem);
BRANCH(elHiggsLeptDau2Fbrem);
BRANCH(elHiggsLeptDau1EcalEn);
BRANCH(elHiggsLeptDau2EcalEn);
std::cout<<"lepton isolation"<<std::endl;
