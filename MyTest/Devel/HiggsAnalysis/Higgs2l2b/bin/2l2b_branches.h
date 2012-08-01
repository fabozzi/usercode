
  // EVENT BRANCHES
  BRANCHFLOAT(met);
  BRANCHFLOAT(t1corrMet);
  BRANCHFLOAT(metPhi);
  BRANCHFLOAT(metSignificance);
  BRANCHINT(numPV);

  std::cout<<"check"<<std::endl;
// -- IN FUTURE NTUPLES --
//  BRANCHFLOAT(rhoRestrictedEta);
// -----------------------

  // MUON BRANCHES
  BRANCHBOOL(passSingleMuTrig);
  BRANCHBOOL(passDoubleMuTrig);
  std::cout<<"check2"<<std::endl;
  BRANCH(CleanJetPt);
  BRANCH(muHiggsLeptDau1Pt);
  BRANCH(muHiggsLeptDau2Pt);
  BRANCH(muHiggsJetDau1Pt);
  BRANCH(muHiggsJetDau1RefitPt);
  BRANCH(muHiggsJetDau2Pt);
  BRANCH(muHiggsJetDau2RefitPt);
  BRANCH(muHiggsLeptDau1Eta);
  BRANCH(muHiggsLeptDau2Eta);
  BRANCH(muHiggsLeptDau1dB);
  BRANCH(muHiggsLeptDau2dB);
  std::cout<<"check2"<<std::endl;
  BRANCHVUINT(muHiggsLeptDau1GlobalMuonBit);
  BRANCHVUINT(muHiggsLeptDau2GlobalMuonBit);
  BRANCHVUINT(muHiggsLeptDau1TrackerMuonBit);
  BRANCHVUINT(muHiggsLeptDau2TrackerMuonBit);
  BRANCH(muHiggsLeptDau1NormChi2);
  BRANCH(muHiggsLeptDau2NormChi2);
  std::cout<<"quality"<<std::endl;


  BRANCH(muHiggsLeptDau1NofPixelHits);
  BRANCH(muHiggsLeptDau2NofPixelHits);
  BRANCH(muHiggsLeptDau1NofMuonHits);
  BRANCH(muHiggsLeptDau2NofMuonHits);
  BRANCH(muHiggsLeptDau1NofMatches);
  BRANCH(muHiggsLeptDau2NofMatches);
  std::cout<<"jet"<<std::endl;
  BRANCH(muHiggsLeptDau1NofTrackerHits);
  BRANCH(muHiggsLeptDau2NofTrackerHits);


  BRANCH(muHiggsJetDau1Eta);
  BRANCH(muHiggsJetDau2Eta);
  BRANCH(muHiggsJetDau1Phi);
  BRANCH(muHiggsJetDau2Phi);
  BRANCH(muHiggsLeptDau1Phi);
  BRANCH(muHiggsLeptDau2Phi);
  BRANCH(muHiggsjjdr);
  BRANCH(muHiggslldr);
  BRANCH(muHiggszzdr);
  BRANCH(muHiggszzdPhi);
  BRANCH(muHiggszzdEta);
  BRANCH(muHiggsjjdPhi);
  BRANCH(muHiggsjjdEta);
  BRANCH(muHiggszllPt);
  BRANCH(muHiggszllPhi);
  BRANCH(muHiggszllCharge);
  BRANCH(muHiggszjjPt);
  BRANCH(muHiggsMass);
  BRANCH(muHiggsPt);
  BRANCH(muHiggsY);
  BRANCH(muHiggstrkMet);
  BRANCH(muHiggszllMass);
  BRANCH(muHiggszjjMass);
  BRANCH(muHiggsJet1TKHE);
  BRANCH(muHiggsJet2TKHE);
  BRANCH(muHiggsJet1CSV);
  BRANCH(muHiggsJet2CSV);
  BRANCH(muHiggsJetDau1PartonFlavour);
  BRANCH(muHiggsJetDau2PartonFlavour);
// FOR FUTURE NTUPLES
//  BRANCH(muHiggsJetDau1ChHadMult);
//  BRANCH(muHiggsJetDau2ChHadMult);
//  BRANCH(muHiggsJetDau1NeuHadMult);
//  BRANCH(muHiggsJetDau2NeuHadMult);
//  BRANCH(muHiggsJetDau1PhotMult);
//  BRANCH(muHiggsJetDau2PhotMult);
//  BRANCH(muHiggsJetDau1PtDJet);
//  BRANCH(muHiggsJetDau2PtDJet);
//-------------------

  BRANCH(muHiggsLeptDau1CombRelIso);
  BRANCH(muHiggsLeptDau2CombRelIso);
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
// ----- FUTURE NTUPLES -------
//  BRANCH(elHiggsLeptDau1EtaSC);
//  BRANCH(elHiggsLeptDau2EtaSC);
// ----------------------------
  BRANCH(elHiggsLeptDau1dB);
  BRANCH(elHiggsLeptDau2dB);
  BRANCH(elHiggsJetDau1Eta);
  BRANCH(elHiggsJetDau2Eta);
  BRANCH(elHiggsJetDau1Phi);
  BRANCH(elHiggsJetDau2Phi);
  BRANCH(elHiggsLeptDau1Phi);
  BRANCH(elHiggsLeptDau2Phi);
  BRANCH(elHiggsjjdr);
  BRANCH(elHiggslldr);
  BRANCH(elHiggszzdr);
  BRANCH(elHiggszzdPhi);
  BRANCH(elHiggszzdEta);
  BRANCH(elHiggsjjdPhi);
  BRANCH(elHiggsjjdEta);
  BRANCH(elHiggszllPt);
  BRANCH(elHiggszllCharge);
  BRANCH(elHiggszllPhi);
  BRANCH(elHiggszjjPt);
  BRANCH(elHiggsMass);
  BRANCH(elHiggsPt);
  BRANCH(elHiggsY);
  BRANCH(elHiggstrkMet);
  BRANCH(elHiggszllMass);
  BRANCH(elHiggszjjMass);
  BRANCH(elHiggsJet1TKHE);
  BRANCH(elHiggsJet2TKHE);
  BRANCH(elHiggsJet1CSV);
  BRANCH(elHiggsJet2CSV);
  BRANCH(elHiggsJetDau1PartonFlavour);
  BRANCH(elHiggsJetDau2PartonFlavour);

// FOR FUTURE NTUPLES
//  BRANCH(elHiggsJetDau1ChHadMult);
//  BRANCH(elHiggsJetDau2ChHadMult);
//  BRANCH(elHiggsJetDau1NeuHadMult);
//  BRANCH(elHiggsJetDau2NeuHadMult);
//  BRANCH(elHiggsJetDau1PhotMult);
//  BRANCH(elHiggsJetDau2PhotMult);
//  BRANCH(elHiggsJetDau1PtDJet);
//  BRANCH(elHiggsJetDau2PtDJet);
// -------------------

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
  BRANCHVUINT(elHiggsLeptDau1isEB);
  BRANCHVUINT(elHiggsLeptDau2isEB);
  BRANCH(elHiggsLeptDau1DeltaPhiAtVtx);
  BRANCH(elHiggsLeptDau2DeltaPhiAtVtx);
  BRANCH(elHiggsLeptDau1HOverE);
  BRANCH(elHiggsLeptDau2HOverE);

