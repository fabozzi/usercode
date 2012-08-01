{


  // Distribution used for Summer2012 MC.

  Double_t Summer2012[60] = {
    2.344E-05,
    2.344E-05,
    2.344E-05,
    2.344E-05,
    4.687E-04,
    4.687E-04,
    7.032E-04,
    9.414E-04,
    1.234E-03,
    1.603E-03,
    2.464E-03,
    3.250E-03,
    5.021E-03,
    6.644E-03,
    8.502E-03,
    1.121E-02,
    1.518E-02,
    2.033E-02,
    2.608E-02,
    3.171E-02,
    3.667E-02,
    4.060E-02,
    4.338E-02,
    4.520E-02,
    4.641E-02,
    4.735E-02,
    4.816E-02,
    4.881E-02,
    4.917E-02,
    4.909E-02,
    4.842E-02,
    4.707E-02,
    4.501E-02,
    4.228E-02,
    3.896E-02,
    3.521E-02,
    3.118E-02,
    2.702E-02,
    2.287E-02,
    1.885E-02,
    1.508E-02,
    1.166E-02,
    8.673E-03,
    6.190E-03,
    4.222E-03,
    2.746E-03,
    1.698E-03,
    9.971E-04,
    5.549E-04,
    2.924E-04,
    1.457E-04,
    6.864E-05,
    3.054E-05,
    1.282E-05,
    5.081E-06,
    1.898E-06,
    6.688E-07,
    2.221E-07,
    6.947E-08,
    2.047E-08
  };  


// Summer11 PU_S4, distribution obtained by averaging the number of 
// interactions in each beam crossing to estimate the true mean.  
// THIS IS THE RECOMMENDED ONE for reweighting.

/*
  Double_t PoissonIntDist[36] = {
    0.104109,
    0.0703573,
    0.0698445,
    0.0698254,
    0.0697054,
    0.0697907,
    0.0696751,
    0.0694486,
    0.0680332,
    0.0651044,
    0.0598036,
    0.0527395,
    0.0439513,
    0.0352202,
    0.0266714,
    0.019411,
    0.0133974,
    0.00898536,
    0.0057516,
    0.00351493,
    0.00212087,
    0.00122891,
    0.00070592,
    0.000384744,
    0.000219377,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.
  };
*/

// Flat10+Tail distribution taken directly from MixingModule input:  
// (Can be used for Spring11 and Summer11 if you don't worry about 
// small shifts in the mean) SHOULD be used for 3-D Reweighting, 
// as this is the "true" input for all Summer11 samples.

/*
  Double_t probdistFlat10[36] = {
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0630151648,
    0.0526654164,
    0.0402754482,
    0.0292988928,
    0.0194384503,
    0.0122016783,
    0.007207042,
    0.004003637,
    0.0020278322,
    0.0010739954,
    0.0004595759,
    0.0002229748,
    0.0001028162,
    4.58337152809607E-05,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.
  };

*/

//  TH1F PUS4_Distr("PUS4_Distr","PUS4_Distr",36,-0.5,35.5);
  TH1F PUS7_Distr("PUS7_Distr","PUS7_Distr",60,0.,60.);
  
  for(int i = 1; i<61; ++i){
    PUS7_Distr.SetBinContent(i,Summer2012[i-1]);
    //    PUS4_Distr.SetBinContent(i,PoissonIntDist[i-1]);
    //    PUS4_Distr.SetBinContent(i,probdistFlat10[i-1]);
  }

  TFile pu_mc("PUDist_Summer12MC_Extended.root", "RECREATE");
  PUS7_Distr.Write();
  pu_mc.Close();
}
