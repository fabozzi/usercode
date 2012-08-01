#include <iostream>
#include <fstream>
#include "HiggsAnalysis/Higgs2l2b/interface/BTagSFUtil.h"
#include "TMath.h"
//#include "HiggsAnalysis/Higgs2l2b/src/SFlightFuncs.C"


/***********************************/
/*           NEW FUNCTION          */
/***********************************/

void BTagSFUtil::SF(const std::string& btagAlgo, const std::string& wp, float pt, float eta){
  int NL = 4; 
  int NM = 3;

  //binning in eta is the same for SFl and mistag, but it changes with the wp

  float etamin_L[4] = {0.0, 0.5, 1.0, 1.5} ;
  float etamax_L[4] = {0.5, 1.0, 1.5, 2.4} ;
  float etamin_M[3] = {0.0, 0.8, 1.6} ;
  float etamax_M[3] = {0.8, 1.6, 2.4} ;
  float etamin(0), etamax(0);
  // SF for b's
  //  float b_SF_Medium = h2_Medium_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
  //  float b_SF_Loose = h2_Loose_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
  //  float SFb_;
  //  float SFlight_, Mistag_;
  SFlightFuncs sfl_func;
  MistagFuncs mt_func;   
  float pt_sfb;
  pt_sfb = pt;

  if(pt>670){ 
    std::cout<<"WARNING: Mistagging rate for pt greater than 670 are not defined, we are going to treat this case as a pt=670 case"<<std::endl;
    pt = 670;}
  if(pt<30) {
    pt_sfb = 30;
  }

  if(wp == "L"){ 
    if(btagAlgo == "CSV")   SFb_ = 1.02658*((1.+(0.0195388*pt_sfb))/(1.+(0.0209145*pt_sfb)));
    if(btagAlgo == "TCHE") SFb_ = 0.603913*((1.+(0.286361*pt_sfb))/(1.+(0.170474*pt_sfb)));
    if(btagAlgo == "JP") SFb_ = 0.969851*((1.+(-6.06362e-05*pt_sfb))/(1.+(-0.000156638*pt_sfb)));
    for(int i=0;i<NL;++i){
      if ((TMath::Abs(eta) >etamin_L[i] || TMath::Abs(eta) == etamin_L[i]) && TMath::Abs(eta) <etamax_L[i]){
	etamin = etamin_L[i];
	etamax = etamax_L[i];
      }
    }
  }// end L 

 if(wp == "M"){ 
   if(btagAlgo == "TCHE") SFb_ = 0.932251*((1.+(0.00335634*pt_sfb))/(1.+(0.00305994*pt_sfb)));
   if(btagAlgo == "CSV")  SFb_ = 0.6981*((1.+(0.414063*pt_sfb))/(1.+(0.300155*pt_sfb)));
   if(btagAlgo == "JP") SFb_ = 0.90806*((1.+(0.000236997*pt_sfb))/(1.+(5.49455e-05*pt_sfb)));
    for(int i=0;i<NM;++i){
      if (( TMath::Abs(eta) >etamin_M[i] || TMath::Abs(eta) == etamin_L[i]) && TMath::Abs(eta) <etamax_M[i]){
	etamin = etamin_M[i];
	etamax = etamax_M[i];
      }
    }
  }// end M 

 // uncomment for debugging
 //  std::cout<<"ETA: "<<eta<<std::endl;
 //  std::cout<<"ETAMIN AND ETAMAX ARE: "<<etamin<<" ,"<<etamax<<std::endl;
  TF1* SFlight_func = sfl_func.GetSFlmean(TString(btagAlgo),TString(wp),etamin, etamax);
  SFlight_ = SFlight_func->Eval(pt);
  TF1* Mistag_func = mt_func.GetMistagmean(TString(btagAlgo),TString(wp),etamin, etamax);
  Mistag_ =  Mistag_func->Eval(pt);
}


BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom3(seed);
  setSFFileName("");

}

/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

void BTagSFUtil::set_fileMedium( TFile* file ) {

  fileMedium_ = file;

  if( fileMedium_==0 ) {
    std::cout << "WARNING!!! File '" << file->GetName() << "' does not exist!! Will not set histograms!!" << std::endl;
    return;
  }

  h2_Medium_BTAGBEFFCORR_ = (TH2D*)fileMedium_->Get("BTAGBEFFCORR");
  h2_Medium_BTAGLEFFCORR_ = (TH2D*)fileMedium_->Get("BTAGLEFFCORR");
  h2_Medium_BTAGLEFF_     = (TH2D*)fileMedium_->Get("BTAGLEFF");

  if( h2_Medium_BTAGBEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGBEFFCORR' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Medium_BTAGLEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFFCORR' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Medium_BTAGLEFF_     == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFF' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }

}


/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

void BTagSFUtil::set_fileLoose( TFile* file ) {

  fileLoose_ = file;

  if( fileLoose_==0 ) {
    std::cout << "WARNING!!! File '" << file->GetName() << "' does not exist!! Will not set histograms!!" << std::endl;
    return;
  }

  h2_Loose_BTAGBEFFCORR_ = (TH2D*)fileLoose_->Get("BTAGBEFFCORR");
  h2_Loose_BTAGLEFFCORR_ = (TH2D*)fileLoose_->Get("BTAGLEFFCORR");
  h2_Loose_BTAGLEFF_     = (TH2D*)fileLoose_->Get("BTAGLEFF");

  if( h2_Loose_BTAGBEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGBEFFCORR' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Loose_BTAGLEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFFCORR' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Loose_BTAGLEFF_     == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFF' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }

}

/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, int pdgIdPart, float Btageff_SF_l, float Btageff_SF_m, float Btagmistag_SF_l,float Btagmistag_SF_m, float Btagmistag_eff_l, float Btagmistag_eff_m) {
  

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > Btageff_SF_m ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > Btageff_SF_l ) {isBTagged_loose=false; }//
    }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    float  Btagmistag_SF = 1.0;
    float  Btagmistag_eff = 1.0;

    if( isBTagged_loose ) {
      Btagmistag_SF = Btagmistag_SF_m;
      Btagmistag_eff = Btagmistag_eff_m;
    } else  {
      Btagmistag_SF = Btagmistag_SF_l;
      Btagmistag_eff = Btagmistag_eff_l;
    }

    float mistagPercent = ( Btagmistag_SF*Btagmistag_eff - Btagmistag_eff ) / ( 1. - Btagmistag_eff );

    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:

    if( coin < mistagPercent ) {
      if( !isBTagged_loose ) {isBTagged_loose = true;}
      else if( !isBTagged_medium ) {isBTagged_medium = true; }

    }
    
  } //if light quark

} //modifyBTagsWithSF






void BTagSFUtil::modifyBTagsWithSF( const std::string& btagAlgo, bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, float sysSF) {
  


  SF(btagAlgo, "M", jetpt, jeteta);
  float  b_SF_Medium = SFb_;
  float light_SF_Medium = SFlight_;
  float light_eff_Medium = Mistag_;
  SF(btagAlgo, "L", jetpt, jeteta);
  float  b_SF_Loose = SFb_;
  float light_SF_Loose = SFlight_;
  float light_eff_Loose = Mistag_;
    
  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    //int iBin_pt = (jetpt<240.) ? h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(jetpt) : h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(239.);
    //int iBin_eta = h2_Medium_BTAGBEFFCORR_->GetYaxis()->FindBin(jeteta);

    // SF for b's

    /* old recipe */
    // float b_SF_Medium = h2_Medium_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
    // float b_SF_Loose = h2_Loose_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );



    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > b_SF_Medium ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > b_SF_Loose ) isBTagged_loose=false; //
    }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    //int iBin_pt = (jetpt<500.) ? h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(jetpt) : h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(499.);
    //int iBin_eta = h2_Medium_BTAGBEFFCORR_->GetYaxis()->FindBin(jeteta);

    // SF for light quarks
    //    float light_SF_Medium = h2_Medium_BTAGLEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
    //    float light_SF_Loose = h2_Loose_BTAGLEFFCORR_->GetBinContent( iBin_pt, iBin_eta );



    //float light_eff_Medium = h2_Medium_BTAGLEFF_->GetBinContent( iBin_pt, iBin_eta );
    //float light_eff_Loose = h2_Loose_BTAGLEFF_->GetBinContent( iBin_pt, iBin_eta );


    float mistagPercent_Loose = ( light_SF_Loose*light_eff_Loose - light_eff_Loose ) / ( 1. - light_eff_Loose );
    float mistagPercent_Medium = ( light_SF_Medium*light_eff_Medium - light_eff_Medium ) / ( 1. - light_eff_Medium );


    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:
    if( !isBTagged_loose ) {
      if( coin < mistagPercent_Loose ) isBTagged_loose = true;
    } else if( !isBTagged_medium ) {
      if( coin < mistagPercent_Medium ) isBTagged_medium = true; 
    }

    
  } //if light quark

} //modifyBTagsWithSF



/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

//This method of getSF uses the functional form of the mistag SF based on 2011 data. We do not need an eta information becasuse the dependence is flat in eta. 
BTagScaleFactor BTagSFUtil::getSF( const std::string& type, float jetpt ) {

  BTagScaleFactor btsf;

  float SFlight = -1;

  if ( type == "medium"){
  
    SFlight = 1.23344 - 0.00033012*jetpt - 0.000000227661*jetpt*jetpt + 0.00000000154075*jetpt*jetpt*jetpt;
    btsf.SF = SFlight;

  }
  else{
    
    SFlight = 1.14062 - 0.000456396*jetpt + 0.0000000838135*jetpt*jetpt;
    btsf.SF = SFlight;

  }

  btsf.SF_err = 0.;
  btsf.eff = 1.;
  btsf.eff_err = 0.;

  return btsf;


}


/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

BTagScaleFactor BTagSFUtil::getSF( const std::string& fileName, float jetpt, float jeteta ) {

   BTagScaleFactor btsf;

   bool foundSF = false;

   ifstream ifs(fileName.c_str());

   if(!ifs.good()){
     std::cout<<"WARNING ! Didn't find file "<<fileName.c_str()<<std::endl;
   }

   while( ifs.good() && !foundSF ) {

     float etaMin, etaMax, ptMin, ptMax, eff, eff_err, SF, SF_err;
     ifs >> etaMin >> etaMax >> ptMin >> ptMax >> eff >> eff_err >> SF >> SF_err;

     if( fabs(jeteta)>=etaMin && fabs(jeteta)<etaMax && jetpt>=ptMin && (jetpt<ptMax||ptMax==999.) ) {
       btsf.SF = SF;
       btsf.SF_err = SF_err;
       btsf.eff = eff;
       btsf.eff_err = eff_err;
       foundSF = true;
     }

   } // while

   ifs.close();

   if( !foundSF ) {
     std::cout << "WARNING! Didn't find SF in file '" << fileName << "'. Setting it to 1." << std::endl;
     btsf.SF = 1.;
     btsf.SF_err = 0.;
     btsf.eff = 1.;
     btsf.eff_err = 0.;
   }

   return btsf;

}
