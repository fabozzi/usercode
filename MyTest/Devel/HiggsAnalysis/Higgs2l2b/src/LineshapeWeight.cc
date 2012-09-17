#include "../interface/LineshapeWeight.h"
#include <iostream>
#include <fstream>
#include "TMath.h"

using namespace std;

LineshapeWeight::LineshapeWeight(){ }

LineshapeWeight::LineshapeWeight( string inputFile ) :
  weightFileName_(inputFile)
{
  double bincenter, initial, pow, powp, powm, out, outp, outm;

  std::ifstream ifs(weightFileName_.c_str());

  while( ifs.good() ) {
    // bincenter, original, shape rew, +unc, -unc, shape+interf rew, +unc, -unc
    ifs >> bincenter >> initial >> pow >> powp >>  powm >> out>> outp >> outm;
    
    bincenters_.push_back(bincenter);
    if(initial > 0){
      weight_.push_back(   TMath::Max(0.,out/initial) );
      weightP_.push_back(  TMath::Max(0.,outp/initial) );
      weightM_.push_back(  TMath::Max(0.,outm/initial) );
    }else{//weights are not defined if initial distribution is 0 => set weight to 0
      weight_.push_back( 0. );
      weightP_.push_back( 0. );
      weightM_.push_back( 0. );
    }
    
  }
  
}

LineshapeWeight::~LineshapeWeight(){}

void LineshapeWeight::getWeight(float m, double & weight, 
				double & weightp, double & weightm) 
{
  weight = 1.;
  weightp = 1.;
  weightm = 1.;
  
  if( m < bincenters_.front() || m >  bincenters_.back() ){ 
    // set weights to 0 if out of range
    weight = 0.;
    weightp = 0.;
    weightm = 0.;
    return;
  }

  std::vector<double>::iterator low;
  low=lower_bound( bincenters_.begin(), bincenters_.end(),m ); 
  int lowindex=(low -  bincenters_.begin());
  if(m == *low ){//exact match
    weight = weight_[lowindex];
    weightp = weightP_[lowindex];
    weightm = weightM_[lowindex];
  } else {
    //linear interpolation
    lowindex--; // lower_bound finds the first element not smaller than X
    weight = weight_[lowindex] +( m - bincenters_[lowindex] )*(weight_[lowindex+1]-weight_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]) ;
    weightp = weightP_[lowindex] +( m - bincenters_[lowindex] )*(weightP_[lowindex+1]-weightP_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex])  ;
    weightm = weightM_[lowindex] +( m - bincenters_[lowindex] )*(weightM_[lowindex+1]-weightM_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex])  ;
  }
  
  cout << m << " " << bincenters_[lowindex] << "  " << bincenters_[lowindex+1] << " " <<  weight << endl;
  
  return;

}


