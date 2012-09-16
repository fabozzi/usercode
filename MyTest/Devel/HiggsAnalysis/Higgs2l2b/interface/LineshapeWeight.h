#ifndef LineshapeWeight_h
#define LineshapeWeight_h

#include <string>
#include <vector>
#include <cmath>

class LineshapeWeight{

 public:

  LineshapeWeight( std::string inputFile );
  ~LineshapeWeight();

  void getWeight(float hmass, double & weight, 
		   double & weightp, double & weightm);

 private: 

  std::string weightFileName_;

  std::vector<double> bincenters_;
  std::vector<double> weight_;
  std::vector<double> weightP_;
  std::vector<double> weightM_;
  
};
#endif
