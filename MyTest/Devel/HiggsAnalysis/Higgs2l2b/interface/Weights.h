#include <string>

class Weights{

 public:
  Weights(std::string filename, std::string name){
    TFile file(filename.c_str());
    wHisto = (TH2F*)file.Get(name.c_str());
  };
  float getEff(float eta, float pt){
    
    int i = wHisto->GetXaxis()->FindBin(eta);
    int j = wHisto->GetYaxis()->FindBin(pt);
    if(i == 0) i = 1;
    if(j == wHisto->GetYaxis()->GetNbins() + 1) j =  wHisto->GetYaxis()->GetNbins();
    float eff = wHisto->GetBinContent(i,j) ;
    return eff;
  }
 
 private:
  TH2F* wHisto;

};
