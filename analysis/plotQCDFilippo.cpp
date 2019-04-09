#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <sstream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TLatex.h"


#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"


void do_fit(float ht_min, float ht_max, int n_bins, TTree* tree);

double ratio(float ht_min, float ht_max, float mt2_low, float mt2_high, float threshold, TTree* tree);

int main( int argc, char* argv[] ) {
    
  if( argc<1 ) {
    std::cout << "USAGE: ./plotQCD [dataFileName] " << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
    
  std::string dataFileName(argv[1]);
  std::cout<< "dataFileName "<< dataFileName << std::endl;

    
  //load data and extract tree
  std::string dataPath("EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/");
  dataPath = dataPath + dataFileName;
  std::cout<< "dataPath "<< dataPath << std::endl;
  std::cout << "Getting data... "; 
  TFile* data = new TFile(dataPath.c_str(), "READ");
  TCanvas* c1 = new TCanvas("c1","First graph");
  TTree* tree = (TTree*)data->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
  std::cout << "Done." << std::endl;


  //now first experiment with first ht region and 60<mt2<100
  std::cout << "Try with first region: 250MeV < ht < 450MeV, 60MeV < mt2 < 100MeV" << std::endl;
  tree->Draw("deltaPhiMin", "ht<450 && ht>250 && mt2>60 && mt2<100");
  //TBranch* DeltaPhiMin = tree->GetBranch("DeltaPhiMin"); 
  
  std::cout << "r_phi (first region)  = "<< ratio(250.0, 450.0, 60.0, 100.0, 0.3, tree) << std::endl;
  c1->SaveAs("plotFilippo/provaqcd.pdf");
 
  do_fit(250.0, 450.0, 8, tree);
  
  delete c1;
  data->Close();
  return 0;
}



double ratio(float ht_min, float ht_max, float mt2_low, float mt2_high, float threshold, TTree* tree){
  //Define conditions for events in given phasespace above/below threshold
  //I need strings which depend on parameters, so I use ostringstream
  std::ostringstream oss_big;
  oss_big << "ht<" << ht_max << " && ht>" << ht_min << " && mt2>" << mt2_low << " && mt2<" << mt2_high << " && deltaPhiMin>" << threshold;
  std::string cond_big = oss_big.str();

  std::ostringstream oss_small;
  oss_small << "ht<" << ht_max << " && ht>" << ht_min << " && mt2>" << mt2_low << " && mt2<" << mt2_high << " && deltaPhiMin<" << threshold;
  std::string cond_small = oss_small.str();
  
  //Find number of entries above and below threshold in the given phasespace region
  Long64_t num_bigphi = tree->GetEntries(cond_big.c_str());
  Long64_t num_smallphi = tree->GetEntries(cond_small.c_str());
  double rat = (double)num_bigphi/(double)num_smallphi;
  
  //Raise warning for low statistics
  if(num_bigphi < 100 || num_smallphi < 100){
    std::ostringstream reg;
    reg << "ht<" << ht_max << " && ht>" << ht_min << " && mt2>" << mt2_low << " && mt2<" << mt2_high;
    std::string region = reg.str();
    std::cout << "WARNING: low statistics in region " << region << std::endl;
   }
  return rat;
}

void do_fit(float ht_min, float ht_max, int n_bins, TTree* tree){
  //ht_min and ht_max are the boundaries of the ht interval to be considered
  //n_bins is the number of bins in mt2
  //tree is the tree containing the data
  double delta_Phi_threshold = 0.3;
  double mt2_min = 60.0; double mt2_max = 100.0;
  double interval_mt2 = (mt2_max - mt2_min)/n_bins;
  
  double mt2_s [n_bins]; double rphi_s [n_bins];
  
  //compute r_phi ratio for every bin:
  for (int i=0; i<n_bins; ++i){
      mt2_s[i] = mt2_min + (double)i*interval_mt2;
      rphi_s[i] = ratio(ht_min, ht_max, mt2_s[i], mt2_s[i] + interval_mt2, delta_Phi_threshold, tree); 
      std::cout << mt2_s[i] << "GeV < mt2 < " << mt2_s[i] + interval_mt2 << " GeV : r_phi = " << rphi_s[i] << std::endl;
      }
};


