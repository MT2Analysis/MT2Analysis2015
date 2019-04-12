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
//double pow_bg(double mt2, double *param);


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
  /* std::cout << "Try with first region: 250MeV < ht < 450MeV, 60MeV < mt2 < 100MeV" << std::endl;
  tree->Draw("deltaPhiMin", "ht<450 && ht>250 && mt2>60 && mt2<100");
  //TBranch* DeltaPhiMin = tree->GetBranch("DeltaPhiMin"); 
  
  std::cout << "r_phi (first region)  = "<< ratio(250.0, 450.0, 60.0, 100.0, 0.3, tree) << std::endl;
  c1->SaveAs("plotFilippo/provaqcd.pdf");
  */

  do_fit(250.0, 450.0, 8, tree);
  
  delete c1;
  data->Close();
  return 0;
}


/*double pow_bg(double mt2, double *param){
    //fitting function for qcd background
    return param[0]*pow(mt2,param[1]);
};*/

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
  
  double mt2_s [n_bins]; 
  
  //trovare nome carino per istogramma
  TH1D* histo = new TH1D("histo", "fit", n_bins, mt2_min, mt2_max); 
   
  //compute r_phi ratio for every bin:
  for (int i=0; i<n_bins; ++i){
      mt2_s[i] = mt2_min + (double)i*interval_mt2;
      double* rphi = new double;
      *rphi = ratio(ht_min, ht_max, mt2_s[i], mt2_s[i] + interval_mt2, delta_Phi_threshold, tree); 
      histo->SetBinContent(i+1, *rphi);
      //std::cout << mt2_s[i] << "GeV < mt2 < " << mt2_s[i] + interval_mt2 << " GeV : r_phi = " << histo->GetBinContent(i) << std::endl;
      std::cout << mt2_s[i] << "GeV < mt2 < " << mt2_s[i] + interval_mt2 << " GeV : r_phi = " << *rphi << std::endl;
      delete rphi;
      }
 
  TF1 *pow_bg = new TF1("pow_bg", "[0]*pow(x,[1])", mt2_min, mt2_max);
  histo->Fit("pow_bg");
 
  TCanvas* cfit = new TCanvas("cfit","Fit with power law"); 
  gPad->SetLogx(); 
  gPad->SetLogy();
  histo->GetXaxis()->SetTitle("MT2");  
  histo->GetYaxis()->SetTitle("r_{#phi}");
  histo->LabelsOption("h","Y"); //IN QUALCHE MODO IL LABEL E' SEMPRE STORTO
  histo->Draw(); //	CAPIRE COME USARE UNO STILE COME NEL FIT DI MASCIOVECCHIO
  cfit->SaveAs("plotFilippo/provafit.pdf");
  
  //CALCOLARE IL RATIO ANCHE PER REGIONI FUORI DA 60-100
  //CAPIRE PERCHE' ERRORE COSI' GRANDE SUI PARAMETRI
  //COPIARE STILE DA CODICE MASCIOVECCHIO
  //PLOTTARE TUTTO, INCLUSA FUNZIONE FITTATA SOPRA, PLOTTATA SU TUTTO IL RANGE DI MT2
  //SCRIVERE NEL PLOT I PARAMETRI OTTENUTI? IN OGNI CASO CAPIRE E SCRIVERE IL CHI2/NDF   
  delete cfit;
  delete histo;
  delete pow_bg;
};


