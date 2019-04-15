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
#include "TLine.h"

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
void set_ratio(float ht_min, float ht_max, TH1D* histogram, float threshold, TTree* tree);

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


  do_fit(250.0, 450.0, 8, tree);
  
  delete c1;
  data->Close();
  return 0;
}


double ratio(float ht_min, float ht_max, float mt2_low, float mt2_high, float threshold, TTree* tree){
  //Define conditions for events in given phasespace above/below threshold
  //I need strings which depend on parameters, so I use ostringstream
  std::ostringstream oss_big;
  oss_big << "150< id && id<160 &&" << "ht<" << ht_max << " && ht>" << ht_min << " && mt2>" << mt2_low << " && mt2<" << mt2_high << " && deltaPhiMin>" << threshold;
  std::string cond_big = oss_big.str();

  std::ostringstream oss_small;
  oss_small << "150< id && id<160 && "<< "ht<" << ht_max << " && ht>" << ht_min << " && mt2>" << mt2_low << " && mt2<" << mt2_high << " && deltaPhiMin<" << threshold;
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

void set_ratio(float ht_min, float ht_max, TH1D* histogram, float threshold, TTree* tree){

  TAxis* Axis = histogram->GetXaxis();
  
  int num_bins = histogram->GetNbinsX(); 
  double rphi;
  double low_lim; double high_lim;

  for (int i = 1; i<=num_bins; i++){
     low_lim = Axis->GetBinLowEdge(i);
     high_lim = Axis->GetBinUpEdge(i);
     rphi = ratio(ht_min, ht_max, low_lim, high_lim, threshold, tree);
     histogram->SetBinContent(i, rphi);
     std::cout << low_lim << "GeV < mt2 < " << high_lim << " GeV : r_phi = " << rphi << std::endl;
     }
  std::cout<<"Ratio was set";
  
  //FOR SOME REASONS MEMEORY LEAK AT END OF FOR LOOP WITH DELETE AXIS
  //delete Axis;

}

void do_fit(float ht_min, float ht_max, int n_bins, TTree* tree){
  //ht_min and ht_max are the boundaries of the ht interval to be considered
  //n_bins is the number of bins in mt2
  //tree is the tree containing the data
  double delta_Phi_threshold = 0.3;
  //limits of fitting region:
  double mt2_min = 60.0; double mt2_max = 100.0;
  
  //max and min mt2 for whole plot (not only region for fitting)
  double mt2_min_global = 50.0; double mt2_max_global = 450.0;
 

  //trovare nome carino per istogramma
 // TH1D* histo = new TH1D("histo", "fit", n_bins, mt2_min, mt2_max); 
  float bins[] = {50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,120.0,200.0,300.0,450.0};
  int  binnum = sizeof(bins)/sizeof(float) - 1;
  TH1D* histo = new TH1D("histo", "CMS simulation, #sqrt{s} = 13 TeV", binnum,  bins);
  //compute r_phi ratio for every bin:
 
  set_ratio(ht_min, ht_max, histo, delta_Phi_threshold, tree);


  TF1 *pow_bg = new TF1("pow_bg", "[0]*pow(x,[1])", mt2_min, mt2_max);
  histo->Fit("pow_bg", "R 0"); //R: fit in the range of the function. don't plot now the fitted function
  //plot it later on whole rangem also beyone 100GeV
 
  TF1 *fitResult = histo->GetFunction("pow_bg");

  double param_a = fitResult->GetParameter(0);
  double param_b = fitResult->GetParameter(1);
  
  //now define fitted function over whole range (i.e. not only fitting region)
  TF1 *fitted_bg = new TF1("pow_bg", "[0]*pow(x,[1])", mt2_min_global, mt2_max_global);
  fitted_bg->SetParameter(0, param_a);
  fitted_bg->SetParameter(1, param_b);

  TCanvas* cfit = new TCanvas("cfit","Fit with power law"); 
  gPad->SetLogx(); 
  gPad->SetLogy();
  histo->GetXaxis()->SetTitle("M_{T2} [GeV]");  
  histo->GetYaxis()->SetTitle("r_{#phi}");
  histo->LabelsOption("h","Y"); //IN QUALCHE MODO IL LABEL E' SEMPRE STORTO
  //histo->SetMarkerStyle(kFullSquare); per qualche motivo non funziona, ottengo una linea continua
  histo->SetStats(0); //don't show information box on graph
  histo->SetMarkerStyle(24);
  histo->SetLineColor(kBlack);
  histo->Draw("P");
  fitted_bg->Draw("SAME");  
  
  //add text with descrption of considered region
  std::ostringstream descr;
  descr << ht_min << " GeV < H_{T} < " << ht_max << " GeV";
  std::string description = descr.str(); 
  TPaveText *pt = new TPaveText(0.4,0.8,0.65,0.85, "NDC");
  pt->AddText(description.c_str());
  pt->Draw("SAME");
  
  //draw vertical lines at 60 and 100 GeV Mt2
  //double x_min = gPad->GetUxmin(); double x_max = gPad->GetUxmax();
  //double y_min = gPad->GetUymin(); double y_max = gPad->GetUymax();
  //double NDC_60 = x_min + (60.0- x_min)/(x_max - x_min);
  //double NDC_100 = x_min + (100.0- x_min)/(x_max - x_min);
  float y_min = 0.1; float y_max = 6.0; //NOT NICE: FIX TO FIND GENERAL WAY TO FIND YMAX ALSO IN LOGPLOTI
  TLine *line1 = new TLine(mt2_min,y_min, mt2_min, y_max);
  TLine *line2 = new TLine(mt2_max,y_min, mt2_max, y_max);
  line1->SetLineStyle(2); line2->SetLineStyle(2);
  line1->Draw("SAME");
  line2->Draw("SAME");

  /*TPaveText *chi2;
  chi2 = new TPaveText(0.22,0.23, 0.5,0.18,"brNDC");
  chi2->SetFillColor(10);   chi2->SetBorderSize(0);
  chi2->AddText(Form("#chi^{2}/ndf = %.1f/%d: %.1f%%", thisFitQCD->GetChisquare(), thisFitQCD->GetNDF(), thisFitQCD->GetProb()*100));
  chi2->Draw("SAME");*/ //SISTEMARE E AGGIUNGERE CHISQUARE

  cfit->SaveAs("plotFilippo/provafit.pdf");
  //problema con getUymax in logascale, da' solo 1 e vegono linee non fino in cima 
  //INSERIRE BARRA ORIZZONTALE PER LARGHEZZA BIN OLTRE AL PUNTO CENTRALE
  //CAPIRE PERCHE' ERRORE COSI' GRANDE SUI PARAMETRI
  //COPIARE STILE DA CODICE MASCIOVECCHIO
  //PLOTTARE ANCHE FIT CON RANGE DELLE INCERTEZZE
  //SCRIVERE NEL PLOT I PARAMETRI OTTENUTI? IN OGNI CASO CAPIRE E SCRIVERE IL CHI2/NDF   
  delete cfit;
  delete histo;
  delete pow_bg;
};


