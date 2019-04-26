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
#include "TMath.h"
#include "TLeaf.h"
#include "TCut.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"


void do_fit(float ht_min, float ht_max, std::string  cond_all, std::string cond_qcd, TTree* tree, bool useMC);
double ratio(float ht_min, float ht_max, float mt2_low, float mt2_high, std::string  condition, float threshold, TTree* tree, bool useMC);
TH1D* set_ratio(float ht_min, float ht_max, std::string  condition, float threshold, TTree* tree);
//double SumWeight(TTree* tree, std::string condition);

int main( int argc, char* argv[] ) {
    
  if( argc<1 ) {
    std::cout << "USAGE: ./plotQCDFilippo [dataFileName] [MC/data] [closureTest=true/false]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

    bool useMC = true;
    bool closureTest = true;  //true or false?? if false metcut. CAPIRE PERCHE'
    std::string dataFileName;
    
  if(argc>1) {
   dataFileName = argv[1];
   std::cout<< "dataFileName "<< dataFileName << std::endl;
   if( argc>2 ) {
       std::string mc_or_data = std::string(argv[2]);
        if( mc_or_data=="data" )  useMC = false;
        }
    
    if( argc>3 ) {
        std::string closurestr(argv[3]);
        if( closurestr=="closureTest" || closurestr=="true" ) {
            closureTest=true;
            std::cout << "-> Running closure test in validation region" << std::endl;
         }
     }
    }
    
  //load data and extract tree
  std::string dataPath("EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/");
  dataPath = dataPath + dataFileName;
  std::cout<< "dataPath "<< dataPath << std::endl;
  std::cout << "Getting data... "; 
  TFile* data = new TFile(dataPath.c_str(), "READ");
  TCanvas* c1 = new TCanvas("c1","First graph");
  TTree* tree = (TTree*)data->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
  std::cout << "Done." << std::endl;

 // don't apply met>250 for validation region 100<MT2<200 from ht-only prescaled triggers
 std::string metcut = closureTest ? "" : "&& (met>250||ht>1000)";
 std::string cond_all = "(id>150&&(id>152||ht<450)&&(id>153||ht<575)&&(id>154||ht<1000)&&(id>155||ht<1500))" + metcut;
 std::string cond_qcd = "(id<160) && (id>=150)" + metcut; //CONTROLLARE se mettere metcut
  
  //std::cout<<SumWeight(tree, "deltaPhiMin>0.3");
  do_fit(250.0, 450.0, cond_all, cond_qcd, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
  do_fit(450.0, 575.0, cond_all, cond_qcd, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
  do_fit(575.0, 1000.0, cond_all, cond_qcd, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
  do_fit(1000.0, 1500.0, cond_all, cond_qcd, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
  do_fit(1500.0, 13000.0, cond_all, cond_qcd, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
  
  delete c1;
  data->Close();
  return 0;
}


TH1D* set_ratio(float ht_min, float ht_max, std::string  condition, float threshold, TTree* tree, std::string qcd_or_all, bool useMC){

    //cuts:
    std::ostringstream oss_big;
    oss_big <<"("<< condition << "&& ht<" << ht_max << " && ht>" << ht_min <<  " && deltaPhiMin>" << threshold<<")";
    std::string cond_big = oss_big.str();
    
    std::ostringstream oss_small;
    oss_small<<"("<< condition << "&& ht<" << ht_max << " && ht>" << ht_min << " && deltaPhiMin<" <<threshold<<")";
    std::string cond_small = oss_small.str();
    
    //name of output histogram:
    std::ostringstream histo_name;
    histo_name <<"h_ratio_"<< qcd_or_all;
    std::string name = histo_name.str();
    
    //set bins depending on all/onlyqcd and ht region FINIRE
    /*int binnum_a = 15; int binnum_b = 16; int binnum_c = 17;
    double bins_a = {50.0, 55.0,  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 120.0, 200.0, 300.0, 500.0};
    double bins_b = {50.0, 55.0,  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 120.0, 200.0, 300.0, 500.0, 800.};
    double bins_c = {50.0, 55.0,  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 120.0, 200.0, 300.0, 500.0, 800., 1100};
    int binnum = binnum_a;
    double bins[binnum];
    if ((ht_min<=500.0 && qcd_or_all == "all") || (ht_min>500.0 && qcd_or_all = "qcd")) {binnum = binnum_b;}
    else if (qcd_or_all == "all" && ht_min>500.0) {binnum = binnum_c;}
     PROVARE CON VECTOR*/
   
    float scaleMC = 0.; // simulate stats for the given lumi [0 means MC stats]
    
    float prescales[3] = {7900., 440.6, 110.2}; // up to run G 27.7 fb-1
    
    float lumiScale = 1.;
    float sfFromSNT[5] = {1.88375, 1.38677, 1.27216, 1.16806, 1.02239};
    float lumiRatioGtoH = 27.261/35.876;
    
    bool onlyUseUpToRunG = false; //SISTEMARE
    
    
    if( !useMC || scaleMC!=0. ) {
        //lumiScale = useMC ? scaleMC : cfg.lumi(); SISTEMARE
        if     ( ht_min < 300. ) lumiScale *= sfFromSNT[0]/prescales[0];
        else if( ht_min< 500. ) lumiScale *= sfFromSNT[1]/prescales[1];
        else if( ht_min < 600. ) lumiScale *= sfFromSNT[2]/prescales[2];
        else if( ht_min < 1100.) lumiScale *= sfFromSNT[3];
        else                          lumiScale *= sfFromSNT[4];
        if (onlyUseUpToRunG)
            lumiScale *= lumiRatioGtoH;
    }
    
    double bins[] = {50.0, 55.0,  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 125.0, 200.0, 300.0, 500.0};
    int  binnum = sizeof(bins)/sizeof(double) - 1;
    
    TH1D* h_small = new TH1D ("h_small", "h_small", binnum, bins);
    TH1D* h_big = new TH1D ("h_big", "h_big", binnum, bins);
    TH1D* h_ratio = new TH1D (name.c_str(), "h_ratio", binnum, bins);
    
    
    TCut selection_big = cond_big.c_str();
    TCut selection_small = cond_small.c_str();
    
    TCut weight = "weight";
    TCut lumi = "lumi";
    TCut selection_big_and_weight = selection_big*weight;
    TCut selection_small_and_weight = selection_small*weight;
    
    tree->Draw("mt2>>h_big", selection_big_and_weight, "goff");
    tree->Draw("mt2>>h_small", selection_small_and_weight, "goff");
    h_ratio->Divide(h_big, h_small);
    
    //std::cout<<"Ratio was set";
    
    return h_ratio;
}

void do_fit(float ht_min, float ht_max, std::string  cond_all, std::string cond_qcd, TTree* tree, bool useMC){
  //ht_min and ht_max are the boundaries of the ht interval to be considered
  //n_bins is the number of bins in mt2
  //tree is the tree containing the data
  double delta_Phi_threshold = 0.3;
  //limits of fitting region:
  double mt2_min = 60.0; double mt2_max = 100.0;
  if(ht_min>=999.0) {mt2_min = 70.0;}  //stay on the safe side for high ht regions
  
  //max and min mt2 for whole plot (not only region for fitting)
  double mt2_min_global = 50.0; double mt2_max_global = 1000.0;
 
  std::cout<<ht_min<<" GeV < H_T < "<<ht_max;

  TH1D* histo_qcd = set_ratio(ht_min, ht_max, cond_qcd, delta_Phi_threshold , tree, "qcd", useMC);
  TH1D* histo_all = set_ratio(ht_min, ht_max, cond_all, delta_Phi_threshold , tree, "all", useMC);

  TF1 *pow_bg = new TF1("pow_bg", "[0]*pow(x,[1])", mt2_min, mt2_max);
  histo_qcd->Fit("pow_bg", "R 0"); //R: fit in the range of the function. don't plot now the fitted function
 
  TF1 *fitResult = histo_qcd->GetFunction("pow_bg");

  double param_a = fitResult->GetParameter(0);
  double param_b = fitResult->GetParameter(1);
   
  double mt2_min_right = mt2_min + 5.;
  double mt2_max_right = mt2_max + 25.;
  double mt2_min_left = mt2_min - 5.;
  double mt2_max_left  = mt2_max;
    
  TF1 *pow_bg_right = new TF1("pow_bg_right", "[0]*pow(x,[1])", mt2_min_right, mt2_max_right);
  TF1 *pow_bg_left = new TF1("pow_bg_left", "[0]*pow(x,[1])", mt2_min_left, mt2_max_left);
  
  std::cout<<"Right variation:"<<std::endl;
  histo_qcd->Fit("pow_bg_right", "R 0");
  TF1 *fitResult_right = histo_qcd->GetFunction("pow_bg_right");
  
  std::cout<<"Left variation:"<<std::endl;
  histo_qcd->Fit("pow_bg_left", "R 0");
  TF1 *fitResult_left = histo_qcd->GetFunction("pow_bg_left");
  
  double param_b_right = fitResult_right->GetParameter(1);
  double param_b_left = fitResult_left->GetParameter(1);
    
  double var_right = param_b_right/param_b;
  double var_left = param_b_left/param_b;
  double var_max = TMath::Max( fabs(var_right-1.0), fabs(var_left-1.0) );
    
  std::cout << "Slope variation: " << var_right << " (right) / " << var_left << "(left)" << std::endl<< "   maximal variation: " << var_max << std::endl;
    
    
  //now define fitted function over whole range (i.e. not only fitting region)
  TF1 *fitted_bg = new TF1("pow_bg_fit", "[0]*pow(x,[1])", mt2_min_global, mt2_max_global);
  fitted_bg->SetParameter(0, param_a);
  fitted_bg->SetParameter(1, param_b);
  fitted_bg->SetParError(0, fitResult->GetParError(0));
  fitted_bg->SetParError(1, fitResult->GetParError(1));
  fitted_bg->GetXaxis()->SetTitle("M_{T2} [GeV]");
  fitted_bg->GetXaxis()->SetNoExponent();
  fitted_bg->GetXaxis()->SetMoreLogLabels();
  fitted_bg->GetYaxis()->SetTitle("r_{#phi}");
  fitted_bg->SetTitle("CMS simulation, #sqrt{s} = 13 TeV");
  fitted_bg->SetFillColor(29);
  fitted_bg->SetFillStyle(3001);
  fitted_bg->SetMinimum(0.01);
  fitted_bg->SetMaximum(10000.0);
    


  TCanvas* cfit = new TCanvas("cfit","Fit with power law");
  //add legend:
  TLegend* legend = new TLegend(0.68,0.72,0.85,0.86);
  gPad->SetLogx(); 
  gPad->SetLogy();
  //gPad->SetTitle( "CMS simulation, #sqrt{s} = 13 TeV");
  
  
  legend->AddEntry( "pow_bg_fit", "Fit", "FL" );
  fitted_bg->Draw("L E3");  //draw fitted function first since it needs more y range
  //fitResult->Draw("E3");
    
    
  histo_qcd->LabelsOption("h","Y"); //IN QUALCHE MODO IL LABEL E' SEMPRE STORTO
  //histo->SetMarkerStyle(kFullSquare); per qualche motivo non funziona, ottengo una linea continua
  histo_qcd->SetStats(0); //don't show information box on graph
  histo_qcd->SetMarkerStyle(24);
  histo_qcd->SetLineColor(kBlack);
  legend->AddEntry("h_ratio_qcd","Simulation (multijet only)","PL");
  histo_qcd->Draw("SAME P");
    
 histo_all->LabelsOption("h","Y"); //IN QUALCHE MODO IL LABEL E' SEMPRE STORTO
    //histo->SetMarkerStyle(kFullSquare); per qualche motivo non funziona, ottengo una linea continua
 histo_all->SetStats(0); //don't show information box on graph
 histo_all->SetMarkerStyle(20);
 histo_all->SetLineColor(kBlack);
 legend->AddEntry("h_ratio_all","Simulation (all)","PL");
 histo_all->Draw("SAME P");


  
  //CAPIRE: I MONOJET ENTRANO IN QCD??
 
  legend->Draw();

  //add text with descrption of considered region
  std::ostringstream descr;
  descr << ht_min << " GeV < H_{T} < " << ht_max << " GeV";
  std::string description = descr.str(); 
  TPaveText *pt = new TPaveText(0.4,0.8,0.65,0.85, "NDC");
  pt->AddText(description.c_str());
  pt->Draw("SAME");
  
  //add lines at 60 and 100 GeV
  TLine *line1 = new TLine(mt2_min,1e-2, mt2_min, 1e4);
  TLine *line2 = new TLine(mt2_max,1e-2, mt2_max, 1e4);
  line1->SetLineStyle(2); line2->SetLineStyle(2);
  line1->Draw("SAME");
  line2->Draw("SAME");
    
  //add chi2 box
  std::ostringstream descr_chi2;
  descr_chi2 <<"#chi^{2}/ndf = "<<fitResult->GetChisquare()<<"/"<< fitResult->GetNDF()<<": "<< fitResult->GetProb()*100.0<<"%";
  std::string description_chi2 = descr_chi2.str();
  TPaveText *chi2;
  chi2 = new TPaveText(0.22,0.23, 0.5,0.18,"brNDC");
  chi2->SetFillColor(10);   chi2->SetBorderSize(0);
  chi2->AddText(description_chi2.c_str());
  chi2->Draw("SAME");
 
  gPad->RedrawAxis();
  
  
    
  //save plot with file name indicating ht area
  std::ostringstream op_nm;
  op_nm <<"plotFilippo/"<< ht_min << "<HT<" << ht_max;
  if (useMC) {op_nm<<"_MC";}
  else {op_nm<<"_data";}
  op_nm<<".pdf";
  std::string output_name = op_nm.str();
  cfit->SaveAs(output_name.c_str());
    
  //aggiungere pesi secondari
  //cominciare a differenziare mc/data
  //PLOTTARE ANCHE FIT CON RANGE DELLE INCERTEZZE
  //SCRIVERE NEL PLOT I PARAMETRI OTTENUTI? IN OGNI CASO CAPIRE E SCRIVERE IL CHI2/NDF
 //CAPIRE COME FUNZIONA MAKEANALYSISFROMINCLUSIVETREE NELLA CLASSE MT2ESTIMATE
    //capire bene perche' i pesi in MC
    //mt2min = 70 per htmin>1000
    //aggiungere il band con incertezza sul fit, capire se è l'incertezza statistica o quella ottenuta variando il range, o entrambe
  //delete cfit;
  //delete pow_bg;
};



/*float lumiRatioGtoH = 27.261/35.876;

float prescales[3] = {7900., 440.6, 110.2}; // up to run G 27.7 fb-1

float lumiScale = 1.;
float sfFromSNT[5] = {1.88375, 1.38677, 1.27216, 1.16806, 1.02239};


if( !useMC || scaleMC!=0. ) {
    lumiScale = useMC ? scaleMC : cfg.lumi();
    if     ( iR->htMin() < 300. ) lumiScale *= sfFromSNT[0]/prescales[0];
    else if( iR->htMin() < 500. ) lumiScale *= sfFromSNT[1]/prescales[1];
    else if( iR->htMin() < 600. ) lumiScale *= sfFromSNT[2]/prescales[2];
    else if( iR->htMin() < 1100.) lumiScale *= sfFromSNT[3];
    else                          lumiScale *= sfFromSNT[4];
    if (onlyUseUpToRunG)
            lumiScale *= lumiRatioGtoH;
    }*/



/*
The measurements
of rb are performed in data using the unprescaled HLT PFHT900 trigger. For the measurements
of fj in the high and extreme HT regions, the unprescaled HLT PFHT900 trigger is used while
for the medium and low HT regions the prescaled HLT PFHT475 and HLT PFHT350 triggers
are used. In the 2016 run the new prescaled trigger path HLT PFHT125 has been included
allowing for this measurment also in the very low HT region. The small contribution from
non-QCD procesess (5–16 %) are subtracted using MC.
 */


/*
 bands[*iR] = new TH1D(Form("band_%s",iR->getName().c_str()), "", 500, matchedEstimate->lDphi->GetXaxis()->GetXmin(), matchedEstimate->lDphi->GetXaxis()->GetXmax());
 (TVirtualFitter::GetFitter())->GetConfidenceIntervals(bands[*iR], 0.68);
 
 */
