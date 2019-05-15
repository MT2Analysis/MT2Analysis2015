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



int main(int argc, char* argv[] ) {
    
    std::cout<<"Importing mc and data ntuples...";
    TFile* data = new TFile("EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/data.root", "READ");
    TTree* tree_data = (TTree*)data->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
    TFile* mc = new TFile("EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/mc.root", "READ");
    TTree* tree_mc = (TTree*)mc->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
    std::cout<<".done"<<std::endl;
    
    
    //first: closure test. check internal consistency of the method.
    
    double params_a_mc[5]    = {5.66652e+03, 6.43300e+03, 2.92651e+03,  3.30363e+03, 1.94651e+02};
    double params_b_mc[5]    = {-1.83833e+00, -1.90208e+00 , -1.74972e+00, -1.78131e+00, -1.11495e+00};
    double params_a_mc_left[5] = {2.74996e+03,  2.69364e+03, 0, 1.78526e+03, 1.12277e+02}; //NOTA:0 AL POSTO DI NAN
    double params_b_mc_left[5] = {-1.67137e+00, -1.70194e+00, 0, -1.64285e+00, -9.91439e-01}; //NOTA:0 AL POSTO DI NAN
    double params_a_mc_right[5] = {1.50507e+04, 1.04282e+04, 1.10357e+04, 4.36061e+03, 4.77744e+02};
    double params_b_mc_right[5] = {-2.06268e+00, -2.00832e+00, -2.05404e+00, -1.84418e+00, -1.31721e+00};
    double errs_a_mc[5] = {2.10810e+03, 3.44801e+03, 7.02880e+02, 1.53675e+03, 8.25640e+01};
    double errs_b_mc[5] = {8.75381e-02, 1.25773e-01, 5.62763e-02, 1.05922e-01, 9.63834e-02};
    double errs_a_mc_left[5] = {7.07395e+02, 9.99335e+02, 0, 5.72578e+02, 3.32497e+01}; //NOTA:0 AL POSTO DI NAN
    double errs_b_mc_left[5] = {6.15207e-02, 8.84721e-02, 0, 7.38730e-02, 6.80277e-02}; //NOTA:0 AL POSTO DI NAN
    double errs_a_mc_right[5] = {6.13545e+03, 6.24749e+03, 2.53558e+03, 1.74272e+03, 1.45854e+02};
    double errs_b_mc_right[5] = {9.40539e-02, 1.37738e-01, 5.26888e-02, 8.92092e-02, 6.75711e-02};
    
    
    
    TCut weight = "weight";
    std::string metcut = "&& (met>250 || ht>1000)";  //cut met for low ht regions
    //std::string metcut = "";
    std::string cut_big_string = "deltaPhiMin>0.3 && id>100 && id<200 && mt2>100 && mt2<200 && nJets>1" + metcut;
    TCut cut_big = cut_big_string.c_str();
    cut_big = cut_big*weight;
    std::string cut_small_string = "deltaPhiMin<0.3 && id>100 && id<200 && mt2>100 && mt2<200 && nJets>1" ; //DOVE METTERE METCUT??
    TCut cut_small = cut_small_string.c_str();
    cut_small = cut_small*weight;
    
    double bins_ht[] = {250.0, 450.0, 575.0, 1000.0, 1500.0, 13000.0};
    int binnum = 5;
    
    std::cout<<"Creating histogram from mc simulation...";
    TH1D* big_mc = new TH1D ("big_mc", "big_mc", binnum, bins_ht);
    tree_mc->Draw("ht>>big_mc", cut_big, "goff");
    TH1D* small_mc = new TH1D ("small_mc", "small_mc", binnum, bins_ht);
    tree_mc->Draw("ht>>small_mc", cut_small, "goff");
    std::cout<<".done"<<std::endl;
    
    TH1D* ratios = new TH1D("ratios", "ratios", binnum, bins_ht);
    for(int i = 0; i<5; ++i) {
        double mt2 = 0.5*(bins_ht[i] + bins_ht[i+1]); //choose bin centre
        double a = params_a_mc[i];
        double b = params_b_mc[i];
        double r = a*TMath::Power(mt2,b);
        ratios->SetBinContent(i+1, r); //CONTROLLARE CHE SIA IL BIN GIUSTO
    }
    
    for(int i = 1; i<=5; ++i) {
        std::cout<<"bin no. "<<i<<" small_mc: "<<small_mc->GetBinContent(i)<<" big_mc: "<<big_mc->GetBinContent(i)<< " ratio: "<< ratios->GetBinContent(i)<<std::endl;
    }
    
    small_mc->Multiply(ratios);
    
    //debug:
    for(int i = 1; i<=5; ++i) {
        std::cout<<"bin no. "<<i<<" small_mc*ratio: "<<small_mc->GetBinContent(i)<<std::endl;
    }
    
    //start plotting; create canvas
    TCanvas* c1 = new TCanvas("closure_test","Closure Test");
    std::cout<<"Canvas was created"<<std::endl;
    //gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0); //don't show stats box
    //get logarithmic graph:
    gPad->SetLogx();
    gPad->SetLogy();
    
    big_mc->SetTitle("Closure Test - MC:  #sqrt{s} = 13 TeV, 100 GeV < M_{T2} < 200 GeV");
    big_mc->GetXaxis()->SetTitle("H_{T} [GeV]");
    big_mc->GetXaxis()->SetNoExponent();
    big_mc->GetXaxis()->SetMoreLogLabels();
    big_mc->SetFillColor(50);
    big_mc->Draw("bar");
    
    small_mc->SetFillColor(8);
    small_mc->Draw("E1 same");
    
    TLegend* legend = new TLegend(0.62,0.70,0.83,0.85);
    legend->AddEntry(big_mc, "MC, high DPhi from simulation", "FL");
    legend->AddEntry(small_mc, "MC, high DPhi from ratio", "FL");
    legend->Draw("same");
    
    
    
    c1->SaveAs("plotFilippo/closure_test.pdf");
    std::cout<<"Closure test histogram was plotted"<<std::endl;
    
    delete c1;
    delete big_mc;
    delete small_mc;
    delete data;
    delete mc;
    //USARE ORA IL RATIO
    //STILE E1 NON FUNZIONA
    //ragionare se ratio up e down sia il modo piu furbo di definire l'errore.
    //fare prima plot senza errore per vedere che tutto funzioni
    return 0;
}


/*
 double ratio(double mt2, int k) {
 //k: index of the corresponding ht range
 double a = params_a_mc[k];
 double b = params_b_mc[k];
 double r = a*TMath::Power(mt2,b);
 return r;
 }
 
 //now upper/lower bounds on ratio. RAGIONARE SU QUALI ERRORI PRENDERE
 double ratio_up(double mt2, int k) {
 //k: index of the corresponding ht range
 double a = params_a_mc[k] + errs_a_mc[k];
 double b = params_b_mc[k] + errs_b_mc[k];
 double r_up = a*TMath::Power(mt2,b);
 return r_up;
 }
 
 double ratio_down(double mt2, int k) {
 //k: index of the corresponding ht range
 double a = params_a_mc[k] - errs_a_mc[k];
 double b = params_b_mc[k] - errs_b_mc[k];
 double r_down = a*TMath::Power(mt2,b);
 return r_down;
 }
 
 */
