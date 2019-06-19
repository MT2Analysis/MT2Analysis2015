 #include <iostream>
#include <map>
#include <string>
#include <sstream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "THStack.h"
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
    //TH1::SetDefaultSumw2(); //get correct errors
    TFile* data = new TFile("EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/data.root", "READ");
    TTree* tree_data = (TTree*)data->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
    TFile* mc = new TFile("EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/mc.root", "READ");
    TTree* tree_mc = (TTree*)mc->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
    std::cout<<".done"<<std::endl;
    
    
    //first: closure test. check internal consistency of the method.
    
    /*double params_a_mc[5]    = {5.66652e+03, 6.43300e+03, 2.92651e+03,  3.30363e+03, 1.94651e+02};
    double params_b_mc[5]    = {-1.83833e+00, -1.90208e+00 , -1.74972e+00, -1.78131e+00, -1.11495e+00};
    double params_a_mc_left[5] = {2.74996e+03,  2.69364e+03,   1.84347e+03, 1.78526e+03, 1.12277e+02};
    double params_b_mc_left[5] = {-1.67137e+00, -1.70194e+00, -1.64348e+00, -1.64285e+00, -9.91439e-01};
    double params_a_mc_right[5] = {1.50507e+04, 1.04282e+04, 1.10357e+04, 4.36061e+03, 4.77744e+02};
    double params_b_mc_right[5] = {-2.06268e+00, -2.00832e+00, -2.05404e+00, -1.84418e+00, -1.31721e+00};
    double errs_a_mc[5] = {2.10810e+03, 3.44801e+03, 7.02880e+02, 1.53675e+03, 8.25640e+01};
    double errs_b_mc[5] = {8.75381e-02, 1.25773e-01, 5.62763e-02, 1.05922e-01, 9.63834e-02};
    double errs_a_mc_left[5] = {7.07395e+02, 9.99335e+02, 3.10767e+02, 5.72578e+02, 3.32497e+01};
    double errs_b_mc_left[5] = {6.15207e-02, 8.84721e-02,  4.01412e-02, 7.38730e-02, 6.80277e-02};
    double errs_a_mc_right[5] = {6.13545e+03, 6.24749e+03, 2.53558e+03, 1.74272e+03, 1.45854e+02};
    double errs_b_mc_right[5] = {9.40539e-02, 1.37738e-01, 5.26888e-02, 8.92092e-02, 6.75711e-02};
    
    double params_a_data[5]            ={ 9.92889e+02, 9.56243e+04, 2.27354e+04, 1.86886e+03, 1.34349e+02};
    double params_b_data[5]            = {-1.47510e+00, -2.57526e+00, -2.25777e+00, -1.65928e+00, -1.03742e+00};
    double params_a_data_left[5]    ={6.99719e+03, 2.62506e+04, 1.17293e+04, 1.07332e+03, 1.32561e+02};
    double params_b_data_left[5]    = {-1.93024e+00, -2.27713e+00, -2.10556e+00, -1.53450e+00, -1.03441e+00};
    double params_a_data_right[5] ={6.73542e+06, 3.76859e+05, 2.05164e+05, 8.23527e+03, 4.14703e+02};
    double params_b_data_right[5] = {-3.51966e+00, -2.88207e+00, -2.76172e+00, -1.99285e+00, -1.29087e+00};
    double errs_a_data[5]                   = {3.19022e+03, 1.30488e+05, 1.45456e+04, 5.39829e+02, 7.98080e+01};
    double errs_b_data[5]                   = {7.61119e-01, 3.21342e-01,  1.50143e-01, 6.58034e-02, 1.35009e-01};
    double errs_a_data_left[5]           = {1.70499e+04, 2.32737e+04 , 5.24764e+03, 2.12914e+02, 5.44539e+01};
    double errs_b_data_left[5]           = {1.70499e+04, 2.12124e-01, 1.06629e-01, 4.57170e-02, 9.43923e-02};
    double errs_a_data_right[5]        = {3.11195e+07, 7.74747e+05, 1.47369e+05, 1.89176e+03, 1.77681e+02};
    double errs_b_data_right[5]        = {1.07970e+00, 4.75440e-01, 1.65481e-01, 5.11765e-02, 9.49688e-02};*/
    
    double params_a_mc[5] = { 4402.027000410028 , 4294.550890222096 , 2869.529313533573 , 2615.214162416537 , 183.71176863956313 };
    double params_b_mc[5] = { -1.741081765408066 , -1.7757538316648727 , -1.7056591920648638 , -1.7120602674860594 , -1.0674835839956383 };
    double err_a_mc_sq[5] = { 24212982.52135981 , 14473352.74981307 , 11664681.575123249 , 1927522.778022502 , 28054.921498601576 };
    double err_b_mc_sq[5] = { 0.029574226373957638 , 0.027358498052237848 , 0.024776716485000974 , 0.012715761381638133 , 0.023593787915964892 };
    double err_ab_mc_sq[5] = { -816.902118956038 , -620.4191332989689 , -522.8320154039699 , -155.92028037602205 , -25.256672080354388 };
    double params_a_data[5] = { 2304.0694758232594 , 45206.11470871173 , 19505.601140819766 , 1816.4281452876933 , 176.6523909711623 };
    double params_b_data[5] = { -1.5284622569141635 , -2.336172386578419 , -2.161617895840904 , -1.6013903195345762 , -1.074963049360108 };
    double err_a_data_sq[5] = { 19891243439.174774 , 11636842552.162848 , 2203686608.6733046 , 7093327.817706462 , 17044.158080991932 };
    double err_b_data_sq[5] = { 0.1052328483130024 , 0.06397564126293889 , 0.03293057413276392 , 0.03149085808415257 , 0.014022881408129819 };
    double err_ab_data_sq[5] = { -7242.474929831113 , -24613.6183971168 , -7703.614213548204 , -456.045884014989 , -15.45953258489528 };

    
    
    TCut weight = "weight";
    //std::string metcut = "&& (met>250 || ht>1000)";  //cut met for low ht regions
    //std::string metcut = "";
    std::string cut_big_string = "deltaPhiMin>0.3 && id>100 && id<200 && mt2>100 && mt2<200 && nJets>1";
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
    
    //make sure we have the right errors (for some strange reason it's not the case...)
   /* for(int i = 1; i<=5; ++i) {
        small_mc->SetBinError(i, TMath::Sqrt(small_mc->GetBinContent(i)));
        big_mc->SetBinError(i, TMath::Sqrt(big_mc->GetBinContent(i)));
    }*/
    
    std::cout<<".done"<<std::endl;
    
    std::cout<<"CLOSURE TEST"<<std::endl;
    TH1D* ratios_mc = new TH1D("ratios_mc", "ratios_mc", binnum, bins_ht);
    for(int i = 0; i<5; ++i) {
        double a = params_a_mc[i];
        double b = params_b_mc[i];
        double r = a*(TMath::Power(2.0, b + 1.0) - 1.0)*TMath::Power(100.0, b)/(b + 1.0);
        
        /*double err_a = fabs(params_a_mc_right[i] - a);
         double err_b = fabs(params_b_mc_right[i] - b);
         if (fabs(params_a_mc_left[i] - a) > err_a  ) err_a = fabs(params_a_mc_left[i] - a);
         if (fabs(params_b_mc_left[i] - a) > err_a  ) err_b = fabs(params_b_mc_left[i] - b);*/
        double der_a = TMath::Power(100.0, b + 1.0)*(-1.0 + TMath::Power(2.0, b + 1.0))/(1.0 + b);
        double der_b = -a*TMath::Power(100.0, b + 1.0)*(-1.0 + TMath::Power(2.0, b + 1.0))/TMath::Power(1.0 + b,2.0) + a*TMath::Power(200.0, 1.0+b)*std::log(2.0)/(1.0 + b)  + a*TMath::Power(100.0, 1.0 + b)*(-1.0 + TMath::Power(2.0, 1+b))*std::log(100.0)/(1.0 + b);
        der_a = der_a/100.0;
        der_b = der_b/100.0;
        //double err_r = TMath::Sqrt( der_a*der_a*err_a*err_a  +  der_b*der_b*err_b*err_b  );
        double err_r = TMath::Sqrt( der_a*der_a*err_a_mc_sq[i]  +  der_b*der_b*err_b_mc_sq[i]  + der_a*der_b*err_ab_mc_sq[i]);
        ratios_mc->SetBinContent(i+1, r);
        ratios_mc->SetBinError(i+1, err_r);
        std::cout<<"bin no. "<<i<<" der_a: "<<der_a<<" der_b: "<<der_b<<std::endl;
        std::cout<<" ratio: "<< ratios_mc->GetBinContent(i+1)<<" err ratio: "<<ratios_mc->GetBinError(i+1)<<std::endl;
    }
    
  
    small_mc->Multiply(ratios_mc);
    
    //debug:
    for(int i = 1; i<=5; ++i) {
        std::cout<<"bin no. "<<i<<" small_mc*ratio: "<<small_mc->GetBinContent(i)<<std::endl;
        std::cout<<"bin no. "<<i<<" relative error: "<<small_mc->GetBinError(i)/small_mc->GetBinContent(i)<<std::endl;
    }
    
    

    
    //start plotting; create canvas
    TCanvas* c1 = new TCanvas("closure_test","Closure Test");
    std::cout<<"CLosure test: canvas was created"<<std::endl;
    //gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0); //don't show stats box
    //get logarithmic graph:
    gPad->SetLogx();
    gPad->SetLogy();
    //c1->SetBottomMargin(0.2);
    
    big_mc->SetTitle("Closure Test - MC:  #sqrt{s} = 13 TeV, 100 GeV < M_{T2} < 200 GeV");
    big_mc->GetXaxis()->SetTitle("H_{T} [GeV]");
    big_mc->GetXaxis()->SetNoExponent();
    big_mc->GetXaxis()->SetMoreLogLabels();
    big_mc->SetFillColor(50);
    big_mc->GetYaxis()->SetRange(0.1, 10000.);
    big_mc->SetMaximum(22000.);
    big_mc->GetYaxis()->SetTitle("N. Events");
    big_mc->GetYaxis()->SetTitleSize(17);
    big_mc->GetYaxis()->SetTitleFont(43);
    big_mc->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    big_mc->GetYaxis()->SetLabelSize(15);
   
    //big_mc->GetHistogram()->SetMaximum(10000.);   // non funziona!!
    //big_mc->GetHistogram()->SetMinimum(0.1);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 0.98);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetTopMargin(0.17);
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    
    gPad->SetLogx();
    gPad->SetLogy();
    
    big_mc->Draw("bar");
    
    small_mc->SetMarkerStyle(20);
    small_mc->SetMarkerColor(1);
    small_mc->SetFillColor(1);
    small_mc->SetLineColor(1);
    small_mc->Draw("EPSAME");
    
    TLegend* legend1 = new TLegend(0.62,0.67,0.83,0.82);
    legend1->AddEntry(big_mc, "MC, high DPhi from simulation", "FL");
    legend1->AddEntry(small_mc, "MC, high DPhi from ratio", "FL");
    legend1->Draw("same");
    
    c1->cd();
    //now plot lower bar with ratio
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.30);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.38);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    gPad->SetLogx();
    
    TH1D* ratio_exp_to_actual = new TH1D ("ratio_exp_to_actual", "ratio_exp_to_actual", binnum, bins_ht);
    ratio_exp_to_actual->Divide(small_mc, big_mc);
    ratio_exp_to_actual->SetLineColor(kBlack);
    ratio_exp_to_actual->SetMinimum(0.0);  // Define Y ..
    ratio_exp_to_actual->SetMaximum(2.5); // .. range
    ratio_exp_to_actual->Sumw2();
    ratio_exp_to_actual->SetStats(0);      // No statistics on lower plot
    ratio_exp_to_actual->SetMarkerStyle(21);
    ratio_exp_to_actual->SetTitle("");  //remove ratio title
    ratio_exp_to_actual->GetYaxis()->SetTitle("ratio exp/mc ");
    ratio_exp_to_actual->GetXaxis()->SetTitle("H_{T} [GeV]");
    ratio_exp_to_actual->Draw("ep");
    
    ratio_exp_to_actual->GetYaxis()->SetNdivisions(505);
    ratio_exp_to_actual->GetYaxis()->SetTitleSize(17);
    ratio_exp_to_actual->GetYaxis()->SetTitleFont(43);
    ratio_exp_to_actual->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_exp_to_actual->GetYaxis()->SetLabelSize(15);
    
    // X axis ratio plot settings
    ratio_exp_to_actual->GetXaxis()->SetTitleSize(17);
    ratio_exp_to_actual->GetXaxis()->SetTitleFont(43);
    ratio_exp_to_actual->GetXaxis()->SetTitleOffset(4.2);
    ratio_exp_to_actual->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_exp_to_actual->GetXaxis()->SetLabelSize(15);


    c1->cd();
    c1->SaveAs("plotFilippo/closure_test.pdf");
    
    std::cout<<"Closure test histogram was plotted"<<std::endl;
    
    
    
    
    
    
    
    
    //NOW PERFORM VALIDATION TEST
    std::cout<<"VALIDATION TEST"<<std::endl;
    
    //std::string metcut = "&& (ht>1000 || met<250)";
    std::string metcut = "";
    
    std::string cond_data = "&& nJets>1 && (id&1==1) && run<=280385"  ;
    
    std::string cut_big_string_nonqcd = "nJets>1 && deltaPhiMin>0.3 && id>300 && mt2>100 && mt2<200" + metcut;  //id>300 o id>200?????
    TCut cut_big_nonqcd = cut_big_string_nonqcd.c_str();
    
    std::string cut_small_string_nonqcd = "nJets>1 && deltaPhiMin<0.3 && id>300 && mt2>100 && mt2<200" + metcut;
    TCut cut_small_nonqcd = cut_small_string_nonqcd.c_str();
    
    std::string cut_small_string_tot = "deltaPhiMin<0.3 && mt2>100 && mt2<200" + metcut + cond_data;
    TCut cut_small_tot = cut_small_string_tot.c_str();  //METCUT???
    
    std::string cut_big_string_tot = "deltaPhiMin>0.3 && mt2>100 && mt2<200" + metcut + cond_data;
    TCut cut_big_tot = cut_big_string_tot.c_str();  //METCUT???
    
    TH1D* big_mc_nonqcd = new TH1D ("big_mc_nonqcd", "big_mc_nonqcd", binnum, bins_ht);
    tree_mc->Draw("ht>>big_mc_nonqcd", cut_big_nonqcd*weight, "goff");
    
    TH1D* small_mc_nonqcd = new TH1D ("small_mc_nonqcd", "small_mc_nonqcd", binnum, bins_ht);
    tree_mc->Draw("ht>>small_mc_nonqcd", cut_small_nonqcd*weight, "goff");
    
    TH1D* small_data_tot = new TH1D("small_data_tot", "small_data_tot", binnum, bins_ht);
    tree_data->Draw("ht>>small_data_tot", cut_small_tot, "goff");
    
    TH1D* big_data_tot = new TH1D("big_data_tot", "big_data_tot", binnum, bins_ht);
    tree_data->Draw("ht>>big_data_tot", cut_big_tot, "goff");
    
    
    //make sure we have the right errors (for some strange reason it's not the case...)
    /*for(int i = 1; i<=5; ++i) {
        big_data_tot->SetBinError(i, TMath::Sqrt(big_data_tot->GetBinContent(i)));
        big_mc_nonqcd->SetBinError(i, TMath::Sqrt(big_mc_nonqcd->GetBinContent(i)));
        small_mc_nonqcd->SetBinError(i, TMath::Sqrt(small_mc_nonqcd->GetBinContent(i)));
        small_data_tot->SetBinError(i, TMath::Sqrt(small_data_tot->GetBinContent(i)));
    }*/
    
    //debug: check errors:
    for(int i = 1; i<=5; ++i) {
        std::cout<<"bin no. "<<i<<" big_data_tot: "<<big_data_tot->GetBinContent(i)<< " err: "<<big_data_tot->GetBinError(i)<<std::endl;
    }
    
    //divide by prescales:
    bool onlyUseUpToRunG = true;
    double lumiRatioGtoH = 27.261/35.876;
    double prescales[3] = {7900., 440.6, 110.2}; // up to run G 27.7 fb-1
    double lumiScale = 35.876;
    if (onlyUseUpToRunG)
        lumiScale *= lumiRatioGtoH;
    double sfFromSNT[5] = {1.88375, 1.38677, 1.27216, 1.16806, 1.02239};
   
    TH1D* prescales_histo = new TH1D("prescales_histo", "prescales_histo", binnum, bins_ht);
    prescales_histo->SetBinContent(1, lumiScale*sfFromSNT[0]/prescales[0]);
    prescales_histo->SetBinContent(2, lumiScale*sfFromSNT[1]/prescales[1]);
    prescales_histo->SetBinContent(3, lumiScale*sfFromSNT[2]/prescales[2]);
    prescales_histo->SetBinContent(4, lumiScale*sfFromSNT[3]);
    prescales_histo->SetBinContent(5, lumiScale*sfFromSNT[4]);
    //no prescale for the two highest ht regions
    
    
    //multiply with prescales:
    big_mc_nonqcd->Multiply(prescales_histo);
    small_mc_nonqcd->Multiply(prescales_histo);
    
    TH1D* small_data_qcd = new TH1D("small_data_qcd", "small_data_qcd", binnum, bins_ht);
    small_data_qcd->Add(small_data_tot, small_mc_nonqcd, 1.0, -1.0);
    
    for(int i = 1; i<=5; ++i) {
        std::cout<<"Bin no. "<<i<<" mc non-qcd estimate DPhi>0.3: "<<big_mc_nonqcd->GetBinContent(i)<<" mc non-qcd estimate DPhi<0.3: "<<small_mc_nonqcd->GetBinContent(i)<<" data tot  Dphi<0.3  "<<small_data_tot->GetBinContent(i)<<" data qcd DPhi<0.3: "<<small_data_qcd->GetBinContent(i)<<std::endl;
    }
    
    
    //compute ratios: see mc case for comments
    TH1D* ratios_data = new TH1D("ratios_data", "ratios_data", binnum, bins_ht);
    for(int i = 0; i<5; ++i) {
        double a = params_a_data[i];
        double b = params_b_data[i];
        double r = a*(TMath::Power(2.0, b + 1.0) - 1.0)*TMath::Power(100.0, b)/(b + 1.0);
        
        /*double err_a = fabs(params_a_data_right[i] - a);
        double err_b = fabs(params_b_data_right[i] - b);
        if (fabs(params_a_data_left[i] - a) > err_a  ) err_a = fabs(params_a_data_left[i] - a);
        if (fabs(params_b_data_left[i] - a) > err_a  ) err_b = fabs(params_b_data_left[i] - b);*/
        double der_a = TMath::Power(100.0, b + 1.0)*(-1.0 + TMath::Power(2.0, b + 1.0))/(1.0 + b);
        double der_b = -a*TMath::Power(100.0, b + 1.0)*(-1.0 + TMath::Power(2.0, b + 1.0))/TMath::Power(1.0 + b,2.0) + a*TMath::Power(200.0, 1.0+b)*std::log(2.0)/(1.0 + b)  + a*TMath::Power(100.0, 1.0 + b)*(-1.0 + TMath::Power(2.0, 1+b))*std::log(100.0)/(1.0 + b);
        der_a = der_a/100.0;
        der_b = der_b/100.0;
        //double err_r = TMath::Sqrt( der_a*der_a*err_a*err_a  +  der_b*der_b*err_b*err_b  );
        double err_r = TMath::Sqrt( der_a*der_a*err_a_data_sq[i]  +  der_b*der_b*err_b_data_sq[i]  + der_a*der_b*err_ab_data_sq[i]);
        ratios_data->SetBinContent(i+1, r);
        ratios_data->SetBinError(i+1, err_r);
        std::cout<<"bin no. "<<i<<" der_a: "<<der_a<<" der_b: "<<der_b<<std::endl;
        std::cout<<" ratio: "<< ratios_data->GetBinContent(i+1)<<" err ratio: "<<ratios_data->GetBinError(i+1)<<std::endl;
        }
    
    //find estimate of number of events in high dphi region using data:
    TH1D* big_data_qcdest = new TH1D("big_data_qcdest", "big_data_qcdest", binnum, bins_ht);
    big_data_qcdest->Multiply(small_data_qcd, ratios_data);
    
    //debug:
    for(int i = 1; i<=5; ++i) {
        std::cout<<"Bin no. "<<i<<" ratio: "<<ratios_data->GetBinContent(i)<<" non-qcd estimate DPhi>0.3: "<<big_mc_nonqcd->GetBinContent(i)<<" qcd estimate DPhi>0.3: "<<big_data_qcdest->GetBinContent(i)<<" tot expected: "<<big_data_tot->GetBinContent(i)<<std::endl;
    }
    
    //now do the plot:
    TCanvas* c2 = new TCanvas("validation_test", "Validation Test");
    std::cout<<"CLosure test: canvas was created"<<std::endl;
    //gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0); //don't show stats box
    //get logarithmic graph:
    gPad->SetLogx();
    gPad->SetLogy();
    
    TPad *pad3 = new TPad("pad1", "pad1", 0, 0.31, 1, 0.98);
    pad3->SetBottomMargin(0); // Upper and lower plot are joined
    pad3->SetTopMargin(0.17);
    pad3->Draw();             // Draw the upper pad: pad3
    pad3->cd();               // pad3 becomes the current pad
    
    
    gPad->SetLogx();
    gPad->SetLogy();
    
    
    big_mc_nonqcd->SetFillColor(9);
    //big_mc_nonqcd->GetYaxis()->SetTitle("N. Events"); //non funziona!!
    //big_mc_nonqcd->GetYaxis()->SetRange(0.1, 25000.);
    //big_mc->GetHistogram()->SetMaximum(10000.);   // non funziona!!
    //big_mc->GetHistogram()->SetMinimum(0.1);
    
    big_data_qcdest->SetFillColor(50);
    
    THStack *hs = new THStack("hs","Validation Test");
    hs->SetTitle("Validation Test:  #sqrt{s} = 13 TeV, 100 GeV < M_{T2} < 200 GeV");
    
    hs->SetMaximum(20000.);
    hs->SetMinimum(70.0);
    
    
    hs->Add(big_mc_nonqcd);
    //big_mc_nonqcd->Draw("bar");
    hs->Add(big_data_qcdest);
    //big_data_qcdest->Draw("bar same");
    
    hs->Draw("bar p");
    hs->GetYaxis()->SetTitle("N. Events");
    //hs->GetXaxis()->SetTitle("H_{T} [GeV]");
    hs->GetXaxis()->SetNoExponent();
    hs->GetXaxis()->SetMoreLogLabels();
    hs->GetYaxis()->SetNdivisions(505);
    hs->GetYaxis()->SetTitleSize(17);
    hs->GetYaxis()->SetTitleFont(43);
    hs->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hs->GetYaxis()->SetLabelSize(15);
    gPad->Modified();
    
    big_data_tot->SetMarkerStyle(20);
    big_data_tot->SetMarkerColor(1);
    big_data_tot->SetMarkerSize(1);
    //big_data_tot->SetFillColor(1);
    big_data_tot->SetLineColor(1);
    
    //make sure we have the right errors (for some strange reason it's not the case...)
    /*for(int i = 1; i<=5; ++i) {
        big_data_tot->SetBinError(i, TMath::Sqrt(big_data_tot->GetBinContent(i)));
    }*/
    big_data_tot->Draw("SAME E1");
   
    
    TLegend* legend2 = new TLegend(0.17,0.67,0.34,0.82);
    legend2->AddEntry(big_mc_nonqcd, "Non multijet simulation", "FL");
    legend2->AddEntry(big_data_qcdest, "Multijet prediction", "FL");
    legend2->AddEntry(big_data_tot, "Data", "FL");

    legend2->Draw("same");
    
    c2->cd();
    //now plot lower bar with ratio
    TPad *pad4 = new TPad("pad4", "pad4", 0, 0.08, 1, 0.30);
    pad4->SetTopMargin(0);
    pad4->SetBottomMargin(0.38);
    pad4->Draw();
    pad4->cd();       // pad2 becomes the current pad
    
    gPad->SetLogx();
    
    TH1D* sim = new TH1D ("sim", "sim", binnum, bins_ht);
    sim->Add(big_data_qcdest, big_mc_nonqcd);
    
    TH1D * data_new = (TH1D*)big_data_tot->Clone("data_new");
    TH1D *sim_new = (TH1D*)sim->Clone("sim_new");
    
    TH1D* ratio_meas_sim = new TH1D ("ratio_meas_sim", "ratio_meas_sim", binnum, bins_ht);
    ratio_meas_sim->Divide(data_new, sim_new);
    ratio_meas_sim->SetLineColor(kBlack);
    ratio_meas_sim->SetMinimum(0.0);  // Define Y ..
    ratio_meas_sim->SetMaximum(2.5); // .. range
    ratio_meas_sim->Sumw2();
    ratio_meas_sim->SetStats(0);      // No statistics on lower plot
    ratio_meas_sim->SetMarkerStyle(21);
    ratio_meas_sim->SetTitle("");  //remove ratio title
    ratio_meas_sim->GetYaxis()->SetTitle("ratio meas/sim");
    ratio_meas_sim->GetXaxis()->SetTitle("H_{T} [GeV]");
    ratio_meas_sim->Draw("ep");
    
    ratio_meas_sim->GetYaxis()->SetNdivisions(505);
    ratio_meas_sim->GetYaxis()->SetTitleSize(17);
    ratio_meas_sim->GetYaxis()->SetTitleFont(43);
    ratio_meas_sim->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_meas_sim->GetYaxis()->SetLabelSize(15);
    
    // X axis ratio plot settings
    ratio_meas_sim->GetXaxis()->SetTitleSize(20);
    ratio_meas_sim->GetXaxis()->SetTitleFont(43);
    ratio_meas_sim->GetXaxis()->SetTitleOffset(4.2);
    ratio_meas_sim->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_meas_sim->GetXaxis()->SetLabelSize(15);
    
    
    c2->cd();
    
    
    c2->SaveAs("plotFilippo/validation_test.pdf");
    std::cout<<"Validation test histogram was plotted"<<std::endl;
    
    
    
    //CONTROLLARE SE MANCHINO DELETE
    delete legend1;
    delete legend2;
    delete pad1;
    delete pad2;
    delete pad3;
    delete pad4;
    delete c1;
    delete c2;
    delete big_data_qcdest;
    delete small_data_qcd;
    delete ratios_data;
    delete ratios_mc;
    delete prescales_histo;
    delete small_mc_nonqcd;
    delete small_data_tot;
    delete big_data_tot;
    delete big_mc;
    delete sim;
    delete sim_new;
    delete data_new;
    delete ratio_meas_sim;
    delete big_mc_nonqcd;
    delete hs;
    delete small_mc;
    delete data;
    delete mc;
    
    //PROBLEMA: big_mc_nonqcd HA VALORI TROPPO ALTI!!!!!
    //nota: risultati attesi leggermente inferiori della stima attesa. ha senso:
    //TRAtioPlot    : AGGIUNGERE BANDA INFERIORE CON RAPPORTO
    //nei plot sistemare ticks dell'asse x in modo che si distinguano chiaramente le cinque aree di ht
    //CONTROLLARE id>300 o id>200?????
    //SetNdivisions
    //aggiungere barra inferiore in validationq
    
    
    return 0;
}







