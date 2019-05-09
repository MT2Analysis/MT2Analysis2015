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


void do_fit(double ht_min, double ht_max, std::string  cond_all, std::string cond_qcd, TTree* tree, bool useMC);
TH1D* set_ratio(double ht_min, double ht_max, std::string  condition, double threshold, TTree* tree);


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
  TTree* tree = (TTree*)data->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
  std::cout << "Done." << std::endl;

 // don't apply met>250 for validation region 100<MT2<200 from ht-only prescaled triggers
 std::string metcut = closureTest ? "" : "&& (met>250||ht>1000)";
 std::string cond_all_mc = "(id>150&&(id>152||ht<450)&&(id>153||ht<575)&&(id>154||ht<1000)&&(id>155||ht<1500))" + metcut;
 std::string cond_qcd_mc = "(id<160) && (id>=150)" + metcut; //CONTROLLARE se mettere metcut

//obviously for the data case I can't just define a condition to get qcd, so the two conditions are formally the same
//(keep them both to keep the same function structure as in mc)
 std::string cond_all_data = "nJets>1" ;  //SISTEMARE CUT!!!
 std::string cond_qcd_data = "nJets>1" ; //CONTROLLARE se mettere metcut
    
  //std::cout<<SumWeight(tree, "deltaPhiMin>0.3");
    if (useMC) {
        do_fit(250.0, 450.0, cond_all_mc, cond_qcd_mc, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(450.0, 575.0, cond_all_mc, cond_qcd_mc, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(575.0, 1000.0, cond_all_mc, cond_qcd_mc, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(1000.0, 1500.0, cond_all_mc, cond_qcd_mc, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(1500.0, 13000.0, cond_all_mc, cond_qcd_mc, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
    }
    
    else if(!useMC) {
        do_fit(250.0, 450.0, cond_all_data, cond_qcd_data, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(450.0, 575.0, cond_all_data, cond_qcd_data, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(575.0, 1000.0, cond_all_data, cond_qcd_data, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(1000.0, 1500.0, cond_all_data, cond_qcd_data, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        do_fit(1500.0, 13000.0, cond_all_data, cond_qcd_data, tree, useMC); std::cout<<""<<std::endl<<""<<std::endl<<""<<std::endl;
        
    }

  data->Close();
  return 0;
}



TH1D* set_ratio(double ht_min, double ht_max, std::string  condition, double threshold, TTree* tree, std::string qcd_or_all, bool useMC){

    //cuts:
    std::ostringstream oss_big;
    oss_big <<"("<< condition.c_str() << "&& ht<" << ht_max << " && ht>" << ht_min <<  " && deltaPhiMin>" << threshold<<")";
    std::string cond_big = oss_big.str();
    
    std::ostringstream oss_small;
    oss_small<<"("<< condition.c_str() << "&& ht<" << ht_max << " && ht>" << ht_min << " && deltaPhiMin<" <<threshold<<")";
    std::string cond_small = oss_small.str();
    
    //name of output histogram:
    std::ostringstream histo_name;
    histo_name <<"h_ratio_"<< qcd_or_all;
    std::string name = histo_name.str();
    
   
    double bins[] = {50.0, 55.0,  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 125.0, 200.0, 300.0, 500.0};
    int  binnum = sizeof(bins)/sizeof(double) - 1;
    
    TH1D* h_small = new TH1D ("h_small", "h_small", binnum, bins);
    TH1D* h_big = new TH1D ("h_big", "h_big", binnum, bins);
    
    //for data subtract nonQCD using montecarlo simulation:
    //if (!useMC) subtract_nonQCD(h_big, h_small, ht_min, ht_max, threshold);

   TH1D* h_ratio = new TH1D (name.c_str(), "h_ratio", binnum, bins);
    
    TCut selection_big = cond_big.c_str();
    TCut selection_small = cond_small.c_str();
    
    TCut weight = "weight";  //for data the variable weight is simply set to 1, so this works formally for both data and mc
    //TCut lumi = "lumi";
    TCut selection_big_and_weight = selection_big*weight;
    TCut selection_small_and_weight = selection_small*weight;
    
    tree->Draw("mt2>>h_big", selection_big_and_weight, "goff");
    tree->Draw("mt2>>h_small", selection_small_and_weight, "goff");
    
    if (!useMC){
        std::cout<<"Data: number of events before subtraction"<<std::endl;
        for (int i = 0; i<=h_big->GetXaxis()->GetNbins(); i++){
            std::cout<<"Bin number: "<<i<<" big dPhi: "<<h_big->GetBinContent(i)<<"  small DPhi: "<<h_small->GetBinContent(i)<<std::endl;
        }
    }
    
    //if data, subtract non qcd using montecarlo simulation
    if (!useMC && qcd_or_all == "qcd") {
        
        //per tagli usati vedi definizione di mc_rest nel codice di Masciovecchio
        std::ostringstream oss_big_mc;
        oss_big_mc <<"( id >= 300" << "&& ht<" << ht_max << " && ht>" << ht_min <<  " && deltaPhiMin>" << threshold<<" && nJets>1)";
        std::string cond_big_mc = oss_big_mc.str();
        
        std::ostringstream oss_small_mc;
        oss_big_mc <<"( id >= 300" << "&& ht<" << ht_max << " && ht>" << ht_min <<  " && deltaPhiMin<" << threshold<<" && nJets>1)";
        std::string cond_small_mc = oss_small_mc.str();
        
        std::cout<<"Subtracting non qcd from mc...";
        
        std::string mcpath = "EventYields_dataETH_SnTMC_35p9ifb/qcdControlRegion/mc.root";
        TFile* mc = new TFile(mcpath.c_str(), "READ");
        TTree* mc_tree = (TTree*)mc->Get("qcdCRtree/HT250toInf_j1toInf_b0toInf/tree_qcdCRtree_HT250toInf_j1toInf_b0toInf");
        
        TCut selection_big_mc = cond_big_mc.c_str();
        TCut selection_small_mc = cond_small_mc.c_str();
        
        TCut selection_big_and_weight_mc = selection_big_mc*weight;  //weight has been defined above
        TCut selection_small_and_weight_mc = selection_small_mc*weight;
        
        TH1D* h_small_mc = new TH1D ("h_small_mc", "h_small_mc", binnum, bins);
        TH1D* h_big_mc = new TH1D ("h_big_mc", "h_big_mc", binnum, bins);
        mc_tree->Draw("mt2>>h_big_mc", selection_big_and_weight_mc, "goff");
        mc_tree->Draw("mt2>>h_small_mc", selection_small_and_weight_mc, "goff");
        //problema: selezione diversa per gli id di montecarlo e data!! non usare selection..., o modificarla
        bool onlyUseUpToRunG = true;
        
        double lumiRatioGtoH = 27.261/35.876;
        double prescales[3] = {7900., 440.6, 110.2}; // up to run G 27.7 fb-1
        double lumiScale = 35.876;
        double sfFromSNT[5] = {1.88375, 1.38677, 1.27216, 1.16806, 1.02239};
        
        //lumiScale = useMC ? scaleMC : cfg.lumi(); SCALEMC???
        //divide mc statistics by prescales depending on trigger used in the given ht region
        if     ( ht_min < 300. ) lumiScale *= sfFromSNT[0]/prescales[0];
        else if( ht_min< 500. ) lumiScale *= sfFromSNT[1]/prescales[1];
        else if( ht_min< 600. ) lumiScale *= sfFromSNT[2]/prescales[2];
        else if( ht_min< 1100.) lumiScale *= sfFromSNT[3];
        else                          lumiScale *= sfFromSNT[4];
        if (onlyUseUpToRunG)
            lumiScale *= lumiRatioGtoH;
        
        lumiScale = - lumiScale; //add minus sign for subtraction
        
        //do the subtraction:
        h_small->Add(h_small_mc, lumiScale);
        h_big->Add(h_big_mc, lumiScale);
        
        std::cout<<" ...done"<<std::endl;
        
        std::cout<<"Data: number of qcd events estimated after subtraction"<<std::endl;
        for (int i = 0; i<=h_big->GetXaxis()->GetNbins(); i++){
            std::cout<<"Bin number: "<<i<<" big dPhi: "<<h_big->GetBinContent(i)<<"  small DPhi: "<<h_small->GetBinContent(i)<<std::endl;
        }
        
        delete h_small_mc;
        delete h_big_mc;
        
       
    }
    
    
    
    h_ratio->Divide(h_big, h_small);
    
    //std::cout<<"Ratio was set";
    delete h_small;
    delete h_big;
    return h_ratio;

}



void do_fit(double ht_min, double ht_max, std::string  cond_all, std::string cond_qcd, TTree* tree, bool useMC){
  //ht_min and ht_max are the boundaries of the ht interval to be considered
  //n_bins is the number of bins in mt2
  //tree is the tree containing the data
  double delta_Phi_threshold = 0.3;
  //limits of fitting region:
  double mt2_min = 60.0; double mt2_max = 100.0;  //fit boundaries
  if(ht_min>=999.0) {mt2_min = 70.0;}  //stay on the safe side for high ht regions
  
  //max and min mt2 for whole plot (not only region for fitting)
  double mt2_min_global = 50.0; double mt2_max_global = 1000.0;
    
    //new names for every energy region (to avoid segfault - still getting segfault, though:()
  std::ostringstream reg_nm;
  reg_nm << ht_min << "<HT<" << ht_max;
  std::string region_name = reg_nm.str();
  
  std::ostringstream pow_bg_nm;
  pow_bg_nm << region_name.c_str()<<"_pow_bg";
  std::string pow_bg_name = pow_bg_nm.str();
   
  std::ostringstream pow_bg_nm_left;
  pow_bg_nm_left << region_name.c_str()<<"_pow_bg_left";
  std::string pow_bg_left_name = pow_bg_nm_left.str();
    
  std::ostringstream pow_bg_nm_right;
  pow_bg_nm_right << region_name.c_str()<<"_pow_bg_right";
  std::string pow_bg_right_name = pow_bg_nm_right.str();
    
  std::ostringstream cfit_nm;
  cfit_nm << region_name.c_str()<<"_cfit";
  std::string cfit_name = cfit_nm.str();
    
  std::ostringstream pow_bg_fit_nm;
  pow_bg_fit_nm << region_name.c_str()<<"_pow_bg_fit";
  std::string pow_bg_fit_name = pow_bg_fit_nm.str();
    
  std::ostringstream band_nm;
  band_nm << region_name.c_str()<<"_band";
  std::string band_name = band_nm.str();
    
    
  std::cout<<ht_min<<" GeV < H_T < "<<ht_max<<std::endl;

  //start the plotting; create canvas
  TCanvas* cfit = new TCanvas(cfit_name.c_str(),"Fit with power law");
  std::cout<<"canvas was created"<<std::endl;
    
  TH1D* histo_qcd = set_ratio(ht_min, ht_max, cond_qcd, delta_Phi_threshold , tree, "qcd", useMC);
  TH1D* histo_all = set_ratio(ht_min, ht_max, cond_all, delta_Phi_threshold , tree, "all", useMC);

  //print ratios
  std::cout<<"Ratios"<<std::endl;
  for (int i = 0; i<=histo_qcd->GetXaxis()->GetNbins(); i++){
      std::cout<<"Bin number: "<<i<<" ratio: "<<histo_qcd->GetBinContent(i)<<"  weight: "<<histo_qcd->GetBinError(i)<<std::endl;
    }
    
  
  //pow_bg for plotting line, band for plotting (68%) confidence interval
  TF1* pow_bg = new TF1(pow_bg_name.c_str(), "[0]*TMath::Power(x,[1])", mt2_min, mt2_max);
  TH1D* band = new TH1D(band_name.c_str(), "Fitted gaussian with .68 conf.band", 500, mt2_min_global, mt2_max_global);

  histo_qcd->Fit(pow_bg_name.c_str(), "R 0"); //R: fit in the range of the function. don't plot now the fitted function
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(band, 0.68);
    
  TF1* fitResult = histo_qcd->GetFunction(pow_bg_name.c_str());

  double param_a = fitResult->GetParameter(0);
  double param_b = fitResult->GetParameter(1);
   
  //compute fit one bin to the left and one bin to the right to look at the fitted paramers' variation
  double mt2_min_right = mt2_min + 5.0;
  double mt2_max_right = mt2_max + 25.0;
  double mt2_min_left = mt2_min - 5.0;
  double mt2_max_left  = mt2_max;
    
  TF1* pow_bg_right = new TF1(pow_bg_right_name.c_str(), "[0]*TMath::Power(x,[1])", mt2_min_right, mt2_max_right);
  
  std::cout<<"Right variation:"<<std::endl;
  histo_qcd->Fit(pow_bg_right_name.c_str(), "R 0");
  TF1* fitResult_right = histo_qcd->GetFunction(pow_bg_right_name.c_str());
  
  TF1* pow_bg_left = new TF1(pow_bg_left_name.c_str(), "[0]*TMath::Power(x,[1])", mt2_min_left, mt2_max_left);
  double param_b_right = fitResult_right->GetParameter(1);
  double var_right = param_b_right/param_b;
    
  std::cout<<"Left variation:"<<std::endl;
  histo_qcd->Fit(pow_bg_left_name.c_str(), "R 0");
  TF1* fitResult_left = histo_qcd->GetFunction(pow_bg_left_name.c_str());
  
  double param_b_left = fitResult_left->GetParameter(1);
  double var_left = param_b_left/param_b;
  //compute now maximal variation:
  double var_max = TMath::Max( fabs(var_right-1.0), fabs(var_left-1.0) );
    
  std::cout<<"Second parameter (b): "<<param_b_right<<" (right) "<<param_b_left<<" (left) "<<std::endl;
  std::cout << "Slope variation: " << var_right << " (right) / " << var_left << " (left) " << std::endl<< "   maximal variation: " << var_max << std::endl;
    
  //get logarithmic graph:
  gPad->SetLogx();
  gPad->SetLogy();
  //gPad->SetTitle( "CMS simulation, #sqrt{s} = 13 TeV");
    
  //set graph title, labels etc.
  band->GetXaxis()->SetTitle("M_{T2} [GeV]");
  band->GetXaxis()->SetNoExponent();
  band->GetXaxis()->SetMoreLogLabels();
  band->GetYaxis()->SetTitle("r_{#phi}");
  if (useMC) band->SetTitle("CMS simulation, #sqrt{s} = 13 TeV");
  else  band->SetTitle("CMS preliminary"); //CAPIRE CROSS SECTION VICINO A CME
  band->SetFillColor(29);
  band->SetFillStyle(3001);
  band->SetMinimum(0.01);
  band->SetMaximum(10000.0);
  band->SetStats(0);
  band->SetLineColor(4);
  band->Draw("C E3");
    

    //now define fitted function over whole range (i.e. not only fitting region)
  TF1* fitted_bg = new TF1(pow_bg_fit_name.c_str(), "[0]*TMath::Power(x,[1])", mt2_min_global, mt2_max_global);
  fitted_bg->SetParameter(0, param_a);
  fitted_bg->SetParameter(1, param_b);
  fitted_bg->SetLineColor(4);
  fitted_bg->SetMinimum(0.01);
  fitted_bg->SetLineWidth(1);
  fitted_bg->SetMaximum(10000.0);
  fitted_bg->Draw("L SAME");
  
    
  std::cout<<"fit was plotted"<<std::endl;
 
  //histo->SetMarkerStyle(kFullSquare); per qualche motivo non funziona, ottengo una linea continua
  //histo_qcd->SetStats(0); //don't show information box on graph
  histo_qcd->SetMarkerStyle(24);
  //histo_qcd->LabelsOption("h","Y"); //IN QUALCHE MODO IL LABEL E' SEMPRE STORTO
  histo_qcd->SetStats(0);
  histo_qcd->SetLineColor(kBlack);
  histo_qcd->Draw("SAME P");
    
 histo_all->LabelsOption("h","Y"); //IN QUALCHE MODO IL LABEL E' SEMPRE STORTO
 //histo->SetMarkerStyle(kFullSquare); per qualche motivo non funziona, ottengo una linea continua
 histo_all->SetStats(0); //don't show information box on graph
 histo_all->SetMarkerStyle(20);
 histo_all->SetLineColor(kBlack);
 histo_all->Draw("SAME P");

  std::cout<<"histos were plotted"<<std::endl;


  //add legend:
  TLegend* legend = new TLegend(0.68,0.72,0.85,0.86);
  legend->AddEntry(band_name.c_str(), "Fit", "FL");
  
  if (useMC){
        legend->AddEntry("h_ratio_qcd","Simulation (multijet only)","PL");
        legend->AddEntry("h_ratio_all","Simulation (all)","PL");
      }
    
  else{
      legend->AddEntry("h_ratio_qcd","Data","PL");
      legend->AddEntry("h_ratio_all","Data after subtraction","PL");
      }
    
  legend->Draw();

  //add text with descrption of considered region
  std::ostringstream descr;
  descr << ht_min << " GeV < H_{T} < " << ht_max << " GeV";
  std::string description = descr.str(); 
  TPaveText* pt = new TPaveText(0.4,0.8,0.65,0.85, "NDC");
  pt->AddText(description.c_str());
  pt->Draw("SAME");
  
  //add lines at 60 and 100 GeV
  TLine* line1 = new TLine(mt2_min,1e-2, mt2_min, 1e4);
  TLine* line2 = new TLine(mt2_max,1e-2, mt2_max, 1e4);
  line1->SetLineStyle(2); line2->SetLineStyle(2);
  line1->Draw("SAME");
  line2->Draw("SAME");
    
  //add chi2 box
  std::ostringstream descr_chi2;
  descr_chi2 <<"#chi^{2}/ndf = "<<fitResult->GetChisquare()<<"/"<< fitResult->GetNDF()<<": "<< fitResult->GetProb()*100.0<<"%";
  std::string description_chi2 = descr_chi2.str();
  TPaveText* chi2;
  chi2 = new TPaveText(0.22,0.23, 0.5,0.18,"brNDC");
  chi2->SetFillColor(10);   chi2->SetBorderSize(0);
  chi2->AddText(description_chi2.c_str());
  chi2->Draw("SAME");
 
  std::cout<<"decos were plotted"<<std::endl;
    
  //save plot with file name indicating ht area
  std::ostringstream op_nm;
  op_nm <<"plotFilippo/"<< ht_min << "<HT<" << ht_max;
  if (useMC) {op_nm<<"_MC";}
  else {op_nm<<"_data";}
  op_nm<<".pdf";
  std::string output_name = op_nm.str();
  cfit->SaveAs(output_name.c_str());
    
    
   //cercare su lui applichi altri tagli particolari, solo nJets>1 ma non cambia nulla
    //FARE SELEZIONI PER DATA. VEDI CODICE DI MASCIOVECCHIO COPIATO IN FONDO CON LE AREE DI PHASE SPACE SELEZIONATE. CAPIRE COSA SIA ID NEL CASO DEI DATI
    
  delete chi2;
  delete pt;
  delete legend;
  delete line1;
  delete line2;
  delete pow_bg;
  delete band;
  delete fitted_bg;
  delete pow_bg_right;
  delete pow_bg_left;
  delete histo_qcd;
  delete histo_all;
  delete cfit;
}




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
 void MT2EstimateQCD::randomizePoisson( float scale, int seed ){
 
 MT2Estimate::randomizePoisson( scale, seed );
 
 TRandom3 rand(seed);
 
 for( int ibin=1; ibin<lDphi->GetXaxis()->GetNbins()+1; ++ibin ) {
 
 int poisson_data = rand.Poisson(scale * lDphi->GetBinContent(ibin));
 lDphi->SetBinContent(ibin, poisson_data);
 //lDphi->SetBinError( ibin, 0. );
 lDphi->SetBinError( ibin, TMath::Sqrt(poisson_data) ); // here i want an approximation of the Poisson error
 
 }
 
 
 */

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


/*
 
 if ( !useMC ) {//= DATA!!!
 // with PFHT900 trigger only (HT only triggered)
 data_4rHat       = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_4rHat"      , regionsSet_rHat  , qcdTree_data, "((id&1)==1 && ht>1000. &&  mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3"+runRange+")" ); // invert deltaPhi; ht>1000 unprescaled triggers
 data_4fJets      = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_4fJets"     , regionsSet_fJets , qcdTree_data, "((id&1)==1 &&              mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3"+runRange+")" ); // invert deltaPhi; HT-only triggers, ps'ed for HT<1000, empty for HT<450
 // without prescaled guys //update to new id stuff ... // it's wrong at the moment
 data_noPS_4fJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_noPS_4fJets", regionsSet_fJets , qcdTree_data, "(((id==1&&ht>100)||(id==2&&ht>450&&ht<1000)||(id==3&&ht<450))&& mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3"+runRange+")" ); // invert deltaPhi; unprescaled triggers, HT-only for ht>100, HTMHT for450< ht<1000, MET for ht<450 (below ht<1000 we live in turnon, better not to subtract an unknown amount of non-QCD [smaller than 18% for VLHT])
 //NON-QCD MC to be subtracted
 rest_4rHat       = MT2EstimateTree::makeAnalysisFromInclusiveTree( "rest_4rHat"      , regionsSet_rHat  , qcdTree_mc  , "(id>=300 && ht>1000. &&  mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3)" ); // invert deltaPhi
 rest_4fJets      = MT2EstimateTree::makeAnalysisFromInclusiveTree( "rest_4fJets"     , regionsSet_fJets , qcdTree_mc  , "(id>=300 &&              mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3)" ); // invert deltaPhi
 }*/
