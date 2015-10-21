#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"


void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawQCDFromDeltaPhi [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  MT2DrawTools::setStyle();


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool useMC = true;

  if( argc>2 ) {

    std::string mc_or_data = std::string(argv[2]); 
    if( mc_or_data=="mc" ) mc_or_data="MC";
    if( mc_or_data=="MC" ) useMC = true;
    else useMC=false;

  } 


  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion/";
  std::string outputdir = qcdCRdir;
  std::string fitsDir = outputdir;
  if( useMC ) fitsDir = fitsDir + "/fitsMC";
  else        fitsDir = fitsDir + "/fitsData";
  system( Form("mkdir -p %s", fitsDir.c_str() ));



  MT2Analysis<MT2EstimateTree>* qcdTree_mc   = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root",   "qcdCRtree" );
  //MT2Analysis<MT2EstimateTree>* qcdTree_data = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/data.root", "qcdCRtree" );
  MT2Analysis<MT2EstimateTree>* mcTruth = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth", cfg.regionsSet(), qcdTree_mc, "id>=153 && id<200 && mt2>200. && deltaPhiMin>0.3" ); // signal region for mcTruth

  std::string mcFile   = qcdCRdir + "qcdEstimateMC.root";
  std::string dataFile = qcdCRdir + "qcdEstimateData.root";

  MT2Analysis<MT2Estimate>* estimateMC     = MT2Analysis<MT2Estimate>::readFromFile( mcFile  , "qcdEstimate" );
  MT2Analysis<MT2Estimate>* estimateData   = MT2Analysis<MT2Estimate>::readFromFile( dataFile, "qcdEstimate" );

  mcTruth->setColor(kQCD);
  estimateMC->setColor(kBlack);
  estimateData->setColor(kBlack);

  std::string plotsDirMC = qcdCRdir + "/plotsMC";
  drawClosure( plotsDirMC, estimateMC, (MT2Analysis<MT2Estimate>*)mcTruth );

  std::string plotsDirData = qcdCRdir + "/plotsData";
  drawClosure( plotsDirData, estimateData, (MT2Analysis<MT2Estimate>*)mcTruth );

  return 0;

} 


//  MT2Analysis<MT2Estimate>* fjets_mc   = MT2Analysis<MT2Estimate>::readFromFile(mcFile  , "f_jets");
//  MT2Analysis<MT2Estimate>* fjets_data = MT2Analysis<MT2Estimate>::readFromFile(dataFile, "f_jets");
//
//
//  std::set<MT2Region> regions = fjets_data->getRegions();
//
//  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
//
//    MT2Estimate* thisFJets_mc   = fjets_mc  ->get(iR);
//    MT2Estimate* thisFJets_data = fjets_data->get(iR);
//
//    TH1D* h1_mc   = thisFJets_mc  ->yield;
//    TH1D* h1_data = thisFJets_data->yield;
//
//
//    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
//    c1->cd();
//
//    float yMax_data = h1_data->GetMaximum()/h1_data->Integral();
//    float yMax = 1.2*yMax_data;
//    TH1D* h2_axes = new TH2D( "axes", "", 10, 1.5, 8.5, 10, 0., yMax );
//    h2_axes->SetXTitle( "Number of Jets" );
//    h2_axes->SetYTitle( "Fraction" );
//    h2_axes->Draw();
//
//    h1_mc->SetLineColor(kRed);
//    h1_mc->SetLineWidth(2);
//
//    h1_data->SetMarkerStyle(20);
//    h1_data->SetMarkerSize(1.6);
//
//    h1_mc  ->Draw("histo norm same");
//    h1_data->Draw("p same norm");
//
//
//    TPaveText* labelTop = MT2DrawTools::getLabelTop( cfg.get_lumi() );
//    labelTop->Draw("same");
//
//    TLegend* legend = new TLegend( 0.55, 0.7, 0.9, 0.9 );
//    legend->SetFillColor(0);
//    legend->SetTextSize(0.038);
//    legend->AddEntry( h1_data, "Data", "P");
//    legend->AddEntry( h1_mc, "MC", "L");
//    legend->Draw("same");
//
//    gPad->RedrawAxis();
//
//    c1->SaveAs( Form("%s/fJets_%s.eps", outputdir.c_str(), iR->getName().c_str()) );
//    c1->SaveAs( Form("%s/fJets_%s.pdf", outputdir.c_str(), iR->getName().c_str()) );
//
//    delete c1;
//    delete h2_axes;
//
//}



void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth ) {


  system(Form("mkdir -p %s", outputdir.c_str()));

  
  std::set<MT2Region> MT2Regions = estimate->getRegions();
  
  TH1D* h_estimate_tot = new TH1D("h_estimate_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_estimate_tot->Sumw2();
  h_estimate_tot->GetYaxis()->SetTitle("Events");
  h_estimate_tot->SetMarkerStyle(20);
  h_estimate_tot->SetMarkerSize(1.6);
  h_estimate_tot->SetLineColor( estimate->getColor() );
  h_estimate_tot->SetMarkerColor( estimate->getColor() );
  
  TH1D* h_mcTruth_tot = new TH1D("h_mcTruth_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_mcTruth_tot->Sumw2();
  h_mcTruth_tot->GetYaxis()->SetTitle("Events");
  h_mcTruth_tot->SetFillColor(0);
  h_mcTruth_tot->SetLineColor( mcTruth->getColor() );
  h_mcTruth_tot->SetMarkerColor( mcTruth->getColor() );
  h_mcTruth_tot->SetMarkerStyle(20);
  h_mcTruth_tot->SetMarkerSize(1.6);
  
  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Data Driven - MC)/#sigma");
  hPull->GetYaxis()->SetTitle("Events");
  
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_estimate = estimate->get(*iMT2)->yield;
      h_estimate->SetMarkerStyle(20);
      h_estimate->SetMarkerSize(1.6);
      h_estimate->SetLineColor( estimate->getColor() );
      h_estimate->SetMarkerColor( estimate->getColor() );


      int nBins = h_estimate->GetXaxis()->GetNbins();
      double err_estimate;
      double int_estimate = h_estimate->IntegralAndError(1, nBins+1, err_estimate);
 
      h_estimate_tot->SetBinContent(iRegion, int_estimate);
      h_estimate_tot->SetBinError(iRegion, err_estimate);
      

      h_estimate_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      TH1D* h_mcTruth = mcTruth->get(*iMT2)->yield;

      h_mcTruth->SetLineColor( mcTruth->getColor() );
      h_mcTruth->SetMarkerColor( mcTruth->getColor() );
      h_mcTruth->SetMarkerStyle(20);
      h_mcTruth->SetMarkerSize(1.6);
      
   
      double err_int;
      double int_mcTruth = h_mcTruth->IntegralAndError(1, nBins+1, err_int);
      h_mcTruth_tot->SetBinContent(iRegion, h_mcTruth->Integral());
      h_mcTruth_tot->SetBinError(iRegion, err_int);
      h_mcTruth_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      if(int_mcTruth>0)
        hPull->Fill((int_estimate-int_mcTruth)/sqrt(err_estimate*err_estimate+err_int*err_int));


      float xMin = h_estimate->GetXaxis()->GetXmin();
      float xMax = h_estimate->GetXaxis()->GetXmax();
      float yMax_1 = h_estimate->GetMaximum()*1.5;
      float yMax_2 = 1.2*(h_estimate->GetMaximum() + h_estimate->GetBinError(h_estimate->GetMaximumBin()));
      float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
      float yMax_3 = h_mcTruth->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_mcTruth->GetMaximum() + h_mcTruth->GetBinError(h_mcTruth->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Events");

      h2_axes->Draw();
 
      //std::vector<std::string> niceNames = iMT2->getNiceNames();
      for( unsigned i=0; i<niceNames.size(); ++i ) {

        float yMax = 0.9-(float)i*0.05;
        float yMin = yMax - 0.05;
        TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
        regionText->SetTextSize(0.035);
        regionText->SetTextFont(42);
        regionText->SetFillColor(0);
        regionText->SetTextAlign(11);
        regionText->AddText( niceNames[i].c_str() );
        regionText->Draw("same");

      }


      TLegend* legend = new TLegend( 0.6, 0.9-2.*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( h_estimate, "data-driven", "P" );
      legend->AddEntry( h_mcTruth, "MC QCD", "P" );

      legend->Draw("same");

      h_estimate->Draw("P same");
      h_mcTruth->Draw("P same");
      //      bgStack.Draw("histoE, same");

      TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();

      std::string thisName = Form("%s_ratio", h_estimate->GetName());
      TH1D* h_ratio = (TH1D*) h_estimate->Clone(thisName.c_str());
      h_ratio->Divide(h_mcTruth);
      h_ratio->SetStats(0);	    
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetLineColor(1);
      //      h_ratio->SetMarkerSize(0.02);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");
      
      

      h_ratio->SetLineWidth(2);
      //h_ratio->Draw("PE");
      
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);
      

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      h_ratio->Draw("pe,same");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/closure_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/closure_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;
      delete h2_axes_ratio;
      //delete h_ratio;
      //delete h_estimate;
      //delete h_mcTruth;
      
      ++iRegion;

  } // for MT2 regions

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = h_estimate_tot->GetMaximum();
  float yMax_2 = h_estimate_tot->GetMaximum() + h_estimate_tot->GetBinError(h_mcTruth_tot->GetMaximumBin());
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = h_mcTruth_tot->GetMaximum();
  float yMax_4 = h_mcTruth_tot->GetMaximum() + h_mcTruth_tot->GetBinError(h_mcTruth_tot->GetMaximumBin());
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  yMax*=20.;
  
  float yMin = 1e-1;
  for( int iBin=1; iBin<h_estimate_tot->GetXaxis()->GetNbins()+1; ++iBin ) {
    if( h_estimate_tot    ->GetBinContent(iBin)>0. && h_estimate_tot    ->GetBinContent(iBin)<yMin ) yMin = h_estimate_tot    ->GetBinContent(iBin);
    if( h_mcTruth_tot->GetBinContent(iBin)>0. && h_mcTruth_tot->GetBinContent(iBin)<yMin ) yMin = h_mcTruth_tot->GetBinContent(iBin);
  }
  yMin /= 3.;
  
  h_mcTruth_tot->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  h_mcTruth_tot->GetYaxis()->SetRangeUser(yMin, yMax);
  h_mcTruth_tot->GetXaxis()->LabelsOption("v");
  h_mcTruth_tot->Draw("PE");


  h_estimate_tot->SetMarkerStyle(20);
  h_estimate_tot->SetMarkerSize(1.6);
  h_estimate_tot->SetLineColor( estimate->getColor() );
  h_estimate_tot->SetMarkerColor( estimate->getColor() );

  h_estimate_tot->Draw("pe,same");

  TLegend* legend = new TLegend( 0.18, 0.7, 0.32, 0.82 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( h_estimate_tot, "data-driven", "PL" );
  legend->AddEntry( h_mcTruth_tot, "QCD MC", "PL" );

  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
  labelTop->Draw("same");
  
  TLine* lHT[3];
  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1] = new TLine(11*iHT, -3., 11*iHT, yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  int nHTRegions = 4;
  std::vector< std::string > htRegions;
  htRegions.push_back("low H_{T}");
  htRegions.push_back("medium H_{T}");
  htRegions.push_back("high H_{T}");
  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[nHTRegions];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.16+0.2*iHT, 0.9-0.06, 0.34+0.2*iHT, 0.9, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );
    
    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.038);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  }

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

  std::string thisName = Form("%s_ratio", h_estimate_tot->GetName());
  TH1D* h_Ratio = (TH1D*) h_estimate_tot->Clone(thisName.c_str());
  for( int iBin=1; iBin<h_Ratio->GetXaxis()->GetNbins()+1; ++iBin ) {
    float mc = h_estimate_tot->GetBinContent(iBin);
    float mc_err = h_estimate_tot->GetBinError(iBin);
    float est = h_mcTruth_tot->GetBinContent(iBin);
    float est_err = h_mcTruth_tot->GetBinError(iBin);
    float denom = sqrt( mc_err*mc_err + est_err*est_err );
    if( denom!=0. && mc>0. ) {
      h_Ratio->SetBinContent( iBin, (mc-est)/denom );
      h_Ratio->SetBinError( iBin, 1. );
    }
  }
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->SetLineWidth(2);
  

  TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(), -3., 3.);
  h2_axes_ratio->SetYTitle("Pull");

  TLine* LineCentral = new TLine(0, 0, MT2Regions.size(), 0);
  LineCentral->SetLineColor(1);


  h2_axes_ratio->Draw("");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");

  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }


  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/closure_allRegions_pull.pdf", outputdir.c_str()) );
  c2->SaveAs( Form("%s/closure_allRegions_pull.eps", outputdir.c_str()) );

  pad2->cd();
  pad2->Clear();

  delete h2_axes_ratio;
  h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(), 0., 2.);
  h2_axes_ratio->SetYTitle("Data / MC");
  h2_axes_ratio->Draw("");

  TLine* lineOne = new TLine(0, 1., MT2Regions.size(), 1.);
  lineOne->SetLineColor(1);
  lineOne->Draw("same");

  delete h_Ratio;
  h_Ratio = (TH1D*) h_estimate_tot->Clone(thisName.c_str());
  h_Ratio->Divide( h_mcTruth_tot );
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->SetLineWidth(2);

  h_Ratio->Draw("pe,same");

  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }


  c2->SaveAs( Form("%s/closure_allRegions_ratio.pdf", outputdir.c_str()) );
  c2->SaveAs( Form("%s/closure_allRegions_ratio.eps", outputdir.c_str()) );


  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  TF1* f1_gaus = new TF1("f1_pull", "gaus", -2., 2.);
  f1_gaus->SetLineColor(kRed);
  hPull->Fit( f1_gaus, "QRL" );
  TPaveText* fitPars = new TPaveText( 0.2, 0.7, 0.5, 0.9, "brNDC" );
  fitPars->SetTextSize(0.03);
  fitPars->SetTextAlign(11);
  fitPars->SetFillColor(0);
  fitPars->AddText("Gaussian Fit:");
  fitPars->AddText(Form("Mean : %.2f +/- %.2f", f1_gaus->GetParameter(1), f1_gaus->GetParError(1) ));
  fitPars->AddText(Form("Sigma: %.2f +/- %.2f", f1_gaus->GetParameter(2), f1_gaus->GetParError(2) ));
  fitPars->Draw("same");
  hPull->Draw("hist same");
  f1_gaus->Draw("l same");
  c3->SaveAs( Form("%s/closure_pull.pdf", outputdir.c_str()) );
  c3->SaveAs( Form("%s/closure_pull.eps", outputdir.c_str()) );


  delete c2;
  delete c3;
  delete h2_axes_ratio;

  delete h_estimate_tot;
  delete h_mcTruth_tot;
  delete hPull;
  
}


