#include <cmath>
#include <iostream>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"


#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateZinvGamma.h"

using namespace RooFit;



struct Purity {

  float purity;
  float purityError;

};



void fitPurity( const std::string& outputdir, TH1D* h1_purityLoose, TH1D* h1_purityTight, RooRealVar* x, std::vector<RooDataSet*> data, TH1D* templPrompt, TH1D* templFake );
void fitSinglePurity( const std::string& outputdir, Purity& loose, Purity& tight, RooRealVar* x, RooDataSet* data, TH1D* h1_templPrompt, TH1D* h1_templFake );





int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v2_Zinv";

  if( argc==1 ) {
    std::cout << "-> You need to pass me the regions set name. Here are some suggestions: " << std::endl;
    std::cout << "  13TeV_CSA14" << std::endl;
    std::cout << "  13TeV_onlyHT" << std::endl;
    std::cout << "  13TeV_onlyJets" << std::endl;
    std::cout << "  13TeV_inclusive" << std::endl;
    exit(101);
  }


  std::string regionsSet = "13TeV_CSA14";
  if( argc>1 ) {
    std::string regionsSet_tmp(argv[1]); 
    regionsSet = regionsSet_tmp;
  }

  std::string mc_or_data = "MC";
  if( argc>2 ) {
    mc_or_data = std::string(argv[2]);
    if( mc_or_data=="data" ) mc_or_data="Data";
    if( mc_or_data=="mc" ) mc_or_data="MC";
  }


  TH1::AddDirectory(kFALSE);


  MT2Analysis<MT2EstimateZinvGamma>* gammaJet_data = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( "GammaControlRegion_" + samples + "_13TeV_CSA14/mc.root", "gammaCR" );

  MT2Analysis<MT2EstimateZinvGamma>* templates_prompt = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( "gammaTemplates" + mc_or_data + "_" + samples + "_" + regionsSet + ".root", "templatesPrompt" );
  MT2Analysis<MT2EstimateZinvGamma>* templates_fake   = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( "gammaTemplates" + mc_or_data + "_" + samples + "_" + regionsSet + ".root", "templatesFake" );


  std::string outputdir = "PurityFits" + mc_or_data + "_" + samples + "_" + regionsSet;
  system( Form( "mkdir -p %s", outputdir.c_str()) );

  std::set<MT2Region> regions = gammaJet_data->getRegions();

  MT2Analysis<MT2Estimate>* purityLoose = new MT2Analysis<MT2Estimate>( "purityLoose", "13TeV_CSA14" );
  MT2Analysis<MT2Estimate>* purityTight = new MT2Analysis<MT2Estimate>( "purity"     , "13TeV_CSA14" );


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    MT2EstimateZinvGamma* thisEstimate = gammaJet_data->get( *iR );

    MT2EstimateZinvGamma* templatePrompt = templates_prompt->get( *(templates_prompt->matchRegion( *iR )) );
    MT2EstimateZinvGamma* templateFake   = templates_fake  ->get( *(templates_fake  ->matchRegion( *iR )) );

    MT2Estimate* thisLoosePurity = purityLoose->get( *iR );
    std::string nameLoose = thisLoosePurity->yield->GetName();

    MT2Estimate* thisTightPurity = purityTight->get( *iR );
    std::string nameTight = thisTightPurity->yield->GetName();

    fitPurity( outputdir, thisLoosePurity->yield, thisTightPurity->yield, thisEstimate->x_, thisEstimate->iso_bins, templatePrompt->iso, templateFake->iso);

    thisLoosePurity->yield->SetName( nameLoose.c_str() );
    thisTightPurity->yield->SetName( nameTight.c_str() );

  }
    

  purityLoose->writeToFile( outputdir + "/purityFit_" + samples + "_" + regionsSet + ".root" );
  purityTight->addToFile( outputdir + "/purityFit_" + samples + "_" + regionsSet + ".root" );

  return 0;

}





void fitPurity( const std::string& outputdir, TH1D* h1_purityLoose, TH1D* h1_purityTight, RooRealVar* x, std::vector<RooDataSet*> data, TH1D* templPrompt, TH1D* templFake ) {


  for( unsigned i=0; i<data.size(); ++i ) {

    int ibin = i+1;

    Purity loose, tight;
    fitSinglePurity( outputdir, loose, tight, x, data[i], templPrompt, templFake );

    h1_purityLoose->SetBinContent( ibin, loose.purity );
    h1_purityLoose->SetBinError  ( ibin, loose.purityError );

    h1_purityTight->SetBinContent( ibin, tight.purity );
    h1_purityTight->SetBinError  ( ibin, tight.purityError );

  }


  return;

}





void fitSinglePurity( const std::string& outputdir, Purity& loose, Purity& tight, RooRealVar* x, RooDataSet* data, TH1D* h1_templPrompt, TH1D* h1_templFake ) {


  float dataIntegral = data->sumEntries();

  if( dataIntegral == 0. ) {
    loose.purity=-1;
    loose.purityError=0.;
    tight.purity=-1;
    tight.purityError=0.;
    return;
  }


  RooDataHist templPrompt("templPrompt", "", *x, h1_templPrompt);
  RooDataHist templFake  ("templFake"  , "", *x, h1_templFake  );

  RooHistPdf pdfPrompt("pdfPrompt", "", *x, templPrompt );
  RooHistPdf pdfFake  ("pdfFake"  , "", *x, templFake   );

  RooRealVar sigFrac("promptFrac","fraction of prompt",0.9,0.,1.) ;
  RooAddPdf  model("model","", RooArgList(pdfPrompt,pdfFake), sigFrac) ;

  int nBins = h1_templPrompt->GetNbinsX();
  float xMin = h1_templPrompt->GetXaxis()->GetXmin();
  float xMax = h1_templPrompt->GetXaxis()->GetXmax();

  float xMaxFit = 0.999*xMax;
  x->setRange( "fittingRange", 0., xMaxFit );
  model.fitTo(*data, SumW2Error(kTRUE), Range("fittingRange")); 

  loose.purity = sigFrac.getVal();
  loose.purityError = sigFrac.getError();

  float sigEff = h1_templPrompt->GetBinContent(1)/h1_templPrompt->Integral(1,nBins);
  float bgEff = h1_templFake->GetBinContent(1)/h1_templFake->Integral(1,nBins);
  float sigFirstBin = sigFrac.getVal()*sigEff;
  float bgFirstBin = (1.-sigFrac.getVal())*bgEff;
  tight.purity = sigFirstBin / (sigFirstBin+bgFirstBin);
  tight.purityError = sigFrac.getError();

  RooPlot* xframe = x->frame();
  data->plotOn(xframe, Binning(nBins, xMin, xMax));
  model.plotOn(xframe);

  // Overlay the background component of model with a dashed line
  model.plotOn(xframe,Components(pdfFake),LineStyle(kDashed)) ;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMaxFit, 10, 0., xframe->GetMaximum()*1.2 );
  h2_axes->Draw();
  xframe->GetYaxis()->SetTitleOffset(1.4); 
  xframe->Draw("same");

  c1->SaveAs(Form("%s/purityFit_%s.eps", outputdir.c_str(), data->GetName()));
  c1->SaveAs(Form("%s/purityFit_%s.png", outputdir.c_str(), data->GetName()));

  delete c1;
  delete xframe;

  return;

}

