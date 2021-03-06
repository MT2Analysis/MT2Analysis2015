#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2DrawTools.h"
#include "../interface/MT2Config.h"





void drawSinglePlot( const MT2Config& cfg, const std::string& name, const MT2Region& region, std::vector<TH1D*> histosData, TH1D* histoMC );
void setHistoTitle( MT2Analysis<MT2EstimateZinvGamma>* analysis, const std::string& title );


int main( int argc, char* argv[] ) {


  std::string regionsSet = "13TeV_inclusive";

  MT2DrawTools::setStyle();


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  MT2Analysis<MT2EstimateZinvGamma>* templatesFakeMC         = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( cfg.getGammaCRdir() + "gammaTemplatesMC_MC.root"  , "templatesFake");
  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptMC       = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( cfg.getGammaCRdir() + "gammaTemplatesMC_MC.root"  , "templatesPrompt");

  MT2Analysis<MT2EstimateZinvGamma>* templatesFakeRC_MC      = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( cfg.getGammaCRdir() + "gammaTemplatesRC_MC.root"  , "templatesFake");
  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRC_MC    = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( cfg.getGammaCRdir() + "gammaTemplatesRC_MC.root"  , "templatesPrompt");

  MT2Analysis<MT2EstimateZinvGamma>* templatesFakeRC_data    = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( cfg.getGammaCRdir() + "gammaTemplatesRC_data.root", "templatesFake");
  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRC_data  = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( cfg.getGammaCRdir() + "gammaTemplatesRC_data.root", "templatesPrompt");


  // setHistoTitle( templatesFake, "#sigma_{i#eta i#eta} Sidebands" );
  setHistoTitle( templatesFakeMC, "MC Fakes" );

  setHistoTitle( templatesFakeRC_data, "#sigma_{i#etai#eta} Sideband Data" );
  setHistoTitle( templatesFakeRC_MC, "#sigma_{i#etai#eta} Sideband MC" );
  //  setHistoTitle( templatesFakeFR_MC, "MC Fakes" );
  // setHistoTitle( templatesPromptRaw, "Data (all)" );
  // setHistoTitle( templatesPrompt, "Data (fake removal)" );
  setHistoTitle( templatesPromptRC_data, "Data (random cone)" );
  setHistoTitle( templatesPromptRC_MC, "MC (random cone)" );
  setHistoTitle( templatesPromptMC, "MC Prompts" );


  std::set<MT2Region> regions = templatesPromptRC_MC->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    std::vector<TH1D*> v_fakes;
    //  v_fakes.push_back( templatesFakeMC->get(*iR)->iso );
    v_fakes.push_back( templatesFakeRC_data->get(*iR)->iso );
    v_fakes.push_back( templatesFakeRC_MC->get(*iR)->iso );
    //  v_fakes.push_back( templatesFakeFR_MC->get(*iR)->iso );
    drawSinglePlot( cfg, "Fake", *iR, v_fakes, templatesFakeMC->get(*iR)->iso );

    
    std::vector<TH1D*> v_prompts;
    //   v_prompts.push_back( templatesPromptRaw->get(*iR)->iso );
    // v_prompts.push_back( templatesPrompt->get(*iR)->iso );
    v_prompts.push_back( templatesPromptRC_data->get(*iR)->iso );
    v_prompts.push_back( templatesPromptRC_MC->get(*iR)->iso );
    //  v_prompts.push_back( templatesPromptMC->get(*iR)->iso );
    drawSinglePlot( cfg, "Prompt", *iR, v_prompts, templatesPromptMC->get(*iR)->iso );

  } 


  return 0;

}



void drawSinglePlot( const MT2Config& cfg, const std::string& name, const MT2Region& region, std::vector<TH1D*> histosData, TH1D* histoMC ) {



  std::string outputdir = cfg.getGammaCRdir() + "/plotsTemplatesData/";
  system( Form("mkdir -p %s", outputdir.c_str()) );



  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
  c1_log->SetLogy();

  float yMinLegend = 0.9 - (histosData.size()+1)*0.065;
  float yMaxScale = (name=="Fake") ? 1.6 : 1.25;
  float xMax = histoMC->GetXaxis()->GetXmax();

  //if( name=="Fake" ) {

  //  histoMC->Scale( histosData[0]->Integral()/histoMC->Integral() );

  //} else { // rescale so that first bin has same content


  if( name=="Fake" ) { 

    // keep first one normalized to lumi, normalize all others to integral of first

    for( unsigned i=1; i<histosData.size(); ++i )
      histosData[i]->Scale( histosData[0]->Integral()/histosData[i]->Integral() );

    histoMC->Scale( histosData[0]->Integral()/histoMC->Integral() );

  } else {

    //// normalize first one to unity
    //histosData[0]->Scale(1./histosData[0]->Integral() );

    // normalize other ones so that first bin has same content:
    for( unsigned i=1; i<histosData.size(); ++i ) 
      histosData[i]->Scale( histosData[0]->GetBinContent(1)/histosData[i]->GetBinContent(1) );

    histoMC->Scale( histosData[0]->GetBinContent(1)/histoMC->GetBinContent(1) );

  }



  float yMax = histosData[0]->GetMaximum()*yMaxScale;


  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, 0., yMax );
  h2_axes->SetXTitle( "Photon Charged Isolation [GeV]" );
  if( name=="Fake" )
    h2_axes->SetYTitle( Form("Events / %.2f GeV", histoMC->GetBinWidth(1)) );
  else
    h2_axes->SetYTitle( "Normalized to First Bin" );
  c1->cd();
  h2_axes->Draw();

  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, 0., xMax, 10, yMax/10000., 3.*yMax );
  h2_axes_log->SetXTitle( "Photon Charged Isolation [GeV]" );
  if( name=="Fake" )
    h2_axes_log->SetYTitle( Form("Events / %.2f GeV", histoMC->GetBinWidth(1)) );
  else
    h2_axes_log->SetYTitle( "Normalized to First Bin" );
  c1_log->cd();
  h2_axes_log->Draw();


  histoMC->SetLineColor( 46 );
  histoMC->SetLineWidth( 2 );
  c1->cd();
  histoMC->Draw("L E same");
  c1_log->cd();
  histoMC->Draw("L E same");

  TLegend* legend = new TLegend( 0.5, yMinLegend, 0.9, 0.9 );
  //TLegend* legend = new TLegend( 0.46, yMinLegend, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);

  std::vector<int> colors;
  colors.push_back( kBlack );
  colors.push_back( kBlack );
  colors.push_back( kBlue );
  colors.push_back( 38 );
  colors.push_back( 46 );

  std::vector<int> markers;
  markers.push_back( 20 );
  markers.push_back( 24 );
  markers.push_back( 25 );

  float markerSize = (name=="Fake") ? 1.6 : 1.3;

  for( unsigned i=0; i<histosData.size(); ++i ) {
    histosData[i]->SetMarkerStyle( markers[i] );
    histosData[i]->SetMarkerColor( colors[i] );
    histosData[i]->SetMarkerSize( markerSize );
    histosData[i]->SetLineColor( colors[i] );
    c1->cd();
    histosData[i]->Draw("p same");
    c1_log->cd();
    histosData[i]->Draw("p same");
    legend->AddEntry( histosData[i], histosData[i]->GetTitle(), "P" );
  }

  legend->AddEntry( histoMC, histoMC->GetTitle(), "L" );


  TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());

  c1->cd();
  legend->Draw("same");
  MT2DrawTools::addLabels( c1, cfg.lumi(), "CMS Preliminary" );
  //  labelTop->Draw("same");
  gPad->RedrawAxis();


  c1_log->cd();
  legend->Draw("same");
  MT2DrawTools::addLabels( c1_log, cfg.lumi(), "CMS Preliminary" );
  //  labelTop->Draw("same");
  gPad->RedrawAxis();


  c1->SaveAs(Form("%s/templateClosure%s_%s.eps", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1->SaveAs(Form("%s/templateClosure%s_%s.png", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1->SaveAs(Form("%s/templateClosure%s_%s.pdf", outputdir.c_str(), name.c_str(), region.getName().c_str()));

  c1_log->SaveAs(Form("%s/templateClosure%s_%s_log.eps", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1_log->SaveAs(Form("%s/templateClosure%s_%s_log.png", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1_log->SaveAs(Form("%s/templateClosure%s_%s_log.pdf", outputdir.c_str(), name.c_str(), region.getName().c_str()));

  delete c1;
  delete h2_axes;
  delete c1_log;
  delete h2_axes_log;

}



void setHistoTitle( MT2Analysis<MT2EstimateZinvGamma>* analysis, const std::string& title ) {

  std::set<MT2Region> regions = analysis->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) 
    analysis->get(*iR)->iso->SetTitle(title.c_str());

}
