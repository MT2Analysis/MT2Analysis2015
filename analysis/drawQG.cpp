#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"

#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"





struct QGHistos {

  TH1D* h1_qFrac;
  TH1D* h1_qgl0;
  TH1D* h1_qgl1;
  TH1D* h1_qgl2;
  TH1D* h1_qgl3;
  TH1D* h1_qglAve;
  TH1D* h1_qglProd;

};



void drawCompare( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* bg, std::vector<MT2Analysis<MT2EstimateTree>* > signals );
QGHistos getHistos( MT2Analysis<MT2EstimateTree>* analysis );
void drawSingleHisto( const std::string& outputdir, const std::string& saveName, TH1D* h1_bg, const std::string& bgName, TH1D* h1_sig0, const std::string& sigName0 );



int main( int argc, char* argv[] ) {


  if( argc==1 ) {
    std::cout << "USAGE: ./drawQG [inputdir]" << std::endl;
    exit(233);
  }


  std::string inputDir = std::string(argv[1]);

  MT2DrawTools::setStyle();

  std::string outputdir = "QGPlots_" + inputDir;
  system( Form("mkdir -p %s", outputdir.c_str()) );

  MT2Analysis<MT2EstimateTree>* analysis_ZJets           = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "ZJets");
  MT2Analysis<MT2EstimateTree>* analysis_T1bbbb_1000_900 = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T1bbbb_2J_mGl1000_mLSP900_post");
  MT2Analysis<MT2EstimateTree>* analysis_T1bbbb_1500_100 = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T1bbbb_2J_mGl1500_mLSP100_post");
  MT2Analysis<MT2EstimateTree>* analysis_T1qqqq_1000_800 = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T1qqqq_2J_mGl1000_mLSP800_post");
  MT2Analysis<MT2EstimateTree>* analysis_T1qqqq_1400_100 = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T1qqqq_2J_mGl1400_mLSP100_post");
  MT2Analysis<MT2EstimateTree>* analysis_T1tttt_1200_800 = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T1tttt_2J_mGl1200_mLSP800_post");
  MT2Analysis<MT2EstimateTree>* analysis_T1tttt_1500_100 = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T1tttt_2J_mGl1500_mLSP100_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2bb_600_580    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2bb_2J_mStop600_mLSP580_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2bb_900_100    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2bb_2J_mStop900_mLSP100_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2qq_1200_100   = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2qq_2J_mStop1200_mLSP100_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2qq_600_550    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2qq_2J_mStop600_mLSP550_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2tt_425_325    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2tt_2J_mStop425_mLSP325_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2tt_500_325    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2tt_2J_mStop500_mLSP325_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2tt_650_325    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2tt_2J_mStop650_mLSP325_post");
  MT2Analysis<MT2EstimateTree>* analysis_T2tt_850_100    = MT2Analysis<MT2EstimateTree>::readFromFile(inputDir + "/analyses.root", "SMS_T2tt_2J_mStop850_mLSP100_post");

  analysis_ZJets           ->setFullName("ZJets");
  analysis_T1bbbb_1000_900 ->setFullName("T1bbbb (1000, 900)");
  analysis_T1bbbb_1500_100 ->setFullName("T1bbbb (1500, 100)");
  analysis_T1qqqq_1000_800 ->setFullName("T1qqqq (1000, 800)");
  analysis_T1qqqq_1400_100 ->setFullName("T1qqqq (1400, 100)");
  analysis_T1tttt_1200_800 ->setFullName("T1tttt (1200, 800)");
  analysis_T1tttt_1500_100 ->setFullName("T1tttt (1500, 100)");
  analysis_T2bb_600_580    ->setFullName("T2bb (600, 580)");
  analysis_T2bb_900_100    ->setFullName("T2bb (900, 100)");
  analysis_T2qq_1200_100   ->setFullName("T2qq (1200, 100)");
  analysis_T2qq_600_550    ->setFullName("T2qq (600, 550)");
  analysis_T2tt_425_325    ->setFullName("T2tt (425, 325)");
  analysis_T2tt_500_325    ->setFullName("T2tt (500, 325)");
  analysis_T2tt_650_325    ->setFullName("T2tt (650, 325)");
  analysis_T2tt_850_100    ->setFullName("T2tt (850, 100)");

  std::vector<MT2Analysis<MT2EstimateTree>*> signals;
  signals.push_back( analysis_T1qqqq_1400_100 );

  drawCompare( outputdir, analysis_ZJets, signals );

  return 0;

}


void drawCompare( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* bg, std::vector<MT2Analysis<MT2EstimateTree>* > signals ) {


  QGHistos vh_bg = getHistos( bg );

  std::vector< QGHistos > vh_sig;
  for( unsigned i=0; i<signals.size(); ++i ) 
    vh_sig.push_back( getHistos( signals[i] ) );

  
  drawSingleHisto( outputdir, "qFrac", vh_bg.h1_qFrac, bg->getFullName(), vh_sig[0].h1_qFrac, signals[0]->getFullName() );  
  drawSingleHisto( outputdir, "qgl0", vh_bg.h1_qgl0, bg->getFullName(), vh_sig[0].h1_qgl0, signals[0]->getFullName() );  
  drawSingleHisto( outputdir, "qgl1", vh_bg.h1_qgl1, bg->getFullName(), vh_sig[0].h1_qgl1, signals[0]->getFullName() );  
  drawSingleHisto( outputdir, "qgl2", vh_bg.h1_qgl2, bg->getFullName(), vh_sig[0].h1_qgl2, signals[0]->getFullName() );  
  drawSingleHisto( outputdir, "qgl3", vh_bg.h1_qgl3, bg->getFullName(), vh_sig[0].h1_qgl3, signals[0]->getFullName() );  
  drawSingleHisto( outputdir, "qglAve", vh_bg.h1_qglAve, bg->getFullName(), vh_sig[0].h1_qglAve, signals[0]->getFullName() );  
  drawSingleHisto( outputdir, "qglProd", vh_bg.h1_qglProd, bg->getFullName(), vh_sig[0].h1_qglProd, signals[0]->getFullName() );  

}



QGHistos getHistos( MT2Analysis<MT2EstimateTree>* analysis ) {

  TTree* tree = analysis->get( *(analysis->getRegions().begin()) )->tree;


  float weight;
  tree->SetBranchAddress( "weight", &weight );
  int nJets;
  tree->SetBranchAddress( "nJets", &nJets );
  int nBJets;
  tree->SetBranchAddress( "nBJets", &nBJets );
  float qgl0;
  tree->SetBranchAddress( "qgl0", &qgl0 );
  float qgl1;
  tree->SetBranchAddress( "qgl1", &qgl1 );
  float qgl2;
  tree->SetBranchAddress( "qgl2", &qgl2 );
  float qgl3;
  tree->SetBranchAddress( "qgl3", &qgl3 );
  float partId0;
  tree->SetBranchAddress( "partId0", &partId0 );
  float partId1;
  tree->SetBranchAddress( "partId1", &partId1 );
  float partId2;
  tree->SetBranchAddress( "partId2", &partId2 );
  float partId3;
  tree->SetBranchAddress( "partId3", &partId3 );
  float qglProd;
  tree->SetBranchAddress( "qglProd", &qglProd );
  float qglAve;
  tree->SetBranchAddress( "qglAve", &qglAve );


  QGHistos histos;
  histos.h1_qFrac = new TH1D( Form("qFrac_%s", analysis->getName().c_str()), "", 4, -0.5, 3.5 );
  histos.h1_qgl0 = new TH1D( Form("qgl0_%s", analysis->getName().c_str()), "", 50, 0., 1.0001 );
  histos.h1_qgl1 = new TH1D( Form("qgl1_%s", analysis->getName().c_str()), "", 50, 0., 1.0001 );
  histos.h1_qgl2 = new TH1D( Form("qgl2_%s", analysis->getName().c_str()), "", 50, 0., 1.0001 );
  histos.h1_qgl3 = new TH1D( Form("qgl3_%s", analysis->getName().c_str()), "", 50, 0., 1.0001 );
  histos.h1_qglProd = new TH1D( Form("qglProd_%s", analysis->getName().c_str()), "", 50, 0., 1.0001 );
  histos.h1_qglAve = new TH1D( Form("qglAve_%s", analysis->getName().c_str()), "", 50, 0., 1.0001 );
  

  float quarks0=0.;
  float quarks1=0.;
  float quarks2=0.;
  float quarks3=0.;
  float all0=0.;
  float all1=0.;
  float all2=0.;
  float all3=0.;


  int nentries = tree->GetEntries();


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 100000 == 0) std::cout << " entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry( iEntry );

    if( nBJets>0 ) continue;

    if( nJets>0 ) {
      all0 += weight;
      if( fabs(partId0)<4.5 ) quarks0+=weight;
      histos.h1_qgl0->Fill( qgl0, weight );
    }

    if( nJets>1 ) {
      all1 += weight;
      if( fabs(partId1)<4.5 ) quarks1+=weight;
      histos.h1_qgl1->Fill( qgl1, weight );
    }

    if( nJets>2 ) {
      all2 += weight;
      if( fabs(partId2)<4.5 ) quarks2+=weight;
      histos.h1_qgl2->Fill( qgl2, weight );
    }

    if( nJets>3 ) {
      all3 += weight;
      if( fabs(partId3)<4.5 ) quarks3+=weight;
      histos.h1_qgl3->Fill( qgl3, weight );
    }

    histos.h1_qglAve->Fill( qglAve, weight );
    histos.h1_qglProd->Fill( qglProd, weight );

  } // for entries


  histos.h1_qFrac->SetBinContent( 1, quarks0/all0 );
  histos.h1_qFrac->SetBinContent( 2, quarks1/all1 );
  histos.h1_qFrac->SetBinContent( 3, quarks2/all2 );
  histos.h1_qFrac->SetBinContent( 4, quarks3/all3 );
  
  histos.h1_qFrac->GetXaxis()->SetBinLabel( 1, "Jet0" );
  histos.h1_qFrac->GetXaxis()->SetBinLabel( 2, "Jet1" );
  histos.h1_qFrac->GetXaxis()->SetBinLabel( 3, "Jet2" );
  histos.h1_qFrac->GetXaxis()->SetBinLabel( 4, "Jet3" );

  histos.h1_qFrac->SetYTitle( "Fraction of Light Quarks" );
  
  histos.h1_qgl0->SetXTitle( "First Jet QGL" );
  histos.h1_qgl1->SetXTitle( "Second Jet QGL" );
  histos.h1_qgl2->SetXTitle( "Third Jet QGL" );
  histos.h1_qgl3->SetXTitle( "Fourth Jet QGL" );

  histos.h1_qglAve->SetXTitle( "Average Jet QGL" );
  histos.h1_qglProd->SetXTitle( "QGL Product" );

  return histos;

}



void drawSingleHisto( const std::string& outputdir, const std::string& saveName, TH1D* h1_bg, const std::string& bgName, TH1D* h1_sig0, const std::string& sigName0 ) {


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  bool drawNormalized = !(saveName=="qFrac");
  float yMax = (drawNormalized ) ? h1_bg->GetMaximum()*1.5/h1_bg->Integral() : 1.;

  TH2D* h2_axes = new TH2D( "axes", "", h1_bg->GetNbinsX(), h1_bg->GetXaxis()->GetXmin(), h1_bg->GetXaxis()->GetXmax(), 10, 0., yMax );
  if( h1_bg->GetXaxis()->GetBinLabel(1)!="" ) {
    for( unsigned iBinx=1; iBinx< h1_bg->GetNbinsX()+1; ++iBinx ) {
      h2_axes->GetXaxis()->SetBinLabel(iBinx, h1_bg->GetXaxis()->GetBinLabel(iBinx));
    }
  } else {
    h2_axes->SetXTitle( h1_bg->GetXaxis()->GetTitle() );
  }
  h2_axes->Draw();

  h1_bg->SetFillColor( 29 );
  h1_bg->SetLineColor( kBlack );

  h1_sig0->SetLineColor( kRed+1 );
  h1_sig0->SetLineWidth( 2 );

  if( drawNormalized ) {
    h1_bg->DrawNormalized("histo same");
    h1_sig0->DrawNormalized("histo same");
  } else {
    h1_bg->Draw("histo same");
    h1_sig0->Draw("histo same");
  }

  TLegend* legend = new TLegend( 0.2, 0.65, 0.55, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_bg, bgName.c_str(), "F" );
  legend->AddEntry( h1_sig0, sigName0.c_str(), "F" );
  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/%s.png", outputdir.c_str(), saveName.c_str()));
  c1->SaveAs(Form("%s/%s.eps", outputdir.c_str(), saveName.c_str()));
  c1->SaveAs(Form("%s/%s.pdf", outputdir.c_str(), saveName.c_str()));

  delete c1;

}
