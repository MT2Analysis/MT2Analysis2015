#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"






int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawSingleLeptonControlRegion [configFileName] [lumi/shape]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  MT2DrawTools::setStyle();

  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  bool shapeNorm=false;
  if( argc>2 ) {
    std::string normType(argv[2]);
    if( normType=="lumi" ) shapeNorm=false;
    else if( normType=="shape" ) shapeNorm=true;
    else {
      std::cout << "-> Only 'lumi' and 'shape' are supported normTypes." << std::endl;
      exit(17);
    }
  }


  std::string slDir = cfg.getEventYieldDir() + "/singleLeptonControlRegion/";
  std::string mcFile = slDir + "/mc.root";
  std::string plotsDir = slDir + "/plots";
  if( shapeNorm ) plotsDir += "_shape";


  MT2DrawTools dt(plotsDir, cfg.lumi());
  dt.set_shapeNorm( shapeNorm );


  MT2Analysis<MT2EstimateTree>* wjet = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "wjet");
  wjet->setFullName("W+jets");
  wjet->setColor(kWJets);

  MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "top");
  top->setFullName("Top");
  top->setColor(kTop);

  MT2Analysis<MT2EstimateTree>* qcd   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "qcd");
  qcd->setFullName("QCD");
  qcd->setColor(kQCD);


  std::string selection = "nJets>=0"; // for now


  MT2Analysis<MT2EstimateTree>* data = 0; //MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "data");
  //data->setFullName("Data");

  std::vector< MT2Analysis<MT2EstimateTree>* > mc;
  mc.push_back(wjet);
  mc.push_back(top);
  mc.push_back(qcd);

  dt.set_data(data);
  dt.set_mc(&mc);

  dt.drawRegionYields_fromTree( "mt2"   , "mt2"   , selection, 60, 0., 600., "M_{T2}", "GeV" );
  dt.drawRegionYields_fromTree( "met"   , "met"   , selection, 50, 0., 500., "Missing E_{T}", "GeV" );
  dt.drawRegionYields_fromTree( "mht"   , "mht"   , selection, 50, 0., 500., "Missing H_{T}", "GeV" );
  dt.drawRegionYields_fromTree( "ht"    , "ht"    , selection, 50, 0., 1000., "H_{T}", "GeV" );
  dt.drawRegionYields_fromTree( "nJets" , "nJets" , selection, 12, 1.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "" );
  dt.drawRegionYields_fromTree( "nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "" );
  dt.drawRegionYields_fromTree( "leptMt", "leptMt", selection, 50, 0., 300., "Lepton+ME_{T} M_{T}", "GeV" );

  selection = "nBJets==1";
  dt.drawRegionYields_fromTree( "leptMt_nB1", "leptMt", selection, 50, 0., 300., "M_{T}", "GeV", "", "#geq2j, b=1" );
  dt.drawRegionYields_fromTree( "topMt", "topMt"  , selection, 50, 0., 300., "Lepton+ME_{T}+b-Jet M_{T}", "GeV", "", "#geq2j, b=1" );
  dt.drawRegionYields_fromTree( "massLeptB", "massLeptB"  , selection, 50, 0., 300., "Lepton+b-Jet Mass", "GeV", "", "#geq2j, b=1" );
  dt.drawRegionYields_fromTree( "deltaR_leptB", "deltaR_leptB"  , selection, 50, 0., 6., "#Delta R( Lepton, b-Jet )", "", "", "#geq2j, b=1" );

  return 0;

}
