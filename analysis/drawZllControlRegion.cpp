#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGaxis.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"


#include "../interface/MT2DrawTools.h"


#include <iostream>
#include "string.h"


#define mt2_cxx
#include "interface/mt2.h"

double lumiErr = 0.12;
bool shapeNorm = false;
bool HFveto = false;


void drawMll( const std::string& outputdir,  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, MT2Analysis<MT2EstimateTree>* data, bool of, float lumi ); 

void randomizePoisson( TH1* histo );

void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data,const MT2Region thisRegion, std::string cut, float lumi);


void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );



int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";


  if( argc>2 ) {
    std::string normType(argv[2]);
    if( normType=="lumi" ) shapeNorm=false;
    else if( normType=="shape" ) shapeNorm=true;
    else {
      std::cout << "-> Only 'lumi' and 'shape' are supported normTypes." << std::endl;
      exit(17);
    }
  }


  if( argc<2 ) {
    std::cout << "USAGE: ./coputeZllGammaRatio [configFileName] regionSet" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string configFileName(argv[1]);
  MT2Config cfg( configFileName);
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat"; 
  std::string samples = cfg.mcSamples();

  regionsSet = cfg.regionsSet();


  std::string outputdir = cfg.getEventYieldDir() + "/zllPurity";
  std::string outputdir_of = cfg.getEventYieldDir() + "/zllPurity_of";

  // std::string outputdir( Form("ZllData_%s", configFileName.c_str() ) );
  // std::string outputdir_of( Form("ZllData_OF_%s", configFileName.c_str()  ) );


  std::cout << "-> Using regions: " << regionsSet << std::endl;

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  double intpart;
  double fracpart = modf(cfg.lumi(), &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  //outputdir += suffix;
  // outputdir_of += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));




  std::string ZllDir = cfg.getEventYieldDir() + "/zllControlRegion";
  std::string ZllDir_of = cfg.getEventYieldDir() + "/zllControlRegion";

  //  std::string ZllDir = "ZllPurity_" + configFileName;
  //  std::string ZllDir_of = "ZllPurity_OF_" + configFileName;


  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run zllPurityTrees first. I need to get the yields from there." << std::endl;    std::cout << "-> Thank you for your cooperation." << std::endl;    exit(197);
  }



  
  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str()  ), "QCD");
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "ZJets");

  
 
  //OPPOSITE FLAVOR TREES
  MT2Analysis<MT2EstimateTree>* Zll_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "DYJets");

  MT2Analysis<MT2EstimateTree>* qcd_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str()  ), "QCD");
  MT2Analysis<MT2EstimateTree>* top_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "Top");
  MT2Analysis<MT2EstimateTree>* wjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "ZJets");
  
  MT2Analysis<MT2EstimateTree>* data_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data_of.root", ZllDir_of.c_str() ) , "zllCR_of");

  Zll_of->setFullName("Z+jets");
  wjets_of->setFullName("W+jets");
  zjets_of->setFullName("Z#nu#nu+jets");
  data_of->setFullName("Data");



  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ZllDir.c_str() ) , "zllCR");
 

  data->setFullName("Data");

  Zll->setFullName("Z+jets");
  wjets->setFullName("W+jets");
  zjets->setFullName("Z#nu#nu+jets");



  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  bgYields.push_back( Zll );
  bgYields.push_back( qcd );
  bgYields.push_back( wjets );
  bgYields.push_back( zjets );
  bgYields.push_back( top );


  // drawMll( outputdir, bgYields, data,  0 , cfg.lumi() );
 
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of; 
  bgYields_of.push_back( Zll_of );
  bgYields_of.push_back( qcd_of );
  bgYields_of.push_back( wjets_of );
  bgYields_of.push_back( zjets_of );
  bgYields_of.push_back( top_of );
  
 
  // drawMll( outputdir_of, bgYields_of, data_of,  1, cfg.lumi() );
  








  MT2DrawTools::setStyle();
 


  std::string    selection = "weight*(abs(Z_mass-91.19)<20 && Z_pt>0)";
  // std::string    selection = "weight*(abs(Z_mass-91.19)<20 && Z_pt>0 && nJetHF30==0)";
  /*
  std::string selection_mass = "weight*(Z_mass>50 && Z_pt>180)";
  std::string      selection_mass_el = "weight*(Z_mass>50 && Z_pt>180 && Z_lepId==11)";
  std::string      selection_mass_mu = "weight*(Z_mass>50 && Z_pt>180 && Z_lepId==13)";


  std::string      selection_mass = "weight";
 
  std::string      selection =  "weight*(Z_mass>20 && Z_pt> 0 )";
  //selection2 = "weight*(Z_mass>80)";
 
  std::string      selection3 = "weight*(Z_mass>70)";
  std::string      selection4 ="weight*(Z_mass>80)";

  std::string selection = "weight*(Z_pt>180 && abs(Z_mass-91.19)<20)";
  */
  drawYields( cfg, data, bgYields, "mt2" , "mt2" , selection, 24, 0, 600, "M_{T2}", "GeV" );
  drawYields( cfg, data, bgYields, "ht" , "ht" , selection, 24, 0, 600, "H_{T}", "GeV" );
  drawYields( cfg, data, bgYields, "met" , "met" , selection, 24, 0, 600, "ME_{T}", "GeV" );

  drawYields( cfg, data, bgYields, "nJets", "nJets", selection, 10, 1.5, 11.5, "Number of Jets (p_{T} > 30 GeV)", "" );
  drawYields( cfg, data,  bgYields, "nBJets", "nBJets", selection, 6, -0.5, 5.5, "Number of b-Jets (p_{T} > 20 GeV)", "" );

  drawYields( cfg, data, bgYields, "Z_pt" , "Z_pt" , selection, 36 , 0, 900, "Z p_{T}", "GeV" );




  drawYields( cfg, data, bgYields, "Z_lepId" , "Z_lepId" , selection, 5, 10,15 , "Lepton Id", "" );


  //drawYields( cfg, data, bgYields, "mt2" , "mt2" , selection, 24, 0, 600, "M_{T2}", "GeV" );

  std::string    selection2 = "weight*(Z_pt>0)";
  drawYields( cfg, data, bgYields, "Z_mass" , "Z_mass" , selection2, 40, 50, 150, "M_{ll}", "GeV" );

  drawYields( cfg, data_of, bgYields_of, "Z_mass_of" , "Z_mass" , selection2, 50, 0, 250, "M_{ll}", "GeV" );






  return 0;
}




















void drawMll( const std::string& outputdir, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields,  MT2Analysis<MT2EstimateTree> * data, bool of , float lumi) {

  MT2DrawTools::setStyle();
  /*
  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other = zll 
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
  }
  */
  TH1F::AddDirectory(kTRUE);

  std::string fullPath = outputdir;  std::string fullPathPlots = outputdir + "/plots";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  std::set<MT2Region> MT2Regions = bgYields[0]->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );



    float bins_nVert[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20,21,22,23,24 ,25, 26,27,28,29,30};
    
    float bins_Z_lepId[] = { 10,11,12,13,14,15};
    
    float bins_nJets[] = {2,3,4,5,6,7,8,9,10,11, 12};
  
    float bins_nJetsHF[] = {0,1,2,3,4,5,6,7,8,9,10,11, 12};
    float bins_nBJets[] = {0,1,2,3,4,5,6};
    //in MT2
    float bins_mt2[] = {0,25,50,75, 100,125, 150,175,200,225,250,275, 300, 350,400,500,600 };
    //float bins_mt2[] = {0,25,50,75, 100,125, 150,175,200,225,250,275, 300,350, 400};
    // float bins_mt2[] = {0,25,50,75, 100,125, 150,175,200,225,250,275, 300,325,350,375, 400,425,450};
    // float bins_mt2[] = {0,25,50,75, 100,125, 150,175,200,225,250,275, 300,325,350,375, 400,425,450,475,500,525,550,575, 600};
    //   float bins_mt2[] = {0,50, 100,150,200,250, 300, 400,500,600,700,800,900,1000,1100,1200,1300,1400, 1500};
    //in HT
    float bins_ht[] =  {0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400, 1500,1600,1700,1800,1900, 2000};

    float bins_mll[] = {50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150};

    float bins_mll_of[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
    float bins_lepPt[] = {0,20,40,60,80,100,120,140,160,180,200,220, 240,260,280,300,320,340,360,380,400};
   
    float bins_lepEta[] = {-3,-2.5,-2.,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3};

    float bins_met[] = {0,50,100,150,200,250,300,350,400,450, 500, 600,700,800,900, 1000};
    float bins_Zpt[] = {0,50,100,150, 200,250,300,350,400,450, 500,550,600,650,700, 750,800,850,900,950, 1000};
    //  float bins_ht[] =  {450,500, 600,700, 800, 900,1000,1100, 1200, 1300, 1500,2000};


    /*    float bins_nJets[] = {2,4,7,12};
    float bins_nBJets[] = {0,1,2,3,6};
    //in MT2
    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 , 1900};
    //in HT
    float bins_ht[] =  {450,575,1000,1500,2000};
    */

    std::string cut;
    std::string cut_HF;
    std::string cut2;
 
    std::string cut3 ;
    std::string cut4;

    std::string cut_mass ;
    std::string cut_mass_el ;
    std::string cut_mass_mu; 

    if(of==0){ 
      // cut =  "weight";
      cut = "weight*(abs(Z_mass-91.19)<20 && Z_pt>0)";
 
      cut_HF = "weight*(abs(Z_mass-91.19)<20 &&Z_pt>0)";
      //    cut_HF = "weight*(abs(Z_mass-91.19)<20 && Z_pt>0)";
 
      cut_mass = "weight*(Z_mass>50 && Z_pt>0)";
      cut_mass_el = "weight*(Z_mass>50 && Z_pt>0 && Z_lepId==11)";
      cut_mass_mu = "weight*(Z_mass>50 && Z_pt>0 && Z_lepId==13)";

      // cut3 = "weight*(abs(Z_mass-91.19)<20&&nBJets<2)";
      //  cut4 ="weight*(abs(Z_mass-91.19)<10&&nBJets<2)";
    }else{

      cut_mass = "weight";
 
      cut =  "weight*(Z_mass>20 && Z_pt> 0 )";
      //cut2 = "weight*(Z_mass>80)";
 
      cut3 = "weight*(Z_mass>70)";
      cut4 ="weight*(Z_mass>80)";
    }
  
    if(HFveto == true){
      cut.pop_back();
      cut_mass.pop_back();
      cut += " && nJetHF30==0)";
      cut_mass += " && nJetHF30==0)";
      if(of==0){
      cut_mass_el.pop_back();      cut_mass_mu.pop_back();
      cut_mass_el += " && nJetHF30==0)";
      cut_mass_mu += " && nJetHF30==0)";
      }
    }

    std::string cut_nJets3 = "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
   
    
    drawStacks( fullPathPlots, bins_nJetsHF,sizeof(bins_nJetsHF)/sizeof(float)-1,  "nJetHF30", bgYields, data , thisRegion, cut_HF , lumi );

    drawStacks( fullPathPlots, bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields, data , thisRegion, cut , lumi );
    drawStacks( fullPathPlots, bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , data, thisRegion , cut , lumi);
    
    drawStacks( fullPathPlots, bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "mt2", bgYields , data, thisRegion  , cut , lumi  );
    
    drawStacks( fullPathPlots, bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "ht", bgYields , data , thisRegion, cut  , lumi);
    drawStacks( fullPathPlots, bins_met,sizeof(bins_met)/sizeof(float)-1,  "met", bgYields , data, thisRegion  , cut , lumi  );     
    drawStacks( fullPathPlots, bins_Zpt,sizeof(bins_Zpt)/sizeof(float)-1,  "Z_pt", bgYields , data, thisRegion  , cut , lumi  );     

    drawStacks( fullPathPlots, bins_nVert,sizeof(bins_nVert)/sizeof(float)-1,  "nVert", bgYields , data, thisRegion  , cut , lumi  );     

    drawStacks( fullPathPlots, bins_Z_lepId,sizeof(bins_Z_lepId)/sizeof(float)-1,  "Z_lepId", bgYields , data, thisRegion  , cut , lumi  );     


    drawStacks( fullPathPlots, bins_lepPt,sizeof(bins_lepPt)/sizeof(float)-1,  "lep_pt0", bgYields , data, thisRegion  , cut , lumi  );   
    drawStacks( fullPathPlots, bins_lepPt,sizeof(bins_lepPt)/sizeof(float)-1,  "lep_pt1", bgYields , data, thisRegion  , cut , lumi  );     

    drawStacks( fullPathPlots, bins_lepEta,sizeof(bins_lepEta)/sizeof(float)-1,  "lep_eta0", bgYields , data, thisRegion  , cut , lumi  );   
    drawStacks( fullPathPlots, bins_lepEta,sizeof(bins_lepEta)/sizeof(float)-1,  "lep_eta1", bgYields , data, thisRegion  , cut , lumi  );     


    if(of==0){
      drawStacks( fullPathPlots,  bins_mll,sizeof(bins_mll)/sizeof(float)-1,  "Z_mass", bgYields , data, thisRegion  , cut_mass , lumi); 

      drawStacks( fullPathPlots,  bins_mll,sizeof(bins_mll)/sizeof(float)-1,  "Z_mass", bgYields , data, thisRegion  , cut_mass_el , lumi); 

      drawStacks( fullPathPlots,  bins_mll,sizeof(bins_mll)/sizeof(float)-1,  "Z_mass", bgYields , data, thisRegion  , cut_mass_mu , lumi);  
    }else{

      drawStacks( fullPathPlots,  bins_mll_of,sizeof(bins_mll_of)/sizeof(float)-1,  "Z_mass", bgYields , data, thisRegion  , cut_mass , lumi); 

    }

       

    /*
      drawStacks( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , data,thisRegion, cut2, lumi );
      drawStacks( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , data, thisRegion , cut2 , lumi);
      drawStacks( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "mt2", bgYields , data, thisRegion  , cut2 , lumi );
      drawStacks( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "ht", bgYields , data, thisRegion, cut2  , lumi );
    */

    /*
    delete c1;
    delete h2_axes;
    delete histo_bg;
    delete h_data; delete histo_data;
    */

    //   }//end of loop over beees


  }// for MT2 regions

}











void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data, const MT2Region thisRegion, std::string cut , float lumi){
 
  std::vector<int> colors;

    colors.push_back(430); // other = zll 
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
  

  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  float xMin = binss[0];
  float xMax = binss[size];

 TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();
 

  TH1D* h_data = new TH1D("h_data","", size, bins);
  TTree *data_Tree = data->get(thisRegion)->tree;
  data_Tree->Project("h_data",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
   
  TH1D* mc_sum = new TH1D("mc_sum", "",size,bins);
  for( unsigned i=0; i<bgYields.size(); ++i ) { 
    if(name=="mt2" && (i>0 && i<4)) continue;
    int index = bgYields.size() - i - 1;
    TH1D* h1_bg_mc = new TH1D("h1_bg_mc","", size  , bins);
    h1_bg_mc->Sumw2();
    TTree *bgTree = bgYields[index]->get(thisRegion)->tree;
    bgTree->Project("h1_bg_mc",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    h1_bg_mc->SetFillColor( colors[index] );
    h1_bg_mc->SetLineColor( kBlack );
    mc_sum->Add(h1_bg_mc);
  }



    std::cout << "Integrals: " << h_data->Integral(0, size+1) << "\t" << mc_sum->Integral(0, size+1) << std::endl;
    float scaleFactor2 = h_data->Integral(0, size+1)/mc_sum->Integral(0, size+1);
    if( shapeNorm ) 
      std::cout << "SF: " << scaleFactor2 << std::endl;

  // std::string cut = "weight*( abs(Z_mass-91.19)<20)";
  TH1D* h_bg = new TH1D("h_bg","",size,bins);

  THStack bgStack("bgStack", "");
  for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
    if(name=="mt2" && (i>0 && i<4)) continue;
    int index = bgYields.size() - i - 1;
    TH1D* h1_bg = new TH1D("h1_bg","", size  , bins);
    h1_bg->Sumw2();
    TTree *bgTree = bgYields[index]->get(thisRegion)->tree;
    bgTree->Project("h1_bg",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    if( shapeNorm ) {
      h1_bg->Scale( scaleFactor2 );
    } 
    h1_bg->SetFillColor( colors[index] );
    h1_bg->SetLineColor( kBlack );
    bgStack.Add(h1_bg);
    h_bg->Add(h1_bg);
  }
    
 
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(kBlack);

  TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_data);
  gr_data->SetMarkerStyle(20);
  gr_data->SetMarkerSize(1.2);
 

 
  // float yMax = 10*(bgStack.GetMaximum());
  //float yMax2 = 10*(h_data->GetMaximum());
  
  float yMax = 1.3*(bgStack.GetMaximum());
  float yMax2 = 1.3*(h_data->GetMaximum());
  if(yMax2>yMax) yMax = yMax2;

   TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0.1, yMax );
  // TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0., yMax );
  if(name  == "ht")  
    h2_axes->SetXTitle("H_{T} [GeV]");
  else if(name == "mt2")  
    h2_axes->SetXTitle("M_{T2} (Z#rightarrowll Removed) [GeV]");
  else if(name == "Z_mass")  
    h2_axes->SetXTitle("M_{ll} [GeV]");
  else if(name == "Z_pt")  
    h2_axes->SetXTitle("Boson p_{T} [GeV]");
  else if(name == "lep_pt0")  
    h2_axes->SetXTitle("Leading Lepton p_{T} [GeV]");
  else if(name == "lep_pt1")  
    h2_axes->SetXTitle("Sub-Leading Lepton p_{T} [GeV]");
  else if(name == "lep_eta0")  
    h2_axes->SetXTitle("Leading Lepton #eta [GeV]");
  else if(name == "lep_eta1")  
    h2_axes->SetXTitle("Sub-Leading Lepton #eta [GeV]");
  else if(name == "met")  
    h2_axes->SetXTitle("ME_{T} [GeV]");
  else if(name == "Z_lepId")  
    h2_axes->SetXTitle("Lepton Id");
  else
    h2_axes->SetXTitle(name.c_str());
  h2_axes->SetYTitle("Events");
  // if(name=="mt2") h2_axes->SetYTitle("Events / (25 GeV)");


  if(name == "Z_mass" && cut == "weight")  
    h2_axes->SetXTitle("M_{e^{#pm}#mu^{#mp}} [GeV]");

  if(cut == "weight*(Z_mass>50 && Z_pt>0 && Z_lepId==11)") h2_axes->SetXTitle("M_{e^{+}e^{-}} [GeV]");
  else if(cut == "weight*(Z_mass>50 && Z_pt>0 && Z_lepId==13)") h2_axes->SetXTitle("M_{#mu^{+}#mu^{-}} [GeV]"); 
  
  //gPad->SetLogy();

  h2_axes->Draw();

  if( shapeNorm ) {
    TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
    normText->SetFillColor(0);
    normText->SetTextSize(0.04);
    normText->AddText( "#splitline{Shape}{Norm.}" );
    pad1->cd();
    //   normText->Draw("same");
  }


  int legendLength = bgYields.size()+1;
  if(name=="mt2") legendLength = bgYields.size()+1-3;

  TLegend* legend = new TLegend( 0.73, 0.9-(legendLength)*0.06, 0.93, 0.9 );
  legend->SetTextSize(0.04);
  //legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry(gr_data, "Data", "p");
  for( unsigned i=0; i<bgYields.size(); ++i ) {  
    if(name=="mt2" && (i>0 && i<4)) continue;
    TH1D* h1_bg1 = bgYields[i]->get(thisRegion)->yield;
    h1_bg1->SetFillColor( colors[i] );
    h1_bg1->SetLineColor( kBlack );
    legend->AddEntry( h1_bg1, bgYields[i]->getFullName().c_str(), "F" );
  }

    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    for( unsigned i=0+1; i<niceNames.size(); ++i ) {
      //    for( unsigned i=0; i<niceNames.size(); ++i ) {
      float yMaxText = 0.9-(float)(i-1)*0.05;
      float yMinText = yMaxText - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMinText, 0.35, yMaxText, "brNDC" );
      regionText->SetTextSize(0.04);
      //    regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );
      regionText->Draw("same");
    }

 
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);

  labelTop->Draw("same");
  legend->Draw("same");
  bgStack.Draw("histo same");
  // h_data->DrawCopy("P same");
  gr_data->Draw("p same");
  gPad->RedrawAxis();




  double error_data;
  double integral_data = h_data->IntegralAndError(0, size+1, error_data);
  double error_mc;
  double integral_mc = h_bg->IntegralAndError(0, size+1, error_mc);
  float scaleFactor = integral_data/integral_mc;
  double error_datamc = scaleFactor*(sqrt( (error_data/integral_mc)*(error_data/integral_mc) + (integral_data*error_mc/(integral_data*integral_data))*(integral_data*error_mc/(integral_data*integral_data)) ));
  // std::cout << error_datamc << std::endl;


  TPaveText* ratioText = new TPaveText( 0.135, -0.049, 0.4, 0.1 , "brNDC" );
  ratioText->SetTextSize(0.04);
  //ratioText->SetTextFont(40);
  ratioText->SetFillColor(0);
  ratioText->SetTextAlign(11);
  ratioText->SetTextColor(2);
  // ratioText->AddText( Form("Data/MC = %.2f", ratio) );
  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
   
  if( !(shapeNorm) ) 
    ratioText->Draw("same");


  TGraphAsymmErrors* gr_ratio = MT2DrawTools::getRatioGraph(  h_data, h_bg );
  gr_ratio->SetMarkerStyle(20);
  gr_ratio->SetMarkerSize(1.2);
 

  canny->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

  TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.0, 2.0  );
  h2_axes_rat->SetYTitle("Data / MC");

  h2_axes_rat->GetXaxis()->SetTitleSize(0.2);
  h2_axes_rat->GetXaxis()->SetTitleOffset(5);
  h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
  h2_axes_rat->GetXaxis()->SetTickLength(0.09);
  h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_rat->GetYaxis()->SetTitleSize(0.2);
  h2_axes_rat->GetYaxis()->SetTitleOffset(0.34);
  h2_axes_rat->GetYaxis()->SetLabelSize(0.17);
  
  TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
  lineSF->SetLineColor(2);

  double x[2]={(double)xMin, (double)xMax};
  double xerr[2]={0., 0.};
  double yerr[2]={error_datamc, error_datamc};
  double y[2]={integral_data/integral_mc, integral_data/integral_mc};
  TGraphErrors* SFband = new TGraphErrors(2, x, y, xerr, yerr);
  SFband->SetLineColor(0);
  SFband->SetFillColor(2);
  SFband->SetFillStyle(3244);


  TLine *line = new TLine(xMin, 1, xMax, 1);
  line->SetLineColor(kBlack);

  h2_axes_rat->Draw();
  line->Draw("same");

  if( !shapeNorm ){
    lineSF->Draw("same");
    SFband->Draw("3,same");}
 

  h_data->Divide(h_bg);

  gr_ratio->Draw("PE same");

  // h_data->Draw("p same");


  gPad->RedrawAxis();

  canny->cd();

  std::string extension= "";
  if(cut == "weight*(Z_mass>50 && Z_pt>0 && Z_lepId==11)") extension = "el";
  else if(cut == "weight*(Z_mass>50 && Z_pt>0 && Z_lepId==13)") extension ="mu";
  //else if(cut == "weight*(abs(Z_mass-91.19)<20 && Z_pt>0)") extension ="HF";
  
  if(HFveto==true)
    extension += "_wHFveto";

  /*
  if(cut == "weight*(abs(Z_mass-91.19)<20)") extension ="mass";
  else if(cut == "weight*(abs(Z_mass-91.19)<10)") extension ="mass10";
  else if(cut == "weight*(abs(Z_mass-91.19)<10&&nJets>2)") extension ="nJets3";
  else if(cut == "weight*(Z_mass>70)") extension = "OF_70GeV";
  else if(cut == "weight*(Z_mass>80)") extension = "OF_80GeV";
  else  extension = "massNbJets";
  // "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
  */

 
  canny->SaveAs( Form("%s/%s_%s_%s.eps", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_%s_%s.png", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_%s_%s.pdf", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );


 
  delete h2_axes;
  delete canny;
  delete h_bg;
  delete h_data;
  delete gr_data;
  delete gr_ratio;
}














void randomizePoisson( TH1* histo ) {

  TRandom3 rand(11);;
  //  TRandom3 rand(13);


  //  std::set<MT2HTRegion> HTRegions = data->getHTRegions();
  //  std::set<MT2SignalRegion> signalRegions = data->getSignalRegions();

  //  for( std::set<MT2HTRegion>::iterator iHT = HTRegions.begin(); iHT!=HTRegions.end(); ++iHT ) {
  //    for( std::set<MT2SignalRegion>::iterator iSR = signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

  for( int ibin=1; ibin<histo->GetXaxis()->GetNbins()+1; ++ibin ) {

    int poisson_data = rand.Poisson(int(histo->GetBinContent(ibin)));
    histo->SetBinContent(ibin, poisson_data);
    histo->SetBinError(ibin, sqrt(poisson_data));
	  
  }  // for bins

  // return histo;
}
















void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units ) {


  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;



  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other=zll
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
    //colors.push_back(); // other
  }

  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsDataMC";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );

    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);


    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
      std::string thisName = "h1_" + bgYields[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if( selection!="" )
	//tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s/puWeight", selection.c_str()) );
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s", selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), "" );

      // h1_mc->Scale( 16.1/20.38 );


      histos_mc.push_back(h1_mc);
    }

    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    TH1D* mc_sum;
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      if( i==0 ) {
        mc_sum = new TH1D( *histos_mc[i] );
        mc_sum->SetName("mc_sum");
      } else {
        mc_sum->Add( histos_mc[i] );
      }
    }

    std::cout << "Integrals: " << h1_data->Integral(0, nBins+1) << "\t" << mc_sum->Integral(0, nBins+1) << std::endl;
    float scaleFactor = h1_data->Integral(0, nBins+1)/mc_sum->Integral(0, nBins+1);   
    if( shapeNorm ) 
      std::cout << "SF: " << scaleFactor << std::endl;
    /*
    double error_data;
    double integral_data = h1_data->IntegralAndError(0, nBins+1, error_data);
    double error_mc;
    double integral_mc = mc_sum->IntegralAndError(0, nBins+1, error_mc);
    double error_datamc = MT2DrawTools::getSFError(integral_data, error_data, integral_mc, error_mc);
    */

    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = bgYields.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      histos_mc[index]->SetLineColor( kBlack );
      if( shapeNorm ) {
        histos_mc[index]->Scale( scaleFactor );
      }

      //  if( shapeNorm ) {
      //   histos_mc[index]->Scale( scaleFactor );
      // } else {
      // 	histos_mc[index]->Scale( 16.1/20.38 );
      //}
      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    TPad *pad1 = MT2DrawTools::getCanvasMainPad();
    TPad *pad2 = MT2DrawTools::getCanvasRatioPad();
        
    TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
    c1_log->cd();
    TPad *pad1_log = MT2DrawTools::getCanvasMainPad( true );
    TPad *pad2_log = MT2DrawTools::getCanvasRatioPad( true );
 
    float yMaxScale = 1.1;
    float yMax1 = h1_data->GetMaximum()*yMaxScale;
    float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = yMaxScale*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( h1_data->GetNbinsX()<2 ) yMax *=3.;

    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    std::string binWidthText;
    if( binWidth>=1. )         binWidthText = (std::string)Form("%.0f", binWidth);
    else if( binWidth>=0.1 )   binWidthText = (std::string)Form("%.1f", binWidth);
    else if( binWidth>=0.01 )  binWidthText = (std::string)Form("%.2f", binWidth);
    else if( binWidth>=0.001 ) binWidthText = (std::string)Form("%.3f", binWidth);
    else                       binWidthText = (std::string)Form("%.4f", binWidth);

   std::string yAxisTitle;
    if( units!="" ) 
      yAxisTitle = (std::string)(Form("Events / (%s %s)", binWidthText.c_str(), units.c_str()));
    else
      yAxisTitle = (std::string)(Form("Events / (%s)", binWidthText.c_str()));


    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
    pad1->Draw();
    pad1->cd();
    h2_axes->Draw();
    
   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.1, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    pad1_log->Draw();
    pad1_log->cd();
    h2_axes_log->Draw();
   

    std::vector<std::string> niceNames = thisRegion.getNiceNames();

    for( unsigned i=0; i<niceNames.size(); ++i ) {
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.04);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );

      pad1->cd();
      regionText->Draw("same");
  
      pad1_log->cd();
      regionText->Draw("same");
    }
    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      pad1->cd();
      //normText->Draw("same");
      pad1_log->cd();
      // normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
    ratioText->SetTextSize(0.04);
    ratioText->SetTextFont(40);
    ratioText->SetTextColor(2);
    ratioText->SetFillColor(0);
    ratioText->SetTextAlign(11);
    ratioText->AddText( Form("Data/MC = %.2f", scaleFactor) );
    //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
     

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
    lineSF->SetLineColor(2);

    float yMinR=0.0;
    float yMaxR=2.0;

  
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
 
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
   
    //    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );
    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    TPaveText* fitText = MT2DrawTools::getFitText( fSF );


    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
     fitText->Draw("same");
    //  ratioText->Draw("same");

    gPad->RedrawAxis();

   /*
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
   
    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */


    c1->cd();
    //   TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();

    h2_axes_ratio->Draw("");
 
    /*  line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same");
    */
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }

    g_ratio->Draw("PE,same");    
    gPad->RedrawAxis();


    c1_log->cd();
    // TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();

    h2_axes_ratio->Draw(""); 
    
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }
    /*
    line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same"); */
    g_ratio->Draw("PE,same");
    gPad->RedrawAxis();


    c1->SaveAs( Form("%s/%s_%s.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    c1_log->SaveAs( Form("%s/%s_%s_log.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h2_axes;

    delete c1_log;
    delete h2_axes_log;
    
    delete h2_axes_ratio;
    
    delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}
