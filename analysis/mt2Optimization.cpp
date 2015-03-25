#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLine.h"
#include "TList.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"

int main( int argc, char* argv[] ){

  if( argc!=2 ) {
    std::cout << "USAGE: ./mt2Optimization_tree [dir]" << std::endl;
    exit(113);
  }

  std::string outputdir = "MT2bins_noCR_v2";
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string dir( argv[1] );
  std::string mc_fileName = dir + "/analyses.root";

  //std::string GJets_fileName = "GammaControlRegion_PHYS14_v2_Zinv_zurich_onlyHT/data.root";
  std::string GJets_fileName = "GammaControlRegion_PHYS14_v2_Zinv_zurich/data.root";

  MT2Analysis<MT2EstimateTree>* data    = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "data" );
  MT2Analysis<MT2EstimateTree>* top     = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "Top" );
  MT2Analysis<MT2EstimateTree>* qcd     = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "QCD" );
  MT2Analysis<MT2EstimateTree>* wjets   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "WJets" );
  MT2Analysis<MT2EstimateTree>* zjets   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "ZJets" );
  MT2Analysis<MT2EstimateTree>* gjets   = MT2Analysis<MT2EstimateTree>::readFromFile( GJets_fileName, "gammaCRtree" );


//  std::set<MT2Region> gammaRegions = gjets->getRegions();
//  int nGammaRegions=(int) gammaRegions.size();
//  
//  TH1F* hmt2g[nGammaRegions];
      

  std::set<MT2Region> regions = data->getRegions();
  int nRegions=(int) regions.size();

  TH1F* hmt2g[nRegions];

  TH1F* hmt2qcd[nRegions];
  TH1F* hmt2z[nRegions];  
  TH1F* hmt2w[nRegions];  
  TH1F* hmt2tt[nRegions]; 

  THStack* hmt2[nRegions];

  TCanvas* c[nRegions];
  TCanvas* cCR[nRegions];
  gStyle->SetOptStat(0);

  int s=0;
  int nBins=140;
  float mt2min=200.;
  float mt2max=3000.;
  float binWidth=(mt2max-mt2min)/nBins;
  
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    float maxMT2 = (float) (*iR).htMax();
    if( maxMT2<0 )
      maxMT2=13e3;

    MT2EstimateTree* gjets_ =  gjets->get( *(gjets  ->matchRegion( *iR )) );
    
    //std::cout << (gjets  ->matchRegion( *iR ))->getName() << std::endl;

    TTree* mt2gjets = gjets_->tree;

    hmt2g[s]   =new TH1F(Form("hmt2g_%d", s), "", nBins, mt2min, mt2max);

    float mt2g, weight_g;
    mt2gjets->SetBranchAddress("mt2", &mt2g);
    mt2gjets->SetBranchAddress("weight", &weight_g);

    for(int e=0; e < mt2gjets->GetEntries(); ++e) {

      mt2gjets->GetEntry(e);
      hmt2g[s]->Fill(mt2g, weight_g);

    }
    

    float mt2CR_threshold=13e3;
    float yCR=0;
     
    //std::cout << std::endl << "First thresholds for region " << iR->getName().c_str() << std::endl;

    for(int i=nBins+1; i>0; --i){
      
      yCR += hmt2g[s]->GetBinContent(i);
 
      if( yCR > 2 && hmt2g[s]->GetBinCenter(i) - binWidth/2 < maxMT2 ){
    	
    	mt2CR_threshold = hmt2g[s]->GetBinCenter(i) - binWidth/2;
    	//std::cout << mt2CR_threshold << std::endl;
    	
    	yCR=0;
     	
	break;      
      }
      
    }  
  
    //Not CR
    mt2CR_threshold=14e3;

    //TTree* this_data  = data  ->get(*iR)->tree;
    TTree* mt2qcd   = qcd   ->get(*iR)->tree;
    TTree* mt2zjets = zjets ->get(*iR)->tree;
    TTree* mt2wjets = wjets ->get(*iR)->tree;
    TTree* mt2top   = top   ->get(*iR)->tree;
    
    hmt2qcd[s] =new TH1F(Form("hmt2qcd_%d", s), "", nBins, mt2min, mt2max); 
    hmt2qcd[s]->SetFillColor(401);
    hmt2qcd[s]->SetLineColor(401);
    
    hmt2z[s]   =new TH1F(Form("hmt2z_%d", s), "", nBins, mt2min, mt2max); 
    hmt2z[s]->SetFillColor(419);
    hmt2z[s]->SetLineColor(419);
    
    hmt2w[s]   =new TH1F(Form("hmt2w_%d", s), "", nBins, mt2min, mt2max); 
    hmt2w[s]->SetFillColor(417);
    hmt2w[s]->SetLineColor(417);
    
    hmt2tt[s]  =new TH1F(Form("hmt2tt_%d", s), "", nBins, mt2min, mt2max); 
    hmt2tt[s]->SetFillColor(855);
    hmt2tt[s]->SetLineColor(855);
    
    hmt2[s]    =new THStack("hmt2", iR->getName().c_str());
    
    cCR[s]= new TCanvas(Form("%s%s.eps", "CR_", iR->getName().c_str()), Form("%s%s.eps", "CR_", iR->getName().c_str()), 600, 600);
    cCR[s]->cd();
    gPad->SetLogy();

    hmt2g[s]->Draw("hist");

    TLine* lCR=new TLine(mt2CR_threshold, 0, mt2CR_threshold, hmt2g[s]->GetMaximum());
    lCR->SetLineColor(2);
    lCR->SetLineWidth(3);
    lCR->Draw("same");

    cCR[s]->SaveAs(Form("%s/%s%s.eps", outputdir.c_str(), "CR_", iR->getName().c_str()));
    
    c[s]= new TCanvas(Form("%s.eps", iR->getName().c_str()), Form("%s.eps", iR->getName().c_str()), 600, 600);
    c[s]->cd();
    gPad->SetLogy();

    float mt2z, mt2w, mt2t, mt2q;
    float weight_z, weight_w, weight_t, weight_q;

    mt2qcd->SetBranchAddress("mt2", &mt2q);
    mt2qcd->SetBranchAddress("weight", &weight_q);
 
    mt2top->SetBranchAddress("mt2", &mt2t);
    mt2top->SetBranchAddress("weight", &weight_t);
 
    mt2wjets->SetBranchAddress("mt2", &mt2w);
    mt2wjets->SetBranchAddress("weight", &weight_w);
    
    mt2zjets->SetBranchAddress("mt2", &mt2z);
    mt2zjets->SetBranchAddress("weight", &weight_z);

    for(int e=0; e < mt2zjets->GetEntries(); ++e) {
      
      mt2zjets->GetEntry(e);
      hmt2z[s]->Fill(mt2z, weight_z);
  
    }
    
    for(int e=0; e < mt2wjets->GetEntries(); ++e) {

      mt2wjets->GetEntry(e);
      hmt2w[s]->Fill(mt2w, weight_w);

    }

    for(int e=0; e < mt2top->GetEntries(); ++e) {

      mt2top->GetEntry(e);
      hmt2tt[s]->Fill(mt2t, weight_t);

    }

    for(int e=0; e < mt2qcd->GetEntries(); ++e) {

      mt2qcd->GetEntry(e);
      if( mt2q < 300. )
	hmt2qcd[s]->Fill(mt2q, weight_q);

    }
 
    hmt2qcd[s]->Sumw2();
    hmt2z[s]->Sumw2();
    hmt2w[s]->Sumw2();
    hmt2tt[s]->Sumw2();

    hmt2[s]->Add(hmt2qcd[s]);
    hmt2[s]->Add(hmt2z[s]);
    hmt2[s]->Add(hmt2w[s]);
    hmt2[s]->Add(hmt2tt[s]);
    
    hmt2[s]->Draw("hist");
    
    hmt2[s]->GetXaxis()->SetTitle("M_{T2} [GeV]");
    hmt2[s]->GetYaxis()->SetTitle("Events");
    hmt2[s]->GetYaxis()->SetTitleOffset(1.5);
    
    hmt2[s]->Draw("hist");
    
    float max = hmt2[s]->GetMaximum();
    TLine* l1=new TLine(mt2CR_threshold, 0, mt2CR_threshold, max);
    l1->SetLineColor(2);
    l1->SetLineWidth(3);
    l1->Draw("same");
    
    float mt2_threshold=13e3;
    float nMT2bins=0;
    float y=0., yNext=0.;
    float yqcd=0., yqcdNext=0;
    int iLast=nBins+1;
    float mt2_lastThreshold = maxMT2;
    float mt2_firstThreshold = mt2min;

    std::cout << std::endl << "Thresholds for region " << iR->getName().c_str() << std::endl;
    if(hmt2qcd[s]->Integral(1, iLast)/(hmt2qcd[s]->Integral(1, iLast)+hmt2z[s]->Integral(1, iLast)+hmt2w[s]->Integral(1, iLast)+hmt2tt[s]->Integral(1, iLast)) > 0.1)
      std::cout << "ATTENTION: QCD is already dominating! Inclusive fraction = " << hmt2qcd[s]->Integral(1, iLast)/(hmt2qcd[s]->Integral(1, iLast)+hmt2z[s]->Integral(1, iLast)+hmt2w[s]->Integral(1, iLast)+hmt2tt[s]->Integral(1, iLast))  << std::endl;

    for(int i=1; i < nBins+1; ++i){

      y += hmt2z[s]->GetBinContent(i);
      y += hmt2w[s]->GetBinContent(i);
      y += hmt2tt[s]->GetBinContent(i);      
      yqcd += hmt2qcd[s]->GetBinContent(i);
      yCR += hmt2g[s]->GetBinContent(i);

      yNext += hmt2z[s]->GetBinContent(i+1);
      yNext += hmt2w[s]->GetBinContent(i+1);
      yNext += hmt2tt[s]->GetBinContent(i+1);
      yqcdNext += hmt2qcd[s]->GetBinContent(i+1);
      
      if( yqcd/(y+yqcd) > 0.1 ){
	
    	y=0;
    	yqcd=0;
	yCR=0;
	yNext=0;
	yqcdNext=0;

	continue;
	
      }
      else {
	
	mt2_threshold = hmt2qcd[s]->GetBinCenter(i) - binWidth/2;
	TLine* l=new TLine(mt2_threshold, 0, mt2_threshold, max);
        l->SetLineColor(1);
        l->SetLineWidth(3);
        l->Draw("same");
        c[s]->Update();

	y=0;
        yqcd=0;
        yCR=0;
	yNext=0;
	yqcdNext=0;

	break;
	
      }
      
    }
    
    if(mt2_threshold < maxMT2)
      mt2_firstThreshold=mt2_threshold;
    else
      mt2_firstThreshold=mt2min;

    mt2_threshold=13e3;
    
    for(int i=nBins+1; i>0; --i){
      
      y += hmt2z[s]->GetBinContent(i);
      y += hmt2w[s]->GetBinContent(i);
      y += hmt2tt[s]->GetBinContent(i);      
      yqcd += hmt2qcd[s]->GetBinContent(i);
      yCR += hmt2g[s]->GetBinContent(i);

      yNext += hmt2z[s]->GetBinContent(i-1);
      yNext += hmt2w[s]->GetBinContent(i-1);
      yNext += hmt2tt[s]->GetBinContent(i-1);
      yqcdNext += hmt2qcd[s]->GetBinContent(i-1);
      
      if( y+yqcd+yNext+yqcdNext > 1. && yqcd/(y+yqcd) < 0.1  && ( hmt2qcd[s]->GetBinCenter(i) - binWidth/2 < mt2_threshold - 99. ) && hmt2qcd[s]->GetBinCenter(i) - binWidth/2 > mt2_firstThreshold + 99. &&  hmt2qcd[s]->GetBinCenter(i) - binWidth/2 < maxMT2 ){
    	
    	mt2_threshold = hmt2qcd[s]->GetBinCenter(i) - binWidth/2;	
    	//std::cout << mt2_threshold << std::endl;

//	if( mt2_threshold < mt2min + 59. )
//	  mt2_threshold = mt2min;
    
	TLine* l=new TLine(mt2_threshold, 0, mt2_threshold, max);
	l->SetLineColor(1);
	l->SetLineWidth(3);
	l->Draw("same");
	c[s]->Update();
	
    	y=0;
    	yqcd=0;
	yCR=0;
    	iLast=i;
	++nMT2bins;
	
	break;      
      }
      
    }

    if( mt2_threshold < 13e3 )
      mt2_lastThreshold=mt2_threshold;
    else
      mt2_lastThreshold=mt2_firstThreshold;

    mt2_threshold=mt2_lastThreshold;

    std::cout << mt2_threshold << std::endl;

    for(int i=nBins+1; i > 0; --i){
      
      y += hmt2z[s]->GetBinContent(i);
      y += hmt2w[s]->GetBinContent(i);
      y += hmt2tt[s]->GetBinContent(i);      
      yqcd += hmt2qcd[s]->GetBinContent(i);
      yCR += hmt2g[s]->GetBinContent(i);

      yNext += hmt2z[s]->GetBinContent(i-1);
      yNext += hmt2w[s]->GetBinContent(i-1);
      yNext += hmt2tt[s]->GetBinContent(i-1);
      yqcdNext += hmt2qcd[s]->GetBinContent(i-1);
      
      if( y+yqcd+yNext+yqcdNext > 1. && hmt2qcd[s]->GetBinCenter(i) - binWidth/2  > mt2_firstThreshold + 99. && hmt2qcd[s]->GetBinCenter(i) - binWidth/2 < mt2_lastThreshold - 99. && hmt2qcd[s]->GetBinCenter(i) - binWidth/2 < mt2_threshold - 99. ){
    	
    	mt2_threshold = hmt2qcd[s]->GetBinCenter(i) - binWidth/2;	
    	std::cout << mt2_threshold << std::endl;
    
	TLine* l=new TLine(mt2_threshold, 0, mt2_threshold, max);
	l->SetLineColor(1);
	l->SetLineWidth(3);
	l->Draw("same");
	c[s]->Update();
	
    	y=0;
    	yqcd=0;
	yCR=0;
    	iLast=i;
	++nMT2bins;
	
    	continue;      
      }
      
    }
    
    std::cout << mt2_firstThreshold << std::endl;
  
//    if(iLast==nBins+1)
//      std::cout << "Total yield = " << hmt2qcd[s]->Integral(1, iLast)+(hmt2qcd[s]->Integral(1, iLast)+hmt2z[s]->Integral(1, iLast)+hmt2w[s]->Integral(1, iLast)+hmt2tt[s]->Integral(1, iLast)) << ", of which QCD is " << hmt2qcd[s]->Integral(1, iLast)/(hmt2qcd[s]->Integral(1, iLast)+hmt2z[s]->Integral(1, iLast)+hmt2w[s]->Integral(1, iLast)+hmt2tt[s]->Integral(1, iLast)) << std::endl;
    
//    if(hmt2qcd[s]->Integral(1, iLast-1)/(hmt2qcd[s]->Integral(1, iLast-1)+hmt2z[s]->Integral(1, iLast-1)+hmt2w[s]->Integral(1, iLast-1)+hmt2tt[s]->Integral(1, iLast-1)) > 0.1)
//      std::cout << "QCD fraction in first bin = " << hmt2qcd[s]->Integral(1, iLast-1)/(hmt2qcd[s]->Integral(1, iLast-1)+hmt2z[s]->Integral(1, iLast-1)+hmt2w[s]->Integral(1, iLast-1)+hmt2tt[s]->Integral(1, iLast-1)) << std::endl;

    c[s]->SaveAs(Form("%s/%s%s.eps", outputdir.c_str(), "",iR->getName().c_str()));

    ++s;
  }
  
}
