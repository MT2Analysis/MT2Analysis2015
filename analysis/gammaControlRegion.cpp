#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateTree.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"


<<<<<<< HEAD
=======
float lumi = 0.1; //fb-1


>>>>>>> e9c4691... CherryPick2, shall conflicts rise?
bool alsoSignals = false;





int round(float d) {
  return (int)(floor(d + 0.5));
}


void computeYield( const MT2Sample& sample, const MT2Config& cfg, 
                   MT2Analysis<MT2EstimateTree>* anaTree,
                   MT2Analysis<MT2EstimateZinvGamma>* prompt=0, MT2Analysis<MT2EstimateZinvGamma>* prompt_pass=0, 
                   MT2Analysis<MT2EstimateZinvGamma>* nip=0, MT2Analysis<MT2EstimateZinvGamma>* nip_pass=0, 
                   MT2Analysis<MT2EstimateZinvGamma>* fake=0, MT2Analysis<MT2EstimateZinvGamma>* fake_pass=0, float isoCut=2.5 );
void roundLikeData( MT2Analysis<MT2EstimateZinvGamma>* data );




int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|            Running gammaControlRegion              |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./gammaControlRegion [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);




  std::string samplesFile = "../samples/samples_" + cfg.mcSamples() + ".dat";
  
  std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "GJets");
  if( samples_gammaJet.size()==0 ) {
    std::cout << "There must be an error: didn't find any gamma+jet files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading QCD samples" << std::endl;

  std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "QCD");
  


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  //std::string outputdir = "GammaControlRegion_%s_%s_%.0ffb", samplesFileName.c_str(), cfg.regionsSet.c_str(), cfg.lumi() ));
  std::string outputdir = cfg.getEventYieldDir() + "/gammaControlRegion"; //"GammaControlRegion_" + cfg.name();
  system(Form("mkdir -p %s", outputdir.c_str()));

  
  MT2Analysis<MT2EstimateZinvGamma>* prompt = new MT2Analysis<MT2EstimateZinvGamma>( "prompt", cfg.regionsSet() );
  MT2Analysis<MT2EstimateZinvGamma>* prompt_pass = new MT2Analysis<MT2EstimateZinvGamma>( "prompt_pass", cfg.regionsSet() );

  MT2Analysis<MT2EstimateZinvGamma>* fake = new MT2Analysis<MT2EstimateZinvGamma>( "fake", cfg.regionsSet() );
  MT2Analysis<MT2EstimateZinvGamma>* fake_pass = new MT2Analysis<MT2EstimateZinvGamma>( "fake_pass", cfg.regionsSet() );

  MT2Analysis<MT2EstimateZinvGamma>* nip = new MT2Analysis<MT2EstimateZinvGamma>( "nip", cfg.regionsSet() );
  MT2Analysis<MT2EstimateZinvGamma>* nip_pass = new MT2Analysis<MT2EstimateZinvGamma>( "nip_pass", cfg.regionsSet() );

  MT2Analysis<MT2EstimateTree>* tree = new MT2Analysis<MT2EstimateTree>( "gammaCRtree", cfg.regionsSet() );
  MT2EstimateTree::addVar( tree, "prompt" );
  MT2EstimateTree::addVar( tree, "iso" );
  


  float isoCut = 2.5; // GeV (absolute iso)

  for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
    computeYield( samples_gammaJet[i], cfg, tree, prompt, prompt_pass, nip, nip_pass, fake, fake_pass, isoCut );
  }

  for( unsigned i=0; i<samples_qcd.size(); ++i ) {
    computeYield( samples_qcd[i], cfg, tree, prompt, prompt_pass, nip, nip_pass, fake, fake_pass, isoCut );
  }


  MT2Analysis<MT2EstimateZinvGamma>* gammaCR_loose = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_loose",  cfg.regionsSet() );
  (*gammaCR_loose) = (*prompt) + (*nip) + (*fake);
 
 
  MT2Analysis<MT2EstimateZinvGamma>* gammaCR         = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR", cfg.regionsSet() );
  MT2Analysis<MT2EstimateZinvGamma>* gammaCR_nipUp   = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_nipUp", cfg.regionsSet() );
  MT2Analysis<MT2EstimateZinvGamma>* gammaCR_nipDown = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_nipDown", cfg.regionsSet() );
  (*gammaCR        ) = (*prompt_pass) + (*fake_pass) +     (*nip_pass);
  (*gammaCR_nipUp  ) = (*prompt_pass) + (*fake_pass) + 2. *(*nip_pass);
  (*gammaCR_nipDown) = (*prompt_pass) + (*fake_pass) + 0.5*(*nip_pass);
  
 
  MT2Analysis<MT2EstimateZinvGamma>* matched = new MT2Analysis<MT2EstimateZinvGamma>( "matched", cfg.regionsSet() ); 
  (*matched) = (*prompt) + (*nip);

  MT2Analysis<MT2EstimateZinvGamma>* matched_pass = new MT2Analysis<MT2EstimateZinvGamma>( "matched_pass", cfg.regionsSet() ); 
  (*matched_pass) = (*prompt_pass) + (*nip_pass);


  MT2Analysis<MT2EstimateSyst>* f = MT2EstimateSyst::makeEfficiencyAnalysis( "f", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt, (MT2Analysis<MT2Estimate>*)matched );
  MT2Analysis<MT2EstimateSyst>* f_pass = MT2EstimateSyst::makeEfficiencyAnalysis( "f_pass", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_pass, (MT2Analysis<MT2Estimate>*)matched_pass );



  MT2Analysis<MT2EstimateSyst>* eff = MT2EstimateSyst::makeEfficiencyAnalysis( "eff", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_pass, (MT2Analysis<MT2Estimate>*)prompt );

  MT2Analysis<MT2EstimateSyst>* purityTight = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)matched_pass, (MT2Analysis<MT2Estimate>*)gammaCR);

  MT2Analysis<MT2EstimateSyst>* purityLoose = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)matched, (MT2Analysis<MT2Estimate>*)gammaCR_loose );


  std::string mcFile = outputdir + "/mc.root";

  gammaCR_loose->writeToFile( mcFile, "RECREATE" );
  gammaCR      ->addToFile( mcFile );
  prompt       ->addToFile( mcFile );
  fake         ->addToFile( mcFile );
  nip          ->addToFile( mcFile );
  prompt_pass  ->addToFile( mcFile );
  fake_pass    ->addToFile( mcFile );
  nip_pass     ->addToFile( mcFile );
  f            ->addToFile( mcFile );
  f_pass       ->addToFile( mcFile );



  std::vector< MT2Analysis< MT2EstimateTree>* > signals;
  if( alsoSignals ) {

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples" << std::endl;
    std::string samplesFile_signals = "samples_signalsSNT.txt";
    std::vector<MT2Sample> samples_sig = MT2Sample::loadSamples(samplesFile_signals,1000 );

    for( unsigned i=0; i<samples_sig.size(); ++i ) {

      MT2Analysis<MT2EstimateTree>* thisSig = new MT2Analysis<MT2EstimateTree>( samples_sig[i].sname, cfg.regionsSet(), samples_sig[i].id );
      MT2EstimateTree::addVar( thisSig, "prompt" );
      MT2EstimateTree::addVar( thisSig, "iso" );

      computeYield( samples_sig[i], cfg, thisSig );

      signals.push_back( thisSig );

    }  // for signals

  }

  for( unsigned i=0; i<signals.size(); ++i )
    signals[i]->addToFile( mcFile);





  purityTight->writeToFile( outputdir + "/purityMC.root" );
  purityLoose->addToFile( outputdir + "/purityMC.root" );
  eff->addToFile( outputdir + "/purityMC.root" );

  // emulate data:
  roundLikeData(gammaCR);
  roundLikeData(gammaCR_loose);
  gammaCR->writeToFile( outputdir + "/data.root" );
  gammaCR_loose->addToFile( outputdir + "/data.root" );
  gammaCR_nipUp->addToFile( outputdir + "/data.root" );
  gammaCR_nipDown->addToFile( outputdir + "/data.root" );
  tree->addToFile( outputdir + "/data.root" );
 

  return 0;

}









void computeYield( const MT2Sample& sample, const MT2Config& cfg,
                   MT2Analysis<MT2EstimateTree>* anaTree,
                   MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* prompt_pass, 
                   MT2Analysis<MT2EstimateZinvGamma>* nip, MT2Analysis<MT2EstimateZinvGamma>* nip_pass, 
                   MT2Analysis<MT2EstimateZinvGamma>* fake, MT2Analysis<MT2EstimateZinvGamma>* fake_pass, float isoCut ) {

  

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;



  bool isQCD  = sample.id>=100 && sample.id<200;

  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);


  int nentries = tree->GetEntries();


  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( !myTree.passSelection("gamma") ) continue;
    if( !myTree.passGammaAdditionalSelection(sample.id) ) continue;

    //if( myTree.gamma_ht>1000. && sample.id==204 ) continue; // remove high-weight spikes (remove GJet_400to600 leaking into HT>1000)

    if( myTree.gamma_idCutBased[0]==0 ) continue;


    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );


    // absolute iso:
    float iso = myTree.gamma_chHadIso[0];
    if( iso>10. ) continue; // preselection anyways in there


    int mcMatchId = myTree.gamma_mcMatchId[0];
    bool isMatched = (mcMatchId==22 || mcMatchId==7);

    bool isPrompt = isMatched && !isQCD;
    bool isNIP    = isMatched && isQCD;
    bool isFake   = !isMatched;

    float ht        = myTree.gamma_ht;
    float met       = myTree.gamma_met_pt;
    float minMTBmet = myTree.gamma_minMTBMet;
    int njets       = myTree.gamma_nJet30;
    int nbjets      = myTree.gamma_nBJet20;    
    float mt2       = (njets>1) ? myTree.gamma_mt2 : myTree.gamma_jet1_pt;




    Double_t weight = myTree.evt_scale1fb*cfg.lumi(); 

    bool passIso = iso<isoCut;

    if( isPrompt && prompt!=0 && prompt_pass!=0 ) {

      MT2EstimateZinvGamma* thisPrompt = prompt->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      if( thisPrompt==0 ) continue;

      thisPrompt->yield->Fill(mt2, weight );
      thisPrompt->fillIso( iso, weight, mt2 );

      if( passIso ) {

        MT2EstimateZinvGamma* thisPrompt_pass = prompt_pass->get( ht, njets, nbjets, met, minMTBmet, mt2 );
        if( thisPrompt_pass==0 ) continue;

        thisPrompt_pass->yield->Fill(mt2, weight );
        thisPrompt_pass->fillIso( iso, weight, mt2 );

      }

    } else if( isNIP && nip!=0 && nip_pass!=0 ) { 
      
      MT2EstimateZinvGamma* thisnip = nip->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      if( thisnip==0 ) continue;

      thisnip->yield->Fill(mt2, weight );
      thisnip->fillIso( iso, weight, mt2 );

      if( passIso ) {

        MT2EstimateZinvGamma* thisnip_pass = nip_pass->get( ht, njets, nbjets, met, minMTBmet, mt2 );
        if( thisnip_pass==0 ) continue;

        thisnip_pass->yield->Fill(mt2, weight );
        thisnip_pass->fillIso( iso, weight, mt2 );

      }

    } else if( isFake && fake!=0 && fake_pass!=0 ) {

      MT2EstimateZinvGamma* thisFake = fake->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      if( thisFake==0 ) continue;

      thisFake->yield->Fill(mt2, weight );
      thisFake->fillIso( iso, weight, mt2 );

      if( passIso ) {

        MT2EstimateZinvGamma* thisFake_pass = fake_pass->get( ht, njets, nbjets, met, minMTBmet, mt2 );
        if( thisFake_pass==0 ) continue;

        thisFake_pass->yield->Fill(mt2, weight );
        thisFake_pass->fillIso( iso, weight, mt2 );

      }

    } // is prompt/nip/fake


    if( passIso ) {
      MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      if( thisTree!=0 ) {
        thisTree->yield->Fill(mt2, weight );
        if( isPrompt )
          thisTree->assignVar( "prompt", 2 );
        else if( isNIP )
          thisTree->assignVar( "prompt", 1 );
        else if( isFake )
          thisTree->assignVar( "prompt", 0 );
        thisTree->assignVar( "iso", iso );
        thisTree->fillTree_gamma(myTree, weight );
      }
    }
       

    
  } // for entries


  if( prompt!=0 ) prompt->finalize();
  if( prompt_pass!=0 ) prompt_pass->finalize();
  if( nip!=0 ) nip->finalize();
  if( nip_pass!=0 ) nip_pass->finalize();
  if( fake!=0 ) fake->finalize();
  if( fake_pass!=0 ) fake_pass->finalize();
  anaTree->finalize();
  

  delete tree;


  file->Close();
  delete file;
  

}



void roundLikeData( MT2Analysis<MT2EstimateZinvGamma>* data ) {



  std::set<MT2Region> regions = data->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    TH1D* thisYield = data->get(*iR)->yield;

    for( int iBin=1; iBin<thisYield->GetNbinsX()+1; ++iBin ) {

      float yield = thisYield->GetBinContent(iBin);
      thisYield->SetBinContent(iBin, round( yield ));
      thisYield->SetBinError(iBin, 0. );

    } // for bins

  } // for regions

}
