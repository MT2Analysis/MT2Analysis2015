#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TList.h"
#include "TObject.h"

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"




void writeToTemplateFile( TFile* file, MT2Analysis<MT2Estimate>* analysis, float err_corr, float err_uncorr );
MT2Analysis<MT2Estimate>* get( const std::string& name, std::vector< MT2Analysis<MT2Estimate>* > analyses, const std::string& name1, const std::string& name2="", const std::string& name3="", const std::string& name4="" );



int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./createDatacards [dir]" << std::endl;
    exit(113);
  } 


  std::string dir( argv[1] );
  std::string fileName = dir + "/analyses.root";

  //std::vector< MT2Analysis<MT2Estimate>* > analyses = MT2Analysis<MT2Estimate>::readAllFromFile( dir + "/analyses.root" );

  //if( analyses.size()==0 ) {
  //  std::cout << "ERROR! No analyses found!" << std::endl;
  //  exit(87);
  //}


  float err_qcd_corr    = 0.0;
  float err_qcd_uncorr  = 1.0;
  float err_llep_corr   = 0.;
  float err_llep_uncorr = 0.075;
  float err_zinv_corr   = 0.;
  float err_zinv_uncorr = 0.;


  MT2Analysis<MT2Estimate>* data  = MT2Analysis<MT2Estimate>::readFromFile( fileName, "data" );
  MT2Analysis<MT2Estimate>* qcd   = MT2Analysis<MT2Estimate>::readFromFile( fileName, "QCD"  );
  qcd->setName("qcd");
  MT2Analysis<MT2Estimate>* zinv  = MT2Analysis<MT2Estimate>::readFromFile( fileName, "ZJets");
  zinv->setName("zinv");
  MT2Analysis<MT2Estimate>* wjets = MT2Analysis<MT2Estimate>::readFromFile( fileName, "WJets");
  MT2Analysis<MT2Estimate>* top   = MT2Analysis<MT2Estimate>::readFromFile( fileName, "Top"  );
  MT2Analysis<MT2Estimate>* llep = new MT2Analysis<MT2Estimate>( *top + *wjets );
  llep->setName( "llep" );
  //*llep += *wjets;
  
  std::cout << "-> Creating BG templates file..." << std::endl;
  TFile* file_templateBG = TFile::Open(Form("%s/bkg_templates.root", dir.c_str()), "recreate");
  writeToTemplateFile( file_templateBG, qcd, err_qcd_corr, err_qcd_uncorr );
  writeToTemplateFile( file_templateBG, llep, err_llep_corr, err_llep_uncorr );
  writeToTemplateFile( file_templateBG, zinv, err_zinv_corr, err_zinv_uncorr );
  file_templateBG->Close();
  std::cout << "-> Created BG templates file: " << file_templateBG->GetName() << std::endl;

  //MT2Analysis<MT2Estimate>* llep = new MT2Analysis<MT2Estimate>( *wjets );
  //*llep += *top;


  //MT2Analysis<MT2Estimate>* data       = get( "data"    , analyses, "data" );
  //MT2Analysis<MT2Estimate>* llep = get( "llep", analyses, "WJets", "Top" );
  //MT2Analysis<MT2Estimate>* qcd        = get( "QCD"     , analyses, "QCD" );
  //MT2Analysis<MT2Estimate>* zinv       = get( "Zinv"    , analyses, "ZJets" );



  std::set<MT2HTRegion> htRegions = data->getHTRegions();
  std::set<MT2SignalRegion> sigRegions = data->getSignalRegions();
  
  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT ) {
    for( std::set<MT2SignalRegion>::iterator iSR=sigRegions.begin(); iSR!=sigRegions.end(); ++iSR ) {

     MT2Region thisRegion(*iHT, *iSR);

     std::string path = dir;
     //std::string path = dir + "/" + iHT->getName() + "/" + iSR->getName();
     std::string histoFileName = path + "/histograms_" + thisRegion.getName() + ".root";
     TFile* histoFile = TFile::Open( histoFileName.c_str() );
     if( histoFile==0 ) {
       std::cout << "ERROR! Didn't find the histogram file for region: " << thisRegion.getName() << std::endl;
       exit(11);
     }
     //TH1D* data = 0;
     //std::vector<TH1D*> processes;
     //TIter next(histoFile->GetListOfKeys());
     //while( TObject *obj = next() ) {
     //  std::string thisName(obj->GetName());
     //  if( thisName == ("yield_data_"+thisRegion.getName()) )
     //    data = (TH1D*)histoFile->Get(thisName.c_str());
     //  else {
     //    TH1D* thisHisto = (TH1D*)histoFile->Get(thisName.c_str());
     //    processes.push_back(thisHisto);
     //  }
     //}

     TH1D* this_data = data->get(thisRegion)->yield;
     TH1D* this_qcd  = qcd ->get(thisRegion)->yield;
     TH1D* this_zinv = zinv->get(thisRegion)->yield;
     TH1D* this_llep = llep->get(thisRegion)->yield;

     std::string datacardName(Form("%s/datacard_%s.txt", path.c_str(), thisRegion.getName().c_str()) );
     ofstream datacard( datacardName.c_str() );


     datacard << "imax 1" << std::endl;
     datacard << "jmax 3" << std::endl;
     datacard << "kmax *" << std::endl;
     datacard << "-------------" << std::endl;
     datacard << std::endl << std::endl;
     
     datacard << "shapes sig "  << thisRegion.getName() << " sig_templates-T2tt-filter/sig_m0-1000_m12-800.root $PROCESS_$CHANNEL_sub $PROCESS_$CHANNEL_sub_$SYSTEMATIC" << std::endl;
     datacard << "shapes qcd "  << thisRegion.getName() << " " << dir << "/bkg_templates.root yield_$PROCESS_$CHANNEL yield_$PROCESS_$CHANNEL" << std::endl;
     datacard << "shapes zinv " << thisRegion.getName() << " " << dir << "/bkg_templates.root yield_$PROCESS_$CHANNEL yield_$PROCESS_$CHANNEL" << std::endl;
     datacard << "shapes llep " << thisRegion.getName() << " " << dir << "/bkg_templates.root yield_$PROCESS_$CHANNEL yield_$PROCESS_$CHANNEL" << std::endl;
     datacard << "-------------" << std::endl;


     datacard << std::endl << std::endl;
     datacard << "bin  " << thisRegion.getName() << std::endl;
     datacard << "observation  " << this_data->Integral(1,this_data->GetXaxis()->GetNbins()) << std::endl;
     datacard << "-------------" << std::endl;
     datacard << std::endl << std::endl;

     // sig qcd zinv llep
     datacard << "bin \t" << thisRegion.getName() << "\t" << thisRegion.getName() << "\t" << thisRegion.getName() << "\t" << thisRegion.getName() << std::endl;
     datacard << "process \t sig \t qcd \t zinv \t llep" << std::endl;
     datacard << "process \t 0 \t 1 \t 2 \t 3" << std::endl;
     datacard << "rate \t XXX \t " << this_qcd->Integral() << " \t " << this_zinv->Integral() << " \t " << this_llep->Integral() << std::endl;
     datacard << "-------------" << std::endl;

     datacard << "syst_sig    lnN \t 1.1 - - -" << std::endl;

     int N_llep = (int)this_llep->Integral();
     float llep_stat_err = (N_llep>0) ? 1./sqrt((float)N_llep) : 1.;
     float llep_tot_err = sqrt( llep_stat_err*llep_stat_err + 0.15*0.15 );
     llep_tot_err+=1.;
     datacard << "syst_ll_corr_" << thisRegion.getName() << "  lnN \t - - - " << llep_tot_err << std::endl;


     if( thisRegion.nBJetsMin()>=2 )
       datacard << "syst_zinv_corr_" << thisRegion.getName() << " lnN \t - - 2.0" << std::endl;
     else {
       datacard << "syst_zinv_corr lnN \t - - 1.2 -" << std::endl;
       if( thisRegion.nBJetsMin()>0 ) {
         datacard << "syst_zinv_Z1b_" << thisRegion.getName() << " lnN \t - - 1.3 -" << std::endl;
       }
     }

     if( this_qcd->Integral()>0. ) {
       //write the bins
     }

      //write bins for zinv

     datacard.close();

     std::cout << "-> Created datacard: " << datacardName << std::endl;

    }
  }

  return 0;

} 



MT2Analysis<MT2Estimate>* get( const std::string& name, std::vector< MT2Analysis<MT2Estimate>* > analyses, const std::string& name1, const std::string& name2, const std::string& name3, const std::string& name4 ) {


  std::cout << "Looking for: " << name << std::endl;
  MT2Analysis<MT2Estimate>* returnAnalysis = new MT2Analysis<MT2Estimate>( name, analyses[0]->getHTRegions(), analyses[0]->getSignalRegions() );

  for( unsigned i=0; i<analyses.size(); ++i ) {

    if( analyses[i]->name == name1 || analyses[i]->name == name2 || analyses[i]->name == name3 || analyses[i]->name == name4 ) {
      std::cout << "  added: " << analyses[i]->name << std::endl;
      (*returnAnalysis) += (*analyses[i]);
    }

  }

  return returnAnalysis;

}



void writeToTemplateFile( TFile* file, MT2Analysis<MT2Estimate>* analysis, float err_corr, float err_uncorr ) {


  file->cd();
  file->mkdir(analysis->name.c_str());
  file->cd(analysis->name.c_str());

  std::set<MT2HTRegion> htRegions = analysis->getHTRegions();
  std::set<MT2SignalRegion> sigRegions = analysis->getSignalRegions();
  
  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT ) {
    for( std::set<MT2SignalRegion>::iterator iSR=sigRegions.begin(); iSR!=sigRegions.end(); ++iSR ) {
  
      MT2Region thisRegion(*iHT, *iSR);

      TH1D* h1 = analysis->get( thisRegion )->yield;

      h1->Write();

      for( unsigned iBin=1; iBin<h1->GetNbinsX()+1; ++iBin ) {

        float binContent = h1->GetBinContent(iBin);

        TH1D* h1_binUp = new TH1D(*h1);
        h1_binUp->SetName(Form("%s_bin_%dUp", h1->GetName(), iBin));
        h1_binUp->SetBinContent( iBin, binContent*( 1. + err_uncorr ) );
        h1_binUp->Write();

        TH1D* h1_binDown = new TH1D(*h1);
        h1_binDown->SetName(Form("%s_bin_%dDown", h1->GetName(), iBin));
        h1_binDown->SetBinContent( iBin, binContent*( 1. - err_uncorr ) );
        h1_binDown->Write();

      } // for bins

    } // for SR
  } // for HT

}

