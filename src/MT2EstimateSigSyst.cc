#include "../interface/MT2EstimateSigSyst.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TEfficiency.h"



MT2EstimateSigSyst::MT2EstimateSigSyst( const std::string& aname, const MT2Region& aregion ) : MT2EstimateSig( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield_systUp = new TH1D( this->getHistoName("yield_systUp").c_str(), "", nBins, bins);
  yield_systUp->Sumw2();
  yield_systDown = new TH1D( this->getHistoName("yield_systDown").c_str(), "", nBins, bins);
  yield_systDown->Sumw2();

}




MT2EstimateSigSyst::MT2EstimateSigSyst( const std::string& aname, const MT2Region& aregion, const MT2EstimateSig& pass, const MT2EstimateSig& tot ) : MT2EstimateSig( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield_systUp = new TH1D( this->getHistoName("yield_systUp").c_str(), "", nBins, bins);
  yield_systUp->Sumw2();
  yield_systDown = new TH1D( this->getHistoName("yield_systDown").c_str(), "", nBins, bins);
  yield_systDown->Sumw2();


  TEfficiency eff( *(pass.yield), *(tot.yield) );

  for( int i=1; i<yield->GetNbinsX()+1; ++i ) {

    yield         ->SetBinContent( i, eff.GetEfficiency(i) );
    yield_systDown->SetBinContent( i, eff.GetEfficiency(i) - eff.GetEfficiencyErrorLow(i) );
    yield_systUp  ->SetBinContent( i, eff.GetEfficiency(i) + eff.GetEfficiencyErrorUp(i) );

  }


}



MT2EstimateSigSyst::MT2EstimateSigSyst( const MT2EstimateSig& rhs ) : MT2EstimateSig(rhs) {

  this->yield_systUp = new TH1D(*(rhs.yield));
  this->yield_systUp   ->SetName(this->getHistoName("yield_systUp").c_str()); 

  this->yield_systDown = new TH1D(*(rhs.yield));
  this->yield_systDown ->SetName(this->getHistoName("yield_systDown").c_str()); 

}




MT2EstimateSigSyst::MT2EstimateSigSyst( const MT2EstimateSigSyst& rhs ) : MT2EstimateSig(rhs) {

  this->yield_systUp = new TH1D(*(rhs.yield_systUp));
  this->yield_systDown = new TH1D(*(rhs.yield_systDown));

}



MT2EstimateSigSyst::~MT2EstimateSigSyst() {

  delete yield_systUp;
  delete yield_systDown;

}







MT2Analysis<MT2EstimateSigSyst>* MT2EstimateSigSyst::makeEfficiencyAnalysis( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateSig>* pass, MT2Analysis<MT2EstimateSig>* all ) {

  std::set<MT2Region> regions = pass->getRegions();

  std::set<MT2EstimateSigSyst*> data;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSig* thisPass = pass->get( *iR );
    MT2EstimateSig* thisAll  = all ->get( *iR );

    MT2EstimateSigSyst* thisEff  = new MT2EstimateSigSyst( aname, *iR, *thisPass, *thisAll );
    data.insert( thisEff );

  } // for regions


  MT2Analysis<MT2EstimateSigSyst>* analysis = new MT2Analysis<MT2EstimateSigSyst>( aname, data );

  return analysis;

}



MT2Analysis<MT2EstimateSigSyst>* MT2EstimateSigSyst::makeAnalysisFromEstimate( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateSig>* estimate ) {

  std::set<MT2Region> regions = estimate->getRegions();

  std::set<MT2EstimateSigSyst*> data;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSig*  thisEstimate = estimate->get( *iR );
    MT2EstimateSigSyst* thisEstimateSyst = new MT2EstimateSigSyst( *thisEstimate );
    data.insert( thisEstimateSyst );

  } // for regions


  MT2Analysis<MT2EstimateSigSyst>* analysis = new MT2Analysis<MT2EstimateSigSyst>( aname, data );

  return analysis;

}



TGraphAsymmErrors* MT2EstimateSigSyst::getGraph() const {

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(0);
  graph->SetName( this->getHistoName("grSyst").c_str() );


  for( int iBin=1; iBin<yield->GetNbinsX()+1; ++iBin ) {

    float x = yield->GetBinCenter(iBin);
    float x_minus = yield->GetBinLowEdge(iBin);
    float x_plus  = yield->GetBinLowEdge(iBin+1);
    float y = yield->GetBinContent(iBin);
    float y_plus  = yield_systUp  ->GetBinContent(iBin);
    float y_minus = yield_systDown->GetBinContent(iBin);
    
    int iPoint = iBin-1;
    graph->SetPoint( iPoint, x, y );
    graph->SetPointError( iPoint, x-x_minus, x_plus-x, y-y_minus, y_plus-y );

  }


  return graph;

}




void MT2EstimateSigSyst::setName( const std::string& newName ) {

  MT2EstimateSig::setName(newName);

  yield_systUp->SetName( this->getHistoName("yield_systUp").c_str() );
  yield_systDown->SetName( this->getHistoName("yield_systDown").c_str() );

}



void MT2EstimateSigSyst::finalize() {

  this->addOverflow();

}


void MT2EstimateSigSyst::addOverflow() {

  MT2EstimateSig::addOverflow();

  MT2EstimateSig::addOverflowSingleHisto( yield_systUp );
  MT2EstimateSig::addOverflowSingleHisto( yield_systDown );

}



void MT2EstimateSigSyst::getShit( TFile* file, const std::string& path ) {

  MT2EstimateSig::getShit(file, path);
  yield_systUp = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield_systUp->GetName()));
  yield_systDown = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield_systDown->GetName()));


}



void MT2EstimateSigSyst::write() const {

  MT2EstimateSig::write();
  yield_systUp->SetLineColor(kRed);
  yield_systDown->SetLineColor(kGreen);
  yield_systUp->Write();
  yield_systDown->Write();

}



void MT2EstimateSigSyst::print(const std::string& ofs){

  Int_t binXmin=1;
  Int_t binXmax=-1;

  Double_t error;
  Double_t integral = yield->Integral(binXmin, binXmax);
  Double_t integral_up = yield_systUp->Integral(binXmin, binXmax);
  Double_t integral_down = yield_systDown->Integral(binXmin, binXmax);
  

  ofstream ofs_file;
  ofs_file.open( ofs, std::ofstream::app );
  if(integral >= 10)
    ofs_file << std::fixed << std::setprecision(1);
  else if(integral < 10)
    ofs_file << std::fixed << std::setprecision(2);
  ofs_file << " & " << integral << " $^{+" << fabs(integral_up-integral) << "}_{-" << fabs(integral_down-integral) << "}$";

}


const MT2EstimateSigSyst& MT2EstimateSigSyst::operator=( const MT2EstimateSigSyst& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield_systUp = new TH1D(*(rhs.yield_systUp));
  this->yield_systDown = new TH1D(*(rhs.yield_systDown));

  this->setName( this->getName() );


  return *this;

}




const MT2EstimateSigSyst& MT2EstimateSigSyst::operator=( const MT2EstimateSig& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield_systUp = new TH1D(*(rhs.yield));
  this->yield_systDown = new TH1D(*(rhs.yield));

  this->setName( this->getName() );


  return *this;

}






MT2EstimateSigSyst MT2EstimateSigSyst::operator+( const MT2EstimateSigSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin+otherBin;
    float newErrUp = sqrt( thisErrUp*thisErrUp + otherErrUp*otherErrUp );
    float newErrDown = sqrt( thisErrDown*thisErrDown + otherErrDown*otherErrDown );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );
    
    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}



MT2EstimateSigSyst MT2EstimateSigSyst::operator-( const MT2EstimateSigSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin-otherBin;
    float newErrUp = sqrt( thisErrUp*thisErrUp + otherErrUp*otherErrUp );
    float newErrDown = sqrt( thisErrDown*thisErrDown + otherErrDown*otherErrDown );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}




MT2EstimateSigSyst MT2EstimateSigSyst::operator/( const MT2EstimateSigSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin/otherBin;
    float newErrUp = sqrt( thisErrUp*thisErrUp/(otherBin*otherBin) + otherErrUp*otherErrUp*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    float newErrDown = sqrt( thisErrDown*thisErrDown/(otherBin*otherBin) + otherErrDown*otherErrDown*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}




MT2EstimateSigSyst MT2EstimateSigSyst::operator*( const MT2EstimateSigSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin*otherBin;
    float newErrUp = sqrt( thisBin*thisBin*otherErrUp*otherErrUp + otherBin*otherBin*thisErrUp*thisErrUp );
    float newErrDown = sqrt( thisBin*thisBin*otherErrDown*otherErrDown + otherBin*otherBin*thisErrDown*thisErrDown );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );
    
  }

  return result;

}



MT2EstimateSigSyst MT2EstimateSigSyst::operator+( const MT2EstimateSig& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin+otherBin;
    float newErrUp   = thisErrUp;
    float newErrDown = thisErrDown;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}


MT2EstimateSigSyst MT2EstimateSigSyst::operator-( const MT2EstimateSig& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin-otherBin;
    float newErrUp   = thisErrUp;
    float newErrDown = thisErrDown;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for ( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}


MT2EstimateSigSyst MT2EstimateSigSyst::operator*( const MT2EstimateSig& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin*otherBin;
    float newErrUp   = thisErrUp*otherBin;
    float newErrDown = thisErrDown*otherBin;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}


MT2EstimateSigSyst MT2EstimateSigSyst::operator/( const MT2EstimateSig& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSigSyst::operator*] ERROR! Can't multiply MT2EstimateSigSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSigSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin/otherBin;
    float newErrUp   = thisErrUp/otherBin;
    float newErrDown = thisErrDown/otherBin;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}





MT2EstimateSigSyst MT2EstimateSigSyst::operator*( float k ) const{

  MT2EstimateSigSyst result(this->getName(), *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);
  
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(k);

  result.yield_systUp = new TH1D(*(this->yield_systUp));
  result.yield_systUp->Scale(k);

  result.yield_systDown = new TH1D(*(this->yield_systDown));
  result.yield_systDown->Scale(k);

  return result;

}



MT2EstimateSigSyst MT2EstimateSigSyst::operator/( float k ) const{

  MT2EstimateSigSyst result(this->getName(), *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(1./k);

  result.yield_systUp = new TH1D(*(this->yield_systUp));
  result.yield_systUp->Scale(1./k);

  result.yield_systDown = new TH1D(*(this->yield_systDown));
  result.yield_systDown->Scale(1./k);

  return result;


}



const MT2EstimateSigSyst& MT2EstimateSigSyst::operator+=( const MT2EstimateSigSyst& rhs ) {

  this->yield->Add(rhs.yield);
  this->yield3d->Add(rhs.yield3d);
  this->yield_systUp->Add(rhs.yield_systUp);
  this->yield_systDown->Add(rhs.yield_systDown);
  return (*this);

}

const MT2EstimateSigSyst& MT2EstimateSigSyst::operator/=( const MT2EstimateSigSyst& rhs ) {


  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin/otherBin;
    float newErrUp = sqrt( thisErrUp*thisErrUp/(otherBin*otherBin) + otherErrUp*otherErrUp*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    float newErrDown = sqrt( thisErrDown*thisErrDown/(otherBin*otherBin) + otherErrDown*otherErrDown*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
        result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return (*this);

}



const MT2EstimateSigSyst& MT2EstimateSigSyst::operator*=( const MT2EstimateSigSyst& rhs ) {


  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin*otherBin;
    float newErrUp = sqrt( thisBin*thisBin*otherErrUp*otherErrUp + otherBin*otherBin*thisErrUp*thisErrUp );
    float newErrDown = sqrt( thisBin*thisBin*otherErrDown*otherErrDown + otherBin*otherBin*thisErrDown*thisErrDown );

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
        result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }


  return (*this);

}



const MT2EstimateSigSyst& MT2EstimateSigSyst::operator*=( const MT2EstimateSig& rhs ) {

  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin = thisBin*otherBin;
    float newErrUp   = thisErrUp*otherBin;
    float newErrDown = thisErrDown*otherBin;

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
        result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }


  return (*this);

}



const MT2EstimateSigSyst& MT2EstimateSigSyst::operator/=( const MT2EstimateSig& rhs ) {

  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin = thisBin/otherBin;
    float newErrUp   = thisErrUp/otherBin;
    float newErrDown = thisErrDown/otherBin;

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
        result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );
    
  }


  return (*this);

}



const MT2EstimateSigSyst& MT2EstimateSigSyst::operator*=( float k ) {

  this->yield->Scale(k);
  this->yield3d->Scale(k);
  this->yield_systUp->Scale(k);
  this->yield_systDown->Scale(k);
  return (*this);

}

const MT2EstimateSigSyst& MT2EstimateSigSyst::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->yield3d->Scale(1./k);
  this->yield_systUp->Scale(1./k);
  this->yield_systDown->Scale(1./k);
  return (*this);

}




// friend functions


MT2EstimateSigSyst operator*( float k, const MT2EstimateSigSyst& rhs ) {

  return rhs*k;

}


MT2EstimateSigSyst operator/( float k, const MT2EstimateSigSyst& rhs ) {

  return rhs/k;

}

