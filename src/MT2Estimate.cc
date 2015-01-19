#include "../interface/MT2Estimate.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>




MT2Estimate::MT2Estimate( const MT2Estimate& rhs ) {

  name = rhs.getName();

  region = new MT2Region(*(rhs.region));

  yield = new TH1D(*(rhs.yield));

}



MT2Estimate::MT2Estimate( const std::string& aname, const MT2Region& aregion ) {

  name = aname;

  region = new MT2Region(aregion);

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield = new TH1D(this->getHistoName("yield").c_str(), "", nBins, bins);
  yield->Sumw2();

}



MT2Estimate::~MT2Estimate() {

  delete region;
  delete yield;

}



std::string MT2Estimate::getHistoName( const std::string& prefix ) const {

  std::string returnName = prefix + "_" + name + "_" + region->getName();
  
  return returnName;

}



void MT2Estimate::setName( const std::string& newName ) {

  name = newName;
  yield->SetName( this->getHistoName("yield").c_str() );

}




void MT2Estimate::addOverflow() {

  MT2Estimate::addOverflowSingleHisto( yield );

}


void MT2Estimate::addOverflowSingleHisto( TH1D* yield ) {

  yield->SetBinContent(yield->GetNbinsX(),
      yield->GetBinContent(yield->GetNbinsX()  )+
      yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
      sqrt(yield->GetBinError(yield->GetNbinsX()  )*
           yield->GetBinError(yield->GetNbinsX()  )+
           yield->GetBinError(yield->GetNbinsX()+1)*
           yield->GetBinError(yield->GetNbinsX()+1)  ));

  yield->SetBinContent(yield->GetNbinsX()+1, 0.);
  yield->SetBinError  (yield->GetNbinsX()+1, 0.);

}



const MT2Estimate& MT2Estimate::operator=( const MT2Estimate& rhs ) {

  if( this->yield == 0 ) { // first time

    this->name = rhs.name;

    this->region = new MT2Region(*(rhs.region));

    this->yield = new TH1D(*(rhs.yield));

  } else { // keep name and histo name, just make histogram identical

    if( this->region!=0 ) delete this->region;
    this->region = new MT2Region(*(rhs.region));

    std::string oldName = this->yield->GetName();
    delete this->yield;
    this->yield = new TH1D(*(rhs.yield));
    this->yield->SetName(oldName.c_str());

  }

  return *this;

}


MT2Estimate MT2Estimate::operator+( const MT2Estimate& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator+] ERROR! Can't add MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  //this->yield->Add(rhs.yield);
  MT2Estimate result(*this);
  result.yield->Add(rhs.yield);

  //return *this;
  return result;

}




MT2Estimate MT2Estimate::operator/( const MT2Estimate& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator/] ERROR! Can't divide MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }


  //this->yield->Divide(rhs.yield);
  MT2Estimate result(name, *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Divide(rhs.yield);

  //return *this;
  return result;

}


MT2Estimate MT2Estimate::operator*( const MT2Estimate& rhs ) const {

  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator*] ERROR! Can't multiply MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  //this->yield->Multiply(rhs.yield);
  MT2Estimate result(name, *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Multiply(rhs.yield);

  //return *this;
  return result;

}



MT2Estimate MT2Estimate::operator/( float k ) const {

  MT2Estimate result(name, *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  return result;

  //this->yield->Scale(1./k);

  //return *this;

}


MT2Estimate MT2Estimate::operator*( float k ) const {

  MT2Estimate result(name, *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);

  return result;


  //this->yield->Scale(k);

  //return *this;

}




MT2Estimate MT2Estimate::operator/=( const MT2Estimate& rhs ) const {

  this->yield->Divide(rhs.yield);
  return (*this);
  //return (*this) / rhs ;

}



MT2Estimate MT2Estimate::operator+=( const MT2Estimate& rhs ) const {

  this->yield->Add(rhs.yield);
  return (*this);
  //return (*this) + rhs ;

}


MT2Estimate MT2Estimate::operator*=( const MT2Estimate& rhs ) const {

  this->yield->Multiply(rhs.yield);
  return (*this);
  //return (*this) * rhs ;

}


MT2Estimate MT2Estimate::operator*=( float k ) const {

  this->yield->Scale(k);
  return (*this);
  //return (*this) * k ;

}


MT2Estimate MT2Estimate::operator/=( float k ) const {

  this->yield->Scale(1./k);
  return (*this);
  //return (*this) / k ;

}






void MT2Estimate::getShit( TFile* file, const std::string& path ) {

  yield = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield->GetName()));

}



void MT2Estimate::print(const std::string& ofs){

  Int_t binXmin=1;
  Int_t binXmax=5;

  Double_t error;
  Double_t integral = yield->IntegralAndError(binXmin, binXmax, error);

  ofstream ofs_file;
  ofs_file.open( ofs, std::ofstream::app );
  if(integral >= 10)
    ofs_file << std::fixed << std::setprecision(1) << " & " << integral << " $\\pm$ " << error;
  else if(integral < 10)
    ofs_file << std::fixed << std::setprecision(2) << " & " << integral << " $\\pm$ " << error;

}
