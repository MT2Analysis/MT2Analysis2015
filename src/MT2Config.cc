#include "../interface/MT2Config.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1F.h"



MT2Config::MT2Config( const std::string& name ) {

  name_ = name;

  std::string configFileName = "cfgs/" + name + ".txt";

  std::cout << std::endl;
  std::cout << "-> Reading config file: " << configFileName << std::endl;
  std::cout << std::endl;

  lumi_ = 0.;
  regionsSet_ = "";
  mcSamples_ = "";
  sigSamples_ = "";
  dataSamples_ = "";
  additionalStuff_ = "";

  gammaTemplateRegions_ = "13TeV_inclusive"; // default
  gammaTemplateType_    = "DataRC"; // default


  std::ifstream IN(configFileName.c_str());
  char buffer[200];
  char StringValue[1000];


  while( IN.getline(buffer, 200, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                                                                                                                                                                                 
    }

    std::cout << buffer << std::endl;

    char this_name_c[200];
    sscanf(buffer, "%s %s", this_name_c, StringValue);
    std::string this_name(this_name_c);

    if( this_name=="lumi" )
      lumi_ = atof(StringValue);
    else if( this_name=="regionsSet" )
      regionsSet_ = std::string(StringValue);
    else if( this_name=="mcSamples" )
      mcSamples_ = std::string(StringValue);
    else if( this_name=="sigSamples" )
      sigSamples_ = std::string(StringValue);
    else if( this_name=="dataSamples" )
      dataSamples_ = std::string(StringValue);
    else if( this_name=="additionalStuff" )
      additionalStuff_ = std::string(StringValue);
    else if( this_name=="gammaTemplateRegions" )
      gammaTemplateRegions_ = std::string(StringValue);
    else if( this_name=="gammaTemplateType" )
      gammaTemplateType_ = std::string(StringValue);

  } // while getline

  std::cout << std::endl;


  if( this->useMC() && lumi_ <=0. ) {
    std::cout << "[MT2Config] ERROR!! If you process MC files you need to set a valid lumi value in your cfg!" << std::endl;
    std::cout << "Exiting..." << std::endl;
    exit(761);
  }

  if( gammaTemplateType_=="dataFR" ) gammaTemplateType_="DataFR"; // data Fake Removal
  if( gammaTemplateType_=="dataRC" ) gammaTemplateType_="DataRC"; // data Random Cone
  if( gammaTemplateType_=="data"   ) {
    std::cout << std::endl;
    std::cout << "[MT2Config::gammaTemplateType] Asking for 'data': will use data Random Cone. (default)" << std::endl;
    std::cout << std::endl;
    gammaTemplateType_="DataRC"; // (default for data)
  }

  if( gammaTemplateType_!="data" && gammaTemplateType_!="DataFR" && gammaTemplateType_!="MC" && gammaTemplateType_!="DataRC" ) {
    std::cout << "[MT2Config::gammaTemplateType] ERROR! gammaTemplateType may only be 'MC' or 'DataFR' or 'DataRC'" << std::endl;
    exit(1111);
  }

     
}


bool MT2Config::useMC() const {

  return mcSamples_!="";

}



bool MT2Config::dummyAnalysis() const {

  return dataSamples_=="datatest";

}


std::string MT2Config::getEventYieldDir() const {

  std::string outputdir = "EventYields_" + name_ + "_";
  if( this->dummyAnalysis() ) outputdir += "dummy_";

  double intpart;
  double fracpart = modf(lumi_, &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("%.0ffb", intpart ) );
  outputdir += suffix;

  return outputdir;

}
