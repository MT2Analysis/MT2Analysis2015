#ifndef MT2Config_h
#define MT2Config_h

#include <string>

class MT2Config {

 public:

  MT2Config( const std::string& configFileName );

  std::string name() const { return name_; };

  float lumi()      const { return lumi_; };
  std::string regionsSet()      const { return regionsSet_; };
  std::string mcSamples()       const { return mcSamples_; };
  std::string sigSamples()      const { return sigSamples_; };
  std::string dataSamples()     const { return dataSamples_; };
  std::string additionalStuff() const { return additionalStuff_; };

  std::string gammaTemplateType() const { return gammaTemplateType_; };
  std::string gammaTemplateRegions() const { return gammaTemplateRegions_; };

  bool useMC() const;

  bool dummyAnalysis() const;

  std::string getEventYieldDir() const;

  void saveAs( const std::string& filename ) const;

 private:

  std::string name_;

  float lumi_;
  std::string regionsSet_;
  std::string mcSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string additionalStuff_;

  std::string gammaTemplateType_;
  std::string gammaTemplateRegions_;

};



#endif
