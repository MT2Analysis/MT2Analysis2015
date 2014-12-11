#ifndef MT2Region_h
#define MT2Region_h

#include <string>
#include <vector>




class MT2HTRegion {

 public:

  MT2HTRegion( const std::string& name );
  MT2HTRegion( const MT2HTRegion& rhs );
  MT2HTRegion( float ahtMin, float ahtMax, float ametMin, const std::string& aHLT_selection="" );
  ~MT2HTRegion() {};


  // univocal identifier:
  std::string getName() const;

  std::vector< std::string > getNiceNames() const;


  float htMin;
  float htMax;
  float metMin;

  std::string HLT_selection;


  bool operator==( const MT2HTRegion& rhs ) const;
  bool operator!=( const MT2HTRegion& rhs ) const;
  bool operator<( const MT2HTRegion& rhs ) const;


 private:

};





class MT2SignalRegion {

 public:

  MT2SignalRegion( const std::string& name );
  MT2SignalRegion( int njmin, int njmax, int nbmin, int nbmax, float mtMaxCut=-1., float mt2MinCut=-1., float mt2MaxCut=-1., bool insideBox=true );
  MT2SignalRegion( const MT2SignalRegion& rhs );

  ~MT2SignalRegion() {};

  // univocal identifier:
  std::string getName() const;

  std::string getNameMt() const;
  std::string getNiceName() const;
  
  int nJetsMin; 
  int nJetsMax;
  int nBJetsMin;
  int nBJetsMax;

  float mtMax;
  float mt2Min;
  float mt2Max;

  bool inBox; 

  bool operator==( const MT2SignalRegion& rhs ) const;
  bool operator!=( const MT2SignalRegion& rhs ) const;
  bool operator<( const MT2SignalRegion& rhs ) const;


 private:

  std::string getNiceJetName( const std::string& pedix, int nmin, int nmax ) const;
  std::string getSingleJetString( const std::string& suffix, int n_min , int n_max=-1 ) const;

};





class MT2Region {

 public:

  MT2Region( const MT2Region& region ) {
    htRegion_ = new MT2HTRegion(*(region.htRegion()));
    sigRegion_ = new MT2SignalRegion(*(region.sigRegion()));
  }

  MT2Region( const MT2HTRegion& htRegion, const MT2SignalRegion& sigRegion ) {
    htRegion_ = new MT2HTRegion(htRegion);
    sigRegion_ = new MT2SignalRegion(sigRegion);
  }
  ~MT2Region() {};


  // this is a univocal identifier (same jet/HT thresholds always correspond to same name)
  std::string getName() const {
    return htRegion_->getName() + "_" + sigRegion_->getName();
  }


  std::vector< std::string > getNiceNames() const;

  void getBins( int& nBins, double*& bins ) const;

  MT2HTRegion* htRegion() const {
    return htRegion_;
  }

  MT2SignalRegion* sigRegion() const {
    return sigRegion_;
  }

  float htMin() { return htRegion_->htMin; };
  float htMax() { return htRegion_->htMax; };
  float metMin() { return htRegion_->metMin; };
  int nJetsMin() { return sigRegion_->nJetsMin; };
  int nJetsMax() { return sigRegion_->nJetsMax; };
  int nBJetsMin() { return sigRegion_->nBJetsMin; };
  int nBJetsMax() { return sigRegion_->nBJetsMax; };


  bool operator==( const MT2Region& rhs ) const {
    return ( *htRegion_==*(rhs.htRegion()) && *sigRegion_==*(rhs.sigRegion()) );
  }

  bool operator!=( const MT2Region& rhs ) const {
    return ( *htRegion_!=*(rhs.htRegion()) || *sigRegion_!=*(rhs.sigRegion()) );
  }

  bool operator<( const MT2Region& rhs ) const {
    if( *htRegion_!=*(rhs.htRegion()) ) {
      return *htRegion_<*(rhs.htRegion());
    } else {
      return *sigRegion_<*(rhs.sigRegion());
    }
    return false;
  }


 private:
  
  MT2HTRegion* htRegion_;
  MT2SignalRegion* sigRegion_;

};




#endif
