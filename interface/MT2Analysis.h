#ifndef MT2Analysis_h
#define MT2Analysis_h


#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "MT2Region.h"
#include "TFile.h"

#include "../interface/mt2.h"




  
template<class T> 
class MT2Analysis {

 public:

  MT2Analysis( const std::string& aname, const std::string& regionsSet="13TeV", int id=-1, const std::string& afullName="" );
  MT2Analysis( const std::string& aname, std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions, int id=-1, const std::string& afullName="" );
  MT2Analysis( const std::string& aname, std::set<MT2Region> regions, int id=-1, const std::string& afullName="" );
  MT2Analysis( const std::string& aname, std::set<T*> data, int id=-1, const std::string& afullName="" );
  MT2Analysis( const MT2Analysis& rhs );
  ~MT2Analysis();


  std::set<MT2Region> multiplyHTandSignal( std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions );
  std::set<MT2Region> getRegions() const { return regions_; };
  std::set<MT2HTRegion> getHTRegions() const; //{ return htRegions_; };
  std::set<MT2SignalRegion> getSignalRegions() const; // { return signalRegions_; };

  MT2Region* getRegion( float ht, int njets, int nbjets, float met=-1., float mt=-1., float mt2=-1. ) const;
  MT2Region* matchRegion( MT2Region region ) const;

  T* get( const MT2Region& r ) const;
  T* get( const MT2Tree& mt2tree ) const;
  T* get( float ht, int njets, int nbjets, float met=-1., float mt=-1., float mt2=-1. ) const;

  std::string getName() const { return name; };
  std::string getFullName() const { return fullName; };

  void setName( const std::string& newName );
  void setFullName( const std::string& newName ) { fullName = newName; };

  const MT2Analysis& operator=( const MT2Analysis& rhs);
  MT2Analysis operator+( const MT2Analysis& rhs) const;
  MT2Analysis operator-( const MT2Analysis& rhs) const;
  void add( const MT2Analysis& rhs );
  MT2Analysis operator/( const MT2Analysis& rhs);
  void divide( const MT2Analysis& rhs );
  MT2Analysis operator*( const MT2Analysis& rhs);
  MT2Analysis operator+=( const MT2Analysis& rhs);
  MT2Analysis operator-=( const MT2Analysis& rhs);
  MT2Analysis operator/=( const MT2Analysis& rhs);
  MT2Analysis operator*=( const MT2Analysis& rhs);

  MT2Analysis operator* ( float k );
  MT2Analysis operator/ ( float k );
  MT2Analysis operator*=( float k );
  MT2Analysis operator/=( float k );

  static MT2Analysis* readFromFile( const std::string& fileName, const std::string& matchName="" );
  static std::vector<MT2Analysis*> readAllFromFile( const std::string& fileName, const std::string& matchName="", bool verbose=true );
  void writeToFile( const std::string& fileName, const std::string& option="RECREATE", bool overwrite=false );
  void addToFile( const std::string& fileName, bool overwrite=false ) {
    return this->writeToFile(fileName,"UPDATE",overwrite);
  }

  static void printFromFile (  const std::string& fileName, const std::string& ofs, const std::string& matchName="" );
  static void print ( const std::vector<MT2Analysis*> analyses, const std::string& ofs, const std::string& matchName="" );
  void print ( const std::string& ofs ) const;

  void printRegions() const;

  void finalize();

  int id;


  std::set<T*> data;



 private:

  void createAnalysisStructure() {

    for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {
      MT2Region thisRegion(*iR);
      T* t = new T(name, thisRegion);
      data.insert(t);
    }

  }

  std::set<MT2Region> regions_;

  std::string name;
  std::string fullName;

};







// constructors

template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, const std::string& regionsSet, int aid, const std::string& afullname ) {


  name = aname;
  fullName = (afullname!="") ? afullname : name;
  id = aid;


  if( regionsSet=="8TeV" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,  750. ));
    htRegions.insert(MT2HTRegion( 750., 1200. ));
    htRegions.insert(MT2HTRegion(1200.,   -1. ));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2, 2, 0, 0));  // 2j0b
    signalRegions.insert(MT2SignalRegion(2, 2, 1, 2));  // 2j1to2b
    signalRegions.insert(MT2SignalRegion(3, 5, 0, 0));  // 3to5j0b
    signalRegions.insert(MT2SignalRegion(3, 5, 1, 1));  // 3to5j1b
    signalRegions.insert(MT2SignalRegion(3, 5, 2, 2));  // 3to5j2b
    signalRegions.insert(MT2SignalRegion(6, -1, 0, 0));  // 6j0b
    signalRegions.insert(MT2SignalRegion(6, -1, 1, 1));  // 6j1b
    signalRegions.insert(MT2SignalRegion(6, -1, 2, 2));  // 6j2b
    signalRegions.insert(MT2SignalRegion(-1, -1, 3, -1));  // 3b

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_inclusive" ) {

    regions_.insert(MT2Region( 450. )); // inclusive 450-inf and no jet requirement

  } else if( regionsSet=="13TeV_inclusive_bjets" ) {

    regions_.insert(MT2Region( 450., -1., 2, -1, 0, 0 )); 
    regions_.insert(MT2Region( 450., -1., 2, -1, 1, 1 )); 
    regions_.insert(MT2Region( 450., -1., 2, -1, 2, 2 )); 
    regions_.insert(MT2Region( 450., -1., 2, -1, 3,-1 )); 

  } else if( regionsSet=="13TeV" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  2)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  2)); 
    signalRegions.insert(MT2SignalRegion(2, -1, 3, -1)); 

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_CSA14" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT" ));
    signalRegions.insert(MT2SignalRegion(4, -1, 2,  2, "loMT" ));
    signalRegions.insert(MT2SignalRegion(2, -1, 3, -1, "loMT" ));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT" ));
    signalRegions.insert(MT2SignalRegion(4, -1, 2,  2, "hiMT" ));
    signalRegions.insert(MT2SignalRegion(2, -1, 3, -1, "hiMT" ));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_PHYS14" ){
    
    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, 8, 0,  0));
    signalRegions.insert(MT2SignalRegion(9, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, 8, 1,  1));
    signalRegions.insert(MT2SignalRegion(9, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "hiMT"));
    
    regions_ = multiplyHTandSignal( htRegions, signalRegions );
    
  } else if( regionsSet=="13TeV_PHYS14_hiHT" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, 8, 0,  0));
    signalRegions.insert(MT2SignalRegion(9, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, 8, 1,  1));
    signalRegions.insert(MT2SignalRegion(9, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "hiMT"));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="13TeV_PHYS14_noMT" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(7, 8, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(9, -1, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(7, 8, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(9, -1, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1 ));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1 ));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="13TeV_PHYS14_hiJHT" ){
    
    regions_.insert(MT2Region(450., 575., 2, 3, 0,  0));
    regions_.insert(MT2Region(575., 1000., 2, 3, 0,  0));
    regions_.insert(MT2Region(1000., -1, 2, 3, 0,  0));

    regions_.insert(MT2Region(450., 575., 4, 6, 0,  0));
    regions_.insert(MT2Region(575., 1000., 4, 6, 0,  0));
    regions_.insert(MT2Region(1000., -1, 4, 6, 0,  0));

    regions_.insert(MT2Region(450., 575., 7, 8, 0,  0));
    regions_.insert(MT2Region(575., 1000., 7, 8, 0,  0));
    regions_.insert(MT2Region(1000., -1, 7, 8, 0,  0));
    
    regions_.insert(MT2Region(450., -1, 9, -1, 0,  0));

    regions_.insert(MT2Region(450., 575., 2, 3, 1,  1));
    regions_.insert(MT2Region(575., 1000., 2, 3, 1,  1));
    regions_.insert(MT2Region(1000., -1, 2, 3, 1,  1));

    regions_.insert(MT2Region(450., 575., 4, 6, 1,  1));
    regions_.insert(MT2Region(575., 1000., 4, 6, 1,  1));
    regions_.insert(MT2Region(1000., -1, 4, 6, 1,  1));

    regions_.insert(MT2Region(450., 575., 7, 8, 1,  1));
    regions_.insert(MT2Region(575., 1000., 7, 8, 1,  1));
    regions_.insert(MT2Region(1000., -1, 7, 8, 1,  1));

    regions_.insert(MT2Region(450., -1, 9, -1, 1,  1));

    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1000., -1, 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1000., -1, 2, 3, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(1000., -1, 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(1000., -1, 4, 6, 2, 2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1000., -1, 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1000., -1, 7, 8, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., -1, 9, -1, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., -1, 9, -1, 2,  2, "hiMT" ));
    
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(1000., -1, 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(1000., -1, 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(1000., -1, 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(1000., -1, 7, -1, 3,  -1, "hiMT"));

  } else if( regionsSet=="13TeV_PHYS14_hiJHT_noMT" ){
    
    regions_.insert(MT2Region(450., 575., 2, 3, 0,  0));
    regions_.insert(MT2Region(575., 1000., 2, 3, 0,  0));
    regions_.insert(MT2Region(1000., -1, 2, 3, 0,  0));

    regions_.insert(MT2Region(450., 575., 4, 6, 0,  0));
    regions_.insert(MT2Region(575., 1000., 4, 6, 0,  0));
    regions_.insert(MT2Region(1000., -1, 4, 6, 0,  0));

    regions_.insert(MT2Region(450., 575., 7, 8, 0,  0));
    regions_.insert(MT2Region(575., 1000., 7, 8, 0,  0));
    regions_.insert(MT2Region(1000., -1, 7, 8, 0,  0));
    
    regions_.insert(MT2Region(450., -1, 9, -1, 0,  0));

    regions_.insert(MT2Region(450., 575., 2, 3, 1,  1));
    regions_.insert(MT2Region(575., 1000., 2, 3, 1,  1));
    regions_.insert(MT2Region(1000., -1, 2, 3, 1,  1));

    regions_.insert(MT2Region(450., 575., 4, 6, 1,  1));
    regions_.insert(MT2Region(575., 1000., 4, 6, 1,  1));
    regions_.insert(MT2Region(1000., -1, 4, 6, 1,  1));

    regions_.insert(MT2Region(450., 575., 7, 8, 1,  1));
    regions_.insert(MT2Region(575., 1000., 7, 8, 1,  1));
    regions_.insert(MT2Region(1000., -1, 7, 8, 1,  1));

    regions_.insert(MT2Region(450., -1, 9, -1, 1,  1));

    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2 ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2 ));
    regions_.insert(MT2Region(1000., -1, 2, 3, 2,  2 ));

    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2 ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2 ));
    regions_.insert(MT2Region(1000., -1, 4, 6, 2, 2 ));

    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2 ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2 ));
    regions_.insert(MT2Region(1000., -1, 7, 8, 2,  2 ));

    regions_.insert(MT2Region(450., -1, 9, -1, 2,  2 ));
        
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1 ));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1 ));

    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1 ));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1 ));

    regions_.insert(MT2Region(1000., -1, 2,  6, 3,  -1 ));
    regions_.insert(MT2Region(1000., -1, 7, -1, 3,  -1 ));

  } else if( regionsSet=="13TeV_CSA14_noMT" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(3, -1, 3, -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_onlyHT" ) {

    regions_.insert(MT2Region( 450.,   575.)); // no cut on jets
    regions_.insert(MT2Region( 575.,  1000.));
    regions_.insert(MT2Region(1000.,    -1.));


  } else if( regionsSet=="13TeV_onlyJets" ) {

    regions_.insert(MT2Region(450., -1., 2,  3, 0,  0)); 
    regions_.insert(MT2Region(450., -1., 4, -1, 0,  0)); 
    regions_.insert(MT2Region(450., -1., 2,  3, 1,  1)); 
    regions_.insert(MT2Region(450., -1., 4, -1, 1,  1)); 
    regions_.insert(MT2Region(450., -1., 2,  3, 2,  2, "loMT" )); 
    regions_.insert(MT2Region(450., -1., 4, -1, 2,  2, "loMT" )); 
    regions_.insert(MT2Region(450., -1., 3, -1, 3, -1, "loMT" )); 
    regions_.insert(MT2Region(450., -1., 2,  3, 2,  2, "hiMT" )); 
    regions_.insert(MT2Region(450., -1., 4, -1, 2,  2, "hiMT" )); 
    regions_.insert(MT2Region(450., -1., 3, -1, 3, -1, "hiMT" )); 

  } else {

    std::cout << "[MT2Analysis::MT2Analysis] Analysis region set '" << regionsSet << "' not implemented yet. Exiting." << std::endl;
    exit(917);

  }

  this->createAnalysisStructure();

}




template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, std::set<MT2Region> regions, int aid, const std::string& afullname ) {

  name = aname;
  fullName = (afullname!="") ? afullname : name;
  id = aid;

  regions_ = regions;

  this->createAnalysisStructure();

}


template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions, int aid, const std::string& afullname ) {

  name = aname;
  fullName = (afullname!="") ? afullname : name;
  id = aid;

  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT )  {
    for( std::set<MT2SignalRegion>::iterator iSR=signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {
      MT2Region thisRegion( *iHT, *iSR );
      regions_.insert(thisRegion);
    }
  }


  this->createAnalysisStructure();

}


template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, std::set<T*> newdata, int aid, const std::string& afullname ) {

  name = aname;
  fullName = (afullname!="") ? afullname : name;
  id = aid;

  for( typename std::set<T*>::iterator idata=newdata.begin(); idata!=newdata.end(); ++idata ) {

    MT2Region* thisRegion = (*idata)->region;
    regions_.insert( *thisRegion );

    this->data.insert( *idata );

  }

}



template<class T> 
MT2Analysis<T>::MT2Analysis( const MT2Analysis& rhs ) {


  regions_ = rhs.getRegions();

  name = rhs.name;
  fullName = rhs.fullName;
  id = rhs.id;

  this->createAnalysisStructure();

  *this = rhs;

}



// destructor

template<class T> 
MT2Analysis<T>::~MT2Analysis() {

  for( typename std::set<T*>::iterator i=data.begin(); i!=data.end(); ++i ) {
    delete *i;
  }

}




// other methods


template<class T>
std::set<MT2Region> MT2Analysis<T>::multiplyHTandSignal( std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions ) {


  std::set<MT2Region> regions;

  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT ) {

    for( std::set<MT2SignalRegion>::iterator iSR=signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

      MT2Region newregion( *iHT, *iSR );
      regions.insert( newregion );

    }

  }

  return regions;

}




template<class T>
std::set<MT2HTRegion> MT2Analysis<T>::getHTRegions() const {

  std::set<MT2HTRegion> regions;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {

    regions.insert( *((*it)->region->htRegion()) );

  }

  return regions;

}



template<class T>
std::set<MT2SignalRegion> MT2Analysis<T>::getSignalRegions() const {

  std::set<MT2SignalRegion> regions;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {

    regions.insert( *((*it)->region->sigRegion()) );

  }

  return regions;

}








template<class T>
void MT2Analysis<T>::printRegions() const {

  std::cout << std::endl;
  std::cout << "-> MT2Analysis '" << name << "' has the following regions: " << std::endl;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) 
    std::cout << "  " << ((*it)->region)->getName() << std::endl;

}



template<class T>
MT2Region* MT2Analysis<T>::getRegion( float ht, int njets, int nbjets, float met, float mt, float mt2 ) const {

  MT2Region* foundRegion = 0;
  
  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {


    float htMin  = (*it)->region->htRegion()->htMin;
    float htMax  = (*it)->region->htRegion()->htMax;
    float metMin = (*it)->region->htRegion()->metMin();

    if( ht<htMin ) continue;
    if( htMax>0. && ht>htMax ) continue;
    if( metMin>0.&& met>0. && met<metMin ) continue;

    int njetsmin  = (*it)->region->sigRegion()->nJetsMin;
    int njetsmax  = (*it)->region->sigRegion()->nJetsMax;
    int nbjetsmin = (*it)->region->sigRegion()->nBJetsMin;
    int nbjetsmax = (*it)->region->sigRegion()->nBJetsMax;

    if( njetsmin >=0 && njets <njetsmin  ) continue;
    if( njetsmax >=0 && njets >njetsmax  ) continue;
    if( nbjetsmin>=0 && nbjets<nbjetsmin ) continue;
    if( nbjetsmax>=0 && nbjets>nbjetsmax ) continue;

    std::string mtCut = (*it)->region->sigRegion()->mtCut;

    if( mtCut == "loMT" ) {

      float mtMin = 0.;
      float mtMax = 200.;
      float mt2Min = 200.;
      float mt2Max = 400.;
      bool insideMT  = ( mt  >= mtMin && mt <= mtMax );
      bool insideMT2 = ( mt2 >= mt2Min && mt2 <= mt2Max );
      bool insideBox = insideMT && insideMT2;
      if( !insideBox ) continue;
 
    } else if( mtCut == "hiMT" ) {

      float mtMin = 0.;
      float mtMax = 200.;
      float mt2Min = 200.;
      float mt2Max = 400.;
      bool insideMT  = ( mt  >= mtMin && mt <= mtMax );
      bool insideMT2 = ( mt2 >= mt2Min && mt2 <= mt2Max );
      bool insideBox = insideMT && insideMT2;
      if( insideBox ) continue; // outside

    } else if( mtCut != "" ) {

      std::cout << "WARNING! Unknown mtCut '" << mtCut << "'! Will ignore." << std::endl;

    }

    foundRegion = (*it)->region;
    break;

  }  // for

  return foundRegion;

}



template<class T>
MT2Region* MT2Analysis<T>::matchRegion( MT2Region region ) const {


  MT2Region* foundRegion = 0;
  
  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {

    if( !( region.isIncluded((*it)->region) ) ) continue;

    foundRegion = (*it)->region;
    break;

  }  // for

  return foundRegion;


}





template<class T>
T* MT2Analysis<T>::get( const MT2Tree& mt2tree ) const {

  return this->get( mt2tree.ht, mt2tree.nJet40, mt2tree.nBJet40, mt2tree.met_pt, mt2tree.minMTBMet, mt2tree.mt2 );

}


template<class T>
T* MT2Analysis<T>::get( float ht, int njets, int nbjets, float met, float mt, float mt2 ) const {

  MT2Region* foundRegion = this->getRegion(ht, njets, nbjets, met, mt, mt2);

  if( foundRegion==0 ) return 0;

  return this->get( *foundRegion );

}


template<class T>
T* MT2Analysis<T>::get( const MT2Region& r ) const {

  bool found = false;

  T* t = 0;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {
    if( *((*it)->region) == r ) {
      t = (*it);
      found = true;
      break;
    }
  }

  if( !found ) return 0;

  return t;

}



template<class T> 
void MT2Analysis<T>::setName( const std::string& newName ) {

  this->name = newName;

  for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t = this->get(thisRegion); 

    t->setName( newName );

  }

}




// operator overloading:

template<class T> 
const MT2Analysis<T>& MT2Analysis<T>::operator=( const MT2Analysis& rhs ) {

  regions_ = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator=] ERROR! Can't equate MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    if( t1==0 ) 
      t1 = new T(*t2);
    else 
      *t1 = *t2;

  }

  return *this;

}



template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator+( const MT2Analysis& rhs ) const {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

      MT2Region thisRegion(*iR);

      T* t1 = this->get(thisRegion); 
      T* t2 = rhs.get(thisRegion); 
      if( t2==0 ) {
        std::cout << "[MT2Analysis::operator+] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
        exit(111);
      }

      T* tnew = new T(*t1 + *t2);
      newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name, newdata);

  return result;

}



template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator-( const MT2Analysis& rhs ) const {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

      MT2Region thisRegion(*iR);

      T* t1 = this->get(thisRegion); 
      T* t2 = rhs.get(thisRegion); 
      if( t2==0 ) {
        std::cout << "[MT2Analysis::operator-] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
        exit(111);
      }

      T* tnew = new T(*t1 - *t2);
      newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name, newdata);

  return result;

}



template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator+=( const MT2Analysis& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator+=] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 += *t2;

  }


  return *this;

}


template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator-=( const MT2Analysis& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator-=] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 -= *t2;

  }


  return *this;

}


template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator/=( const MT2Analysis& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator/= ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 /= *t2;

  }


  return *this;

}



template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator*=( const MT2Analysis& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator*= ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 *= *t2;

  }


  return *this;

}




template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator/=( float k ) {

  std::set<MT2Region> regions = this->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    *t1 /= k;

  }


  return *this;

}



template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator*=( float k ) {


  std::set<MT2Region> regions = this->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    *t1 /= k;

  }


  return *this;

}





//template<class T>
//void MT2Analysis<T>::add( const MT2Analysis& rhs ) {
//
//  (*this) = (*this) + rhs;
//
//}




template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator/( const MT2Analysis& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator/] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    T* tnew = new T(*t1 / *t2);
    newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name, newdata);

  return result;

}





template<class T>
void MT2Analysis<T>::divide( const MT2Analysis& rhs ) {

  (*this) = (*this) / rhs;

}




template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator*( const MT2Analysis& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator*] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    T* tnew = new T(*t1 * (*t2) );
    newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name, newdata);

  return result;

}





template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator*( float k ) {


  std::set<MT2Region> regions = this->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    *t1 *= k;

  }


  return *this;

}




template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator/( float k ) {


  std::set<MT2Region> regions = this->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    *t1 /= k;

  }


  return *this;

}



template<class T> 
void MT2Analysis<T>::writeToFile( const std::string& fileName, const std::string& option, bool overwrite ) {

  TFile* file = TFile::Open(fileName.c_str(), option.c_str() );
  file->cd();

  if( file->GetDirectory(this->name.c_str()) ) {
    file->cd();
    if( overwrite ) {
      file->rmdir(this->name.c_str());
    } else {
      std::cout << "[MT2Analysis::writeToFile] Directory '" << this->name << "' already exists in file '" << fileName << "'. Will not overwrite." << std::endl;
      return;
    }
  }

  file->mkdir(this->name.c_str());
  file->cd(this->name.c_str());

  std::set<MT2Region> regions = this->getRegions();
  for( std::set<MT2Region>::iterator it=regions.begin(); it!=regions.end(); ++it ) {
    file->cd();
    file->mkdir(Form("%s/%s", this->name.c_str(), it->getName().c_str()) );
  }

  
  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {
    file->cd();
    file->cd(Form("%s/%s", this->name.c_str(), (*it)->region->getName().c_str()) );
    (*it)->write();
  }

  file->Close();

  std::cout << "-> Wrote '" << this->name << "' to file: " << fileName << std::endl;

}



template<class T>
void MT2Analysis<T>::printFromFile( const std::string& fileName, const std::string& ofs, const std::string& matchName ) {

  std::vector<MT2Analysis<T>*> analyses = readAllFromFile(fileName, false);

  if( analyses.size()==0 ) {
    std::cout << "[MT2Analysis::printFromFile] WARNING!!! Didn't find any MT2Analysis in file " << fileName << std::endl;
    return 0;
  }
  else
    print( analyses, ofs, matchName );

}



template<class T>
void MT2Analysis<T>::print( std::vector<MT2Analysis<T>*> analyses, const std::string& ofs, const std::string& matchName ) {

  MT2Analysis<T>* analysis = 0;
  if( matchName=="" ) {
    analysis = *(analyses.begin());
  } else {
    for( typename std::vector<MT2Analysis<T>*>::iterator iAn=analyses.begin(); iAn!=analyses.end(); ++iAn ) {
      if( (*iAn)->name == matchName ) {
        analysis = new MT2Analysis<T>(*(*iAn));
        break;
      }
    }
    if( analysis==0 ) {
      std::cout << "[MT2Analysis::print] WARNING!!! Didn't find any MT2Analysis named '" << matchName << std::endl;
      return 0;
    }
  }

  if( analyses.size()>1 && matchName=="" ) {
    std::cout << "[MT2Analysis::print] WARNING!!! Multiple analyses found, but reading only one ('" << analysis->name << "')" << std::endl;
  } else {
    std::cout << "[MT2Analysis::print] Grabbed MT2Analysis '" << analysis->name << std::endl;
  }

  analysis->print( ofs );

}


template<class T>
void MT2Analysis<T>::print( const std::string& ofs ) const {

  ifstream isExist(ofs);
  if(isExist) {
    
    system( Form("rm %s", ofs.c_str()) );
    //std::cout << "File " << ofs.c_str() << " already exist!" << std::endl;
    //exit(1);

  }

  ofstream ofs_file;
  if (ofs_file)
    ofs_file.open( ofs, std::ofstream::app );
  
  std::set<MT2HTRegion> htRegions = this->getHTRegions();
  std::set<MT2SignalRegion> sigRegions = this->getSignalRegions();
  std::set<MT2Region> mt2Regions = this->getRegions();

  std::string oldName="";

  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT ) {
    
    if( iHT->getNiceName() == oldName ) continue;
    
    std::string htRegionName = iHT->getNiceName();
    ofs_file << htRegionName << std::endl;

    for ( std::set<MT2SignalRegion>::iterator iSR=sigRegions.begin(); iSR!=sigRegions.end(); ++iSR ){
      
      std::string sigRegionName = iSR->getNiceName();
      ofs_file << " & "  << sigRegionName;

    }
    
    ofs_file << std::endl;

    for( std::set<MT2Region>::iterator imt2=mt2Regions.begin(); imt2!=mt2Regions.end(); ++imt2 ) {
             	
      if( *(imt2->htRegion()) != (*iHT) )
	continue;

      T* thisT = this->get(*imt2);
      thisT->print( ofs );

    } // for mt2 regions                                                                                                                                                          

    ofs_file << std::endl << std::endl;
    
    oldName = iHT->getNiceName();
 
  } // for ht regions

  std::cout << "-> Printed analysis '" << name << "' to: " << ofs << std::endl;
  
}



template<class T> 
std::vector<MT2Analysis<T>*> MT2Analysis<T>::readAllFromFile( const std::string& fileName, const std::string& matchName, bool verbose ) {

  TFile* file = TFile::Open(fileName.c_str());
 
  if( file==0 ) {
    std::cout << "[MT2Analysis::readAllFromFile] ERROR! Can't open file: " << fileName << std::endl;
    exit(1357);
  }

  if( verbose ) std::cout << "[MT2Analysis] Reading analyses from file: " << file->GetName() << std::endl;

  std::vector<MT2Analysis<T>*> analyses;


  TIter next(file->GetListOfKeys());

  // these are the uppermost dirs in the file, 
  // so one dir per analysis name (tyipically only one anyways)
  while(TObject *obj = next()) { 

    TString thisdirname( obj->GetName() );
    if( thisdirname.BeginsWith("ProcessID") ) continue;

    std::string analysisName(obj->GetName());
    file->cd(analysisName.c_str()); 

    std::set<MT2Region> regions;

    TIter next2(gDirectory->GetListOfKeys());

    // these are the directiories inside the analysis dir
    // there will be one directory per region
    // and the dir name is the region identifying name
    while(TObject *obj2 = next2()) { 

      std::string regionName(obj2->GetName());
      MT2Region region(regionName);
      regions.insert(region);

      
      file->cd(analysisName.c_str());

    } // while regions


    //TString analysisName_tstr(analysisName);
    //if( matchExpression!="" && !(analysisName_tstr.Contains(matchExpression)) ) continue;
    if( matchName!=analysisName ) continue;

    // now that we know name and region structure we can istantiate an MT2Analysis:
    MT2Analysis<T>* analysis = new MT2Analysis<T>( analysisName, regions );

    // second loop to set everything
    file->cd(analysisName.c_str()); 
    TIter nextAgain(gDirectory->GetListOfKeys());

    while(TObject *obj2 = nextAgain()) { // loop on regions

      std::string regionName(obj2->GetName());
      MT2Region region(regionName);
      
      std::string path = analysisName + "/" + regionName;
      file->cd(path.c_str());

      T* thisT = analysis->get(region);
      thisT->getShit( file, path );

    } // while regions

    analyses.push_back( analysis );

    if( verbose ) std::cout << "  -> added: " << analysis->name << std::endl;

  } // while analysis names


  return analyses;
  
}



template<class T> 
MT2Analysis<T>* MT2Analysis<T>::readFromFile( const std::string& fileName, const std::string& matchName ) {


  std::vector<MT2Analysis<T>*> analyses = readAllFromFile(fileName, matchName, false);

  if( analyses.size()==0 ) {
    if( matchName!="" )
      std::cout << "[MT2Analysis::readFromFile] WARNING!!! Didn't find any MT2Analysis matching '" << matchName << "' in file " << fileName << std::endl;
    else
      std::cout << "[MT2Analysis::readFromFile] WARNING!!! Didn't find any MT2Analysis in file " << fileName << std::endl;
    return 0;
  }


  MT2Analysis<T>* analysis = *(analyses.begin());

  if( analyses.size()>1 && matchName=="" ) {
    std::cout << "[MT2Analysis::readFromFile] WARNING!!! Multiple analyses found in file: " << fileName << std::endl;
    std::cout << "[MT2Analysis::readFromFile] but reading only one ('" << analysis->name << "')" << std::endl;
    std::cout << "[MT2Analysis::readFromFile] (if you want to read all of them you should use readAllFromFile)" << std::endl;
  } else {
    std::cout << "[MT2Analysis::readFromFile] Grabbed MT2Analysis '" << analysis->name << "' from file " << fileName << std::endl;
  }

  return analysis;

}




template<class T> 
void MT2Analysis<T>::finalize() {

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) 
    (*it)->finalize();


}


#endif
