#ifndef PixelNtuplizer_h
#define PixelNtuplizer_h

/** \class PixelNtuplizer
 *
 * PixelNtuplizer is the EDProducer subclass which finds seeds
 *
 * \author Oliver Gutsche, Fermilab
 *
 * \version   1st Version Aug. 01, 2005  

 *
 ************************************************************/

#include "DataFormats/SiPixelCluster/interface/SiPixelClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TTree;
class TFile;

class PixelNtuplizer : public edm::EDAnalyzer
{
 public:
  
  explicit PixelNtuplizer(const edm::ParameterSet& conf);
  virtual ~PixelNtuplizer();
  virtual void beginJob(const edm::EventSetup& es);
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& es);

 protected:
  void fillClust (SiPixelClusterCollection::ContainerIterator & matchIt);
  //void fillRecHit( SiPixelRecHitCollection::ContainerIterator & rmatchIt);

  
 private:
  edm::ParameterSet conf_;
  
  //--- Structures for ntupling:
  struct Det 
  {
    int subdet;
    int layer;
    int ladder;
    float z;
    float r;
    float thickness;
    int cols;
    int rows;
    
    void init();
  } det_;

  struct Sim 
  {
    float x;
    float y;
    float xIn;
    float xOut;
    float yIn;
    float yOut;
    void init();
  } sim_;

  struct Clust 
  {
    float x;
    float y;
    float ch;
    int size;
    int sizeX;
    int sizeY;
    int maxPixelCol;
    int maxPixelRow;
    int minPixelCol;
    int minPixelRow;
    unsigned int geoId;
    bool edgeHitX;
    bool edgeHitY;
    void init();
  } clust_;


  struct RecHit 
  {
    float x;
    float y;
    float xErr;
    float yErr;
    void init();
  } recHit_;


  TFile * tfile_;
  TTree * t_;
};


#endif
