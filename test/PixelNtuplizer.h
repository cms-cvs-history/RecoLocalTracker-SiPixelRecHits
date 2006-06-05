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

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

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
  void fillClust (SiPixelClusterCollection::ContainerIterator & );
  void fillRecHit( LocalPoint, LocalError );
  void fillDet(DetId &, int );
  void fillSim(std::vector<PSimHit>::iterator); 
 
 private:
  edm::ParameterSet conf_;
  
  //--- Structures for ntupling:
  struct Det 
  {
    float thickness;
    int cols;
    int rows;
    int layer;
    int ladder;
    int module;
    int disk;
    int blade;
    int panel;
    int plaquette;

    void init();
  } det_;

  struct vertex
  {
    int num;
    float r;
    float z;

    void init();
  } vertex_;

  struct track
  {
    float eta;
    float phi;

    void init();
  } track_;

  struct sim 
  {
    float x;
    float y;
    float px;
    float py;
    float pz;
    float eloss;
    float phi;
    float theta;
    int subdetid;
    int isflipped;
    // alpha and beta are related to phi and theta, but
    // more standard:
    float alpha;   
    float beta;
    int PID;
    unsigned int TID;



    void init();
  } sim_;

  struct clust 
  {
    float x;
    float y;
    float charge;
    int size;
    int size_x;
    int size_y;
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
    float xx;
    float xy;
    float yy;

    void init();
  } recHit_;

  TFile * tfile_;
  TTree * t_;
};


#endif
