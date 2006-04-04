// File: PixelNtuplizer.cc
// Description:  see PixelNtuplizer.h
// Authors:  Petar Maksimovic, Jason Shaev, JHU
// History: 4/4/06   creation
//--------------------------------------------

#include "RecoLocalTracker/SiPixelRecHits/test/PixelNtuplizer.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// Old
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"



#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


// For ROOT
#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

// STD
#include <memory>
#include <string>
#include <iostream>

using namespace std;


PixelNtuplizer::PixelNtuplizer(edm::ParameterSet const& conf) : 
  conf_(conf), tfile_(0), t_(0)
{
}


// Virtual destructor needed.
PixelNtuplizer::~PixelNtuplizer() { }  

// End job: write and close the ntuple file
void PixelNtuplizer::endJob() 
{
  std::cout << " PixelNtuplizer::endJob" << std::endl;
  tfile_->Write();
  tfile_->Close();
}



void PixelNtuplizer::beginJob(const edm::EventSetup& es)
{
  std::cout << " PixelNtuplizer::beginJob" << std::endl;

  // put here whatever you want to do at the beginning of the job
  tfile_ = new TFile ( "pix_ntuple.root", "RECREATE" );


  t_ = new TTree("PixNtuple","Pixel hit analyzer ntuple");
  int bufsize = 64000;

  // Create one branch. If splitlevel is set, event is a superbranch
  // creating a sub branch for each data member of the Event object.
  
  //tree->Branch("event", "Event", &event, bsize,split);
  
  std::cout << "Making det branch:" << std::endl;
  t_->Branch("det",   &det_,   "subdet/I:layer/I:ladder/I:z/F:r/F:thickness/F:cols/I:rows:I", bufsize);

  std::cout << "Making sim hit branch:" << std::endl;
  t_->Branch("sim",   &sim_,   "x/F:y/F:xIn/F:xOut/F:yIn/F:yOut/F", bufsize);

  std::cout << "Making cluster branch:" << std::endl;
  t_->Branch("clust", &clust_, "x/F:y/F:ch/F:size/I:sizeX/I:sizeY/I:maxPixelCol/I:maxPixelRow/I:minPixelCol/I:minPixelRow/I:geoId/I", bufsize);

  std::cout << "Making recHit branch:" << std::endl;
  t_->Branch("recHit", &recHit_, "x/F:y/F:xErr/F:yErr/F", bufsize);

  std::cout << "Made all branches." << std::endl;
}



// Functions that gets called by framework every event
void PixelNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  std::string rechitProducer = conf_.getParameter<std::string>("RecHitProducer");

  
  // Get event setup (to get global transformation)
  edm::ESHandle<TrackerGeometry> geom;
  es.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);

  
  // PSim hits
  std::vector<PSimHit> pixelPSimHits_;
  pixelPSimHits_.clear();      // once pixelPSimHits_ gets moved into a class variable

  edm::Handle<edm::PSimHitContainer> PixelBarrelHitsLowTof;
  edm::Handle<edm::PSimHitContainer> PixelBarrelHitsHighTof;
  edm::Handle<edm::PSimHitContainer> PixelEndcapHitsLowTof;
  edm::Handle<edm::PSimHitContainer> PixelEndcapHitsHighTof;

  e.getByLabel("SimG4Object","TrackerHitsPixelBarrelLowTof", PixelBarrelHitsLowTof);
  e.getByLabel("SimG4Object","TrackerHitsPixelBarrelHighTof", PixelBarrelHitsHighTof);
  e.getByLabel("SimG4Object","TrackerHitsPixelEndcapLowTof", PixelEndcapHitsLowTof);
  e.getByLabel("SimG4Object","TrackerHitsPixelEndcapHighTof", PixelEndcapHitsHighTof);

  pixelPSimHits_.insert(pixelPSimHits_.end(), PixelBarrelHitsLowTof->begin(), PixelBarrelHitsLowTof->end()); 
  pixelPSimHits_.insert(pixelPSimHits_.end(), PixelBarrelHitsHighTof->begin(), PixelBarrelHitsHighTof->end());
  pixelPSimHits_.insert(pixelPSimHits_.end(), PixelEndcapHitsLowTof->begin(), PixelEndcapHitsLowTof->end()); 
  pixelPSimHits_.insert(pixelPSimHits_.end(), PixelEndcapHitsHighTof->begin(), PixelEndcapHitsHighTof->end());
  
  std::cout 
    <<" FOUND " << pixelPSimHits_.size() 
    <<" Pixel PSim Hits"<<std::endl;


  //--- Fetch Pixel Clusters
  edm::Handle<SiPixelClusterCollection> clustColl;
  e.getByType(clustColl);

  std::cout 
    <<" FOUND " 
    << const_cast<SiPixelClusterCollection*>(clustColl.product())->size()
    << " Pixel Clusters" << std::endl;


  //--- Fetch Pixel RecHits
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByType(recHitColl);
  
  std::cout 
    <<" FOUND " 
    << const_cast<SiPixelRecHitCollection*>(recHitColl.product())->size()
    << " Pixel RecHits" << std::endl;



  //--- Iterate over PSim hit
  bool PRINT = true;

  for (std::vector<PSimHit>::iterator isim = pixelPSimHits_.begin(),
	 isimEnd = pixelPSimHits_.end(); 
       isim != isimEnd; ++isim){
    //DetId detid((*isim).detUnitId());
    DetId detIdObj((*isim).detUnitId());
    unsigned int subid=detIdObj.subdetId();
    if (! ((subid==PixelSubdetector::PixelBarrel) || 
	   (subid== PixelSubdetector::PixelEndcap))) {
      std::cout << "Huh?  Not a pixel PSimHit.  Weird... Skipped it." 
		<< std::endl;
      continue;
    }

    //--- Since TTree::Fill() will simply take a snapshot of what's 
    //--- in the memory, we need to re-initialize the cache used for 
    //--- each of the branches.
    det_.init();
    sim_.init();
    clust_.init();
    recHit_.init();


    //--- Get the information for this PSimHit.
    sim_.xIn  = (*isim).entryPoint().x(); // width (row index, in col direction)
    sim_.yIn  = (*isim).entryPoint().y(); // length (col index, in row direction) 
    sim_.xOut = (*isim).exitPoint().x();
    sim_.yOut = (*isim).exitPoint().y();
    sim_.x    = 0.5*(sim_.xIn + sim_.xOut);
    sim_.y    = 0.5*(sim_.yIn + sim_.yOut);

    det_.subdet = detIdObj.subdetId();
    det_.layer = 9999;
    det_.ladder = 9999;


    // Get the geom-detector
    const PixelGeomDetUnit * theGeomDet =
      dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detIdObj) );
    det_.z = theGeomDet->surface().position().z();
    det_.r = theGeomDet->surface().position().perp();

    const BoundPlane& plane = theGeomDet->surface(); //for transf.
    
    det_.thickness = theGeomDet->specificSurface().bounds().thickness();
    det_.cols = theGeomDet->specificTopology().ncolumns();
    det_.rows = theGeomDet->specificTopology().nrows();


    // Get topology
    const RectangularPixelTopology * topol = 
      dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));
    

    //------------------------------------------------------------------------
    //--- Find the matching PixelCluster by iterating over the Clusters for
    //--- this DetUnit and finding the closest one.
    //------------------------------------------------------------------------

    float dr_min=99999., dx=99999., dy=99999.;
    const SiPixelClusterCollection::Range clustRange = 
      clustColl->get(detIdObj());

    SiPixelClusterCollection::ContainerIterator clustIt;
    SiPixelClusterCollection::ContainerIterator matchIt = clustRange.second; //end 

    for (clustIt = clustRange.first; clustIt != clustRange.second; 
	 ++clustIt ) {
      //
      float x = clustIt->x();
      float y = clustIt->y();
      //
      LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));
      if(PRINT)cout<<lp.x()<<" "<<lp.y()<<endl;
      //
      float dr = 
	(sim_.x-lp.x())*(sim_.x-lp.x()) + (sim_.y-lp.y())*(sim_.y-lp.y());
      if (dr < dr_min){
	// &&& Instead, let's store the iterator and then fill in all the
	// &&& vars for that hit
        matchIt = clustIt;
	dr_min=dr;
	clust_.x = lp.x();
	clust_.y = lp.y();
      }  
      if (PRINT) cout<<"simhit "<<sim_.x<<" "<<sim_.y<<" "<<dr<<endl;
    }
    

    //--- If there was a match, capture the cluster quantities
    if (dr_min < 99999.) { 
      fillClust(matchIt);           // copy relevant stuff from matchIt into clust_
    }
    else {
      clust_.init();            // wipe it all out
    }
      
      

    //------------------------------------------------------------------------
    //--- Find the matching PixelRecHit by iterating over the RecHits for
    //--- this DetUnit and finding the closest one.
    //------------------------------------------------------------------------

#if 0
    const SiPixelRecHitCollection::Range recHitRange = 
      recHitColl->get(detIdObj());
    

    SiPixelRecHitCollection::ContainerIterator recHitIt;
    SiPixelRecHitCollection::ContainerIterator rmatchIt = recHitRange.second; //end 

    float dr0=99999., dx=99999., dy=99999.;
    for (recHitIt = recHitRange.first; recHitIt != recHitRange.second; 
	 ++recHitIt ) {
      //
      float x = recHitIt->x();
      float y = recHitIt->y();
      //
      LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));
      if(PRINT)cout<<lp.x()<<" "<<lp.y()<<endl;
      //
      float dr = 
	(sim_.x-lp.x())*(sim_.x-lp.x()) + (sim_.y-lp.y())*(sim_.y-lp.y());
      if (dr < dr_min){
	// &&& Instead, let's store the iterator and then fill in all the
	// &&& vars for that hit
        rmatchIt = recHitIt;
	dr_min=dr;
	recHit_.x = lp.x();
	recHit_.y = lp.y();
      }  
      if (PRINT) cout<<"simhit "<<sim_.x<<" "<<sim_.y<<" "<<dr<<endl;
    }
    

    //--- If there was a match, capture the recHiter quantities
    if (dr_min < 99999.) { 
      fillRecHit(rmatchIt);           // copy relevant stuff from rmatchIt into recHit_
    }
    else {
      recHit_.init();            // wipe it all out
      //if(PRINT) cout<<" no match for this hit "<<endl;
    }

#endif



    //++++++++++
    t_->Fill();
    //++++++++++

  }
}




void 
PixelNtuplizer::fillClust(SiPixelClusterCollection::ContainerIterator & matchIt)
{

  //       if(PRINT) cout<<"clus "<<ch<<" "<<size<<" "<<sizeX<<" "<<sizeY<<" "
  // 		    <<x<<" "<<y<<" "<<geoId<<" "<<edgeHitX<<" "
  // 		    <<edgeHitY<<endl;
  //       if(PRINT) cout<<"match "<<dr0<<" "<<dx<<" "<<dy<<endl;
  
  //const vector<Pixel>  = clustIt->pixels();
  clust_.ch = (matchIt->charge())/1000.; // convert ke to electrons
  clust_.size = matchIt->size();
  clust_.sizeX = matchIt->sizeX();
  clust_.sizeY = matchIt->sizeY();
  clust_.x = matchIt->x();
  clust_.y = matchIt->y();
  clust_.maxPixelCol = matchIt->maxPixelCol();
  clust_.maxPixelRow = matchIt->maxPixelRow();
  clust_.minPixelCol = matchIt->minPixelCol();
  clust_.minPixelRow = matchIt->minPixelRow();
  
  clust_.geoId = matchIt->geographicalId();
  clust_.edgeHitX = matchIt->edgeHitX();
  clust_.edgeHitY = matchIt->edgeHitY();
}



void PixelNtuplizer::Det::init()
{
  static float dummy_float = 9999.0;
  static int dummy_int = 9999;
  static unsigned int dummy_uint = 9999;
  static  bool dummy_bool = false;
  
  subdet = dummy_int;
  layer  = dummy_int;
  ladder  = dummy_int;
  z = dummy_float;
  r = dummy_float;
  thickness = dummy_float;
  cols  = dummy_int;
  rows  = dummy_int;
}


void PixelNtuplizer::Sim::init()
{
  static float dummy_float = 9999.0;
  static int dummy_int = 9999;
  static unsigned int dummy_uint = 9999;
  static  bool dummy_bool = false;
  
  x = dummy_float;
  y = dummy_float;
  xIn = dummy_float;
  xOut = dummy_float;
  yIn = dummy_float;
  yOut = dummy_float;
}



void PixelNtuplizer::Clust::init()
{
  static float dummy_float = 9999.0;
  static int dummy_int = 9999;
  static unsigned int dummy_uint = 9999;
  static  bool dummy_bool = false;
  
  ch = dummy_float; // convert ke to electrons
  size = dummy_int;
  sizeX = dummy_int;
  sizeY = dummy_int;
  x = dummy_float;
  y = dummy_float;
  maxPixelCol = dummy_int;
  maxPixelRow = dummy_int;
  minPixelCol = dummy_int;
  minPixelRow = dummy_int;
  
  geoId = dummy_uint;
  edgeHitX = dummy_bool;
  edgeHitY = dummy_bool;
}




void PixelNtuplizer::RecHit::init()
{
  static float dummy_float = 9999.0;
  static int dummy_int = 9999;
  static unsigned int dummy_uint = 9999;
  static  bool dummy_bool = false;

  x = dummy_float;
  y = dummy_float;
  xErr = dummy_float;
  yErr = dummy_float; 
}


