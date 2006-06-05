// File: PixelNtuplizer.cc
// Description:  see PixelNtuplizer.h
// Authors:  Petar Maksimovic, Jason Shaev, JHU...Vincenzo Chiochia, CERN
// History: 4/4/06   creation
//--------------------------------------------

#include "RecoLocalTracker/SiPixelRecHits/test/PixelNtuplizer.h"


// DataFormats
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelClusterCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

// Old
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


// SimDataFormats
//#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/EmbdSimTrack.h"
#include "SimDataFormats/Vertex/interface/EmbdSimVertex.h"
//
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
  
   t_->Branch("det", &det_, "thickness/F:cols/I:rows/I:layer/I:ladder/I:module/I:disk/I:blade/I:panel/I:plaquette/I", bufsize);

    std::cout << "Making vertex branch:" << std::endl;
  t_->Branch("vertex",   &vertex_,   "num/I:r/F:z/F", bufsize);

  std::cout << "Making track branch:" << std::endl;
   t_->Branch("track", &track_, "eta/F:phi/F", bufsize);

  std::cout << "Making sim hit branch:" << std::endl;
  t_->Branch("sim",   &sim_,   "x/F:y/F:px/F:py/F:pz/F:eloss/F:phi/F:theta/F:subdetid/I:isflipped/I:alpha/F:beta/F:PID/I:TID/i", bufsize);

  std::cout << "Making cluster branch:" << std::endl;
  t_->Branch("clust", &clust_, "x/F:y/F:charge/F:size/I:size_x/I:size_y/I:maxPixelCol/I:maxPixelRow/I:minPixelCol/I:minPixelRow/I:geoId/i:edgeHitX/O:edgeHitY/O", bufsize);

  std::cout << "Making recHit branch:" << std::endl;
  t_->Branch("recHit", &recHit_, "x/F:y/F:xx/F:xy/F:yy/F", bufsize);

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
  std::string clusterCollLabel = conf_.getUntrackedParameter<std::string>("ClusterCollLabel","pixClust"); 
  edm::Handle< edm::DetSetVector<SiPixelCluster> > clustColl;
  e.getByLabel(clusterCollLabel, clustColl);
  edm::DetSetVector<SiPixelCluster>::const_iterator DSViter = clustColl->begin();

/*
  std::cout 
    <<" FOUND " 
    << const_cast<SiPixelClusterCollection*>(clustColl.product())->size()
    << " Pixel Clusters" << std::endl;
*/

  //--- Fetch Pixel RecHits
  std::string recHitCollLabel = conf_.getUntrackedParameter<std::string>("RecHitCollLabel","pixRecHitConverter");
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel(recHitCollLabel, recHitColl);
  
  std::cout 
    <<" FOUND " 
    << const_cast<SiPixelRecHitCollection*>(recHitColl.product())->size()
    << " Pixel RecHits" << std::endl;



  //--- Iterate over PSim hit
  bool PRINT = true;

  std::vector<PSimHit>::iterator isim;
  std::vector<PSimHit>::iterator isimEnd;

  for (isim = pixelPSimHits_.begin(),
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
    vertex_.init();
    track_.init();
    sim_.init();
    clust_.init();
    recHit_.init();

    fillDet(detIdObj, subid);

    fillSim(isim);
    sim_.subdetid = subid;


    // Get the geom-detector
    const PixelGeomDetUnit * theGeomDet =
      dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detIdObj) );
    vertex_.z = theGeomDet->surface().position().z();
    vertex_.r = theGeomDet->surface().position().perp();

    const BoundPlane& plane = theGeomDet->surface(); //for transf.
    
    det_.thickness = theGeomDet->specificSurface().bounds().thickness();
    det_.cols = theGeomDet->specificTopology().ncolumns();
    det_.rows = theGeomDet->specificTopology().nrows();

 // Is flipped ?
    float tmp1 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
    float tmp2 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
    //cout << " 1: " << tmp1 << " 2: " << tmp2 << endl;
    if ( tmp2<tmp1 ) sim_.isflipped = 1;
    else sim_.isflipped = 0;
 
    // Get topology
    const RectangularPixelTopology * topol = 
      dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));
    

    //------------------------------------------------------------------------
    //--- Find the matching PixelCluster by iterating over the Clusters for
    //--- this DetUnit and finding the closest one.
    //------------------------------------------------------------------------

    float dr_min=99999., dx=99999., dy=99999.;
    //const SiPixelClusterCollection::Range clustRange = 
     // clustColl->get(detIdObj());

    edm::DetSet<SiPixelCluster>::const_iterator clustIt;
    edm::DetSet<SiPixelCluster>::const_iterator matchIt = DSViter->data.end(); //end 

    for (clustIt = DSViter->data.begin(); clustIt != DSViter->data.end(); 
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


    float dr_rh0=99999., dx_rh=99999., dy_rh=99999.;

    SiPixelRecHitCollection::range rechitRange = recHitColl->get(detIdObj);
    SiPixelRecHitCollection::const_iterator rechitIt; 
    for (rechitIt = rechitRange.first; rechitIt != rechitRange.second; 
         ++rechitIt ) {
      LocalPoint lp = rechitIt->localPosition();
      LocalError le = rechitIt->localPositionError();
      float dr_rh = 
        (sim_.x-lp.x())*(sim_.x-lp.x()) + (sim_.y-lp.y())*(sim_.y-lp.y());

      if(dr_rh<dr_rh0){
        dr_rh0=dr_rh;
	fillRecHit(lp,le);
      }
      if(PRINT)cout<<"---> rechit found: x " << lp.x() <<"  y "<<lp.y() << " distance " << dr_rh0 << endl;
    }

    if(dr_rh0==99999.) {
      recHit_.x = 9999.;
      recHit_.y = 9999.;
    }
    



//     const SiPixelRecHitCollection::Range recHitRange = 
//       recHitColl->get(detIdObj());
    
//     SiPixelRecHitCollection::ContainerIterator recHitIt;
//     SiPixelRecHitCollection::ContainerIterator rmatchIt = recHitRange.second; //end 

//     float dr0=99999., dx=99999., dy=99999.;
//     for (recHitIt = recHitRange.first; recHitIt != recHitRange.second; 
// 	 ++recHitIt ) {
//       //
//       float x = recHitIt->x();
//       float y = recHitIt->y();
//       //
//       LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));
//       if(PRINT)cout<<lp.x()<<" "<<lp.y()<<endl;
//       //
//       float dr = 
// 	(sim_.x-lp.x())*(sim_.x-lp.x()) + (sim_.y-lp.y())*(sim_.y-lp.y());
//       if (dr < dr_min){
// 	// &&& Instead, let's store the iterator and then fill in all the
// 	// &&& vars for that hit
//         rmatchIt = recHitIt;
// 	dr_min=dr;
// 	recHit_.x = lp.x();
// 	recHit_.y = lp.y();
//       }  
//       if (PRINT) cout<<"simhit "<<sim_.x<<" "<<sim_.y<<" "<<dr<<endl;
//     }
    

//     //--- If there was a match, capture the recHiter quantities
//     if (dr_min < 99999.) { 
//       fillRecHit(rmatchIt);           // copy relevant stuff from rmatchIt into recHit_
//     }
//     else {
//       recHit_.init();            // wipe it all out
//       //if(PRINT) cout<<" no match for this hit "<<endl;
//     }



    //++++++++++
    t_->Fill();
    //++++++++++

  }
}

// void PixelNtuplizer::fillRecHit( SiPixelRecHitCollection::ContainerIterator recHitIt)
// {
//       LocalPoint lp = rechitIt->localPosition();
//       LocalError le = rechitIt->localPositionError();
//       float dr_rh = 
// 	(sim_xpos-lp.x())*(sim_xpos-lp.x()) + (sim_ypos-lp.y())*(sim_ypos-lp.y());
//
//       if(dr_rh<dr_rh0){
// 	dr_rh0=dr_rh;
// 	recHit_.x = lp.x();
// 	recHit_.y = lp.y();
// 	recHit_.xx = le.xx();
// 	recHit_.xy = le.xy();
// 	recHit_.yy = le.yy();
//       }
//       if(PRINT)cout<<"---> rechit found: x " << lp.x() <<"  y "<<lp.y() << " distance " << dr_rh0 << endl;
//     }
//
//     if(dr_rh0==99999.) {
//       recHit_.x = 9999.;
//       recHit_.y = 9999.;
//     }
//    
//     if(dr0<99999.) { // some match was found
//       if(PRINT) cout<<"match "<<dr0<<" "<<dx<<" "<<dy<<endl;
//       //hdr->Fill(sqrt(dr0)*10000.);
//       //hresX1->Fill(dx);
//       //hresY1->Fill(dy);
//     } else {
//       //hdr->Fill(999.);
//       if(PRINT) cout<<" no match for this hit "<<endl;
//       clust_.x = 9999.;
//       clust_.y = 9999.;
//     }
//
// }

void PixelNtuplizer::fillRecHit(LocalPoint lp, LocalError le) {
        recHit_.x = lp.x();
        recHit_.y = lp.y();
        recHit_.xx = le.xx();
        recHit_.xy = le.xy();
        recHit_.yy = le.yy();
}


void PixelNtuplizer::fillSim(std::vector<PSimHit>::iterator isim) {
    float sim_x1 = (*isim).entryPoint().x(); // width (row index, in col direction)
    float sim_y1 = (*isim).entryPoint().y(); // length (col index, in row direction) 
    float sim_x2 = (*isim).exitPoint().x();
    float sim_y2 = (*isim).exitPoint().y();
    float sim_xpos = 0.5*(sim_x1+sim_x2);
    float sim_ypos = 0.5*(sim_y1+sim_y2);         
    sim_.x     = sim_xpos;
    sim_.y     = sim_ypos;
   
    sim_.px  = (*isim).momentumAtEntry().x();
    sim_.py  = (*isim).momentumAtEntry().y();
    sim_.pz  = (*isim).momentumAtEntry().z();
    sim_.eloss = (*isim).energyLoss();
    sim_.phi   = (*isim).phiAtEntry();
    sim_.theta = (*isim).thetaAtEntry();
    sim_.PID = (*isim).particleType();
    sim_.TID = (*isim).trackId();

    //--- Fill alpha and beta -- more useful for exploring the residuals...
    sim_.beta  = atan2(sim_.pz, sim_.py);
    sim_.alpha = atan2(sim_.pz, sim_.px);
 }


void PixelNtuplizer::fillDet(DetId &tofill, int subdetid)
{
  if(subdetid==1) {
      det_.layer  = PXBDetId::PXBDetId(tofill).layer();
      det_.ladder = PXBDetId::PXBDetId(tofill).ladder();
      det_.module = PXBDetId::PXBDetId(tofill).module();
    } else {
      det_.disk      =  PXFDetId::PXFDetId(tofill).disk();
      det_.blade     =  PXFDetId::PXFDetId(tofill).blade();
      det_.panel     =  PXFDetId::PXFDetId(tofill).panel();
      det_.plaquette =  PXFDetId::PXFDetId(tofill).module();
      
      //Following Danek's advice...
      unsigned int side = PXFDetId::PXFDetId(tofill).side();
      if (side==1) det_.disk = - det_.disk; 
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
  clust_.charge = (matchIt->charge())/1000.; // convert ke to electrons
  clust_.size = matchIt->size();
  clust_.size_x = matchIt->sizeX();
  clust_.size_y = matchIt->sizeY();
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
  float dummy_float = 9999.0;
  int dummy_int = 9999;
  //unsigned int dummy_uint = 9999;
  //bool dummy_bool = false;
  
  thickness = dummy_float;
  cols  = dummy_int;
  rows  = dummy_int;
  layer = dummy_int;
  ladder = dummy_int;
  module = dummy_int;
  disk = dummy_int;
  blade = dummy_int;
  panel = dummy_int;
  plaquette = dummy_int;
}

void PixelNtuplizer::vertex::init()
{
  float dummy_float = 9999.0;
  int dummy_int = 9999;
  //unsigned int dummy_uint = 9999;
  //bool dummy_bool = false;

  num = dummy_int;
  r = dummy_float;
  z = dummy_float;
 }

void PixelNtuplizer::track::init()
{
  float dummy_float = 9999.0;
  //int dummy_int = 9999;
  //unsigned int dummy_uint = 9999;
  //bool dummy_bool = false;

  eta = dummy_float;
  phi = dummy_float;
}

void PixelNtuplizer::sim::init()
{
  float dummy_float = 9999.0;
  int dummy_int = 9999;
  //unsigned int dummy_uint = 9999;
  //bool dummy_bool = false;
  
  x = dummy_float;
  y = dummy_float;
  px = dummy_float;
  py = dummy_float;
  pz = dummy_float;
  eloss = dummy_float;
  phi = dummy_float;
  theta = dummy_float;
  isflipped = dummy_int;
}



void PixelNtuplizer::clust::init()
{
  float dummy_float = 9999.0;
  int dummy_int = 9999;
  unsigned int dummy_uint = 9999;
  bool dummy_bool = false;
  
  charge = dummy_float; // convert ke to electrons
  size = dummy_int;
  size_x = dummy_int;
  size_y = dummy_int;
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
  float dummy_float = 9999.0;
  //int dummy_int = 9999;
  //unsigned int dummy_uint = 9999;
  //bool dummy_bool = false;

  x = dummy_float;
  y = dummy_float;
  xx = dummy_float;
  xy = dummy_float; 
  yy = dummy_float;
}


