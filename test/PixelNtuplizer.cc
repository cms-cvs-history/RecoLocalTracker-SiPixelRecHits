// File: PixelNtuplizer.cc
// Description:  see PixelNtuplizer.h
// Authors:  Petar Maksimovic, Jason Shaev, JHU...Vincenzo Chiochia, CERN
// History: 4/4/06   creation, new version: 6/6/06
//--------------------------------------------
/*
   The new version works by iterating over all detids. It then grabs rechits from the collection based upon this detid. This rechit can then be associated to a sim hit using the TrackerAssociator tool. It grabs the cluster using the cluster function of the rec hit.
*/
	
#include "RecoLocalTracker/SiPixelRecHits/test/PixelNtuplizer.h"


// DataFormats
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/Ref.h"

// Old
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


// SimDataFormats
//#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/EmbdSimTrack.h"
#include "SimDataFormats/Vertex/interface/EmbdSimVertex.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

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
using namespace edm;

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
  
  t_->Branch("evt", &evt_, "run/I:evtnum/I", bufsize);
  
   t_->Branch("det", &det_, "thickness/F:cols/I:rows/I:layer/I:ladder/I:module/I:disk/I:blade/I:panel/I:plaquette/I", bufsize);

    std::cout << "Making vertex branch:" << std::endl;
  t_->Branch("vertex",   &vertex_,   "num/I:r/F:z/F", bufsize);

  std::cout << "Making track branch:" << std::endl;
   t_->Branch("track", &track_, "eta/F:phi/F", bufsize);

  std::cout << "Making sim hit branch:" << std::endl;
  t_->Branch("sim",   &sim_,  
 
"x/F:y/F:px/F:py/F:pz/F:eloss/F:phi/F:theta/F:subdetid/I:isflipped/I:alpha/F:beta/F:PID/I:TID/i:x1/F:x2:y1:y2:z1:z2:row1:row2:col1:col2:gx1:gx2:gy1:gy2:gz1:gz2", bufsize);

  std::cout << "Making cluster branch:" << std::endl;
  t_->Branch("clust", &clust_, "row/F:col/F:x/F:y/F:charge/F:size/I:size_x/I:size_y/I:maxPixelCol/I:maxPixelRow/I:minPixelCol/I:minPixelRow/I:geoId/i:edgeHitX/O:edgeHitY/O", bufsize);
  
  std::cout << "Making pixinfo branch:" << std::endl;
  t_->Branch("npix", &pixinfo_.npix, "npix/I", bufsize);
  t_->Branch("rowpix", pixinfo_.row, "row[npix]/F", bufsize);
  t_->Branch("colpix", pixinfo_.col, "col[npix]/F", bufsize);
  t_->Branch("adc", pixinfo_.adc, "adc[npix]/F", bufsize);
  t_->Branch("xpix", pixinfo_.x, "x[npix]/F", bufsize);
  t_->Branch("ypix", pixinfo_.y, "y[npix]/F", bufsize);
  t_->Branch("gxpix", pixinfo_.gx, "gx[npix]/F", bufsize);
  t_->Branch("gypix", pixinfo_.gy, "gy[npix]/F", bufsize);
  t_->Branch("gzpix", pixinfo_.gz, "gz[npix]/F", bufsize);
  

  std::cout << "Making recHit branch:" << std::endl;
  t_->Branch("recHit", &recHit_, "x/F:y:xx:xy:yy:row:col:gx:gy:gz", bufsize);

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

  //--- Fetch Pixel RecHits
  std::string recHitCollLabel = conf_.getUntrackedParameter<std::string>("RecHitCollLabel","pixRecHitConverter");
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel(recHitCollLabel, recHitColl);
  
  std::cout 
    <<" FOUND " 
    << const_cast<SiPixelRecHitCollection*>(recHitColl.product())->size()
    << " Pixel RecHits" << std::endl;

  //--- Get the simtracks for matching
  Handle<edm::EmbdSimTrackContainer> simtracks;
  e.getByLabel("SimG4Object",simtracks);

  bool PRINT = true;

  TrackerHitAssociator associate(e);

  //-----Iterate over detunits
   for (TrackerGeometry::DetContainer::const_iterator it = geom->dets().begin(); it != geom->dets().end(); it++) {

    //--- Since TTree::Fill() will simply take a snapshot of what's 
    //--- in the memory, we need to re-initialize the cache used for 
    //--- each of the branches.
    init();

    DetId detId = ((*it)->geographicalId());
    unsigned int subid=detId.subdetId();
    if (! ((subid==1) || 
	   (subid==2))) {
      continue;
    } // end subid if

    // Get the geom-detector
    const PixelGeomDetUnit * theGeomDet =
      dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );
    vertex_.z = theGeomDet->surface().position().z();
    vertex_.r = theGeomDet->surface().position().perp();

    const BoundPlane& plane = theGeomDet->surface(); //for transf.
    
    det_.thickness = theGeomDet->specificSurface().bounds().thickness();
    det_.cols = theGeomDet->specificTopology().ncolumns();
    det_.rows = theGeomDet->specificTopology().nrows();

    SiPixelRecHitCollection::range pixelrechitRange = (recHitColl.product())->get(detId);
    SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorBegin = pixelrechitRange.first;
    SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorEnd = pixelrechitRange.second;
    SiPixelRecHitCollection::const_iterator pixeliter = pixelrechitRangeIteratorBegin;
    std::vector<PSimHit> matched;

    //----Loop over rechits for this detId
    for ( ; pixeliter != pixelrechitRangeIteratorEnd; ++pixeliter) {
 
        matched.clear();
        matched = associate.associateHit(*pixeliter);

	 // Get the geom-detector
         const PixelGeomDetUnit * theGeomDet =
      		dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );
    	 vertex_.z = theGeomDet->surface().position().z();
    	 vertex_.r = theGeomDet->surface().position().perp();
	 
	  const RectangularPixelTopology * topol = 
         dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));

	  fillRecHit(pixeliter, topol, theGeomDet);

        if (!matched.empty()) {

           //---Loop over sim hits, fill closest
	   float closest=9999.;
	   std::vector<PSimHit>::const_iterator old = matched.begin();

	   for (std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++) {
        	fillSim(m, subid, theGeomDet, topol);
		float x_res = sim_.x - recHit_.x;

		if (x_res < closest) {
		   closest = x_res;
		} else {
		   fillSim(old, subid, theGeomDet, topol);
		}
	
	   } // end sim hit loop

	   edm::Ref<edm::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = pixeliter->cluster();

	  fillEvt(e);
	  
	  fillClust(*clust, topol, theGeomDet);
       	   
	  fillDet(detId, subid);

	  fillTrack(*(simtracks.product()));
    
  	 det_.thickness = theGeomDet->specificSurface().bounds().thickness();
    	 det_.cols = theGeomDet->specificTopology().ncolumns();
	 det_.rows = theGeomDet->specificTopology().nrows();
	 
	 fillPix(*clust, topol, theGeomDet);


	   //++++++++++
	   t_->Fill();
	   //++++++++++
	   init();
	   
	} // end matched if
	
    } // end rechit loop	

   } //end det loop
   
} // end analyze function

void PixelNtuplizer::fillRecHit(SiPixelRecHitCollection::const_iterator pixeliter,
				const RectangularPixelTopology * topol, const PixelGeomDetUnit * PixGeom) {
	LocalPoint lp = pixeliter->localPosition();
	LocalError le = pixeliter->localPositionError();

        recHit_.x = lp.x();
        recHit_.y = lp.y();
        recHit_.xx = le.xx();
        recHit_.xy = le.xy();
        recHit_.yy = le.yy();
	MeasurementPoint mp = topol->measurementPosition(LocalPoint(recHit_.x, recHit_.y));
	recHit_.row = mp.x();
	recHit_.col = mp.y();
	GlobalPoint GP = PixGeom->surface().toGlobal(Local3DPoint(lp.x(),lp.y(),lp.z()));
	recHit_.gx = GP.x();
	recHit_.gy = GP.y();
	recHit_.gz = GP.z();
}


void PixelNtuplizer::fillSim(std::vector<PSimHit>::const_iterator isim, unsigned int subid, const PixelGeomDetUnit * theGeomDet, 
                             const RectangularPixelTopology * topol) {
    sim_.x1 = (*isim).entryPoint().x(); // width (row index, in col direction)
    sim_.y1 = (*isim).entryPoint().y(); // length (col index, in row direction) 
    sim_.z1 = (*isim).entryPoint().z(); 
    sim_.x2 = (*isim).exitPoint().x();
    sim_.y2 = (*isim).exitPoint().y();
    sim_.z2 = (*isim).exitPoint().z();
  GlobalPoint GP1 = theGeomDet->surface().toGlobal(Local3DPoint((*isim).entryPoint().x(),(*isim).entryPoint().y(),(*isim).entryPoint().z()));
    GlobalPoint GP2 = theGeomDet->surface().toGlobal(Local3DPoint((*isim).exitPoint().x(),(*isim).exitPoint().y(),(*isim).exitPoint().z()));
    sim_.gx1 = GP1.x();
    sim_.gx2 = GP2.x();
    sim_.gy1 = GP1.y();
    sim_.gy2 = GP2.y();
    sim_.gz1 = GP1.z();
    sim_.gz2 = GP2.z();
    
    MeasurementPoint mp1 =
    topol->measurementPosition(LocalPoint((*isim).entryPoint().x(),(*isim).entryPoint().y(),(*isim).entryPoint().z()));
    MeasurementPoint mp2 =
    topol->measurementPosition(LocalPoint((*isim).exitPoint().x(),(*isim).exitPoint().y(), (*isim).exitPoint().z()));
    sim_.row1 = mp1.x();
    sim_.col1 = mp1.y();
    sim_.row2 = mp2.x();
    sim_.col2 = mp2.y();
    
    float sim_xpos = 0.5*(sim_.x1+sim_.x2);
    float sim_ypos = 0.5*(sim_.y1+sim_.y2);         
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
    sim_.subdetid = subid;

    //--- Fill alpha and beta -- more useful for exploring the residuals...
    sim_.beta  = atan2(sim_.pz, sim_.py);
    sim_.alpha = atan2(sim_.pz, sim_.px);

 // Is flipped ?
    float tmp1 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
    float tmp2 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
    //cout << " 1: " << tmp1 << " 2: " << tmp2 << endl;
    if ( tmp2<tmp1 ) sim_.isflipped = 1;
    else sim_.isflipped = 0;
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
PixelNtuplizer::fillClust(const SiPixelCluster & matchIt, const RectangularPixelTopology * topol, const PixelGeomDetUnit * PixGeom) 
{

  //       if(PRINT) cout<<"clus "<<ch<<" "<<size<<" "<<sizeX<<" "<<sizeY<<" "
  // 		    <<x<<" "<<y<<" "<<geoId<<" "<<edgeHitX<<" "
  // 		    <<edgeHitY<<endl;
  //       if(PRINT) cout<<"match "<<dr0<<" "<<dx<<" "<<dy<<endl;
  
  //const vector<Pixel>  = clustIt->pixels();
  clust_.charge = (matchIt.charge())/1000.; // convert ke to electrons
  clust_.size = matchIt.size();
  clust_.size_x = matchIt.sizeX();
  clust_.size_y = matchIt.sizeY();
  clust_.row = matchIt.x();
  clust_.col = matchIt.y();
  
  LocalPoint lp = topol->localPosition(MeasurementPoint(clust_.row,clust_.col));
  clust_.x = lp.x();
  clust_.y = lp.y();
  
  clust_.maxPixelCol = matchIt.maxPixelCol();
  clust_.maxPixelRow = matchIt.maxPixelRow();
  clust_.minPixelCol = matchIt.minPixelCol();
  clust_.minPixelRow = matchIt.minPixelRow();
  
  clust_.geoId = matchIt.geographicalId();
  clust_.edgeHitX = matchIt.edgeHitX();
  clust_.edgeHitY = matchIt.edgeHitY();
  
}

void
PixelNtuplizer::fillPix(const SiPixelCluster & LocPix, 
			const RectangularPixelTopology * topol, const PixelGeomDetUnit * PixGeom)
{
   const std::vector<SiPixelCluster::Pixel>& pixvector = LocPix.pixels();
   for( ; pixinfo_.npix < pixvector.size(); ++pixinfo_.npix)
   {
      SiPixelCluster::Pixel holdpix = pixvector[pixinfo_.npix];
      pixinfo_.row[pixinfo_.npix] = holdpix.x;
      pixinfo_.col[pixinfo_.npix] = holdpix.y;
      pixinfo_.adc[pixinfo_.npix] = holdpix.adc;
      LocalPoint lp = topol->localPosition(MeasurementPoint(holdpix.x, holdpix.y));
      pixinfo_.x[pixinfo_.npix] = lp.x();
      pixinfo_.y[pixinfo_.npix]= lp.y();
      GlobalPoint GP =  PixGeom->surface().toGlobal(Local3DPoint(lp.x(),lp.y(),lp.z()));
  pixinfo_.gx[pixinfo_.npix] = GP.x();	
  pixinfo_.gy[pixinfo_.npix]= GP.y();
  pixinfo_.gz[pixinfo_.npix]= GP.z();
   }
}

void 
PixelNtuplizer::fillTrack(const edm::EmbdSimTrackContainer& trks) {
//  std::cout << "Looking for simtrack number " << sim_.TID << std::endl;
//  for (unsigned int i=0; i<trks.size(); ++i) {
//    std::cout << i << ": " << trks[i] << std::endl;
//  }
  track_.eta = trks[sim_.TID].momentum().eta();
  track_.phi = trks[sim_.TID].momentum().phi();
}

void
PixelNtuplizer::fillEvt(const edm::Event& E)
{
   evt_.run = E.id().run();
   evt_.evtnum = E.id().event();
}

void PixelNtuplizer::init() {
	evt_.init();
	det_.init();
	vertex_.init();
	clust_.init();
	sim_.init();
	track_.init();
	pixinfo_.init();
	recHit_.init();
}

void PixelNtuplizer::evt::init()
{
  int dummy_int = 9999;
  run = dummy_int;
  evtnum = dummy_int;
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
  x1 = x2 = y1 = y2 = z1 = z2 = dummy_float;
  row1 = dummy_float;
  row2 = dummy_float;
  col1 = dummy_float;
  col2 = dummy_float;
  gx1 = dummy_float;
  gx2 = dummy_float;
  gy1 = dummy_float;
  gy2 = dummy_float;
  gz1 = dummy_float;
  gz2 = dummy_float;
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
  row = dummy_float;
  col = dummy_float;
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

void PixelNtuplizer::pixinfo::init()
{
   //float dummy_float = 9999.0;
   npix = 0;
  /* for(int i = 0; i != maxpix; ++i)
   {
      row[i] = dummy_float;
      col[i] = dummy_float;
      adc[i] = dummy_float;
      x[i] = dummy_float;
      y[i] = dummy_float;
      gx[i] = dummy_float;
      gy[i] = dummy_float;
      gz[i] = dummy_float;
   }*/
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
  row = dummy_float;
  col = dummy_float;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
}


