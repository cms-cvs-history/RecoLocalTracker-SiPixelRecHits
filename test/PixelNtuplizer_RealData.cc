// File: PixelNtuplizer_RealData.cc
// Description: Adapted tree structure from PixelNtuplizer package
//   		Adapted method for finding residual information from TrackerValidationVariables package
// Authors: Andrew York, Tennessee
//
//
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "DQM/SiStripCommon/interface/SiStripHistoId.h"
#include "DQM/TrackerMonitorTrack/interface/MonitorTrackResiduals.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"  
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"

//----------------------------------------------------------------------

#include "RecoLocalTracker/SiPixelRecHits/test/PixelNtuplizer_RealData.h"

// DataFormats
#include "DataFormats/Math/interface/deltaPhi.h"
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
#include "DQM/SiPixelCommon/interface/SiPixelFolderOrganizer.h"
#include "DQM/SiPixelMonitorTrack/interface/SiPixelMonitorTrackResiduals.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"

// SimDataFormats
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
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

PixelNtuplizer_RD::PixelNtuplizer_RD(edm::ParameterSet const& ps) : 
  conf_(ps),
  tfile_(0), 
  t_(0)
{}

// Virtual destructor needed.
PixelNtuplizer_RD::~PixelNtuplizer_RD() { }  

// End job: write and close the ntuple file
void PixelNtuplizer_RD::endJob() 
{
  counters_.trackCounter = TrackCounter;
  counters_.pixelTrackCounter = PixelTrackCounter;
  //++++++++++
  tc_->Fill();
  //++++++++++
  std::string outputFileName = conf_.getParameter<std::string>("OutputFile");
  std::cout << " PixelNtuplizer_RealData::endJob" << std::endl;
  tfile_->Write();
  tfile_->Close();

}


void PixelNtuplizer_RD::beginJob(const edm::EventSetup& es)
{
  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  es.get<IdealMagneticFieldRecord>().get(magneticField_);

  // put here whatever you want to do at the beginning of the job
  std::string outputFile = conf_.getParameter<std::string>("OutputFile");
  tfile_ = new TFile ( outputFile.c_str() , "RECREATE" );

  t_ = new TTree("PixNtuple", "Pixel hit analyzer ntuple");
  ts_ = new TTree("StripNtuple", "Strip hit analyzer ntuple");
  tc_ = new TTree("CounterNtuple", "Counters filled every event");
  int bufsize = 64000;

  // Create one branch. If splitlevel is set, event is a superbranch
  // creating a sub branch for each data member of the Event object.
  t_->Branch("evt", &evt_, "run/I:evtnum", bufsize);
  
  t_->Branch("det", &det_, "thickness/F:cols/I:rows/I:layer/I:ladder/I:module/I:disk/I:blade/I:panel/I:plaquette/I", bufsize);

  std::cout << "Making vertex branch:" << std::endl;
  t_->Branch("vertex",   &vertex_,   "r/F:z", bufsize);

  std::cout << "Making cluster branch:" << std::endl;
  t_->Branch("Cluster", &clust_, "row/F:col:x:y:charge:size/I:size_x:size_y:maxPixelCol:maxPixelRow:minPixelCol:minPixelRow:geoId/i:edgeHitX/I:edgeHitY:clust_alpha/F:clust_beta", bufsize);
  
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

  std::cout << "Making rechit branch:" << std::endl;
  t_->Branch("RecHit", &rechit_, "localX/F:localY:globalX:globalY:globalZ:residualX:residualY:resErrX:resErrY:resXprime:resXprimeErr", bufsize);

  std::cout << "Making track branch:" << std::endl;
  t_->Branch("track", &track_, "pt/F:px:py:pz:globalEta:globalPhi:localEta:localPhi", bufsize);

  std::cout << "Making tracker hit branch:" << std::endl;
  ts_->Branch("TrackerHit", &trackerhits_, "globalX/F:globalY:globalZ", bufsize);

  std::cout << "Making counter branch:" << std::endl;
  tc_->Branch("Counters", &counters_, "trackCounter/I:pixelTrackCounter", bufsize);
  
  std::cout << "Made all branches." << std::endl;
  TrackCounter = 0;
  PixelTrackCounter = 0;

}


// Functions that get called by framework every event
void PixelNtuplizer_RD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  counters_.init();
  edm::Handle<std::vector<Trajectory> > trajCollectionHandle;
  iEvent.getByLabel(conf_.getParameter<std::string>("trajectoryInput"),trajCollectionHandle);
  
  TrajectoryStateCombiner tsoscomb;
  edm::LogVerbatim("TrackerValidationVariables") << "trajColl->size(): " << trajCollectionHandle->size() ;
  for(std::vector<Trajectory>::const_iterator it = trajCollectionHandle->begin(), itEnd = trajCollectionHandle->end(); 
      it!=itEnd;++it){

    TrackCounter++;
    bool trajContainsPixelHit = false;

    std::vector<TrajectoryMeasurement> checkColl = it->measurements();
    for(std::vector<TrajectoryMeasurement>::const_iterator checkTraj = checkColl.begin(), checkTrajEnd = checkColl.end(); 
	checkTraj != checkTrajEnd; ++checkTraj) {

      if(! checkTraj->updatedState().isValid()) continue;
      TransientTrackingRecHit::ConstRecHitPointer testhit = checkTraj->recHit();
      if(! testhit->isValid() || testhit->geographicalId().det() != DetId::Tracker ) continue;
      uint testSubDetID = (testhit->geographicalId().subdetId());
      if(testSubDetID == PixelSubdetector::PixelBarrel || testSubDetID == PixelSubdetector::PixelEndcap) 
	trajContainsPixelHit = true;
    }

    if (!trajContainsPixelHit) continue;

    PixelTrackCounter++;

    std::vector<TrajectoryMeasurement> tmColl = it->measurements();
    for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = tmColl.begin(), itTrajEnd = tmColl.end(); 
	itTraj != itTrajEnd; ++itTraj) {

      if(! itTraj->updatedState().isValid()) continue;
      
      
      TrajectoryStateOnSurface tsos = tsoscomb( itTraj->forwardPredictedState(), itTraj->backwardPredictedState() );
      TransientTrackingRecHit::ConstRecHitPointer hit = itTraj->recHit();
      if(! hit->isValid() || hit->geographicalId().det() != DetId::Tracker ) {
	continue; 
      } else {
	const DetId & hit_detId = hit->geographicalId();
	uint IntRawDetID = (hit_detId.rawId());
	uint IntSubDetID = (hit_detId.subdetId());
	
	if(IntSubDetID == 0 )
	  continue;
	
	align::LocalVector res = tsos.localPosition() - hit->localPosition();

	LocalError err1 = tsos.localError().positionError();
	LocalError err2 = hit->localPositionError();
	
	float errX = std::sqrt( err1.xx() + err2.xx() );
	float errY = std::sqrt( err1.yy() + err2.yy() );
	
	LogDebug("PixelNtuplizer_RealData") << "Residual x/y " << res.x() << '/' << res.y() 
					       << ", Error x/y " << errX << '/' << errY;		

	// begin partly copied from Tifanalyser 

	const GeomDetUnit* detUnit = hit->detUnit();
	double dPhi = -999, dR = -999, dZ = -999, phiorientation = -999;
	double R = 0.;
	double origintointersect = 0.;	


	if(detUnit) {
	  const Surface& surface = hit->detUnit()->surface();
	  LocalPoint lPModule(0.,0.,0.), lPhiDirection(1.,0.,0.), lROrZDirection(0.,1.,0.);
	  GlobalPoint gPModule       = surface.toGlobal(lPModule),
	    gPhiDirection  = surface.toGlobal(lPhiDirection),
	    gROrZDirection = surface.toGlobal(lROrZDirection);
	  phiorientation = deltaPhi(gPhiDirection.phi(),gPModule.phi()) >= 0 ? +1. : -1.;

	  dPhi = tsos.globalPosition().phi() - hit->globalPosition().phi();
	  
	  if(IntSubDetID == PixelSubdetector::PixelBarrel || IntSubDetID == PixelSubdetector::PixelEndcap) {

            const TrackerGeometry& theTracker(*tkGeom_);

	    const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(hit_detId) );

            const RectangularPixelTopology * topol = dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));

            // get the enclosed persistent hit
            const TrackingRecHit *persistentHit = hit->hit();
            // check if it's not null, and if it's a valid pixel hit
            if ((persistentHit != 0) && (typeid(*persistentHit) == typeid(SiPixelRecHit))) {
              // tell the C++ compiler that the hit is a pixel hit
              const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>( hit->hit() );
              // get the edm::Ref to the cluster
              edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = (*pixhit).cluster();
              //  check if the ref is not null
              if (clust.isNonnull()) {

		init();
		rechit_.localX = hit->localPosition().x();
		rechit_.localY = hit->localPosition().y();
		rechit_.globalX = hit->globalPosition().x();
		rechit_.globalY = hit->globalPosition().y();
		rechit_.globalZ = hit->globalPosition().z();
		rechit_.residualX = res.x();
		rechit_.resErrX = errX;
		rechit_.resErrY = errY;
	        rechit_.resXprime = (res.x())*(phiorientation );
	        rechit_.resXprimeErr = errX;
		dZ = gROrZDirection.z() - gPModule.z();
		if(dR != -999) rechit_.residualY = dR;
		else if(dZ != -999) rechit_.residualY = res.y() * (dZ >=0.? +1 : -1) ;
		else rechit_.residualY = res.y();

                // get the contents
                fillEvt(iEvent);
                fillDet(hit_detId, IntSubDetID, theGeomDet);
		fillVertex(theGeomDet);
                fillClust(*clust, topol, theGeomDet);
                fillPix(*clust, topol, theGeomDet);
                fillTrack( tsos );
                //fillRecHit(pixeliter, topol, theGeomDet);

      	        //++++++++++
	        t_->Fill();
	        //++++++++++

              } // end if(cluster is valid)
            } // end if(hit is pixel hit)

	  } else if (IntSubDetID == StripSubdetector::TIB || IntSubDetID == StripSubdetector::TOB ||
		     IntSubDetID == StripSubdetector::TID || IntSubDetID == StripSubdetector::TEC) {
	    /*	    const RadialStripTopology* theTopol = dynamic_cast<const RadialStripTopology*>(&(detUnit->topology()));
	    origintointersect =  static_cast<float>(theTopol->originToIntersection());
	    	    
	    MeasurementPoint theMeasHitPos = theTopol->measurementPosition(hit->localPosition());
	    MeasurementPoint theMeasStatePos = theTopol->measurementPosition(tsos.localPosition());
	    Measurement2DVector residual =  theMeasStatePos - theMeasHitPos;
	    
	    MeasurementError theMeasHitErr = theTopol->measurementError(hit->localPosition(),err2);
	    MeasurementError theMeasStateErr = theTopol->measurementError(tsos.localPosition(),err1);

	    double localPitch = theTopol->localPitch(hit->localPosition());
	    float xPrime = residual.x()*localPitch ;
	    float measErr = std::sqrt( theMeasHitErr.uu()+theMeasStateErr.uu())*localPitch;

	    R = origintointersect;
	    dR = theTopol->yDistanceToIntersection( tsos.localPosition().y()) - 
	      theTopol->yDistanceToIntersection( hit->localPosition().y());
	    */
	    trackerhits_.init();
	    trackerhits_.globalX = hit->globalPosition().x();
	    trackerhits_.globalY = hit->globalPosition().y();
	    trackerhits_.globalZ = hit->globalPosition().z();

	    //++++++++++
	    ts_->Fill();
	    //++++++++++

	  } else {
	    edm::LogWarning("PixelNtuplizer_RealData") << "@SUB=PixelNtuplizer_RealData::fillHitQuantities" 
							  << "No valid tracker subdetector " << IntSubDetID;
	    rechit_.resXprime = -999;
	  }  // if-else to differentiate pixel hits vs tracker hits	  
	}  // end if(good detUnit)
      }  // end else (if-else to screen out invalid hits)
    }  // end loop over trajectory measurements (rec hits) 
  }  // end loop over trajectories
}  // end analyze function


void PixelNtuplizer_RD::fillEvt(const edm::Event& E)
{
  evt_.run = E.id().run();
  evt_.evtnum = E.id().event();
}

void PixelNtuplizer_RD::fillDet(const DetId &tofill, uint subdetid, const PixelGeomDetUnit* PixGeom)
{
  if (subdetid==1) 
    {
      det_.layer  = PXBDetId::PXBDetId(tofill).layer();
      det_.ladder = PXBDetId::PXBDetId(tofill).ladder();
      det_.module = PXBDetId::PXBDetId(tofill).module();
    } 
  else
    {
      det_.disk      =  PXFDetId::PXFDetId(tofill).disk();
      det_.blade     =  PXFDetId::PXFDetId(tofill).blade();
      det_.panel     =  PXFDetId::PXFDetId(tofill).panel();
      det_.plaquette =  PXFDetId::PXFDetId(tofill).module();
      
      //Following Danek's advice...
      //      unsigned int side = PXFDetId::PXFDetId(tofill).side();
      //      if (side==1) det_.disk = - det_.disk; 
    }
  det_.thickness = PixGeom->specificSurface().bounds().thickness();
  det_.cols = PixGeom->specificTopology().ncolumns();
  det_.rows = PixGeom->specificTopology().nrows();
}

void PixelNtuplizer_RD::fillVertex(const PixelGeomDetUnit* PixGeom)
{
  vertex_.z = PixGeom->surface().position().z();
  vertex_.r = PixGeom->surface().position().perp();
}

void PixelNtuplizer_RD::fillClust(const SiPixelCluster& matchIt, const RectangularPixelTopology* topol, const PixelGeomDetUnit* PixGeom) 
{
  clust_.charge = (matchIt.charge())/1000.0; // convert electrons to kilo-electrons
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
  
  clust_.geoId = PixGeom->geographicalId().rawId();

  // Replace with the topology methods
  // edge method moved to topologi class
  clust_.edgeHitX = (int) ( topol->isItEdgePixelInX( clust_.minPixelRow ) || topol->isItEdgePixelInX( clust_.maxPixelRow ) );
  clust_.edgeHitY = (int) ( topol->isItEdgePixelInY( clust_.minPixelCol ) || topol->isItEdgePixelInY( clust_.maxPixelCol ) );

  // calculate alpha and beta from cluster position

  // get cluster center of gravity (of charge)
  float xcenter = matchIt.x();
  float ycenter = matchIt.y();
  
  // get the cluster position in local coordinates (cm) 
  LocalPoint mylp = topol->localPosition( MeasurementPoint(xcenter, ycenter) );

  // get the cluster position in global coordinates (cm)
  GlobalPoint gp = PixGeom->surface().toGlobal( mylp );
  float gp_mod = sqrt( gp.x()*gp.x() + gp.y()*gp.y() + gp.z()*gp.z() );

  // normalize
  float gpx = gp.x()/gp_mod;
  float gpy = gp.y()/gp_mod;
  float gpz = gp.z()/gp_mod;

  // make a global vector out of the global point; this vector will point from the 
  // origin of the detector to the cluster
  GlobalVector gv(gpx, gpy, gpz);

  // make local unit vector along local X axis
  const Local3DVector lvx(1.0, 0.0, 0.0);

  // get the unit X vector in global coordinates
  GlobalVector gvx = PixGeom->surface().toGlobal( lvx );

  // make local unit vector along local Y axis
  const Local3DVector lvy(0.0, 1.0, 0.0);

  // get the unit Y vector in global coordinates
  GlobalVector gvy = PixGeom->surface().toGlobal( lvy );
   
  // make local unit vector along local Z axis
  const Local3DVector lvz(0.0, 0.0, 1.0);

  // get the unit Z vector in global coordinates
  GlobalVector gvz = PixGeom->surface().toGlobal( lvz );
    
  // calculate the components of gv (the unit vector pointing to the cluster) 
  // in the local coordinate system given by the basis {gvx, gvy, gvz}
  // note that both gv and the basis {gvx, gvy, gvz} are given in global coordinates
  float gv_dot_gvx = gv.x()*gvx.x() + gv.y()*gvx.y() + gv.z()*gvx.z();
  float gv_dot_gvy = gv.x()*gvy.x() + gv.y()*gvy.y() + gv.z()*gvy.z();
  float gv_dot_gvz = gv.x()*gvz.x() + gv.y()*gvz.y() + gv.z()*gvz.z();

  // calculate angles
  clust_.clust_alpha = atan2( gv_dot_gvz, gv_dot_gvx );
  clust_.clust_beta  = atan2( gv_dot_gvz, gv_dot_gvy );

}

void PixelNtuplizer_RD::fillPix(const SiPixelCluster & LocPix, const RectangularPixelTopology * topol, const PixelGeomDetUnit * PixGeom)
{
  const std::vector<SiPixelCluster::Pixel>& pixvector = LocPix.pixels();
  for ( ; pixinfo_.npix < (int)pixvector.size(); ++pixinfo_.npix)
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

void PixelNtuplizer_RD::fillTrack(TrajectoryStateOnSurface& tsos) 
{
  track_.pt = tsos.globalMomentum().transverse();
  track_.px = tsos.globalMomentum().x();
  track_.py = tsos.globalMomentum().y();
  track_.pz = tsos.globalMomentum().z();
  track_.globalPhi = tsos.globalDirection().phi();
  track_.globalEta = tsos.globalDirection().eta();
  track_.localPhi = tsos.localDirection().phi();
  track_.localEta = tsos.localDirection().eta();
}


void PixelNtuplizer_RD::init() 
{
  evt_.init();
  det_.init();
  vertex_.init();
  clust_.init();
  pixinfo_.init();
  rechit_.init();
  track_.init();
}

void PixelNtuplizer_RD::EvtStruct::init()
{
  int dummy_int = -9999;

  run = dummy_int;
  evtnum = dummy_int;
}

void PixelNtuplizer_RD::DetStruct::init()
{
  int dummy_int = -9999;
  float dummy_float = -9999.0;
 
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

void PixelNtuplizer_RD::VertexStruct::init()
{
  float dummy_float = -9999.0;

  r = dummy_float;
  z = dummy_float;
 }

void PixelNtuplizer_RD::ClusterStruct::init()
{
  float dummy_float = -9999.0;
  int dummy_int = -9999;
  unsigned int dummy_uint = 9999;
  
  row = dummy_float;
  col = dummy_float;
  x = dummy_float;
  y = dummy_float;
  charge = dummy_float;
  size = dummy_int;
  size_x = dummy_int;
  size_y = dummy_int;
  maxPixelCol = dummy_int;
  maxPixelRow = dummy_int;
  minPixelCol = dummy_int;
  minPixelRow = dummy_int;
  geoId = dummy_uint;
  edgeHitX = dummy_int;
  edgeHitY = dummy_int;    
  clust_alpha = dummy_float;
  clust_beta = dummy_float;
}

void PixelNtuplizer_RD::PixInfoStruct::init()
{
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
     }
  */
} 

void PixelNtuplizer_RD::RecHitStruct::init()
{
  float dummy_float = -9999.0;

  localX = dummy_float;
  localY = dummy_float;
  globalX = dummy_float;
  globalY = dummy_float;
  globalZ = dummy_float;
  residualX = dummy_float;
  residualY = dummy_float;
  resErrX = dummy_float;
  resErrY = dummy_float;
  resXprime = dummy_float;
  resXprimeErr = dummy_float;  
}

void PixelNtuplizer_RD::TrackStruct::init()
{
  float dummy_float = -9999.0;

  pt = dummy_float; 
  px = dummy_float;
  py = dummy_float;
  pz = dummy_float;   
  globalEta = dummy_float;
  globalPhi = dummy_float;
  localEta = dummy_float;
  localPhi = dummy_float;
}

void PixelNtuplizer_RD::TrackerHitStruct::init()
{
  float dummy_float = -9999.0;

  globalX = dummy_float;
  globalY = dummy_float;
  globalZ = dummy_float;
}

void PixelNtuplizer_RD::CounterStruct::init()
{
  trackCounter = 0;
  pixelTrackCounter = 0;
}

// define this as a plug-in
//
DEFINE_FWK_MODULE(PixelNtuplizer_RD);
