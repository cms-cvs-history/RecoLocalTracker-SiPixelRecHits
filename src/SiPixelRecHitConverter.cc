/** SiPixelRecHitConverter.cc
 * -------------------------------------------- 
 * Description:  see SiPixelRecHitConverter.h
 * Author:  P. Maksimovic, JHU
 * History: Feb 27, 2006, initial version
 * -------------------------------------------- 
 */


// Our own stuff
#include "RecoLocalTracker/SiPixelRecHits/interface/SiPixelRecHitConverter.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/CPEFromDetPosition.h"

// Geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
class GeometricDet;   // hack in 0.2.0pre5, OK for pre6 -- still needed?
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

// Data Formats
#include "DataFormats/SiPixelCluster/interface/SiPixelClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/DetId/interface/DetId.h"

// Framework
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

// STL
#include <vector>
#include <memory>
#include <string>
#include <iostream>

// MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Magnetic Field
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

namespace cms
{

  //---------------------------------------------------------------------------
  //!  Constructor: set the ParameterSet and defer all thinking to setupCPE().
  //---------------------------------------------------------------------------
  SiPixelRecHitConverter::SiPixelRecHitConverter(edm::ParameterSet const& conf) 
    : 
    conf_(conf),
    cpeName_("None"),     // bogus
    cpe_(0),              // the default, in case we fail to make one
    ready_(false)         // since we obviously aren't
  {
    //--- Declare to the EDM what kind of collections we will be making.
    produces<SiPixelRecHitCollection>();

    //--- Make the algorithm(s) according to what the user specified
    //--- in the ParameterSet.
    //setupCPE();
  }


  // Virtual destructor needed, just in case.
  SiPixelRecHitConverter::~SiPixelRecHitConverter() { }  

  //---------------------------------------------------------------------------
  // Begin job: get magnetic field
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::beginJob(const edm::EventSetup& c) {
    edm::ESHandle<MagneticField> magfield;
    c.get<IdealMagneticFieldRecord>().get(magfield);
    setupCPE(magfield.product());
  }

  //---------------------------------------------------------------------------
  //! The "Event" entrypoint: gets called by framework for every event
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::produce(edm::Event& e, const edm::EventSetup& es)
  {
    // retrieve producer name of input SiPixelClusterCollection
    std::string clusterCollLabel = conf_.getUntrackedParameter<std::string>("ClusterCollLabel","pixClust");

    // Step A.1: get input data
    edm::Handle<SiPixelClusterCollection> clusterColl;
    e.getByLabel(clusterCollLabel, clusterColl);

    // Step A.2: get event setup
    edm::ESHandle<TrackerGeometry> geom;
    es.get<TrackerDigiGeometryRecord>().get( geom );


    // Step B: create empty output collection
    std::auto_ptr<SiPixelRecHitCollection> output(new SiPixelRecHitCollection);

    // Step C: Iterate over DetIds and invoke the strip CPE algorithm
    // on each DetUnit
    run( clusterColl.product(), *output, geom );

    // Step D: write output to file
    e.put(output);

  }



  //---------------------------------------------------------------------------
  //!  Set up the specific algorithm we are going to use.  
  //!  TO DO: in the future, we should allow for a different algorithm for 
  //!  each detector subset (e.g. barrel vs forward, per layer, etc).
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::setupCPE(const MagneticField* mag) 
  {
    cpeName_ = conf_.getParameter<std::string>("CPE");

    if ( cpeName_ == "FromDetPosition" ) {
      cpe_ = new CPEFromDetPosition(conf_,mag);
      ready_ = true;
    } 
    else {
      edm::LogError("SiPixelRecHitConverter") 
		<<" Cluster parameter estimator " << cpeName_ << " is invalid.\n"
		<< "Possible choices:\n" 
		<< "    - FromDetPosition";
      ready_ = false;
    }
  }


  //---------------------------------------------------------------------------
  //!  Iterate over DetUnits, then over Clusters and invoke the CPE on each,
  //!  and make a RecHit to store the result.
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::run(const SiPixelClusterCollection* input, 
				   SiPixelRecHitCollection &output,
				   edm::ESHandle<TrackerGeometry> & geom)
  {
    if ( ! ready_ ) {
      edm::LogError("SiPixelRecHitConverter") << " at least one CPE is not ready -- can't run!";
      // TO DO: throw an exception here?  The user may want to know...
      assert(0);
      return;   // clusterizer is invalid, bail out
    }

    int numberOfDetUnits = 0;
    int numberOfClusters = 0;
    
    // get vector of detunit ids
    const std::vector<unsigned int> detIDs = input->detIDs();
    

    //--- Loop over detunits.
    std::vector<unsigned int>::const_iterator 
      detunit_it  = detIDs.begin(),
      detunit_end = detIDs.end();
    for ( ; detunit_it != detunit_end; ++detunit_it ) {
      //
      unsigned int detid = *detunit_it;
      ++numberOfDetUnits;
      const SiPixelClusterCollection::Range clustRange = input->get(detid);
      SiPixelClusterCollection::ContainerIterator 
	   clustIt = clustRange.first;
      SiPixelClusterCollection::ContainerIterator 
	endClustIt = clustRange.second;

      // TO DO: if we were to allow concurrent clusterizers, we would
      // be picking one from the map and running it.

      // geometry information for this DetUnit. TrackerGeom:DetContainer is
      // a vector<GeomDetUnit*>.  We need to call 
      //     const GeomDet&  TrackerGeom::idToDet(DetId) const;
      // to get a GeomDet (a GeomDetUnit*) from a DetId.

      // convert detid (unsigned int) to DetId
      DetId detIdObject( detid );      

      const GeomDetUnit * genericDet = geom->idToDetUnit( detIdObject );

      const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(genericDet);
      assert(pixDet);   // ensures that geoUnit != 0 too
      // &&& Need to throw an exception instead!

      //--- Vector for temporary storage.  Need to do it this way
      //--- since put() method will map the range to DetId.
      //old SiPixelRecHitCollection::Container recHitsOnDetUnit;  
       edm::OwnVector<SiPixelRecHit> recHitsOnDetUnit;  

      // &&& We really should preallocate the size of this container since
      // &&& we *know* how many PixelClusters there are on this det unit!

       //
       // Not needed anymore: setTheDet is called directly from localPosition
       // and localError
       //      
       //if (dynamic_cast<CPEFromDetPosition*>(cpe_))
       // (dynamic_cast<CPEFromDetPosition*>(cpe_))->setTheDet( *pixDet );  // &&& not in the base class,

      for ( ; clustIt != endClustIt; ++clustIt ) {
	std::pair<LocalPoint,LocalError> lv = 
	  cpe_->localParameters( *clustIt, *genericDet );
	LocalPoint lp( lv.first );
	LocalError le( lv.second );
	
	// Make a RecHit and add it to a temporary OwnVector.
	recHitsOnDetUnit.push_back( new SiPixelRecHit( lp, le, detIdObject, &*clustIt) );
      }

      //--- At the end of this det unit, move all hits to the
      //--- real PixelRecHitColl.
      //old SiPixelRecHitCollection::Range inputRange;
      //old inputRange.first = recHitsOnDetUnit.begin();
      //old inputRange.second = recHitsOnDetUnit.end();
      //old output.put(inputRange, detid);

      if (recHitsOnDetUnit.size() > 0) {
	output.put(detIdObject, 
		   recHitsOnDetUnit.begin(), recHitsOnDetUnit.end());
	LogDebug("SiPixelRecHitConverter") << " Found " 
					   << recHitsOnDetUnit.size() << " RecHits on" << detid;	
      }
      // numberOfRecHits += recHitsOnDetUnit.size();
    }
    // end of the loop over detunits
    
    LogDebug ("SiPixelRecHitConverter") 
	      << cpeName_ << " converted " << numberOfClusters 
	      << " SiPixelClusters into SiPixelRecHits, in " << numberOfDetUnits << " DetUnits."; 
	   
  }


  


}  // end of namespace cms
