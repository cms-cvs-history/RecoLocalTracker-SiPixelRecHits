#ifndef PixelNtuplizer_RD_h
#define PixelNtuplizer_RD_h

/** \class PixelNtuplizer_RealData
 *
 *
 ************************************************************/
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
#include <memory>
#include <vector>


#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h"

class MagneticField;
class TrackerGeometry;


//___________________________________________________________

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"


class TTree;
class TFile;
class RectangularPixelTopology;

class PixelNtuplizer_RD : public edm::EDAnalyzer
{
 public:
  
  explicit PixelNtuplizer_RD(const edm::ParameterSet& ps);
  virtual ~PixelNtuplizer_RD();
  virtual void beginJob(const edm::EventSetup& iSetup);
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

 protected:
  void fillEvt(const edm::Event&);
  void fillDet(const DetId &, uint, const PixelGeomDetUnit*);
  void fillVertex(const PixelGeomDetUnit*);
  void fillClust(const SiPixelCluster&, const RectangularPixelTopology*, const PixelGeomDetUnit*);
  void fillPix(const SiPixelCluster&, const RectangularPixelTopology*, const PixelGeomDetUnit*);
  void fillTrack(TrajectoryStateOnSurface&);
  
 private:
  edm::ParameterSet conf_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  edm::ESHandle<MagneticField> magneticField_;

  int TrackCounter;
  int PixelTrackCounter;

  TFile* tfile_;
  TTree* t_;  // tree filled on every pixel rec hit
  TTree* ts_; // tree filled on every strip rec hit
  TTree* tc_; // tree filled every job (eventually every event)
  
  void init();
  
  //--- Structures for ntupling:

  struct EvtStruct {

    int run;
    int evtnum;

    void init();
  } evt_;
  
  struct DetStruct {

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

  struct VertexStruct {

    float r;
    float z;

    void init();
  } vertex_;

  struct ClusterStruct {

    float row;
    float col;
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
    uint32_t geoId;
    int edgeHitX;
    int edgeHitY;    
    float clust_alpha; // alpha from cluster position w.r.t detector center
    float clust_beta;  // beta from cluster position w.r.t detector center

    void init();
  } clust_;

  struct PixInfoStruct {

    int npix;
    float row[100];
    float col[100];
    float adc[100];
    // Just added
    float x[100];
    float y[100];
    float gx[100];
    float gy[100];
    float gz[100];

    void init();
  } pixinfo_;

  struct RecHitStruct{

    float localX;
    float localY;
    float globalX;
    float globalY;
    float globalZ;
    float residualX;
    float residualY;
    float resErrX;
    float resErrY;
    float resXprime;
    float resXprimeErr;

    void init();
    } rechit_;

  struct TrackStruct{

    float pt;
    float px;
    float py;
    float pz;
    float globalEta;
    float globalPhi;
    float localEta;
    float localPhi;
    /*    float kappa;
    float chi2;
    float normchi2;
    float d0;
    float dz; */

    void init();
    } track_;

  struct TrackerHitStruct{

    float globalX;
    float globalY;
    float globalZ;

    void init();
    } trackerhits_;

  struct CounterStruct{

    int trackCounter;
    int pixelTrackCounter;

    void init();
    } counters_;

};

#endif
