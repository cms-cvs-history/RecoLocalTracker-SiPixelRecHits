import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# DQM services
process.load("DQMServices.Core.DQM_cfg")

# Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

#fake conditions
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
#process.GlobalTag.connect ="sqlite_file:/afs/cern.ch/user/m/malgeri/public/globtag/CRZT210_V1.db"
process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "CRUZET4_V2::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

##
## Load and Configure track selection for alignment
##
#process.load("Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi")
#process.AlignmentTrackSelector.ptMin = 3.0


##
## Load and Configure TrackRefitter
##
process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff")
process.TrackRefitter.src = 'ALCARECOTkAlCosmicsCTF0T'
process.TrackRefitter.TrajectoryInEvent = True

process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")


##
## Load and Configure OfflineValidation
##
process.load("Alignment.OfflineValidation.TrackerOfflineValidation_cfi")


process.load("RecoLocalTracker.SiPixelRecHits.pixelNtuplizer_cfi")
process.PixelNtuplizer.OutputFile = 'AlCa1TTree.root'


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_0_0_pre2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0001/9C36DFC0-1CB4-DD11-9D8B-001617DBD288.root'
    ) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) )

process.p = cms.Path(process.offlineBeamSpot*process.TrackRefitter*process.PixelNtuplizer)
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = 'Info'
process.TrackerDigiGeometryESModule.applyAlignment = True
