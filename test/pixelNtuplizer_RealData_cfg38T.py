import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# DQM services
process.load("DQMServices.Core.DQM_cfg")

# Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

#fake conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
#process.GlobalTag.connect ="sqlite_file:/afs/cern.ch/user/m/malgeri/public/globtag/CRZT210_V1.db"
process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "CRAFT_V3P::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

##
## Load and Configure track selection for alignment
##
#process.load("Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi")
#process.AlignmentTrackSelector.ptMin = 3.0


##
## Load and Configure TrackRefitter
##
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitterP5.src = 'ALCARECOTkAlCosmicsCosmicTF0T'
process.TrackRefitterP5.TrajectoryInEvent = True

process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")


##
## Load and Configure OfflineValidation
##
process.load("Alignment.OfflineValidation.TrackerOfflineValidation_cfi")


process.load("RecoLocalTracker.SiPixelRecHits.pixelNtuplizerRealData_cfi")
process.PixelNtuplizer_RD.OutputFile = 'W43ALCARECOTkAlCosmicsCosmicTF0TTTree.root'
process.PixelNtuplizer_RD.trajectoryInput = 'TrackRefitterP5'


process.source = cms.Source("PoolSource",
    # replace with your files
    #lastRun = cms.untracked.uint32(64789),
    #timetype = cms.string('runnumber'),
    #firstRun = cms.untracked.uint32(64108),
    #interval = cms.uint32(1),
    fileNames = cms.untracked.vstring(

        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/061080D1-2AA2-DD11-BD61-001617DBD5AC.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/0A3A6473-02A3-DD11-8EE6-001617C3B65A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/0E028847-C1A2-DD11-8A2B-001617E30CE8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/12783584-02A3-DD11-8CF3-001D09F28F25.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/16903375-E3A2-DD11-8794-000423D94E1C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/44A2CE9D-F3A2-DD11-AEF1-001D09F23D1D.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/5C625580-02A3-DD11-8144-001D09F2924F.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/6C68553B-A0A2-DD11-8F7A-001617DC1F70.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/6E38AE71-02A3-DD11-AE98-001617E30D4A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/8227ED9D-F3A2-DD11-9586-001D09F23A20.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/A4968345-B1A2-DD11-B4EE-000423D99160.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/AA178A3C-D2A2-DD11-AE10-000423D98750.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/AAE6F486-02A3-DD11-BAAB-00304879FBB2.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/C8BEB9EE-0CA3-DD11-8FD6-0019DB29C614.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/CE024F3B-A0A2-DD11-9CC7-000423D99B3E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/CE4FBA9D-F3A2-DD11-9DB2-001D09F24D67.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/D47CB39E-F3A2-DD11-8639-0019B9F704D1.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/DA81A547-C1A2-DD11-B449-000423D6B42C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/DAAFC130-B1A2-DD11-AF83-000423D98844.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0016/F48A209F-F3A2-DD11-AB49-0019B9F72F97.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/0E486601-19A4-DD11-9B95-000423D9870C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/1AF96CC4-E7A3-DD11-9D20-001D09F290D8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/2A2AB5C2-E7A3-DD11-B3B2-001D09F25456.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/34103149-1EA3-DD11-A88D-000423D991F0.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/3A89E283-2AA4-DD11-B43E-000423D99394.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/3E104082-2AA4-DD11-ACBE-000423D94534.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/4A386415-72A3-DD11-9994-001D09F23F2A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/4C93B04C-50A3-DD11-8902-001617E30D4A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/4EA51201-F8A3-DD11-B1D6-000423D9863C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/523306BC-F8A3-DD11-8BCD-000423D98EA8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/5A0F49C9-E7A3-DD11-B6AD-001D09F24489.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/5C480B27-1FA3-DD11-8FCC-000423D996C8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/6874CEB4-3AA4-DD11-BF2E-000423D99614.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/6A0ABBD5-08A4-DD11-BFA0-000423D952C0.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/741EBDD3-08A4-DD11-8488-000423D98750.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/80781BBF-3AA4-DD11-ABF4-001617C3B6DE.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/96F0695A-50A3-DD11-8F20-00304879FA4A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/A8E714FB-18A4-DD11-9268-001617C3B5D8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/AE9B16C3-3AA4-DD11-924E-001617E30E2C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/B8157604-F8A3-DD11-BB95-000423D98BC4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/D8FD6E8D-31A4-DD11-B270-000423D991F0.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/DEFC6145-1EA3-DD11-9966-000423D95220.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/E02E6AD2-08A4-DD11-A351-000423D98E30.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/EEACA44A-61A3-DD11-B008-0019B9F730D2.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/F45B2731-61A3-DD11-B06B-0019B9F72BFF.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0017/FEB82AFA-18A4-DD11-A201-001617C3B70E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/08D27392-4BA4-DD11-B79F-0019DB29C5FC.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/145B21F9-E2A4-DD11-8DE1-000423D944DC.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/32CC0F4D-F3A4-DD11-A473-000423D98920.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/68A0CAC4-ABA5-DD11-9080-001617E30CA4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/70B4DBFC-E2A4-DD11-BC89-001617DBCF6A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/7212D29C-46A5-DD11-9CDC-001D09F251BD.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/74F49628-58A5-DD11-863D-001D09F252F3.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/96541FA0-BCA5-DD11-8879-000423D944F0.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/96F9AD7A-36A5-DD11-A632-000423D98834.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/9E45C322-F3A4-DD11-B6B9-000423D9A212.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/9EA3101C-8AA5-DD11-A579-000423D9970C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/A68D5810-5DA4-DD11-BCAE-001617C3B6C6.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/BAD90288-34A5-DD11-B969-000423D6A6F4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/C2D1E5B0-9BA5-DD11-8408-000423D99660.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/CA9C1A7B-10A5-DD11-BCD9-000423D6CA6E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/D85CEF9A-25A5-DD11-99E4-001D09F24EE3.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/D89C6219-8AA5-DD11-8A9A-000423D99AAA.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/D8E531C6-ABA5-DD11-942F-0019DB29C614.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/E2968F9B-4BA4-DD11-83BB-000423D6B5C4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0018/EE465C9E-46A5-DD11-8F1C-001D09F28C1E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/0658CC8E-EEA5-DD11-86E6-001D09F24303.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/08D7D91B-00A6-DD11-83CE-000423D6BA18.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/1001A15E-96A6-DD11-A4DB-000423D99AAE.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/14A05CD1-42A6-DD11-A429-000423D991D4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/1864E202-32A6-DD11-B35A-000423D9517C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/1A793687-A6A6-DD11-A948-000423D99AAE.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/1C7F89A3-BCA5-DD11-B571-000423D951D4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/28E8E8AF-53A6-DD11-97E3-000423D9997E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/4C99B63D-C8A6-DD11-9E17-000423D996B4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/52391DCA-D9A6-DD11-A238-0019B9F7310E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/568DBF3E-A7A6-DD11-8A43-001D09F25325.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/6ED5B5D4-42A6-DD11-B93A-0016177CA778.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/7451A83D-10A6-DD11-988F-000423D944F0.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/90FB2EA9-8EA6-DD11-ACAE-001D09F290BF.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/9207C4D6-63A6-DD11-B062-000423D6B2D8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/9AD981C8-D9A6-DD11-8EFB-0030487A18A4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/A4E16164-DEA5-DD11-8A46-00304879FA4C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/B2C45C84-CDA5-DD11-B4DA-000423D98DB4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/BC79F68A-CDA5-DD11-B311-001617E30D40.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/BE3E09C1-D5A6-DD11-8159-000423D6B2D8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/CEBF6218-21A6-DD11-991B-001D09F2532F.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/DA210EDA-63A6-DD11-AF86-000423D98B08.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/DC51D21B-B8A6-DD11-B3E6-001D09F2906A.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0019/EAB7F355-75A6-DD11-9A2C-000423D6CA02.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/160104AD-0BA7-DD11-BAA4-000423D99896.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/22215EB1-1CA7-DD11-B283-001D09F276CF.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/2C1B455A-2EA7-DD11-B1E7-000423D94AA8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/32FC868C-FBA6-DD11-89B9-001D09F2503C.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/4246AFF2-1CA7-DD11-9CCF-0019B9F704D6.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/6C39508D-FBA6-DD11-8920-001D09F23A4D.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/C45A27AF-0BA7-DD11-80CA-001617E30CE8.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/D4948CFE-E9A6-DD11-8060-001617C3B76E.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/DC0B9FFE-E9A6-DD11-98C9-000423D98DC4.root',
        '/store/data/Commissioning08/Cosmics/ALCARECO/CRAFT_V3P_StreamALCARECOTkAlCosmics0T_v11/0020/E488569C-28A7-DD11-8110-000423D6CAF2.root'

    ) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) )

process.p = cms.Path(process.offlineBeamSpot*process.TrackRefitterP5*process.PixelNtuplizer_RD)
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = 'Info'
process.TrackerDigiGeometryESModule.applyAlignment = True
