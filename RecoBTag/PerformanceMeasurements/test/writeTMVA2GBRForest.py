import FWCore.ParameterSet.Config as cms

process = cms.Process("writeGBRForests")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1) # NB: needs to be set to 1 so that GBRForestWriter::analyze method gets called exactly once
)

process.source = cms.Source("EmptySource")

process.load('Configuration/StandardSequences/Services_cff')

process.gbrForestWriter = cms.EDAnalyzer("GBRForestWriter",
    jobs = cms.VPSet(
        cms.PSet(
            inputFileName = cms.FileInPath('RecoBTag/PerformanceMeasurements/test/TMVA_weights.xml'),
            inputFileType = cms.string("XML"),
            inputVariables = cms.vstring(['TagVarCSV_jetNTracks', 'TagVarCSV_trackSip3dSig_0', 'TagVarCSV_trackSip3dSig_1', 'TagVarCSV_trackSip3dSig_2', 'TagVarCSV_trackSip3dSig_3', 'TagVarCSV_trackSip3dSigAboveCharm', 'TagVarCSV_trackPtRel_0', 'TagVarCSV_trackPtRel_1', 'TagVarCSV_trackPtRel_2', 'TagVarCSV_trackPtRel_3', 'TagVarCSV_trackEtaRel_0', 'TagVarCSV_trackEtaRel_1', 'TagVarCSV_trackEtaRel_2', 'TagVarCSV_trackEtaRel_3', 'TagVarCSV_trackDeltaR_0', 'TagVarCSV_trackDeltaR_1', 'TagVarCSV_trackDeltaR_2', 'TagVarCSV_trackDeltaR_3', 'TagVarCSV_trackPtRatio_0', 'TagVarCSV_trackPtRatio_1', 'TagVarCSV_trackPtRatio_2', 'TagVarCSV_trackPtRatio_3', 'TagVarCSV_trackJetDist_0', 'TagVarCSV_trackJetDist_1', 'TagVarCSV_trackJetDist_2', 'TagVarCSV_trackJetDist_3', 'TagVarCSV_trackDecayLenVal_0', 'TagVarCSV_trackDecayLenVal_1', 'TagVarCSV_trackDecayLenVal_2', 'TagVarCSV_trackDecayLenVal_3', 'TagVarCSV_trackSumJetEtRatio', 'TagVarCSV_trackSumJetDeltaR', 'TagVarCSV_vertexMass', 'TagVarCSV_vertexNTracks', 'TagVarCSV_vertexEnergyRatio', 'TagVarCSV_vertexJetDeltaR', 'TagVarCSV_flightDistance2dSig', 'TagVarCSV_jetNSecondaryVertices']),
            spectatorVariables = cms.vstring(),
            methodName = cms.string("BDT"),
            outputFileType = cms.string("SQLLite"),
            outputRecord = cms.string("CSVHI_BDT_v1")
        )
    )
)

process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'sqlite_file:CSVHI_BDT_v1.db'#btag_CombinedMVAv2_BDT_TMVAv420_GBRForest_74X_v1.db'

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string('CSVHI_BDT_v1'),
            tag = cms.string('CSVHI_BDT_v1'),
            label = cms.untracked.string('CSVHI_BDT_v1'),
        )
    )
)

process.p = cms.Path(process.gbrForestWriter)
