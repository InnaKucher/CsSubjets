### HiForest Configuration
# Collisions: pp
# Type: Data
# Input: AOD

import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

#####################################################################################
# HiForest labelling info
#####################################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V3",)
import subprocess
version = subprocess.Popen(["(cd $CMSSW_BASE/src && git describe --tags)"], stdout=subprocess.PIPE, shell=True).stdout.read()
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #"/store/group/phys_heavyions/velicanu/reco/HIPhysicsMinBiasUPC/v0/000/262/548/recoExpress_84.root"
                            '/store/hidata/HIRun2015/HIHardProbes/AOD/PromptReco-v1/000/262/548/00000/9091A3D3-AD9A-E511-8B73-02163E011A3E.root'
				)
)


# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100))


#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v13', '')  #for now track GT manually, since centrality tables updated ex post facto
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_PbPb5020
process = overrideJEC_PbPb5020(process)

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("HiForestAOD.root"))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

####################################################################################

#############################
# Jets
#############################

#require the pu algo to use a certain threshold of towers for bg subtraction
process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_puLimitedDataPbPb")
#or don't do that
#process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_puUnlimitedDataPbPb")

#####################################################################################

############################
# Event Analysis
############################
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_PbPb_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import addHLTdummybranches
addHLTdummybranches(process)

process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzerCS_cfi")
process.pfcandAnalyzerCS.skipCharged = False
process.pfcandAnalyzerCS.pfPtMin = 0
process.load("HeavyIonsAnalysis.JetAnalysis.hcalNoise_cff")

#####################################################################################

#########################
# Track Analyzer
#########################
process.load('HeavyIonsAnalysis.JetAnalysis.ExtraTrackReco_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_cff')
# process.load("HeavyIonsAnalysis.TrackAnalysis.METAnalyzer_cff")


####################################################################################

#####################
# Photons
#####################
process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doGenParticles = False
process.ggHiNtuplizerGED = process.ggHiNtuplizer.clone(recoPhotonSrc = cms.InputTag('gedPhotonsTmp'),
                                                       recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerGED')
)


####################################################################################

#####################
# tupel and necessary PAT sequences
#####################

process.load("HeavyIonsAnalysis.VectorBosonAnalysis.tupelSequence_PbPb_cff")

#####################################################################################

#replace pp CSVv2 with PbPb CSVv2 (positive and negative taggers unchanged!)
process.load('RecoBTag.CSVscikit.csvscikitTagJetTags_cfi')
process.load('RecoBTag.CSVscikit.csvscikitTaggerProducer_cfi')
process.akPu4PFCombinedSecondaryVertexV2BJetTags = process.pfCSVscikitJetTags.clone()
process.akPu4PFCombinedSecondaryVertexV2BJetTags.tagInfos=cms.VInputTag(cms.InputTag("akPu4PFImpactParameterTagInfos"), cms.InputTag("akPu4PFSecondaryVertexTagInfos"))

process.akCs4PFCombinedSecondaryVertexV2BJetTags = process.pfCSVscikitJetTags.clone()
process.akCs4PFCombinedSecondaryVertexV2BJetTags.tagInfos=cms.VInputTag(cms.InputTag("akCs4PFImpactParameterTagInfos"), cms.InputTag("akCs4PFSecondaryVertexTagInfos"))

process.akPu4CaloCombinedSecondaryVertexV2BJetTags = process.pfCSVscikitJetTags.clone()
process.akPu4CaloCombinedSecondaryVertexV2BJetTags.tagInfos=cms.VInputTag(cms.InputTag("akPu4CaloImpactParameterTagInfos"), cms.InputTag("akPu4CaloSecondaryVertexTagInfos"))
process.CSVscikitTags.weightFile=cms.FileInPath('HeavyIonsAnalysis/JetAnalysis/data/bTagCSVv2PbPb_758p3_Jan2017_BDTG_weights.xml')

#########################
# Main analysis list
#########################

process.ana_step = cms.Path(
#			    process.hltanalysis *
#			    process.hltobject *
                            process.centralityBin *
                            process.hiEvtAnalyzer*
                            process.jetSequences +
#                            process.ggHiNtuplizer +
#                            process.ggHiNtuplizerGED +
                            process.pfcandAnalyzer +
                            process.pfcandAnalyzerCS +
                            process.HiForest +
#                            process.trackSequencesPbPb +
                            process.hcalNoise #+
                            #process.tupelPatSequence
                            )

#####################################################################################

#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pcollisionEventSelection = cms.Path(process.collisionEventSelectionAOD)
process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter )

process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
process.phfCoincFilter1 = cms.Path(process.hfCoincFilter)
process.phfCoincFilter2 = cms.Path(process.hfCoincFilter2)
process.phfCoincFilter3 = cms.Path(process.hfCoincFilter3)
process.phfCoincFilter4 = cms.Path(process.hfCoincFilter4)
process.phfCoincFilter5 = cms.Path(process.hfCoincFilter5)

process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)

process.pAna = cms.EndPath(process.skimanalysis)

# Customization
##########################################UE##########################################
process.ak4PFJets.jetPtMin = cms.double(0.0)
process.akSoftDrop4PFJets.src = cms.InputTag("particleFlowTmp")
process.akCsSoftDrop4PFJets.jetPtMin = cms.double(0.0)

process.akCsSoftDrop4PFpatJetsWithBtagging.getJetMCFlavour = cms.bool(False)
process.akCsSoftDrop4PFJetAnalyzer.doExtendedFlavorTagging = cms.untracked.bool(False)
process.akCsSoftDrop4PFJetAnalyzer.jetFlavourInfos    = cms.InputTag("akCsSoftDrop4PFPatJetFlavourAssociation")
process.akCsSoftDrop4PFJetAnalyzer.subjetFlavourInfos = cms.InputTag("akCsSoftDrop4PFPatJetFlavourAssociation","SubJets")
process.akCsSoftDrop4PFJetAnalyzer.groomedJets        = cms.InputTag("akCsSoftDrop4PFJets")
process.akCsSoftDrop4PFJetAnalyzer.isPythia6 = cms.untracked.bool(False)

#process.akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex = process.akCsSoftDrop4PFJetTracksAssociatorAtVertex.clone()
#process.akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex.jets = cms.InputTag('akCsSoftDrop4PFJets','SubJets')
#process.akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex.coneSize = cms.double(0.25)

process.akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorExplicit",
        jets = cms.InputTag('akCsSoftDrop4PFJets','SubJets'),
        tracks = cms.InputTag('highPurityTracks')
)

process.akCsSoftDrop4PFSubjetImpactParameterTagInfos = process.akCsSoftDrop4PFImpactParameterTagInfos.clone()
process.akCsSoftDrop4PFSubjetImpactParameterTagInfos.jetTracks = cms.InputTag("akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos = process.akCsSoftDrop4PFSecondaryVertexTagInfos.clone()
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag('akCsSoftDrop4PFSubjetImpactParameterTagInfos')

#doing ghost-vertex reclustering for subjets
#starting with IVF vertexing
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
process.inclusiveVertexFinder.primaryVertices = cms.InputTag("hiSelectedVertex")
process.inclusiveVertexFinder.tracks = cms.InputTag("highPurityTracks")
process.trackVertexArbitrator.primaryVertices = cms.InputTag("hiSelectedVertex")
process.trackVertexArbitrator.tracks = cms.InputTag("highPurityTracks")
process.inclusiveSecondaryVertices.primaryVertices= cms.InputTag("hiSelectedVertex")
process.inclusiveSecondaryVertices.tracks = cms.InputTag("highPurityTracks")

process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.useExternalSV = cms.bool(True)
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.extSVCollection = cms.InputTag("inclusiveSecondaryVertices")
#process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.extSVCollection = cms.InputTag("hiSelectedVertex")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.extSVDeltaRToJet = cms.double(1.0) #make this big to make sure ghost vertexing works 
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.useSVClustering = cms.bool(True)
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.useSVMomentum = cms.bool(True)
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.fatJets = cms.InputTag("ak4PFJets")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.groomedFatJets = cms.InputTag("akCsSoftDrop4PFJets")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.jetAlgorithm = cms.string("AntiKt")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.rParam = cms.double(0.4)
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = cms.double(0.2)

process.akCsSoftDrop4PFJetAnalyzer.trackSelection = process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.trackSelection
process.akCsSoftDrop4PFJetAnalyzer.trackPairV0Filter = process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.vertexCuts.v0Filter

process.akCsSoftDrop4PFCombinedSubjetSecondaryVertexBJetTags = process.akCsSoftDrop4PFCombinedSecondaryVertexBJetTags.clone(
        tagInfos = cms.VInputTag(cms.InputTag("akCsSoftDrop4PFSubjetImpactParameterTagInfos"),
                cms.InputTag("akCsSoftDrop4PFSubjetSecondaryVertexTagInfos"))
)
process.akCsSoftDrop4PFJetBtaggingSV *= process.akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex+process.akCsSoftDrop4PFSubjetImpactParameterTagInfos+process.akCsSoftDrop4PFSubjetJetProbabilityBJetTags+process.inclusiveVertexing+process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos+process.akCsSoftDrop4PFCombinedSubjetSecondaryVertexBJetTags

######################################################################################
from CondCore.DBCommon.CondDBSetup_cfi import *
process.uetable = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
          cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("UETableCompatibilityFormat_PF_v02_offline"),
                   label = cms.untracked.string("UETable_PF")
          ),
          cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("UETableCompatibilityFormat_Calo_v02_offline"),
                   label = cms.untracked.string("UETable_Calo")
          )
      ), 
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
)
process.es_prefer_uetable = cms.ESPrefer('PoolDBESSource','uetable')
##########################################UE##########################################
