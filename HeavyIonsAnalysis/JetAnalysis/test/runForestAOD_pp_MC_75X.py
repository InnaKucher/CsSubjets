### HiForest Configuration
# Collisions: pp
# Type: MC
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
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
        #"root://eoscms.cern.ch//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch2/TTbar/RECO/Events_1.root"
                            #    'file:/afs/cern.ch/work/k/kjung/holdingCell/Pythia6_pp_HLTReco.root'
                           #'/store/himc/HINppWinter16DR/Pythia8_bJet120_pp502_TuneCUETP8M1/AODSIM/75X_mcRun2_asymptotic_ppAt5TeV_v3-v1/80000/02C10805-D6E4-E511-9C81-0CC47A4D9A10.root'), 
			   '/store/himc/HINppWinter16DR/Pythia8_Dijet120_pp502_TuneCUETP8M1/AODSIM/75X_mcRun2_asymptotic_ppAt5TeV_v3-v1/60000/1AC3D4FC-2DF0-E511-96CA-0CC47A04D082.root') 
#			   eventsToProcess = cms.untracked.VEventRange('1:42101','1:52848')
)

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000))


#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_asymptotic_ppAt5TeV_v3', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag


# Customization
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pp5020
process = overrideJEC_pp5020(process)

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

process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_nominalPP")

# Include this to turn on storing the jet constituents and new jet variables for q/g separation
#process.ak4PFJetAnalyzer.doJetConstituents = cms.untracked.bool(True)
#process.ak4PFJetAnalyzer.doNewJetVars = cms.untracked.bool(True)
# Use this version for JEC
#process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_JECPP")

#####################################################################################

############################
# Event Analysis
############################
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi') #use data version to avoid PbPb MC
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.hiEvtAnalyzer.doMC = cms.bool(True) #general MC info
process.hiEvtAnalyzer.doHiMC = cms.bool(False) #HI specific MC info

process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
process.HiGenParticleAna.genParticleSrc = cms.untracked.InputTag("genParticles")
process.HiGenParticleAna.doHI = False
process.HiGenParticleAna.stableOnly = False
process.HiGenParticleAna.etaMax = cms.untracked.double(999)
process.HiGenParticleAna.ptMin = cms.untracked.double(0)
process.load('HeavyIonsAnalysis.EventAnalysis.runanalyzer_cff')
process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_pp_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.pfcandAnalyzer.pfCandidateLabel = cms.InputTag("particleFlow")
process.pfcandAnalyzer.doVS = cms.untracked.bool(False)
process.pfcandAnalyzer.doUEraw_ = cms.untracked.bool(False)
process.pfcandAnalyzer.genLabel = cms.InputTag("genParticles")

#####################################################################################

#########################
# Track Analyzer
#########################
process.load('HeavyIonsAnalysis.JetAnalysis.ExtraTrackReco_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_cff')

# Use this instead for track corrections
## process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_Corr_cff')

#####################################################################################

#####################
# photons
######################
process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.gsfElectronLabel   = cms.InputTag("gedGsfElectrons")
process.ggHiNtuplizer.recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerpp')
process.ggHiNtuplizer.VtxLabel           = cms.InputTag("offlinePrimaryVertices")
process.ggHiNtuplizer.particleFlowCollection = cms.InputTag("particleFlow")
process.ggHiNtuplizer.doVsIso            = cms.bool(False)
process.ggHiNtuplizer.doElectronVID      = cms.bool(True)
process.ggHiNtuplizerGED = process.ggHiNtuplizer.clone(recoPhotonSrc = cms.InputTag('gedPhotons'),
                                                       recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerppGED'))

####################################################################################
#####################
# Electron ID
#####################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be processed
# DataFormat.AOD or DataFormat.MiniAOD
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce. Check here https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_7_4
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
####################################################################################


#####################
# tupel and necessary PAT sequences
#####################

process.load("HeavyIonsAnalysis.VectorBosonAnalysis.tupelSequence_pp_mc_cff")

#####################################################################################

#########################
# Main analysis list
#########################
process.ana_step = cms.Path(#process.hltanalysis *
                            process.hiEvtAnalyzer *
                            process.HiGenParticleAna*
                            process.jetSequences +
                            #process.egmGsfElectronIDSequence + #Should be added in the path for VID module
                            #process.ggHiNtuplizer +
                            #process.ggHiNtuplizerGED +
                            process.pfcandAnalyzer +
                            process.HiForest +
			    process.trackSequencesPP +
                            process.runAnalyzer 
  #                          process.tupelPatSequence
)

#####################################################################################

#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

process.NoScraping = cms.EDFilter("FilterOutScraping",
 applyfilter = cms.untracked.bool(True),
 debugOn = cms.untracked.bool(False),
 numtrack = cms.untracked.uint32(10),
 thresh = cms.untracked.double(0.25)
)

process.pPAprimaryVertexFilter = cms.Path(process.PAprimaryVertexFilter)
process.pBeamScrapingFilter=cms.Path(process.NoScraping)

process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")

process.pVertexFilterCutG = cms.Path(process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.skimanalysis)

# Customization
process.akSoftDrop4PFPatJetFlavourAssociation.jets="ak4PFJets"
process.akSoftDrop4PFPatJetFlavourAssociation.groomedJets=cms.InputTag("akSoftDrop4PFJets")
process.akSoftDrop4PFPatJetFlavourAssociation.subjets= cms.InputTag('akSoftDrop4PFJets','SubJets')
process.akSoftDrop4PFJets.useSoftDrop = True
process.akSoftDrop4PFpatJetsWithBtagging.getJetMCFlavour = cms.bool(False)
process.akSoftDrop4PFJetAnalyzer.doExtendedFlavorTagging = cms.untracked.bool(True)
process.akSoftDrop4PFJetAnalyzer.jetFlavourInfos    = cms.InputTag("akSoftDrop4PFPatJetFlavourAssociation")
process.akSoftDrop4PFJetAnalyzer.subjetFlavourInfos = cms.InputTag("akSoftDrop4PFPatJetFlavourAssociation","SubJets")
process.akSoftDrop4PFJetAnalyzer.groomedJets        = cms.InputTag("akSoftDrop4PFJets")
process.akSoftDrop4PFJetAnalyzer.isPythia6 = cms.untracked.bool(False)

#old-school dR association is no good for subjets... get lots of tracks associated to both subjets
#process.akSoftDrop4PFSubjetJetTracksAssociatorAtVertex = process.akSoftDrop4PFJetTracksAssociatorAtVertex.clone()
#process.akSoftDrop4PFSubjetJetTracksAssociatorAtVertex.jets = cms.InputTag('akSoftDrop4PFJets','SubJets')
#process.akSoftDrop4PFSubjetJetTracksAssociatorAtVertex.coneSize = cms.double(0.25)

#turn on explicit subjet-track association
process.akSoftDrop4PFSubjetJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorExplicit",
	jets = cms.InputTag('akSoftDrop4PFJets','SubJets'),
	tracks = cms.InputTag('highPurityTracks')
)

process.akSoftDrop4PFSubjetImpactParameterTagInfos = process.akSoftDrop4PFImpactParameterTagInfos.clone()
process.akSoftDrop4PFSubjetImpactParameterTagInfos.jetTracks = cms.InputTag("akSoftDrop4PFSubjetJetTracksAssociatorAtVertex")
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos = process.akSoftDrop4PFSecondaryVertexTagInfos.clone()
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag('akSoftDrop4PFSubjetImpactParameterTagInfos')

process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.fatJets = cms.InputTag('ak4PFJets')
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.groomedFatJets = cms.InputTag('akSoftDrop4PFJets')
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.useSVClustering = cms.bool(False)
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.useExternalSV = cms.bool(True)
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.jetAlgorithm = cms.string('AntiKt')
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.useSVMomentum = cms.bool(True)
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.rParam = cms.double(0.4)
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.extSVCollection = cms.InputTag('inclusiveSecondaryVertices')
process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = cms.double(0.1)

process.akSoftDrop4PFJetAnalyzer.trackSelection = process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.trackSelection
process.akSoftDrop4PFJetAnalyzer.trackPairV0Filter = process.akSoftDrop4PFSubjetSecondaryVertexTagInfos.vertexCuts.v0Filter

process.akSoftDrop4PFCombinedSubjetSecondaryVertexBJetTags = process.akSoftDrop4PFCombinedSecondaryVertexBJetTags.clone(
	tagInfos = cms.VInputTag(cms.InputTag("akSoftDrop4PFSubjetImpactParameterTagInfos"),
                cms.InputTag("akSoftDrop4PFSubjetSecondaryVertexTagInfos"))
)
#process.akSoftDrop4PFCombinedSubjetSecondaryVertexV2BJetTags = process.akSoftDrop4PFCombinedSecondaryVertexV2BJetTags.clone(
#        tagInfos = cms.VInputTag(cms.InputTag("akSoftDrop4PFSubjetImpactParameterTagInfos"),
#                cms.InputTag("akSoftDrop4PFSubjetSecondaryVertexTagInfos"))
#)

process.akSoftDrop4PFJetBtaggingSV *= process.akSoftDrop4PFSubjetJetTracksAssociatorAtVertex+process.akSoftDrop4PFSubjetImpactParameterTagInfos+process.akSoftDrop4PFSubjetSecondaryVertexTagInfos+process.akSoftDrop4PFCombinedSubjetSecondaryVertexBJetTags
##+process.akSoftDrop4PFCombinedSubjetSecondaryVertexV2BJetTags

#process.printEventAKSoftDrop4PFJets = cms.EDAnalyzer("printJetFlavourInfo",
#                                                     jetFlavourInfos    = cms.InputTag("akSoftDrop4PFPatJetFlavourAssociation"),
#                                                     subjetFlavourInfos = cms.InputTag("akSoftDrop4PFPatJetFlavourAssociation","SubJets"),
#                                                     groomedJets        = cms.InputTag("akSoftDrop4PFJets"),
#                                                     )


#process.ana_step *= process.printEventAKSoftDrop4PFJets

