### HiForest Configuration
# Collisions: PbPb
# Type: MonteCarlo
# Input: AOD

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')

options.setDefault('maxEvents', 1)
options.setDefault('outputFile', 'JetTree_qcd_PuJets')
options.setDefault('inputFiles',['file:qcd_99.root'])

options.parseArguments()

print('Running on files ',options.inputFiles)
print('Output file ',options.outputFile)
print('N events ',options.maxEvents)



process = cms.Process('HiForest')
process.options = cms.untracked.PSet()


process.Timing = cms.Service("Timing")

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
                            fileNames = cms.untracked.vstring(options.inputFiles
                               # "file:fcr_99.root"#samples/PbPb_MC_RECODEBUG.root"
                                )
                            )

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v13', '') #for now track GT manually, since centrality tables updated ex post facto
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_PbPb5020
process = overrideJEC_PbPb5020(process)

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
# process.centralityBin.Centrality = cms.InputTag("hiCentrality")
# process.centralityBin.centralityVariable = cms.string("HFtowers")
#process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

####################################################################################

#############################
# Jets
#############################
# full gen jets followed by filters to select signal-only genjets
process.load('HeavyIonsAnalysis.JetAnalysis.GenJetSequence')
process.load('HeavyIonsAnalysis.JetAnalysis.hiSignalGenFilters')


# nominal jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.FullJetSequence_nominalPbPb')
# replace above with this one for JEC:
#process.load('HeavyIonsAnalysis.JetAnalysis.FullJetSequence_JECPbPb')

#rho analyzer
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJRhoAnalyzer_cff')

####################################################################################

#############################
# Gen Analyzer
#############################
process.load('HeavyIonsAnalysis.EventAnalysis.HiMixAnalyzerRECO_cff')
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.runanalyzer_cff')
process.HiGenParticleAna.genParticleSrc = cms.untracked.InputTag("genParticles")
# Temporary disactivation - until we have DIGI & RECO in CMSSW_7_5_7_patch4
process.HiGenParticleAna.doHI = False


#####################################################################################

############################
# Event Analysis
############################
# process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
# process.hiEvtAnalyzer.doMC = cms.bool(True) #general MC info
# process.hiEvtAnalyzer.doHiMC = cms.bool(True) #HI specific MC info
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzerCS_cfi")
process.pfcandAnalyzerCS.skipCharged = False
process.pfcandAnalyzerCS.pfPtMin = 0

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
# Photons
#####################

process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizerGED = process.ggHiNtuplizer.clone(recoPhotonSrc = cms.InputTag('gedPhotonsTmp'),
                                                       recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerGED')
)

#####################################################################################


#####################
# tupel and necessary PAT sequences
#####################

process.load("HeavyIonsAnalysis.VectorBosonAnalysis.tupelSequence_PbPb_mc_cff")

#####################################################################################

#########################
# Main analysis list
#########################

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("outforest.root"),
    outputCommands = cms.untracked.vstring(['keep *'])
)

#process.pAna = cms.EndPath(process.out)#process.skimanalysis)

## Postfix                                                                                                                   
postfix = "PFlow"
## Various collection names                                                                                                  
genParticles = 'genParticles'
jetSource = 'pfJetsPFBRECO'+postfix
genJetCollection = 'ak4GenJetsNoNu'+postfix
pfCandidates = 'particleFlowTmp'
pvSource = 'offlinePrimaryVertices'
svSource = 'inclusiveCandidateSecondaryVertices'
muSource = ''#muons'                                                                                                         
elSource = 'gedGsfElectrons'
patMuons = ''#selectedPatMuons'                                                                                              
trackSource = 'generalTracks'

## Add GenParticlePruner for boosted b-tagging studies
#if not options.runOnData:
process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
                                                     src = cms.InputTag(genParticles),
                                                     select = cms.vstring("drop  *  ", #by default
                                                                          "keep ( status = 3 || (status>=21 && status<=29) )", #keep hard process particles
                                                                          "keep abs(pdgId) = 13 || abs(pdgId) = 15" #keep muons and taus
                                                                          )
                                                     )

#-------------------------------------
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *

from RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff import *
process.btagana = bTagAnalyzer.clone()

process.btagana.tracksColl            = cms.InputTag(trackSource)
process.btagana.useSelectedTracks     = True  ## False if you want to run on all tracks : for commissioning studies          
process.btagana.useTrackHistory       = False ## Can only be used with GEN-SIM-RECODEBUG files                               
process.btagana.fillsvTagInfo         = False ## True if you want to store information relative to the svTagInfos, set to Fa\lse if produceJetTrackTree is set to False                                                                                   
process.btagana.produceJetTrackTree   = False ## True if you want to keep info for tracks associated to jets : for commissio\ning studies                                                                                                                 
process.btagana.produceAllTrackTree   = False ## True if you want to keep info for all tracks : for commissioning studies    
process.btagana.producePtRelTemplate  = False ## options.producePtRelTemplate  ## True for performance studies                        
#------------------                                                                                                          
process.btagana.storeTagVariables     = False ## True if you want to keep TagInfo TaggingVariables                           
process.btagana.storeCSVTagVariables  = True  ## True if you want to keep CSV TaggingVariables                               
process.btagana.primaryVertexColl     = cms.InputTag(pvSource)
process.btagana.Jets                  = cms.InputTag('akPu4PFpatJetsWithBtagging')#selectedPatJets')#+postfix)
process.btagana.muonCollectionName    = cms.InputTag(muSource)
process.btagana.patMuonCollectionName = cms.InputTag(patMuons)
process.btagana.use_ttbar_filter      = cms.bool(False)#options.useTTbarFilter)
process.btagana.triggerTable          = cms.InputTag('TriggerResults::HLT') # Data and MC                                    
process.btagana.genParticles          = cms.InputTag(genParticles)
process.btagana.candidates            = cms.InputTag(pfCandidates)
process.btagana.storePatMuons         = False




process.btagana.hitAssociator.associateStrip = cms.bool(False)
process.btagana.hitAssociator.associatePixel = cms.bool(False)

#process.selectedPatJets.src=cms.InputTag("akPu4PFpatJetsWithBtagging")

process.btagana.tracksColl=cms.InputTag("hiGeneralTracks")
process.btagana.trackProducer = cms.untracked.InputTag("hiGeneralTracks")

bTagInfos = [
    'akPu4PFImpactParameterTagInfos'
    ,'akPu4PFSecondaryVertexTagInfos'
    ,'akPu4PFSecondaryVertexNegativeTagInfos']

process.akPu4PFpatJetsWithBtagging.tagInfoSources=cms.VInputTag(bTagInfos)

#process.btagana.fillsvTagInfo = cms.string('')
process.btagana.ipTagInfos = cms.string('akPu4PFImpactParameter')
process.btagana.ipTagInfosCTag = cms.string('akPu4PFImpactParameter')
process.btagana.softPFElectronTagInfos = cms.string('')
process.btagana.softPFElectronTagInfosCTag = cms.string('')
process.btagana.softPFMuonTagInfos = cms.string('')
process.btagana.softPFMuonTagInfosCTag = cms.string('')
process.btagana.svNegTagInfos = cms.string('akPu4PFSecondaryVertexNegative')
process.btagana.svNegTagInfosCTag = cms.string('akPu4PFSecondaryVertexNegative')
process.btagana.svTagInfos = cms.string('akPu4PFSecondaryVertex')
process.btagana.svTagInfosCTag = cms.string('akPu4PFSecondaryVertex')

process.btagana.svComputer=cms.string('combinedSecondaryVertexComputer')



process.btagana.combinedSVBJetTags = cms.string('akPu4PFCombinedSecondaryVertexBJetTags')
process.btagana.combinedSVNegBJetTags = cms.string('akPu4PFNegativeCombinedSecondaryVertexBJetTags')
process.btagana.combinedSVPosBJetTags = cms.string('akPu4PFPositiveCombinedSecondaryVertexBJetTags')
process.btagana.jetBPBJetTags = cms.string('akPu4PFJetBProbabilityBJetTags')
process.btagana.jetPBJetTags = cms.string('akPu4PFJetProbabilityBJetTags')
process.btagana.simpleSVHighEffBJetTags = cms.string('akPu4PFSimpleSecondaryVertexHighEffBJetTags')
process.btagana.simpleSVHighPurBJetTags = cms.string('akPu4PFSimpleSecondaryVertexHighPurBJetTags')
process.btagana.simpleSVNegHighEffBJetTags = cms.string('akPu4PFNegativeSimpleSecondaryVertexHighEffBJetTags')
process.btagana.simpleSVNegHighPurBJetTags = cms.string('akPu4PFNegativeSimpleSecondaryVertexHighPurBJetTags')
process.btagana.trackCHEBJetTags = cms.string('akPu4PFTrackCountingHighEffBJetTags')
process.btagana.trackCHPBJetTags = cms.string('akPu4PFTrackCountingHighPurBJetTags')

process.btagana.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")


process.load('RecoBTag.CSVscikit.csvscikitTagJetTags_cfi')
process.load('RecoBTag.CSVscikit.csvscikitTaggerProducer_cfi')

process.akPu4PFCSVscikitJetTags = process.pfCSVscikitJetTags.clone()

process.akPu4PFpatJetsWithBtagging.discriminatorSources.append(cms.InputTag("akPu4PFCSVscikitJetTags"))
process.akPu4PFCSVscikitJetTags.tagInfos=cms.VInputTag(cms.InputTag("akPu4PFImpactParameterTagInfos"), cms.InputTag("akPu4PFSecondaryVertexTagInfos"))
process.CSVscikitTags.weightFile=cms.FileInPath('RecoBTag/PerformanceMeasurements/data/TMVA_weights2.xml')

process.jetSequences.replace(process.akPu4PFCombinedSecondaryVertexV2BJetTags,cms.Sequence(process.akPu4PFCombinedSecondaryVertexV2BJetTags+process.akPu4PFCSVscikitJetTags))


process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')
process.offlinePrimaryVertices.TrackLabel = cms.InputTag('hiGeneralTracks')


process.ana_step = cms.Path(
# Temporary disactivation - until we have DIGI & RECO in CMSSW_7_5_7_patch4
# process.mixAnalyzer *
                            process.runAnalyzer *
#                            process.hltanalysis *
                            process.centralityBin *
#                            process.hiEvtAnalyzer*
#                            process.HiGenParticleAna*
                            process.akHiGenJets +
                            process.hiSignalGenFilters + 
                            process.jetSequences 
                         # +  process.selectedPatJets*
#                           + process.offlinePrimaryVertices * 
                            + process.prunedGenParticlesBoost*process.btagana 


#                            process.hiFJRhoAnalyzer 
#                            process.ggHiNtuplizer +
#                            process.ggHiNtuplizerGED +
#                            process.pfcandAnalyzer +
#                            process.pfcandAnalyzerCS +
#                            process.HiForest +
#                            process.trackSequencesPbPb #+
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


# Customization
##########################################UE##########################################
from CondCore.DBCommon.CondDBSetup_cfi import *
process.uetable = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
          cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("UETableCompatibilityFormat_PF_HYDJET_5020GeV_754_38T_v02_mc"),
                   label = cms.untracked.string("UETable_PF")
          ),
          cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("UETableCompatibilityFormat_Calo_HYDJET_5020GeV_754_38T_v02_mc"),
                   label = cms.untracked.string("UETable_Calo")
          )
      ),
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
)
process.es_prefer_uetable = cms.ESPrefer('PoolDBESSource','uetable')
##########################################UE##########################################


