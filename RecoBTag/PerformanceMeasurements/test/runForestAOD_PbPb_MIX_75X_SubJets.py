### HiForest Configuration
# Collisions: PbPb
# Type: MonteCarlo
# Input: AOD

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')

options.setDefault('maxEvents', 10)
options.setDefault('outputFile', 'JetTree_qcd_CsJets')
#options.setDefault('inputFiles',['file:qcd_99.root'])
#options.setDefault('inputFiles',['/store/himc/HINPbPbWinter16DR/Pythia6_bJetFCR30_pp502_Hydjet_Cymbal_MB/AODSIM/75X_mcRun2_HeavyIon_v14-v1/120000/0421C390-6BF6-E611-8856-008CFA0021D4.root'])
options.setDefault('inputFiles',['/store/himc/HINPbPbWinter16DR/Pythia6_Dijet120_pp502_Hydjet_MB/AODSIM/75X_mcRun2_HeavyIon_v13-v1/00000/04A04E01-A80D-E611-835A-02163E012AD1.root'])

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
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
      tag = cms.string("CentralityTable_HFtowers200_HydjetCymbal5Ev8_v758x03_mc"),
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
      label = cms.untracked.string("HFtowersCymbal5Ev8")
   ),
])


from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_PbPb5020
process = overrideJEC_PbPb5020(process)

#process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
# process.centralityBin.Centrality = cms.InputTag("hiCentrality")
# process.centralityBin.centralityVariable = cms.string("HFtowers")
#process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("Cymbal5Ev8")

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
#process.load('HeavyIonsAnalysis.JetAnalysis.FullJetSequence_nominalPbPb')
#jet sequence as in Kurts' ex
process.load('HeavyIonsAnalysis.JetAnalysis.FullJetSequence_puLimitedPbPb')
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

process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
process.hiEvtAnalyzer.doCentrality = cms.bool(True)
process.hiEvtAnalyzer.CentralitySrc = cms.InputTag("hiCentrality")
process.hiEvtAnalyzer.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
process.hiEvtAnalyzer.doMC = cms.bool(True) #general MC info
process.hiEvtAnalyzer.doHiMC = cms.bool(True) #HI specific MC info
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
#process.btagana.Jets                  = cms.InputTag('akCs4PFpatJetsWithBtagging')
process.btagana.Jets                  = cms.InputTag('akCsSoftDrop4PFpatJetsWithBtagging')
process.btagana.runSubJets            = cms.bool(True)
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
    'akCs4PFImpactParameterTagInfos'
    ,'akCs4PFSecondaryVertexTagInfos'
    ,'akCs4PFSecondaryVertexNegativeTagInfos']

process.akCs4PFpatJetsWithBtagging.tagInfoSources=cms.VInputTag(bTagInfos)

#process.btagana.fillsvTagInfo = cms.string('')

process.btagana.ipTagInfos = cms.string('akCs4PFImpactParameter')
process.btagana.ipTagInfosCTag = cms.string('akCs4PFImpactParameter')
process.btagana.softPFElectronTagInfos = cms.string('')
process.btagana.softPFElectronTagInfosCTag = cms.string('')
process.btagana.softPFMuonTagInfos = cms.string('')
process.btagana.softPFMuonTagInfosCTag = cms.string('')
process.btagana.svNegTagInfos = cms.string('akCs4PFSecondaryVertexNegative')
process.btagana.svNegTagInfosCTag = cms.string('akCs4PFSecondaryVertexNegative')
process.btagana.svTagInfos = cms.string('akCs4PFSecondaryVertex')
process.btagana.svTagInfosCTag = cms.string('akCs4PFSecondaryVertex')

process.btagana.svComputer=cms.string('combinedSecondaryVertexComputer')

process.btagana.combinedSVBJetTags = cms.string('akCs4PFCombinedSecondaryVertexBJetTags')
process.btagana.combinedSVNegBJetTags = cms.string('akCs4PFNegativeCombinedSecondaryVertexBJetTags')
process.btagana.combinedSVPosBJetTags = cms.string('akCs4PFPositiveCombinedSecondaryVertexBJetTags')
process.btagana.jetBPBJetTags = cms.string('akCs4PFJetBProbabilityBJetTags')
process.btagana.jetPBJetTags = cms.string('akCs4PFJetProbabilityBJetTags')
process.btagana.simpleSVHighEffBJetTags = cms.string('akCs4PFSimpleSecondaryVertexHighEffBJetTags')
process.btagana.simpleSVHighPurBJetTags = cms.string('akCs4PFSimpleSecondaryVertexHighPurBJetTags')
process.btagana.simpleSVNegHighEffBJetTags = cms.string('akCs4PFNegativeSimpleSecondaryVertexHighEffBJetTags')
process.btagana.simpleSVNegHighPurBJetTags = cms.string('akCs4PFNegativeSimpleSecondaryVertexHighPurBJetTags')
process.btagana.trackCHEBJetTags = cms.string('akCs4PFTrackCountingHighEffBJetTags')
process.btagana.trackCHPBJetTags = cms.string('akCs4PFTrackCountingHighPurBJetTags')

process.btagana.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")


process.load('RecoBTag.CSVscikit.csvscikitTagJetTags_cfi')
process.load('RecoBTag.CSVscikit.csvscikitTaggerProducer_cfi')

process.akCs4PFCSVscikitJetTags = process.pfCSVscikitJetTags.clone()

process.akCs4PFpatJetsWithBtagging.discriminatorSources.append(cms.InputTag("akCs4PFCSVscikitJetTags"))
process.akCs4PFCSVscikitJetTags.tagInfos=cms.VInputTag(cms.InputTag("akCs4PFImpactParameterTagInfos"), cms.InputTag("akCs4PFSecondaryVertexTagInfos"))
process.CSVscikitTags.weightFile=cms.FileInPath('RecoBTag/PerformanceMeasurements/data/TMVA_weights2.xml')

process.jetSequences.replace(process.akCs4PFCombinedSecondaryVertexV2BJetTags,cms.Sequence(process.akCs4PFCombinedSecondaryVertexV2BJetTags+process.akCs4PFCSVscikitJetTags))


process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')
process.offlinePrimaryVertices.TrackLabel = cms.InputTag('hiGeneralTracks')


process.ana_step = cms.Path(
# Temporary disactivation - until we have DIGI & RECO in CMSSW_7_5_7_patch4
# process.mixAnalyzer *
                            process.runAnalyzer *
#                            process.hltanalysis *
                            process.centralityBin *
                            process.hiEvtAnalyzer*
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


#Jets customization
##########################################UE##########################################
process.ak4PFJets.jetPtMin = cms.double(0.0)
process.akSoftDrop4PFJets.src = cms.InputTag("particleFlowTmp")
process.akCsSoftDrop4PFJets.jetPtMin = cms.double(0.0)

process.akCsSoftDrop4PFPatJetFlavourAssociation.jets="akCs4PFJets"
process.akCsSoftDrop4PFPatJetFlavourAssociation.unsubtractedJets = cms.InputTag("ak4PFJets")
process.akCsSoftDrop4PFPatJetFlavourAssociation.redoUESubtraction = cms.untracked.bool(True)
process.akCsSoftDrop4PFPatJetFlavourAssociation.groomedJets=cms.InputTag("akCsSoftDrop4PFJets")
process.akCsSoftDrop4PFPatJetFlavourAssociation.subjets= cms.InputTag('akCsSoftDrop4PFJets','SubJets')
process.akCsSoftDrop4PFPatJetFlavourAssociation.etaMap = cms.InputTag("hiFJRhoProducer","mapEtaEdges")
process.akCsSoftDrop4PFPatJetFlavourAssociation.rho = cms.InputTag("hiFJGridEmptyAreaCalculator","mapToRhoCorr1Bin")
process.akCsSoftDrop4PFPatJetFlavourAssociation.rhom = cms.InputTag("hiFJGridEmptyAreaCalculator","mapToRhoMCorr1Bin")

process.akCsSoftDrop4PFpatJetsWithBtagging.getJetMCFlavour = cms.bool(False)
process.akCsSoftDrop4PFJetAnalyzer.doExtendedFlavorTagging = cms.untracked.bool(True)
process.akCsSoftDrop4PFJetAnalyzer.jetFlavourInfos    = cms.InputTag("akCsSoftDrop4PFPatJetFlavourAssociation")
process.akCsSoftDrop4PFJetAnalyzer.subjetFlavourInfos = cms.InputTag("akCsSoftDrop4PFPatJetFlavourAssociation","SubJets")
process.akCsSoftDrop4PFJetAnalyzer.groomedJets        = cms.InputTag("akCsSoftDrop4PFJets")
process.akCsSoftDrop4PFJetAnalyzer.isPythia6 = cms.untracked.bool(True)

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


process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.extSVDeltaRToJet = cms.double(1.0) #make this big to make sure ghost vertexing works                                                                    
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.useSVClustering = cms.bool(True)
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.fatJets = cms.InputTag("ak4PFJets")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.groomedFatJets = cms.InputTag("akCsSoftDrop4PFJets")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.jetAlgorithm = cms.string("AntiKt")
process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.rParam = cms.double(0.4)

process.akCsSoftDrop4PFJetAnalyzer.trackSelection = process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.trackSelection
process.akCsSoftDrop4PFJetAnalyzer.trackPairV0Filter = process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.vertexCuts.v0Filter

process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = cms.double(0.2)
process.akCsSoftDrop4PFCombinedSubjetSecondaryVertexBJetTags = process.akCsSoftDrop4PFCombinedSecondaryVertexBJetTags.clone(
        tagInfos = cms.VInputTag(cms.InputTag("akCsSoftDrop4PFSubjetImpactParameterTagInfos"),
                cms.InputTag("akCsSoftDrop4PFSubjetSecondaryVertexTagInfos"))
)
process.akCsSoftDrop4PFJetBtaggingSV *= process.akCsSoftDrop4PFSubjetJetTracksAssociatorAtVertex+process.akCsSoftDrop4PFSubjetImpactParameterTagInfos+process.akCsSoftDrop4PFSubjetJetProbabilityBJetTags+process.inclusiveVertexing+process.akCsSoftDrop4PFSubjetSecondaryVertexTagInfos+process.akCsSoftDrop4PFCombinedSubjetSecondaryVertexBJetTags



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


