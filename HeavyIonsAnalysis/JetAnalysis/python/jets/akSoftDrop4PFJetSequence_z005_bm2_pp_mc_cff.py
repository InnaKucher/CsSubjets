

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDrop4PFz005bm2match = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDrop4PFz005bm2Jets"),
    matched = cms.InputTag("ak4GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akSoftDrop4PFz005bm2matchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDrop4GenJets"),
    matched = cms.InputTag("ak4GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akSoftDrop4PFz005bm2parton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDrop4PFz005bm2Jets")
                                                        )

akSoftDrop4PFz005bm2corr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDrop4PFz005bm2Jets"),
    payload = "AK4PF_offline"
    )

akSoftDrop4PFz005bm2JetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akSoftDrop4CaloJets'))

#akSoftDrop4PFz005bm2clean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4GenJets'))

akSoftDrop4PFz005bm2bTagger = bTaggers("akSoftDrop4PFz005bm2",0.4,1,1)

#create objects locally since they dont load properly otherwise
#akSoftDrop4PFz005bm2match = akSoftDrop4PFz005bm2bTagger.match
akSoftDrop4PFz005bm2parton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDrop4PFz005bm2Jets"), matched = cms.InputTag("genParticles"))
akSoftDrop4PFz005bm2PatJetFlavourAssociationLegacy = akSoftDrop4PFz005bm2bTagger.PatJetFlavourAssociationLegacy
akSoftDrop4PFz005bm2PatJetPartons = akSoftDrop4PFz005bm2bTagger.PatJetPartons
akSoftDrop4PFz005bm2JetTracksAssociatorAtVertex = akSoftDrop4PFz005bm2bTagger.JetTracksAssociatorAtVertex
akSoftDrop4PFz005bm2JetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akSoftDrop4PFz005bm2SimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz005bm2bTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz005bm2SimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz005bm2bTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz005bm2CombinedSecondaryVertexBJetTags = akSoftDrop4PFz005bm2bTagger.CombinedSecondaryVertexBJetTags
akSoftDrop4PFz005bm2CombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz005bm2bTagger.CombinedSecondaryVertexV2BJetTags
akSoftDrop4PFz005bm2JetBProbabilityBJetTags = akSoftDrop4PFz005bm2bTagger.JetBProbabilityBJetTags
akSoftDrop4PFz005bm2SoftPFMuonByPtBJetTags = akSoftDrop4PFz005bm2bTagger.SoftPFMuonByPtBJetTags
akSoftDrop4PFz005bm2SoftPFMuonByIP3dBJetTags = akSoftDrop4PFz005bm2bTagger.SoftPFMuonByIP3dBJetTags
akSoftDrop4PFz005bm2TrackCountingHighEffBJetTags = akSoftDrop4PFz005bm2bTagger.TrackCountingHighEffBJetTags
akSoftDrop4PFz005bm2TrackCountingHighPurBJetTags = akSoftDrop4PFz005bm2bTagger.TrackCountingHighPurBJetTags
akSoftDrop4PFz005bm2PatJetPartonAssociationLegacy = akSoftDrop4PFz005bm2bTagger.PatJetPartonAssociationLegacy

akSoftDrop4PFz005bm2ImpactParameterTagInfos = akSoftDrop4PFz005bm2bTagger.ImpactParameterTagInfos
akSoftDrop4PFz005bm2ImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDrop4PFz005bm2JetProbabilityBJetTags = akSoftDrop4PFz005bm2bTagger.JetProbabilityBJetTags

akSoftDrop4PFz005bm2SecondaryVertexTagInfos = akSoftDrop4PFz005bm2bTagger.SecondaryVertexTagInfos
akSoftDrop4PFz005bm2SimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz005bm2bTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz005bm2SimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz005bm2bTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz005bm2CombinedSecondaryVertexBJetTags = akSoftDrop4PFz005bm2bTagger.CombinedSecondaryVertexBJetTags
akSoftDrop4PFz005bm2CombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz005bm2bTagger.CombinedSecondaryVertexV2BJetTags

akSoftDrop4PFz005bm2SecondaryVertexNegativeTagInfos = akSoftDrop4PFz005bm2bTagger.SecondaryVertexNegativeTagInfos
akSoftDrop4PFz005bm2NegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz005bm2bTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz005bm2NegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz005bm2bTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz005bm2NegativeCombinedSecondaryVertexBJetTags = akSoftDrop4PFz005bm2bTagger.NegativeCombinedSecondaryVertexBJetTags
akSoftDrop4PFz005bm2PositiveCombinedSecondaryVertexBJetTags = akSoftDrop4PFz005bm2bTagger.PositiveCombinedSecondaryVertexBJetTags
akSoftDrop4PFz005bm2NegativeCombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz005bm2bTagger.NegativeCombinedSecondaryVertexV2BJetTags
akSoftDrop4PFz005bm2PositiveCombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz005bm2bTagger.PositiveCombinedSecondaryVertexV2BJetTags

akSoftDrop4PFz005bm2SoftPFMuonsTagInfos = akSoftDrop4PFz005bm2bTagger.SoftPFMuonsTagInfos
akSoftDrop4PFz005bm2SoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDrop4PFz005bm2SoftPFMuonBJetTags = akSoftDrop4PFz005bm2bTagger.SoftPFMuonBJetTags
akSoftDrop4PFz005bm2SoftPFMuonByIP3dBJetTags = akSoftDrop4PFz005bm2bTagger.SoftPFMuonByIP3dBJetTags
akSoftDrop4PFz005bm2SoftPFMuonByPtBJetTags = akSoftDrop4PFz005bm2bTagger.SoftPFMuonByPtBJetTags
akSoftDrop4PFz005bm2NegativeSoftPFMuonByPtBJetTags = akSoftDrop4PFz005bm2bTagger.NegativeSoftPFMuonByPtBJetTags
akSoftDrop4PFz005bm2PositiveSoftPFMuonByPtBJetTags = akSoftDrop4PFz005bm2bTagger.PositiveSoftPFMuonByPtBJetTags
akSoftDrop4PFz005bm2PatJetFlavourIdLegacy = cms.Sequence(akSoftDrop4PFz005bm2PatJetPartonAssociationLegacy*akSoftDrop4PFz005bm2PatJetFlavourAssociationLegacy)
#Not working with our PU sub
akSoftDrop4PFz005bm2PatJetFlavourAssociation = akSoftDrop4PFz005bm2bTagger.PatJetFlavourAssociation
akSoftDrop4PFz005bm2PatJetFlavourId = cms.Sequence(akSoftDrop4PFz005bm2PatJetPartons*akSoftDrop4PFz005bm2PatJetFlavourAssociation)

#adding the subjet taggers
akSoftDrop4PFz005bm2SubjetImpactParameterTagInfos = akSoftDrop4PFz005bm2bTagger.SubjetImpactParameterTagInfos
akSoftDrop4PFz005bm2SubjetSecondaryVertexTagInfos = akSoftDrop4PFz005bm2bTagger.SubjetSecondaryVertexTagInfos
akSoftDrop4PFz005bm2SubjetJetTracksAssociatorAtVertex = akSoftDrop4PFz005bm2bTagger.SubjetJetTracksAssociatorAtVertex
akSoftDrop4PFz005bm2CombinedSubjetSecondaryVertexBJetTags = akSoftDrop4PFz005bm2bTagger.CombinedSubjetSecondaryVertexBJetTags
akSoftDrop4PFz005bm2CombinedSubjetSecondaryVertexV2BJetTags = akSoftDrop4PFz005bm2bTagger.CombinedSubjetSecondaryVertexV2BJetTags

akSoftDrop4PFz005bm2JetBtaggingIP       = cms.Sequence(akSoftDrop4PFz005bm2ImpactParameterTagInfos *
            (akSoftDrop4PFz005bm2TrackCountingHighEffBJetTags +
             akSoftDrop4PFz005bm2TrackCountingHighPurBJetTags +
             akSoftDrop4PFz005bm2JetProbabilityBJetTags +
             akSoftDrop4PFz005bm2JetBProbabilityBJetTags 
            )
            )

akSoftDrop4PFz005bm2JetBtaggingSV = cms.Sequence(akSoftDrop4PFz005bm2ImpactParameterTagInfos
            *
            akSoftDrop4PFz005bm2SecondaryVertexTagInfos
            * (akSoftDrop4PFz005bm2SimpleSecondaryVertexHighEffBJetTags+
                akSoftDrop4PFz005bm2SimpleSecondaryVertexHighPurBJetTags+
                akSoftDrop4PFz005bm2CombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz005bm2CombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDrop4PFz005bm2JetBtaggingNegSV = cms.Sequence(akSoftDrop4PFz005bm2ImpactParameterTagInfos
            *
            akSoftDrop4PFz005bm2SecondaryVertexNegativeTagInfos
            * (akSoftDrop4PFz005bm2NegativeSimpleSecondaryVertexHighEffBJetTags+
                akSoftDrop4PFz005bm2NegativeSimpleSecondaryVertexHighPurBJetTags+
                akSoftDrop4PFz005bm2NegativeCombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz005bm2PositiveCombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz005bm2NegativeCombinedSecondaryVertexV2BJetTags+
                akSoftDrop4PFz005bm2PositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDrop4PFz005bm2JetBtaggingMu = cms.Sequence(akSoftDrop4PFz005bm2SoftPFMuonsTagInfos * (akSoftDrop4PFz005bm2SoftPFMuonBJetTags
                +
                akSoftDrop4PFz005bm2SoftPFMuonByIP3dBJetTags
                +
                akSoftDrop4PFz005bm2SoftPFMuonByPtBJetTags
                +
                akSoftDrop4PFz005bm2NegativeSoftPFMuonByPtBJetTags
                +
                akSoftDrop4PFz005bm2PositiveSoftPFMuonByPtBJetTags
              )
            )

akSoftDrop4PFz005bm2JetBtagging = cms.Sequence(akSoftDrop4PFz005bm2JetBtaggingIP
            *akSoftDrop4PFz005bm2JetBtaggingSV
            *akSoftDrop4PFz005bm2JetBtaggingNegSV
#            *akSoftDrop4PFz005bm2JetBtaggingMu
            )

akSoftDrop4PFz005bm2patJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akSoftDrop4PFz005bm2Jets"),
        genJetMatch          = cms.InputTag("akSoftDrop4PFz005bm2match"),
        genPartonMatch       = cms.InputTag("akSoftDrop4PFz005bm2parton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akSoftDrop4PFz005bm2corr")),
        #JetPartonMapSource   = cms.InputTag("akSoftDrop4PFz005bm2PatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akSoftDrop4PFz005bm2PatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akSoftDrop4PFz005bm2PatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akSoftDrop4PFz005bm2JetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akSoftDrop4PFz005bm2SimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2SimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2CombinedSecondaryVertexBJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2CombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2JetBProbabilityBJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2JetProbabilityBJetTags"),
            #cms.InputTag("akSoftDrop4PFz005bm2SoftPFMuonByPtBJetTags"),
            #cms.InputTag("akSoftDrop4PFz005bm2SoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2TrackCountingHighEffBJetTags"),
            cms.InputTag("akSoftDrop4PFz005bm2TrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akSoftDrop4PFz005bm2JetID"),
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = True,
        addGenPartonMatch = True,
        addGenJetMatch = True,
        embedGenJetMatch = True,
        embedGenPartonMatch = True,
        # embedCaloTowers = False,
        # embedPFCandidates = True
        )

akSoftDrop4PFz005bm2Njettiness = Njettiness.clone(
		    src = cms.InputTag("akSoftDrop4PFz005bm2Jets"),
           	    R0  = cms.double( 0.4)
)
akSoftDrop4PFz005bm2patJetsWithBtagging.userData.userFloats.src += ['akSoftDrop4PFz005bm2Njettiness:tau1','akSoftDrop4PFz005bm2Njettiness:tau2','akSoftDrop4PFz005bm2Njettiness:tau3']

akSoftDrop4PFz005bm2JetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akSoftDrop4PFz005bm2patJetsWithBtagging"),
                                                             genjetTag = 'ak4GenJets',
                                                             rParam = 0.4,
                                                             matchJets = cms.untracked.bool(False),
                                                             matchTag = 'patJetsWithBtagging',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlow'),
                                                             trackTag = cms.InputTag("generalTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
							     doSubEvent = True,
                                                             useHepMC = cms.untracked.bool(False),
							     genParticles = cms.untracked.InputTag("genParticles"),
							     eventInfoTag = cms.InputTag("generator"),
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akSoftDrop4PFz005bm2"),
                                                             jetName = cms.untracked.string("akSoftDrop4PFz005bm2"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(False),
							     doSubJets = cms.untracked.bool(True),
                                                             doGenSubJets = cms.untracked.bool(True),     
                                                             subjetGenTag = cms.untracked.InputTag("akSoftDrop4GenJets"),
                                                             doGenTaus = True
                                                            )

akSoftDrop4PFz005bm2JetSequence_mc = cms.Sequence(
                                                  #akSoftDrop4PFz005bm2clean
                                                  #*
                                                  akSoftDrop4PFz005bm2match
                                                  #*
                                                  #akSoftDrop4PFz005bm2matchGroomed
                                                  *
                                                  akSoftDrop4PFz005bm2parton
                                                  *
                                                  akSoftDrop4PFz005bm2corr
                                                  *
                                                  #akSoftDrop4PFz005bm2JetID
                                                  #*
                                                  #akSoftDrop4PFz005bm2PatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akSoftDrop4PFz005bm2PatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akSoftDrop4PFz005bm2JetTracksAssociatorAtVertex
                                                  *
                                                  akSoftDrop4PFz005bm2JetBtagging
                                                  *
                                                  akSoftDrop4PFz005bm2Njettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akSoftDrop4PFz005bm2patJetsWithBtagging
                                                  *
                                                  akSoftDrop4PFz005bm2JetAnalyzer
                                                  )

akSoftDrop4PFz005bm2JetSequence_data = cms.Sequence(akSoftDrop4PFz005bm2corr
                                                    *
                                                    #akSoftDrop4PFz005bm2JetID
                                                    #*
                                                    akSoftDrop4PFz005bm2JetTracksAssociatorAtVertex
                                                    *
                                                    akSoftDrop4PFz005bm2JetBtagging
                                                    *
                                                    akSoftDrop4PFz005bm2Njettiness 
                                                    *
                                                    akSoftDrop4PFz005bm2patJetsWithBtagging
                                                    *
                                                    akSoftDrop4PFz005bm2JetAnalyzer
                                                    )

akSoftDrop4PFz005bm2JetSequence_jec = cms.Sequence(akSoftDrop4PFz005bm2JetSequence_mc)
akSoftDrop4PFz005bm2JetSequence_mb = cms.Sequence(akSoftDrop4PFz005bm2JetSequence_mc)

akSoftDrop4PFz005bm2JetSequence = cms.Sequence(akSoftDrop4PFz005bm2JetSequence_mc)
