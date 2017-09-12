

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDrop4PFz01bm1match = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDrop4PFz01bm1Jets"),
    matched = cms.InputTag("ak4GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akSoftDrop4PFz01bm1matchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDrop4GenJets"),
    matched = cms.InputTag("ak4GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akSoftDrop4PFz01bm1parton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDrop4PFz01bm1Jets")
                                                        )

akSoftDrop4PFz01bm1corr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDrop4PFz01bm1Jets"),
    payload = "AK4PF_offline"
    )

akSoftDrop4PFz01bm1JetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akSoftDrop4CaloJets'))

#akSoftDrop4PFz01bm1clean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4GenJets'))

akSoftDrop4PFz01bm1bTagger = bTaggers("akSoftDrop4PFz01bm1",0.4,1,1)

#create objects locally since they dont load properly otherwise
#akSoftDrop4PFz01bm1match = akSoftDrop4PFz01bm1bTagger.match
akSoftDrop4PFz01bm1parton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDrop4PFz01bm1Jets"), matched = cms.InputTag("genParticles"))
akSoftDrop4PFz01bm1PatJetFlavourAssociationLegacy = akSoftDrop4PFz01bm1bTagger.PatJetFlavourAssociationLegacy
akSoftDrop4PFz01bm1PatJetPartons = akSoftDrop4PFz01bm1bTagger.PatJetPartons
akSoftDrop4PFz01bm1JetTracksAssociatorAtVertex = akSoftDrop4PFz01bm1bTagger.JetTracksAssociatorAtVertex
akSoftDrop4PFz01bm1JetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akSoftDrop4PFz01bm1SimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz01bm1bTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz01bm1SimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz01bm1bTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz01bm1CombinedSecondaryVertexBJetTags = akSoftDrop4PFz01bm1bTagger.CombinedSecondaryVertexBJetTags
akSoftDrop4PFz01bm1CombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01bm1bTagger.CombinedSecondaryVertexV2BJetTags
akSoftDrop4PFz01bm1JetBProbabilityBJetTags = akSoftDrop4PFz01bm1bTagger.JetBProbabilityBJetTags
akSoftDrop4PFz01bm1SoftPFMuonByPtBJetTags = akSoftDrop4PFz01bm1bTagger.SoftPFMuonByPtBJetTags
akSoftDrop4PFz01bm1SoftPFMuonByIP3dBJetTags = akSoftDrop4PFz01bm1bTagger.SoftPFMuonByIP3dBJetTags
akSoftDrop4PFz01bm1TrackCountingHighEffBJetTags = akSoftDrop4PFz01bm1bTagger.TrackCountingHighEffBJetTags
akSoftDrop4PFz01bm1TrackCountingHighPurBJetTags = akSoftDrop4PFz01bm1bTagger.TrackCountingHighPurBJetTags
akSoftDrop4PFz01bm1PatJetPartonAssociationLegacy = akSoftDrop4PFz01bm1bTagger.PatJetPartonAssociationLegacy

akSoftDrop4PFz01bm1ImpactParameterTagInfos = akSoftDrop4PFz01bm1bTagger.ImpactParameterTagInfos
akSoftDrop4PFz01bm1ImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDrop4PFz01bm1JetProbabilityBJetTags = akSoftDrop4PFz01bm1bTagger.JetProbabilityBJetTags

akSoftDrop4PFz01bm1SecondaryVertexTagInfos = akSoftDrop4PFz01bm1bTagger.SecondaryVertexTagInfos
akSoftDrop4PFz01bm1SimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz01bm1bTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz01bm1SimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz01bm1bTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz01bm1CombinedSecondaryVertexBJetTags = akSoftDrop4PFz01bm1bTagger.CombinedSecondaryVertexBJetTags
akSoftDrop4PFz01bm1CombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01bm1bTagger.CombinedSecondaryVertexV2BJetTags

akSoftDrop4PFz01bm1SecondaryVertexNegativeTagInfos = akSoftDrop4PFz01bm1bTagger.SecondaryVertexNegativeTagInfos
akSoftDrop4PFz01bm1NegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz01bm1bTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz01bm1NegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz01bm1bTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz01bm1NegativeCombinedSecondaryVertexBJetTags = akSoftDrop4PFz01bm1bTagger.NegativeCombinedSecondaryVertexBJetTags
akSoftDrop4PFz01bm1PositiveCombinedSecondaryVertexBJetTags = akSoftDrop4PFz01bm1bTagger.PositiveCombinedSecondaryVertexBJetTags
akSoftDrop4PFz01bm1NegativeCombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01bm1bTagger.NegativeCombinedSecondaryVertexV2BJetTags
akSoftDrop4PFz01bm1PositiveCombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01bm1bTagger.PositiveCombinedSecondaryVertexV2BJetTags

akSoftDrop4PFz01bm1SoftPFMuonsTagInfos = akSoftDrop4PFz01bm1bTagger.SoftPFMuonsTagInfos
akSoftDrop4PFz01bm1SoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDrop4PFz01bm1SoftPFMuonBJetTags = akSoftDrop4PFz01bm1bTagger.SoftPFMuonBJetTags
akSoftDrop4PFz01bm1SoftPFMuonByIP3dBJetTags = akSoftDrop4PFz01bm1bTagger.SoftPFMuonByIP3dBJetTags
akSoftDrop4PFz01bm1SoftPFMuonByPtBJetTags = akSoftDrop4PFz01bm1bTagger.SoftPFMuonByPtBJetTags
akSoftDrop4PFz01bm1NegativeSoftPFMuonByPtBJetTags = akSoftDrop4PFz01bm1bTagger.NegativeSoftPFMuonByPtBJetTags
akSoftDrop4PFz01bm1PositiveSoftPFMuonByPtBJetTags = akSoftDrop4PFz01bm1bTagger.PositiveSoftPFMuonByPtBJetTags
akSoftDrop4PFz01bm1PatJetFlavourIdLegacy = cms.Sequence(akSoftDrop4PFz01bm1PatJetPartonAssociationLegacy*akSoftDrop4PFz01bm1PatJetFlavourAssociationLegacy)
#Not working with our PU sub
akSoftDrop4PFz01bm1PatJetFlavourAssociation = akSoftDrop4PFz01bm1bTagger.PatJetFlavourAssociation
akSoftDrop4PFz01bm1PatJetFlavourId = cms.Sequence(akSoftDrop4PFz01bm1PatJetPartons*akSoftDrop4PFz01bm1PatJetFlavourAssociation)

#adding the subjet taggers
akSoftDrop4PFz01bm1SubjetImpactParameterTagInfos = akSoftDrop4PFz01bm1bTagger.SubjetImpactParameterTagInfos
akSoftDrop4PFz01bm1SubjetSecondaryVertexTagInfos = akSoftDrop4PFz01bm1bTagger.SubjetSecondaryVertexTagInfos
akSoftDrop4PFz01bm1SubjetJetTracksAssociatorAtVertex = akSoftDrop4PFz01bm1bTagger.SubjetJetTracksAssociatorAtVertex
akSoftDrop4PFz01bm1CombinedSubjetSecondaryVertexBJetTags = akSoftDrop4PFz01bm1bTagger.CombinedSubjetSecondaryVertexBJetTags
akSoftDrop4PFz01bm1CombinedSubjetSecondaryVertexV2BJetTags = akSoftDrop4PFz01bm1bTagger.CombinedSubjetSecondaryVertexV2BJetTags

akSoftDrop4PFz01bm1PatJetFlavourId = cms.Sequence(akSoftDrop4PFz01bm1PatJetPartons*akSoftDrop4PFz01bm1PatJetFlavourAssociation)

akSoftDrop4PFz01bm1JetBtaggingIP       = cms.Sequence(akSoftDrop4PFz01bm1ImpactParameterTagInfos *
            (akSoftDrop4PFz01bm1TrackCountingHighEffBJetTags +
             akSoftDrop4PFz01bm1TrackCountingHighPurBJetTags +
             akSoftDrop4PFz01bm1JetProbabilityBJetTags +
             akSoftDrop4PFz01bm1JetBProbabilityBJetTags 
            )
            )

akSoftDrop4PFz01bm1JetBtaggingSV = cms.Sequence(akSoftDrop4PFz01bm1ImpactParameterTagInfos
            *
            akSoftDrop4PFz01bm1SecondaryVertexTagInfos
            * (akSoftDrop4PFz01bm1SimpleSecondaryVertexHighEffBJetTags+
                akSoftDrop4PFz01bm1SimpleSecondaryVertexHighPurBJetTags+
                akSoftDrop4PFz01bm1CombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz01bm1CombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDrop4PFz01bm1JetBtaggingNegSV = cms.Sequence(akSoftDrop4PFz01bm1ImpactParameterTagInfos
            *
            akSoftDrop4PFz01bm1SecondaryVertexNegativeTagInfos
            * (akSoftDrop4PFz01bm1NegativeSimpleSecondaryVertexHighEffBJetTags+
                akSoftDrop4PFz01bm1NegativeSimpleSecondaryVertexHighPurBJetTags+
                akSoftDrop4PFz01bm1NegativeCombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz01bm1PositiveCombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz01bm1NegativeCombinedSecondaryVertexV2BJetTags+
                akSoftDrop4PFz01bm1PositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDrop4PFz01bm1JetBtaggingMu = cms.Sequence(akSoftDrop4PFz01bm1SoftPFMuonsTagInfos * (akSoftDrop4PFz01bm1SoftPFMuonBJetTags
                +
                akSoftDrop4PFz01bm1SoftPFMuonByIP3dBJetTags
                +
                akSoftDrop4PFz01bm1SoftPFMuonByPtBJetTags
                +
                akSoftDrop4PFz01bm1NegativeSoftPFMuonByPtBJetTags
                +
                akSoftDrop4PFz01bm1PositiveSoftPFMuonByPtBJetTags
              )
            )

akSoftDrop4PFz01bm1JetBtagging = cms.Sequence(akSoftDrop4PFz01bm1JetBtaggingIP
            *akSoftDrop4PFz01bm1JetBtaggingSV
            *akSoftDrop4PFz01bm1JetBtaggingNegSV
#            *akSoftDrop4PFz01bm1JetBtaggingMu
            )

akSoftDrop4PFz01bm1patJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akSoftDrop4PFz01bm1Jets"),
        genJetMatch          = cms.InputTag("akSoftDrop4PFz01bm1match"),
        genPartonMatch       = cms.InputTag("akSoftDrop4PFz01bm1parton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akSoftDrop4PFz01bm1corr")),
        #JetPartonMapSource   = cms.InputTag("akSoftDrop4PFz01bm1PatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akSoftDrop4PFz01bm1PatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akSoftDrop4PFz01bm1PatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akSoftDrop4PFz01bm1JetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akSoftDrop4PFz01bm1SimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1SimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1CombinedSecondaryVertexBJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1CombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1JetBProbabilityBJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1JetProbabilityBJetTags"),
            #cms.InputTag("akSoftDrop4PFz01bm1SoftPFMuonByPtBJetTags"),
            #cms.InputTag("akSoftDrop4PFz01bm1SoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1TrackCountingHighEffBJetTags"),
            cms.InputTag("akSoftDrop4PFz01bm1TrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akSoftDrop4PFz01bm1JetID"),
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

akSoftDrop4PFz01bm1Njettiness = Njettiness.clone(
		    src = cms.InputTag("akSoftDrop4PFz01bm1Jets"),
           	    R0  = cms.double( 0.4)
)
akSoftDrop4PFz01bm1patJetsWithBtagging.userData.userFloats.src += ['akSoftDrop4PFz01bm1Njettiness:tau1','akSoftDrop4PFz01bm1Njettiness:tau2','akSoftDrop4PFz01bm1Njettiness:tau3']

akSoftDrop4PFz01bm1JetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akSoftDrop4PFz01bm1patJetsWithBtagging"),
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
                                                             bTagJetName = cms.untracked.string("akSoftDrop4PFz01bm1"),
                                                             jetName = cms.untracked.string("akSoftDrop4PFz01bm1"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(False),
							     doSubJets = cms.untracked.bool(True),
                                                             doGenSubJets = cms.untracked.bool(True),     
                                                             subjetGenTag = cms.untracked.InputTag("akSoftDrop4GenJets"),
                                                             doGenTaus = True
                                                            )

akSoftDrop4PFz01bm1JetSequence_mc = cms.Sequence(
                                                  #akSoftDrop4PFz01bm1clean
                                                  #*
                                                  akSoftDrop4PFz01bm1match
                                                  #*
                                                  #akSoftDrop4PFz01bm1matchGroomed
                                                  *
                                                  akSoftDrop4PFz01bm1parton
                                                  *
                                                  akSoftDrop4PFz01bm1corr
                                                  *
                                                  #akSoftDrop4PFz01bm1JetID
                                                  #*
                                                  #akSoftDrop4PFz01bm1PatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akSoftDrop4PFz01bm1PatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akSoftDrop4PFz01bm1JetTracksAssociatorAtVertex
                                                  *
                                                  akSoftDrop4PFz01bm1JetBtagging
                                                  *
                                                  akSoftDrop4PFz01bm1Njettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akSoftDrop4PFz01bm1patJetsWithBtagging
                                                  *
                                                  akSoftDrop4PFz01bm1JetAnalyzer
                                                  )

akSoftDrop4PFz01bm1JetSequence_data = cms.Sequence(akSoftDrop4PFz01bm1corr
                                                    *
                                                    #akSoftDrop4PFz01bm1JetID
                                                    #*
                                                    akSoftDrop4PFz01bm1JetTracksAssociatorAtVertex
                                                    *
                                                    akSoftDrop4PFz01bm1JetBtagging
                                                    *
                                                    akSoftDrop4PFz01bm1Njettiness 
                                                    *
                                                    akSoftDrop4PFz01bm1patJetsWithBtagging
                                                    *
                                                    akSoftDrop4PFz01bm1JetAnalyzer
                                                    )

akSoftDrop4PFz01bm1JetSequence_jec = cms.Sequence(akSoftDrop4PFz01bm1JetSequence_mc)
akSoftDrop4PFz01bm1JetSequence_mb = cms.Sequence(akSoftDrop4PFz01bm1JetSequence_mc)

akSoftDrop4PFz01bm1JetSequence = cms.Sequence(akSoftDrop4PFz01bm1JetSequence_mc)
