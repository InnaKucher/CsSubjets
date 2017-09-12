

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDrop4PFz01b1match = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDrop4PFz01b1Jets"),
    matched = cms.InputTag("ak4GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akSoftDrop4PFz01b1matchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDrop4GenJets"),
    matched = cms.InputTag("ak4GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akSoftDrop4PFz01b1parton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDrop4PFz01b1Jets")
                                                        )

akSoftDrop4PFz01b1corr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDrop4PFz01b1Jets"),
    payload = "AK4PF_offline"
    )

akSoftDrop4PFz01b1JetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akSoftDrop4CaloJets'))

#akSoftDrop4PFz01b1clean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4GenJets'))

akSoftDrop4PFz01b1bTagger = bTaggers("akSoftDrop4PFz01b1",0.4,1,1)

#create objects locally since they dont load properly otherwise
#akSoftDrop4PFz01b1match = akSoftDrop4PFz01b1bTagger.match
akSoftDrop4PFz01b1parton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDrop4PFz01b1Jets"), matched = cms.InputTag("genParticles"))
akSoftDrop4PFz01b1PatJetFlavourAssociationLegacy = akSoftDrop4PFz01b1bTagger.PatJetFlavourAssociationLegacy
akSoftDrop4PFz01b1PatJetPartons = akSoftDrop4PFz01b1bTagger.PatJetPartons
akSoftDrop4PFz01b1JetTracksAssociatorAtVertex = akSoftDrop4PFz01b1bTagger.JetTracksAssociatorAtVertex
akSoftDrop4PFz01b1JetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akSoftDrop4PFz01b1SimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz01b1bTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz01b1SimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz01b1bTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz01b1CombinedSecondaryVertexBJetTags = akSoftDrop4PFz01b1bTagger.CombinedSecondaryVertexBJetTags
akSoftDrop4PFz01b1CombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01b1bTagger.CombinedSecondaryVertexV2BJetTags
akSoftDrop4PFz01b1JetBProbabilityBJetTags = akSoftDrop4PFz01b1bTagger.JetBProbabilityBJetTags
akSoftDrop4PFz01b1SoftPFMuonByPtBJetTags = akSoftDrop4PFz01b1bTagger.SoftPFMuonByPtBJetTags
akSoftDrop4PFz01b1SoftPFMuonByIP3dBJetTags = akSoftDrop4PFz01b1bTagger.SoftPFMuonByIP3dBJetTags
akSoftDrop4PFz01b1TrackCountingHighEffBJetTags = akSoftDrop4PFz01b1bTagger.TrackCountingHighEffBJetTags
akSoftDrop4PFz01b1TrackCountingHighPurBJetTags = akSoftDrop4PFz01b1bTagger.TrackCountingHighPurBJetTags
akSoftDrop4PFz01b1PatJetPartonAssociationLegacy = akSoftDrop4PFz01b1bTagger.PatJetPartonAssociationLegacy

akSoftDrop4PFz01b1ImpactParameterTagInfos = akSoftDrop4PFz01b1bTagger.ImpactParameterTagInfos
akSoftDrop4PFz01b1ImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDrop4PFz01b1JetProbabilityBJetTags = akSoftDrop4PFz01b1bTagger.JetProbabilityBJetTags

akSoftDrop4PFz01b1SecondaryVertexTagInfos = akSoftDrop4PFz01b1bTagger.SecondaryVertexTagInfos
akSoftDrop4PFz01b1SimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz01b1bTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz01b1SimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz01b1bTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz01b1CombinedSecondaryVertexBJetTags = akSoftDrop4PFz01b1bTagger.CombinedSecondaryVertexBJetTags
akSoftDrop4PFz01b1CombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01b1bTagger.CombinedSecondaryVertexV2BJetTags

akSoftDrop4PFz01b1SecondaryVertexNegativeTagInfos = akSoftDrop4PFz01b1bTagger.SecondaryVertexNegativeTagInfos
akSoftDrop4PFz01b1NegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDrop4PFz01b1bTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akSoftDrop4PFz01b1NegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDrop4PFz01b1bTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akSoftDrop4PFz01b1NegativeCombinedSecondaryVertexBJetTags = akSoftDrop4PFz01b1bTagger.NegativeCombinedSecondaryVertexBJetTags
akSoftDrop4PFz01b1PositiveCombinedSecondaryVertexBJetTags = akSoftDrop4PFz01b1bTagger.PositiveCombinedSecondaryVertexBJetTags
akSoftDrop4PFz01b1NegativeCombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01b1bTagger.NegativeCombinedSecondaryVertexV2BJetTags
akSoftDrop4PFz01b1PositiveCombinedSecondaryVertexV2BJetTags = akSoftDrop4PFz01b1bTagger.PositiveCombinedSecondaryVertexV2BJetTags

akSoftDrop4PFz01b1SoftPFMuonsTagInfos = akSoftDrop4PFz01b1bTagger.SoftPFMuonsTagInfos
akSoftDrop4PFz01b1SoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDrop4PFz01b1SoftPFMuonBJetTags = akSoftDrop4PFz01b1bTagger.SoftPFMuonBJetTags
akSoftDrop4PFz01b1SoftPFMuonByIP3dBJetTags = akSoftDrop4PFz01b1bTagger.SoftPFMuonByIP3dBJetTags
akSoftDrop4PFz01b1SoftPFMuonByPtBJetTags = akSoftDrop4PFz01b1bTagger.SoftPFMuonByPtBJetTags
akSoftDrop4PFz01b1NegativeSoftPFMuonByPtBJetTags = akSoftDrop4PFz01b1bTagger.NegativeSoftPFMuonByPtBJetTags
akSoftDrop4PFz01b1PositiveSoftPFMuonByPtBJetTags = akSoftDrop4PFz01b1bTagger.PositiveSoftPFMuonByPtBJetTags
akSoftDrop4PFz01b1PatJetFlavourIdLegacy = cms.Sequence(akSoftDrop4PFz01b1PatJetPartonAssociationLegacy*akSoftDrop4PFz01b1PatJetFlavourAssociationLegacy)
#Not working with our PU sub
akSoftDrop4PFz01b1PatJetFlavourAssociation = akSoftDrop4PFz01b1bTagger.PatJetFlavourAssociation
akSoftDrop4PFz01b1PatJetFlavourId = cms.Sequence(akSoftDrop4PFz01b1PatJetPartons*akSoftDrop4PFz01b1PatJetFlavourAssociation)

#adding the subjet taggers
akSoftDrop4PFz01b1SubjetImpactParameterTagInfos = akSoftDrop4PFz01b1bTagger.SubjetImpactParameterTagInfos
akSoftDrop4PFz01b1SubjetSecondaryVertexTagInfos = akSoftDrop4PFz01b1bTagger.SubjetSecondaryVertexTagInfos
akSoftDrop4PFz01b1SubjetJetTracksAssociatorAtVertex = akSoftDrop4PFz01b1bTagger.SubjetJetTracksAssociatorAtVertex
akSoftDrop4PFz01b1CombinedSubjetSecondaryVertexBJetTags = akSoftDrop4PFz01b1bTagger.CombinedSubjetSecondaryVertexBJetTags
akSoftDrop4PFz01b1CombinedSubjetSecondaryVertexV2BJetTags = akSoftDrop4PFz01b1bTagger.CombinedSubjetSecondaryVertexV2BJetTags

akSoftDrop4PFz01b1JetBtaggingIP       = cms.Sequence(akSoftDrop4PFz01b1ImpactParameterTagInfos *
            (akSoftDrop4PFz01b1TrackCountingHighEffBJetTags +
             akSoftDrop4PFz01b1TrackCountingHighPurBJetTags +
             akSoftDrop4PFz01b1JetProbabilityBJetTags +
             akSoftDrop4PFz01b1JetBProbabilityBJetTags 
            )
            )

akSoftDrop4PFz01b1JetBtaggingSV = cms.Sequence(akSoftDrop4PFz01b1ImpactParameterTagInfos
            *
            akSoftDrop4PFz01b1SecondaryVertexTagInfos
            * (akSoftDrop4PFz01b1SimpleSecondaryVertexHighEffBJetTags+
                akSoftDrop4PFz01b1SimpleSecondaryVertexHighPurBJetTags+
                akSoftDrop4PFz01b1CombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz01b1CombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDrop4PFz01b1JetBtaggingNegSV = cms.Sequence(akSoftDrop4PFz01b1ImpactParameterTagInfos
            *
            akSoftDrop4PFz01b1SecondaryVertexNegativeTagInfos
            * (akSoftDrop4PFz01b1NegativeSimpleSecondaryVertexHighEffBJetTags+
                akSoftDrop4PFz01b1NegativeSimpleSecondaryVertexHighPurBJetTags+
                akSoftDrop4PFz01b1NegativeCombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz01b1PositiveCombinedSecondaryVertexBJetTags+
                akSoftDrop4PFz01b1NegativeCombinedSecondaryVertexV2BJetTags+
                akSoftDrop4PFz01b1PositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDrop4PFz01b1JetBtaggingMu = cms.Sequence(akSoftDrop4PFz01b1SoftPFMuonsTagInfos * (akSoftDrop4PFz01b1SoftPFMuonBJetTags
                +
                akSoftDrop4PFz01b1SoftPFMuonByIP3dBJetTags
                +
                akSoftDrop4PFz01b1SoftPFMuonByPtBJetTags
                +
                akSoftDrop4PFz01b1NegativeSoftPFMuonByPtBJetTags
                +
                akSoftDrop4PFz01b1PositiveSoftPFMuonByPtBJetTags
              )
            )

akSoftDrop4PFz01b1JetBtagging = cms.Sequence(akSoftDrop4PFz01b1JetBtaggingIP
            *akSoftDrop4PFz01b1JetBtaggingSV
            *akSoftDrop4PFz01b1JetBtaggingNegSV
#            *akSoftDrop4PFz01b1JetBtaggingMu
            )

akSoftDrop4PFz01b1patJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akSoftDrop4PFz01b1Jets"),
        genJetMatch          = cms.InputTag("akSoftDrop4PFz01b1match"),
        genPartonMatch       = cms.InputTag("akSoftDrop4PFz01b1parton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akSoftDrop4PFz01b1corr")),
        #JetPartonMapSource   = cms.InputTag("akSoftDrop4PFz01b1PatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akSoftDrop4PFz01b1PatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akSoftDrop4PFz01b1PatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akSoftDrop4PFz01b1JetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akSoftDrop4PFz01b1SimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1SimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1CombinedSecondaryVertexBJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1CombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1JetBProbabilityBJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1JetProbabilityBJetTags"),
            #cms.InputTag("akSoftDrop4PFz01b1SoftPFMuonByPtBJetTags"),
            #cms.InputTag("akSoftDrop4PFz01b1SoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1TrackCountingHighEffBJetTags"),
            cms.InputTag("akSoftDrop4PFz01b1TrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akSoftDrop4PFz01b1JetID"),
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

akSoftDrop4PFz01b1Njettiness = Njettiness.clone(
		    src = cms.InputTag("akSoftDrop4PFz01b1Jets"),
           	    R0  = cms.double( 0.4)
)
akSoftDrop4PFz01b1patJetsWithBtagging.userData.userFloats.src += ['akSoftDrop4PFz01b1Njettiness:tau1','akSoftDrop4PFz01b1Njettiness:tau2','akSoftDrop4PFz01b1Njettiness:tau3']

akSoftDrop4PFz01b1JetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akSoftDrop4PFz01b1patJetsWithBtagging"),
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
                                                             bTagJetName = cms.untracked.string("akSoftDrop4PFz01b1"),
                                                             jetName = cms.untracked.string("akSoftDrop4PFz01b1"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(False),
							     doSubJets = cms.untracked.bool(True),
                                                             doGenSubJets = cms.untracked.bool(True),     
                                                             subjetGenTag = cms.untracked.InputTag("akSoftDrop4GenJets"),
                                                             doGenTaus = True
                                                            )

akSoftDrop4PFz01b1JetSequence_mc = cms.Sequence(
                                                  #akSoftDrop4PFz01b1clean
                                                  #*
                                                  akSoftDrop4PFz01b1match
                                                  #*
                                                  #akSoftDrop4PFz01b1matchGroomed
                                                  *
                                                  akSoftDrop4PFz01b1parton
                                                  *
                                                  akSoftDrop4PFz01b1corr
                                                  *
                                                  #akSoftDrop4PFz01b1JetID
                                                  #*
                                                  #akSoftDrop4PFz01b1PatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akSoftDrop4PFz01b1PatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akSoftDrop4PFz01b1JetTracksAssociatorAtVertex
                                                  *
                                                  akSoftDrop4PFz01b1JetBtagging
                                                  *
                                                  akSoftDrop4PFz01b1Njettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akSoftDrop4PFz01b1patJetsWithBtagging
                                                  *
                                                  akSoftDrop4PFz01b1JetAnalyzer
                                                  )

akSoftDrop4PFz01b1JetSequence_data = cms.Sequence(akSoftDrop4PFz01b1corr
                                                    *
                                                    #akSoftDrop4PFz01b1JetID
                                                    #*
                                                    akSoftDrop4PFz01b1JetTracksAssociatorAtVertex
                                                    *
                                                    akSoftDrop4PFz01b1JetBtagging
                                                    *
                                                    akSoftDrop4PFz01b1Njettiness 
                                                    *
                                                    akSoftDrop4PFz01b1patJetsWithBtagging
                                                    *
                                                    akSoftDrop4PFz01b1JetAnalyzer
                                                    )

akSoftDrop4PFz01b1JetSequence_jec = cms.Sequence(akSoftDrop4PFz01b1JetSequence_mc)
akSoftDrop4PFz01b1JetSequence_mb = cms.Sequence(akSoftDrop4PFz01b1JetSequence_mc)

akSoftDrop4PFz01b1JetSequence = cms.Sequence(akSoftDrop4PFz01b1JetSequence_mc)
