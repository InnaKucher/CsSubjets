

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akVs6Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akVs6CaloJets"),
    matched = cms.InputTag("ak6HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.6
    )

akVs6CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("ak6HiGenJets"),
    matched = cms.InputTag("ak6HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.6
    )

akVs6Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akVs6CaloJets")
                                                        )

akVs6Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akVs6CaloJets"),
    payload = "AK6Calo_offline"
    )

akVs6CaloJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akVs6CaloJets'))

#akVs6Caloclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak6HiSignalGenJets'))

akVs6CalobTagger = bTaggers("akVs6Calo",0.6,False,False)

#create objects locally since they dont load properly otherwise
#akVs6Calomatch = akVs6CalobTagger.match
akVs6Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akVs6CaloJets"), matched = cms.InputTag("hiSignalGenParticles"))
akVs6CaloPatJetFlavourAssociationLegacy = akVs6CalobTagger.PatJetFlavourAssociationLegacy
akVs6CaloPatJetPartons = akVs6CalobTagger.PatJetPartons
akVs6CaloJetTracksAssociatorAtVertex = akVs6CalobTagger.JetTracksAssociatorAtVertex
akVs6CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akVs6CaloSimpleSecondaryVertexHighEffBJetTags = akVs6CalobTagger.SimpleSecondaryVertexHighEffBJetTags
akVs6CaloSimpleSecondaryVertexHighPurBJetTags = akVs6CalobTagger.SimpleSecondaryVertexHighPurBJetTags
akVs6CaloCombinedSecondaryVertexBJetTags = akVs6CalobTagger.CombinedSecondaryVertexBJetTags
akVs6CaloCombinedSecondaryVertexV2BJetTags = akVs6CalobTagger.CombinedSecondaryVertexV2BJetTags
akVs6CaloJetBProbabilityBJetTags = akVs6CalobTagger.JetBProbabilityBJetTags
akVs6CaloSoftPFMuonByPtBJetTags = akVs6CalobTagger.SoftPFMuonByPtBJetTags
akVs6CaloSoftPFMuonByIP3dBJetTags = akVs6CalobTagger.SoftPFMuonByIP3dBJetTags
akVs6CaloTrackCountingHighEffBJetTags = akVs6CalobTagger.TrackCountingHighEffBJetTags
akVs6CaloTrackCountingHighPurBJetTags = akVs6CalobTagger.TrackCountingHighPurBJetTags
akVs6CaloPatJetPartonAssociationLegacy = akVs6CalobTagger.PatJetPartonAssociationLegacy

akVs6CaloImpactParameterTagInfos = akVs6CalobTagger.ImpactParameterTagInfos
akVs6CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akVs6CaloJetProbabilityBJetTags = akVs6CalobTagger.JetProbabilityBJetTags

akVs6CaloSecondaryVertexTagInfos = akVs6CalobTagger.SecondaryVertexTagInfos
akVs6CaloSimpleSecondaryVertexHighEffBJetTags = akVs6CalobTagger.SimpleSecondaryVertexHighEffBJetTags
akVs6CaloSimpleSecondaryVertexHighPurBJetTags = akVs6CalobTagger.SimpleSecondaryVertexHighPurBJetTags
akVs6CaloCombinedSecondaryVertexBJetTags = akVs6CalobTagger.CombinedSecondaryVertexBJetTags
akVs6CaloCombinedSecondaryVertexV2BJetTags = akVs6CalobTagger.CombinedSecondaryVertexV2BJetTags

akVs6CaloSecondaryVertexNegativeTagInfos = akVs6CalobTagger.SecondaryVertexNegativeTagInfos
akVs6CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akVs6CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akVs6CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akVs6CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akVs6CaloNegativeCombinedSecondaryVertexBJetTags = akVs6CalobTagger.NegativeCombinedSecondaryVertexBJetTags
akVs6CaloPositiveCombinedSecondaryVertexBJetTags = akVs6CalobTagger.PositiveCombinedSecondaryVertexBJetTags
akVs6CaloNegativeCombinedSecondaryVertexV2BJetTags = akVs6CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
akVs6CaloPositiveCombinedSecondaryVertexV2BJetTags = akVs6CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags

akVs6CaloSoftPFMuonsTagInfos = akVs6CalobTagger.SoftPFMuonsTagInfos
akVs6CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akVs6CaloSoftPFMuonBJetTags = akVs6CalobTagger.SoftPFMuonBJetTags
akVs6CaloSoftPFMuonByIP3dBJetTags = akVs6CalobTagger.SoftPFMuonByIP3dBJetTags
akVs6CaloSoftPFMuonByPtBJetTags = akVs6CalobTagger.SoftPFMuonByPtBJetTags
akVs6CaloNegativeSoftPFMuonByPtBJetTags = akVs6CalobTagger.NegativeSoftPFMuonByPtBJetTags
akVs6CaloPositiveSoftPFMuonByPtBJetTags = akVs6CalobTagger.PositiveSoftPFMuonByPtBJetTags
akVs6CaloPatJetFlavourIdLegacy = cms.Sequence(akVs6CaloPatJetPartonAssociationLegacy*akVs6CaloPatJetFlavourAssociationLegacy)
#Not working with our PU sub
akVs6CaloPatJetFlavourAssociation = akVs6CalobTagger.PatJetFlavourAssociation
akVs6CaloPatJetFlavourId = cms.Sequence(akVs6CaloPatJetPartons*akVs6CaloPatJetFlavourAssociation)

#adding the subjet taggers
#SUBJETDUMMY_akVs6CaloSubjetImpactParameterTagInfos = akVs6CalobTagger.SubjetImpactParameterTagInfos
#SUBJETDUMMY_akVs6CaloSubjetJetProbabilityBJetTags = akVs6CalobTagger.SubjetJetProbabilityBJetTags
#SUBJETDUMMY_akVs6CaloSubjetSecondaryVertexTagInfos = akVs6CalobTagger.SubjetSecondaryVertexTagInfos
#SUBJETDUMMY_akVs6CaloSubjetJetTracksAssociatorAtVertex = akVs6CalobTagger.SubjetJetTracksAssociatorAtVertex
#SUBJETDUMMY_akVs6CaloCombinedSubjetSecondaryVertexBJetTags = akVs6CalobTagger.CombinedSubjetSecondaryVertexBJetTags
#SUBJETDUMMY_akVs6CaloCombinedSubjetSecondaryVertexV2BJetTags = akVs6CalobTagger.CombinedSubjetSecondaryVertexV2BJetTags

akVs6CaloJetBtaggingIP       = cms.Sequence(akVs6CaloImpactParameterTagInfos *
            (akVs6CaloTrackCountingHighEffBJetTags +
             akVs6CaloTrackCountingHighPurBJetTags +
             akVs6CaloJetProbabilityBJetTags +
             akVs6CaloJetBProbabilityBJetTags 
            )
            )

akVs6CaloJetBtaggingSV = cms.Sequence(akVs6CaloImpactParameterTagInfos
            *
            akVs6CaloSecondaryVertexTagInfos
            * (akVs6CaloSimpleSecondaryVertexHighEffBJetTags+
                akVs6CaloSimpleSecondaryVertexHighPurBJetTags+
                akVs6CaloCombinedSecondaryVertexBJetTags+
                akVs6CaloCombinedSecondaryVertexV2BJetTags
              )
            )

akVs6CaloJetBtaggingNegSV = cms.Sequence(akVs6CaloImpactParameterTagInfos
            *
            akVs6CaloSecondaryVertexNegativeTagInfos
            * (akVs6CaloNegativeSimpleSecondaryVertexHighEffBJetTags+
                akVs6CaloNegativeSimpleSecondaryVertexHighPurBJetTags+
                akVs6CaloNegativeCombinedSecondaryVertexBJetTags+
                akVs6CaloPositiveCombinedSecondaryVertexBJetTags+
                akVs6CaloNegativeCombinedSecondaryVertexV2BJetTags+
                akVs6CaloPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akVs6CaloJetBtaggingMu = cms.Sequence(akVs6CaloSoftPFMuonsTagInfos * (akVs6CaloSoftPFMuonBJetTags
                +
                akVs6CaloSoftPFMuonByIP3dBJetTags
                +
                akVs6CaloSoftPFMuonByPtBJetTags
                +
                akVs6CaloNegativeSoftPFMuonByPtBJetTags
                +
                akVs6CaloPositiveSoftPFMuonByPtBJetTags
              )
            )

akVs6CaloJetBtagging = cms.Sequence(akVs6CaloJetBtaggingIP
            *akVs6CaloJetBtaggingSV
            *akVs6CaloJetBtaggingNegSV
#            *akVs6CaloJetBtaggingMu
            )

akVs6CalopatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akVs6CaloJets"),
        genJetMatch          = cms.InputTag("akVs6Calomatch"),
        genPartonMatch       = cms.InputTag("akVs6Caloparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akVs6Calocorr")),
        #JetPartonMapSource   = cms.InputTag("akVs6CaloPatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akVs6CaloPatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akVs6CaloPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akVs6CaloJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akVs6CaloSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akVs6CaloSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akVs6CaloCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akVs6CaloCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akVs6CaloJetBProbabilityBJetTags"),
            cms.InputTag("akVs6CaloJetProbabilityBJetTags"),
            #cms.InputTag("akVs6CaloSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akVs6CaloSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akVs6CaloTrackCountingHighEffBJetTags"),
            cms.InputTag("akVs6CaloTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akVs6CaloJetID"),
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

akVs6CaloNjettiness = Njettiness.clone(
		    src = cms.InputTag("akVs6CaloJets"),
           	    R0  = cms.double( 0.6)
)
akVs6CalopatJetsWithBtagging.userData.userFloats.src += ['akVs6CaloNjettiness:tau1','akVs6CaloNjettiness:tau2','akVs6CaloNjettiness:tau3']

akVs6CaloJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akVs6CalopatJetsWithBtagging"),
                                                             genjetTag = 'ak6HiGenJets',
                                                             rParam = 0.6,
                                                             matchJets = cms.untracked.bool(False),
                                                             matchTag = 'patJetsWithBtagging',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
                                                             trackTag = cms.InputTag("hiGeneralTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
							     doSubEvent = True,
                                                             useHepMC = cms.untracked.bool(False),
							     genParticles = cms.untracked.InputTag("genParticles"),
							     eventInfoTag = cms.InputTag("generator"),
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akVs6Calo"),
                                                             jetName = cms.untracked.string("akVs6Calo"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(True),
							     doSubJets = cms.untracked.bool(False),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("ak6GenJets"),
							     doExtendedFlavorTagging = cms.untracked.bool(True),
							     jetFlavourInfos = cms.InputTag("akVs6CaloPatJetFlavourAssociation"),
							     subjetFlavourInfos = cms.InputTag("akVs6CaloPatJetFlavourAssociation","SubJets"),
							     groomedJets = cms.InputTag("akVs6CaloJets"),
							     isPythia6 = cms.untracked.bool(False),
                                                             doGenTaus = True
                                                            )

akVs6CaloJetSequence_mc = cms.Sequence(
                                                  #akVs6Caloclean
                                                  #*
                                                  akVs6Calomatch
                                                  #*
                                                  #akVs6CalomatchGroomed
                                                  *
                                                  akVs6Caloparton
                                                  *
                                                  akVs6Calocorr
                                                  *
                                                  #akVs6CaloJetID
                                                  #*
                                                  #akVs6CaloPatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akVs6CaloPatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akVs6CaloJetTracksAssociatorAtVertex
                                                  *
                                                  akVs6CaloJetBtagging
                                                  *
                                                  akVs6CaloNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akVs6CalopatJetsWithBtagging
                                                  *
                                                  akVs6CaloJetAnalyzer
                                                  )

akVs6CaloJetSequence_data = cms.Sequence(akVs6Calocorr
                                                    *
                                                    #akVs6CaloJetID
                                                    #*
                                                    akVs6CaloJetTracksAssociatorAtVertex
                                                    *
                                                    akVs6CaloJetBtagging
                                                    *
                                                    akVs6CaloNjettiness 
                                                    *
                                                    akVs6CalopatJetsWithBtagging
                                                    *
                                                    akVs6CaloJetAnalyzer
                                                    )

akVs6CaloJetSequence_jec = cms.Sequence(akVs6CaloJetSequence_mc)
akVs6CaloJetSequence_mb = cms.Sequence(akVs6CaloJetSequence_mc)

akVs6CaloJetSequence = cms.Sequence(akVs6CaloJetSequence_mc)
