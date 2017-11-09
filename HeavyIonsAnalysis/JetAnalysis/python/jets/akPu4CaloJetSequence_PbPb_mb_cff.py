

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akPu4Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPu4CaloJets"),
    matched = cms.InputTag("ak4HiCleanedGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akPu4CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("ak4HiGenJets"),
    matched = cms.InputTag("ak4HiCleanedGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akPu4Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akPu4CaloJets")
                                                        )

akPu4Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akPu4CaloJets"),
    payload = "AKPu4Calo_offline"
    )

akPu4CaloJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akPu4CaloJets'))

#akPu4Caloclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiCleanedGenJets'))

akPu4CalobTagger = bTaggers("akPu4Calo",0.4,False,False)

#create objects locally since they dont load properly otherwise
#akPu4Calomatch = akPu4CalobTagger.match
akPu4Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akPu4CaloJets"), matched = cms.InputTag("selectedPartons"))
akPu4CaloPatJetFlavourAssociationLegacy = akPu4CalobTagger.PatJetFlavourAssociationLegacy
akPu4CaloPatJetPartons = akPu4CalobTagger.PatJetPartons
akPu4CaloJetTracksAssociatorAtVertex = akPu4CalobTagger.JetTracksAssociatorAtVertex
akPu4CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akPu4CaloSimpleSecondaryVertexHighEffBJetTags = akPu4CalobTagger.SimpleSecondaryVertexHighEffBJetTags
akPu4CaloSimpleSecondaryVertexHighPurBJetTags = akPu4CalobTagger.SimpleSecondaryVertexHighPurBJetTags
akPu4CaloCombinedSecondaryVertexBJetTags = akPu4CalobTagger.CombinedSecondaryVertexBJetTags
akPu4CaloCombinedSecondaryVertexV2BJetTags = akPu4CalobTagger.CombinedSecondaryVertexV2BJetTags
akPu4CaloJetBProbabilityBJetTags = akPu4CalobTagger.JetBProbabilityBJetTags
akPu4CaloSoftPFMuonByPtBJetTags = akPu4CalobTagger.SoftPFMuonByPtBJetTags
akPu4CaloSoftPFMuonByIP3dBJetTags = akPu4CalobTagger.SoftPFMuonByIP3dBJetTags
akPu4CaloTrackCountingHighEffBJetTags = akPu4CalobTagger.TrackCountingHighEffBJetTags
akPu4CaloTrackCountingHighPurBJetTags = akPu4CalobTagger.TrackCountingHighPurBJetTags
akPu4CaloPatJetPartonAssociationLegacy = akPu4CalobTagger.PatJetPartonAssociationLegacy

akPu4CaloImpactParameterTagInfos = akPu4CalobTagger.ImpactParameterTagInfos
akPu4CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akPu4CaloJetProbabilityBJetTags = akPu4CalobTagger.JetProbabilityBJetTags

akPu4CaloSecondaryVertexTagInfos = akPu4CalobTagger.SecondaryVertexTagInfos
akPu4CaloSimpleSecondaryVertexHighEffBJetTags = akPu4CalobTagger.SimpleSecondaryVertexHighEffBJetTags
akPu4CaloSimpleSecondaryVertexHighPurBJetTags = akPu4CalobTagger.SimpleSecondaryVertexHighPurBJetTags
akPu4CaloCombinedSecondaryVertexBJetTags = akPu4CalobTagger.CombinedSecondaryVertexBJetTags
akPu4CaloCombinedSecondaryVertexV2BJetTags = akPu4CalobTagger.CombinedSecondaryVertexV2BJetTags

akPu4CaloSecondaryVertexNegativeTagInfos = akPu4CalobTagger.SecondaryVertexNegativeTagInfos
akPu4CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akPu4CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akPu4CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akPu4CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akPu4CaloNegativeCombinedSecondaryVertexBJetTags = akPu4CalobTagger.NegativeCombinedSecondaryVertexBJetTags
akPu4CaloPositiveCombinedSecondaryVertexBJetTags = akPu4CalobTagger.PositiveCombinedSecondaryVertexBJetTags
akPu4CaloNegativeCombinedSecondaryVertexV2BJetTags = akPu4CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
akPu4CaloPositiveCombinedSecondaryVertexV2BJetTags = akPu4CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags

akPu4CaloSoftPFMuonsTagInfos = akPu4CalobTagger.SoftPFMuonsTagInfos
akPu4CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akPu4CaloSoftPFMuonBJetTags = akPu4CalobTagger.SoftPFMuonBJetTags
akPu4CaloSoftPFMuonByIP3dBJetTags = akPu4CalobTagger.SoftPFMuonByIP3dBJetTags
akPu4CaloSoftPFMuonByPtBJetTags = akPu4CalobTagger.SoftPFMuonByPtBJetTags
akPu4CaloNegativeSoftPFMuonByPtBJetTags = akPu4CalobTagger.NegativeSoftPFMuonByPtBJetTags
akPu4CaloPositiveSoftPFMuonByPtBJetTags = akPu4CalobTagger.PositiveSoftPFMuonByPtBJetTags
akPu4CaloPatJetFlavourIdLegacy = cms.Sequence(akPu4CaloPatJetPartonAssociationLegacy*akPu4CaloPatJetFlavourAssociationLegacy)
#Not working with our PU sub
akPu4CaloPatJetFlavourAssociation = akPu4CalobTagger.PatJetFlavourAssociation
akPu4CaloPatJetFlavourId = cms.Sequence(akPu4CaloPatJetPartons*akPu4CaloPatJetFlavourAssociation)

#adding the subjet taggers
#SUBJETDUMMY_akPu4CaloSubjetImpactParameterTagInfos = akPu4CalobTagger.SubjetImpactParameterTagInfos
#SUBJETDUMMY_akPu4CaloSubjetJetProbabilityBJetTags = akPu4CalobTagger.SubjetJetProbabilityBJetTags
#SUBJETDUMMY_akPu4CaloSubjetSecondaryVertexTagInfos = akPu4CalobTagger.SubjetSecondaryVertexTagInfos
#SUBJETDUMMY_akPu4CaloSubjetJetTracksAssociatorAtVertex = akPu4CalobTagger.SubjetJetTracksAssociatorAtVertex
#SUBJETDUMMY_akPu4CaloCombinedSubjetSecondaryVertexBJetTags = akPu4CalobTagger.CombinedSubjetSecondaryVertexBJetTags
#SUBJETDUMMY_akPu4CaloCombinedSubjetSecondaryVertexV2BJetTags = akPu4CalobTagger.CombinedSubjetSecondaryVertexV2BJetTags

akPu4CaloJetBtaggingIP       = cms.Sequence(akPu4CaloImpactParameterTagInfos *
            (akPu4CaloTrackCountingHighEffBJetTags +
             akPu4CaloTrackCountingHighPurBJetTags +
             akPu4CaloJetProbabilityBJetTags +
             akPu4CaloJetBProbabilityBJetTags 
            )
            )

akPu4CaloJetBtaggingSV = cms.Sequence(akPu4CaloImpactParameterTagInfos
            *
            akPu4CaloSecondaryVertexTagInfos
            * (akPu4CaloSimpleSecondaryVertexHighEffBJetTags+
                akPu4CaloSimpleSecondaryVertexHighPurBJetTags+
                akPu4CaloCombinedSecondaryVertexBJetTags+
                akPu4CaloCombinedSecondaryVertexV2BJetTags
              )
            )

akPu4CaloJetBtaggingNegSV = cms.Sequence(akPu4CaloImpactParameterTagInfos
            *
            akPu4CaloSecondaryVertexNegativeTagInfos
            * (akPu4CaloNegativeSimpleSecondaryVertexHighEffBJetTags+
                akPu4CaloNegativeSimpleSecondaryVertexHighPurBJetTags+
                akPu4CaloNegativeCombinedSecondaryVertexBJetTags+
                akPu4CaloPositiveCombinedSecondaryVertexBJetTags+
                akPu4CaloNegativeCombinedSecondaryVertexV2BJetTags+
                akPu4CaloPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akPu4CaloJetBtaggingMu = cms.Sequence(akPu4CaloSoftPFMuonsTagInfos * (akPu4CaloSoftPFMuonBJetTags
                +
                akPu4CaloSoftPFMuonByIP3dBJetTags
                +
                akPu4CaloSoftPFMuonByPtBJetTags
                +
                akPu4CaloNegativeSoftPFMuonByPtBJetTags
                +
                akPu4CaloPositiveSoftPFMuonByPtBJetTags
              )
            )

akPu4CaloJetBtagging = cms.Sequence(akPu4CaloJetBtaggingIP
            *akPu4CaloJetBtaggingSV
            *akPu4CaloJetBtaggingNegSV
#            *akPu4CaloJetBtaggingMu
            )

akPu4CalopatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akPu4CaloJets"),
        genJetMatch          = cms.InputTag("akPu4Calomatch"),
        genPartonMatch       = cms.InputTag("akPu4Caloparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu4Calocorr")),
        #JetPartonMapSource   = cms.InputTag("akPu4CaloPatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akPu4CaloPatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akPu4CaloPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akPu4CaloJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akPu4CaloSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akPu4CaloSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akPu4CaloCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akPu4CaloCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akPu4CaloJetBProbabilityBJetTags"),
            cms.InputTag("akPu4CaloJetProbabilityBJetTags"),
            #cms.InputTag("akPu4CaloSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akPu4CaloSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akPu4CaloTrackCountingHighEffBJetTags"),
            cms.InputTag("akPu4CaloTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akPu4CaloJetID"),
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

akPu4CaloNjettiness = Njettiness.clone(
		    src = cms.InputTag("akPu4CaloJets"),
           	    R0  = cms.double( 0.4)
)
akPu4CalopatJetsWithBtagging.userData.userFloats.src += ['akPu4CaloNjettiness:tau1','akPu4CaloNjettiness:tau2','akPu4CaloNjettiness:tau3']

akPu4CaloJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akPu4CalopatJetsWithBtagging"),
                                                             genjetTag = 'ak4HiGenJets',
                                                             rParam = 0.4,
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
                                                             bTagJetName = cms.untracked.string("akPu4Calo"),
                                                             jetName = cms.untracked.string("akPu4Calo"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(True),
							     doSubJets = cms.untracked.bool(False),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("ak4GenJets"),
							     doExtendedFlavorTagging = cms.untracked.bool(True),
							     jetFlavourInfos = cms.InputTag("akPu4CaloPatJetFlavourAssociation"),
							     subjetFlavourInfos = cms.InputTag("akPu4CaloPatJetFlavourAssociation","SubJets"),
							     groomedJets = cms.InputTag("akPu4CaloJets"),
							     isPythia6 = cms.untracked.bool(False),
                                                             doGenTaus = True
                                                            )

akPu4CaloJetSequence_mc = cms.Sequence(
                                                  #akPu4Caloclean
                                                  #*
                                                  akPu4Calomatch
                                                  #*
                                                  #akPu4CalomatchGroomed
                                                  *
                                                  akPu4Caloparton
                                                  *
                                                  akPu4Calocorr
                                                  *
                                                  #akPu4CaloJetID
                                                  #*
                                                  #akPu4CaloPatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akPu4CaloPatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akPu4CaloJetTracksAssociatorAtVertex
                                                  *
                                                  akPu4CaloJetBtagging
                                                  *
                                                  akPu4CaloNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akPu4CalopatJetsWithBtagging
                                                  *
                                                  akPu4CaloJetAnalyzer
                                                  )

akPu4CaloJetSequence_data = cms.Sequence(akPu4Calocorr
                                                    *
                                                    #akPu4CaloJetID
                                                    #*
                                                    akPu4CaloJetTracksAssociatorAtVertex
                                                    *
                                                    akPu4CaloJetBtagging
                                                    *
                                                    akPu4CaloNjettiness 
                                                    *
                                                    akPu4CalopatJetsWithBtagging
                                                    *
                                                    akPu4CaloJetAnalyzer
                                                    )

akPu4CaloJetSequence_jec = cms.Sequence(akPu4CaloJetSequence_mc)
akPu4CaloJetSequence_mb = cms.Sequence(akPu4CaloJetSequence_mc)

akPu4CaloJetSequence = cms.Sequence(akPu4CaloJetSequence_mb)
