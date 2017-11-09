

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akCsFilter4PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akCsFilter4PFJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.4
    )

akCsFilter4PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akFilter4HiGenJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.4
    )

akCsFilter4PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsFilter4PFJets")
                                                        )

akCsFilter4PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akCsFilter4PFJets"),
    payload = "AK4PF_offline"
    )

akCsFilter4PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akCsFilter4CaloJets'))

#akCsFilter4PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiSignalGenJets'))

akCsFilter4PFbTagger = bTaggers("akCsFilter4PF",0.4,False,True)

#create objects locally since they dont load properly otherwise
#akCsFilter4PFmatch = akCsFilter4PFbTagger.match
akCsFilter4PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsFilter4PFJets"), matched = cms.InputTag("hiSignalGenParticles"))
akCsFilter4PFPatJetFlavourAssociationLegacy = akCsFilter4PFbTagger.PatJetFlavourAssociationLegacy
akCsFilter4PFPatJetPartons = akCsFilter4PFbTagger.PatJetPartons
akCsFilter4PFJetTracksAssociatorAtVertex = akCsFilter4PFbTagger.JetTracksAssociatorAtVertex
akCsFilter4PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akCsFilter4PFSimpleSecondaryVertexHighEffBJetTags = akCsFilter4PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsFilter4PFSimpleSecondaryVertexHighPurBJetTags = akCsFilter4PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsFilter4PFCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.CombinedSecondaryVertexBJetTags
akCsFilter4PFCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.CombinedSecondaryVertexV2BJetTags
akCsFilter4PFJetBProbabilityBJetTags = akCsFilter4PFbTagger.JetBProbabilityBJetTags
akCsFilter4PFSoftPFMuonByPtBJetTags = akCsFilter4PFbTagger.SoftPFMuonByPtBJetTags
akCsFilter4PFSoftPFMuonByIP3dBJetTags = akCsFilter4PFbTagger.SoftPFMuonByIP3dBJetTags
akCsFilter4PFTrackCountingHighEffBJetTags = akCsFilter4PFbTagger.TrackCountingHighEffBJetTags
akCsFilter4PFTrackCountingHighPurBJetTags = akCsFilter4PFbTagger.TrackCountingHighPurBJetTags
akCsFilter4PFPatJetPartonAssociationLegacy = akCsFilter4PFbTagger.PatJetPartonAssociationLegacy

akCsFilter4PFImpactParameterTagInfos = akCsFilter4PFbTagger.ImpactParameterTagInfos
akCsFilter4PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsFilter4PFJetProbabilityBJetTags = akCsFilter4PFbTagger.JetProbabilityBJetTags

akCsFilter4PFSecondaryVertexTagInfos = akCsFilter4PFbTagger.SecondaryVertexTagInfos
akCsFilter4PFSimpleSecondaryVertexHighEffBJetTags = akCsFilter4PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsFilter4PFSimpleSecondaryVertexHighPurBJetTags = akCsFilter4PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsFilter4PFCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.CombinedSecondaryVertexBJetTags
akCsFilter4PFCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.CombinedSecondaryVertexV2BJetTags

akCsFilter4PFSecondaryVertexNegativeTagInfos = akCsFilter4PFbTagger.SecondaryVertexNegativeTagInfos
akCsFilter4PFNegativeSimpleSecondaryVertexHighEffBJetTags = akCsFilter4PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akCsFilter4PFNegativeSimpleSecondaryVertexHighPurBJetTags = akCsFilter4PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akCsFilter4PFNegativeCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akCsFilter4PFPositiveCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akCsFilter4PFNegativeCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akCsFilter4PFPositiveCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akCsFilter4PFSoftPFMuonsTagInfos = akCsFilter4PFbTagger.SoftPFMuonsTagInfos
akCsFilter4PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsFilter4PFSoftPFMuonBJetTags = akCsFilter4PFbTagger.SoftPFMuonBJetTags
akCsFilter4PFSoftPFMuonByIP3dBJetTags = akCsFilter4PFbTagger.SoftPFMuonByIP3dBJetTags
akCsFilter4PFSoftPFMuonByPtBJetTags = akCsFilter4PFbTagger.SoftPFMuonByPtBJetTags
akCsFilter4PFNegativeSoftPFMuonByPtBJetTags = akCsFilter4PFbTagger.NegativeSoftPFMuonByPtBJetTags
akCsFilter4PFPositiveSoftPFMuonByPtBJetTags = akCsFilter4PFbTagger.PositiveSoftPFMuonByPtBJetTags
akCsFilter4PFPatJetFlavourIdLegacy = cms.Sequence(akCsFilter4PFPatJetPartonAssociationLegacy*akCsFilter4PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub
akCsFilter4PFPatJetFlavourAssociation = akCsFilter4PFbTagger.PatJetFlavourAssociation
akCsFilter4PFPatJetFlavourId = cms.Sequence(akCsFilter4PFPatJetPartons*akCsFilter4PFPatJetFlavourAssociation)

#adding the subjet taggers
akCsFilter4PFSubjetImpactParameterTagInfos = akCsFilter4PFbTagger.SubjetImpactParameterTagInfos
akCsFilter4PFSubjetJetProbabilityBJetTags = akCsFilter4PFbTagger.SubjetJetProbabilityBJetTags
akCsFilter4PFSubjetSecondaryVertexTagInfos = akCsFilter4PFbTagger.SubjetSecondaryVertexTagInfos
akCsFilter4PFSubjetJetTracksAssociatorAtVertex = akCsFilter4PFbTagger.SubjetJetTracksAssociatorAtVertex
#akCsFilter4PFCombinedSubjetSecondaryVertexBJetTags = akCsFilter4PFbTagger.CombinedSubjetSecondaryVertexBJetTags
#akCsFilter4PFCombinedSubjetSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.CombinedSubjetSecondaryVertexV2BJetTags
akCsFilter4PFSubjetCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.SubjetCombinedSecondaryVertexBJetTags
akCsFilter4PFSubjetCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.SubjetCombinedSecondaryVertexV2BJetTags

akCsFilter4PFSubjetNegativeSimpleSecondaryVertexHighEffBJetTags = akCsFilter4PFbTagger.SubjetNegativeSimpleSecondaryVertexHighEffBJetTags
akCsFilter4PFSubjetSimpleSecondaryVertexHighPurBJetTags = akCsFilter4PFbTagger.SubjetSimpleSecondaryVertexHighPurBJetTags
akCsFilter4PFSubjetSecondaryVertexNegativeTagInfos = akCsFilter4PFbTagger.SubjetSecondaryVertexNegativeTagInfos
akCsFilter4PFSubjetSimpleSecondaryVertexHighEffBJetTags = akCsFilter4PFbTagger.SubjetSimpleSecondaryVertexHighEffBJetTags
akCsFilter4PFSubjetNegativeSimpleSecondaryVertexHighPurBJetTags = akCsFilter4PFbTagger.SubjetNegativeSimpleSecondaryVertexHighPurBJetTags
akCsFilter4PFSubjetNegativeCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.SubjetNegativeCombinedSecondaryVertexBJetTags
akCsFilter4PFSubjetPositiveCombinedSecondaryVertexBJetTags = akCsFilter4PFbTagger.SubjetPositiveCombinedSecondaryVertexBJetTags
akCsFilter4PFSubjetNegativeCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.SubjetNegativeCombinedSecondaryVertexV2BJetTags
akCsFilter4PFSubjetPositiveCombinedSecondaryVertexV2BJetTags = akCsFilter4PFbTagger.SubjetPositiveCombinedSecondaryVertexV2BJetTags
akCsFilter4PFSubjetTrackCountingHighEffBJetTags = akCsFilter4PFbTagger.SubjetTrackCountingHighEffBJetTags
akCsFilter4PFSubjetTrackCountingHighPurBJetTags = akCsFilter4PFbTagger.SubjetTrackCountingHighPurBJetTags
akCsFilter4PFSubjetJetBProbabilityBJetTags = akCsFilter4PFbTagger.SubjetJetBProbabilityBJetTags


akCsFilter4PFJetBtaggingIP       = cms.Sequence(akCsFilter4PFImpactParameterTagInfos *
            (akCsFilter4PFTrackCountingHighEffBJetTags +
             akCsFilter4PFTrackCountingHighPurBJetTags +
             akCsFilter4PFJetProbabilityBJetTags +
             akCsFilter4PFJetBProbabilityBJetTags 
            )
            )

akCsFilter4PFSubjetBtaggingIP       = cms.Sequence(akCsFilter4PFSubjetImpactParameterTagInfos *
            (akCsFilter4PFSubjetTrackCountingHighEffBJetTags +
             akCsFilter4PFSubjetTrackCountingHighPurBJetTags +
             akCsFilter4PFSubjetJetProbabilityBJetTags +
             akCsFilter4PFSubjetJetBProbabilityBJetTags
            )
            )


akCsFilter4PFJetBtaggingSV = cms.Sequence(akCsFilter4PFImpactParameterTagInfos
            *
            akCsFilter4PFSecondaryVertexTagInfos
            * (akCsFilter4PFSimpleSecondaryVertexHighEffBJetTags+
                akCsFilter4PFSimpleSecondaryVertexHighPurBJetTags+
                akCsFilter4PFCombinedSecondaryVertexBJetTags+
                akCsFilter4PFCombinedSecondaryVertexV2BJetTags
              )
            )

akCsFilter4PFSubjetBtaggingSV = cms.Sequence(akCsFilter4PFSubjetImpactParameterTagInfos
            *
            akCsFilter4PFSubjetSecondaryVertexTagInfos
            * (akCsFilter4PFSubjetSimpleSecondaryVertexHighEffBJetTags+
                akCsFilter4PFSubjetSimpleSecondaryVertexHighPurBJetTags+
                akCsFilter4PFSubjetCombinedSecondaryVertexBJetTags+
                akCsFilter4PFSubjetCombinedSecondaryVertexV2BJetTags
              )
            )

akCsFilter4PFJetBtaggingNegSV = cms.Sequence(akCsFilter4PFImpactParameterTagInfos
            *
            akCsFilter4PFSecondaryVertexNegativeTagInfos
            * (akCsFilter4PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akCsFilter4PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akCsFilter4PFNegativeCombinedSecondaryVertexBJetTags+
                akCsFilter4PFPositiveCombinedSecondaryVertexBJetTags+
                akCsFilter4PFNegativeCombinedSecondaryVertexV2BJetTags+
                akCsFilter4PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akCsFilter4PFSubjetBtaggingNegSV = cms.Sequence(akCsFilter4PFSubjetImpactParameterTagInfos
            *
            akCsFilter4PFSubjetSecondaryVertexNegativeTagInfos
            * (akCsFilter4PFSubjetNegativeSimpleSecondaryVertexHighEffBJetTags+
                akCsFilter4PFSubjetNegativeSimpleSecondaryVertexHighPurBJetTags+
                akCsFilter4PFSubjetNegativeCombinedSecondaryVertexBJetTags+
                akCsFilter4PFSubjetPositiveCombinedSecondaryVertexBJetTags+
                akCsFilter4PFSubjetNegativeCombinedSecondaryVertexV2BJetTags+
                akCsFilter4PFSubjetPositiveCombinedSecondaryVertexV2BJetTags
              )
            )


akCsFilter4PFJetBtaggingMu = cms.Sequence(akCsFilter4PFSoftPFMuonsTagInfos * (akCsFilter4PFSoftPFMuonBJetTags
                +
                akCsFilter4PFSoftPFMuonByIP3dBJetTags
                +
                akCsFilter4PFSoftPFMuonByPtBJetTags
                +
                akCsFilter4PFNegativeSoftPFMuonByPtBJetTags
                +
                akCsFilter4PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akCsFilter4PFJetBtagging = cms.Sequence(akCsFilter4PFJetBtaggingIP
            *akCsFilter4PFJetBtaggingSV
            *akCsFilter4PFJetBtaggingNegSV
#            *akCsFilter4PFJetBtaggingMu
            )

akCsFilter4PFSubjetBtagging = cms.Sequence(akCsFilter4PFSubjetBtaggingIP
            *akCsFilter4PFSubjetBtaggingSV
            *akCsFilter4PFSubjetBtaggingNegSV
)


akCsFilter4PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akCsFilter4PFJets"),
        genJetMatch          = cms.InputTag("akCsFilter4PFmatch"),
        genPartonMatch       = cms.InputTag("akCsFilter4PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akCsFilter4PFcorr")),
        #JetPartonMapSource   = cms.InputTag("akCsFilter4PFPatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akCsFilter4PFPatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akCsFilter4PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akCsFilter4PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akCsFilter4PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akCsFilter4PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akCsFilter4PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akCsFilter4PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akCsFilter4PFJetBProbabilityBJetTags"),
            cms.InputTag("akCsFilter4PFJetProbabilityBJetTags"),
            #cms.InputTag("akCsFilter4PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akCsFilter4PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akCsFilter4PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akCsFilter4PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akCsFilter4PFJetID"),
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

akCsFilter4PFpatSubjetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akCsFilter4PFJets","Subjets"),
        genJetMatch          = cms.InputTag("akCsFilter4PFmatch"),
        genPartonMatch       = cms.InputTag("akCsFilter4PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("")),
        JetPartonMapSource   = cms.InputTag("akCsFilter4PFPatJetFlavourAssociation"),
        JetFlavourInfoSource   = cms.InputTag("akCsFilter4PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akCsFilter4PFSubjetJetTracksAssociatorAtVertex"),
        useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akCsFilter4PFSubjetSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akCsFilter4PFSubjetSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akCsFilter4PFSubjetCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akCsFilter4PFSubjetCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akCsFilter4PFSubjetJetBProbabilityBJetTags"),
            cms.InputTag("akCsFilter4PFSubjetJetProbabilityBJetTags"),
            cms.InputTag("akCsFilter4PFSubjetTrackCountingHighEffBJetTags"),
            cms.InputTag("akCsFilter4PFSubjetTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akCsFilter4PFJetID"),
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = True,
        addGenPartonMatch = False,
        addGenJetMatch = False,
        embedGenJetMatch = False,
        embedGenPartonMatch = False,
        )


akCsFilter4PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akCsFilter4PFJets"),
           	    R0  = cms.double( 0.4)
)
akCsFilter4PFpatJetsWithBtagging.userData.userFloats.src += ['akCsFilter4PFNjettiness:tau1','akCsFilter4PFNjettiness:tau2','akCsFilter4PFNjettiness:tau3']

akCsFilter4PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akCsFilter4PFpatJetsWithBtagging"),
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
                                                             bTagJetName = cms.untracked.string("akCsFilter4PF"),
                                                             jetName = cms.untracked.string("akCsFilter4PF"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(True),
							     doSubJets = cms.untracked.bool(True),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("akFilter4GenJets"),
							     doExtendedFlavorTagging = cms.untracked.bool(True),
							     jetFlavourInfos = cms.InputTag("akCsFilter4PFPatJetFlavourAssociation"),
							     subjetFlavourInfos = cms.InputTag("akCsFilter4PFPatJetFlavourAssociation","SubJets"),
							     groomedJets = cms.InputTag("akCsFilter4PFJets"),
							     isPythia6 = cms.untracked.bool(False),
                                                             doGenTaus = True
                                                            )

akCsFilter4PFJetSequence_mc = cms.Sequence(
                                                  #akCsFilter4PFclean
                                                  #*
                                                  akCsFilter4PFmatch
                                                  #*
                                                  #akCsFilter4PFmatchGroomed
                                                  *
                                                  akCsFilter4PFparton
                                                  *
                                                  akCsFilter4PFcorr
                                                  *
                                                  #akCsFilter4PFJetID
                                                  #*
                                                  #akCsFilter4PFPatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akCsFilter4PFPatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akCsFilter4PFJetTracksAssociatorAtVertex
                                                  *
                                                  akCsFilter4PFJetBtagging
                                                  *
                                                  akCsFilter4PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akCsFilter4PFpatJetsWithBtagging
                                                  *
                                                  akCsFilter4PFJetAnalyzer
                                                  *
                                                  akCsFilter4PFSubjetJetTracksAssociatorAtVertex
                                                  *
                                                  akCsFilter4PFSubjetBtagging
                                                  *
                                                  akCsFilter4PFpatSubjetsWithBtagging
                                                  )

akCsFilter4PFJetSequence_data = cms.Sequence(akCsFilter4PFcorr
                                                    *
                                                    #akCsFilter4PFJetID
                                                    #*
                                                    akCsFilter4PFJetTracksAssociatorAtVertex
                                                    *
                                                    akCsFilter4PFJetBtagging
                                                    *
                                                    akCsFilter4PFNjettiness 
                                                    *
                                                    akCsFilter4PFpatJetsWithBtagging
                                                    *
                                                    akCsFilter4PFJetAnalyzer
                                                    )

akCsFilter4PFJetSequence_jec = cms.Sequence(akCsFilter4PFJetSequence_mc)
akCsFilter4PFJetSequence_mb = cms.Sequence(akCsFilter4PFJetSequence_mc)

akCsFilter4PFJetSequence = cms.Sequence(akCsFilter4PFJetSequence_mc)
