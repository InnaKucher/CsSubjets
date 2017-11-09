import FWCore.ParameterSet.Config as cms

pfCSVscikitJetTags  = cms.EDProducer(
#pfCombinedCvsLJetTags  = cms.EDProducer(
   "JetTagProducer",
   #jetTagComputer = cms.string('charmTagsComputerCvsL'),
   jetTagComputer = cms.string('CSVscikitTags'),
   tagInfos = cms.VInputTag(
      cms.InputTag('pfImpactParameterTagInfos'),
      #cms.InputTag('pfInclusiveSecondaryVertexFinderCvsLTagInfos'),
      cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfos'),
      #cms.InputTag('softPFMuonsTagInfos'),
      #cms.InputTag('softPFElectronsTagInfos'),
      )
)

pfCSVscikitSubjetTags  = cms.EDProducer(
                                     "JetTagProducer",
                                     jetTagComputer = cms.string('CSVscikitSubjetTags'),
                                     tagInfos = cms.VInputTag(
                                       cms.InputTag('pfSubjetImpactParameterTagInfos'),
                                       cms.InputTag('pfSubjetInclusiveSecondaryVertexFinderTagInfos'),
                                     ) 
                                  )

#pfCSVscikitDummy = pfCSVscikitJetTags.clone(
##pfCombinedCvsBJetTags = pfCombinedCvsLJetTags.clone(
#   #jetTagComputer = cms.string('charmTagsComputerCvsB')
#   jetTagComputer = cms.string('CSVscikitTagsdummy')
#   )
