from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'btag_csJets_fcr170v1'

config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runForestAOD_PbPb_MIX_75X_CsJets.py'
config.JobType.maxMemoryMB = 2400


config.section_('Data')
config.Data.inputDataset ='/Pythia6_bJetFCR170_pp502_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'

config.Data.inputDBS = 'global'#'phys03'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 200
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/%s' % (getUsernameFromSiteDB(),config.General.requestName)
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/ikucher/btag_csJets_fcr170v1'
config.Data.outputDatasetTag = config.General.requestName

config.section_('Site')
config.Site.storageSite = 'T2_FR_CCIN2P3'
config.Site.whitelist = ['T2_FR_CCIN2P3']
config.Site.whitelist = ['T2_ES_CIEMAT']
config.Site.whitelist = ['T2_FR_GRIF_LLR']
config.Site.whitelist = ['T2_US_Purdue']

#config.section_('Site')
#config.Site.storageSite = 'T2_FR_GRIF_LLR'

