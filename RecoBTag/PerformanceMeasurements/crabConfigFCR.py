from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'btagana_fcr170'

config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runForestAOD_PbPb_MIX_75X.py'
config.JobType.maxMemoryMB = 2400


config.section_('Data')
config.Data.inputDataset ='/Pythia6_bJetFCR170_pp502_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'

config.Data.inputDBS = 'global'#'phys03'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 200
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/%s' % (getUsernameFromSiteDB(),config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

config.section_('Site')
config.Site.storageSite = 'T2_FR_GRIF_LLR'
#config.Site.whitelist = ['T2_US_MIT']
