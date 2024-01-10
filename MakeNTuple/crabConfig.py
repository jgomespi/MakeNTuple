#if __name__ == '__main__':

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from CRABClient.UserUtilities import config 
#from httplib import HTTPException
from http.client import HTTPException
from multiprocessing import Process
#from WMCore.Configuration import Configuration
config = config()

config.section_("General")
config.General.workArea = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Analysis/test/Ntuple_build/CMSSW_10_6_26/src/MakeNtuple/MakeNTuple/crab_projects/'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'
#config.JobType.pyCfgParams=["tag=MC_preVFP","maxEvents=-1"] 
config.JobType.inputFiles = ['xangle.root','PreliminaryEfficiencies_July132020_1D2DMultiTrack.root','pixelEfficiencies_multiRP_reMiniAOD.root','pixelEfficiencies_singleRP.root','pixelEfficiencies_multiRP.root','pixelEfficiencies_radiation_reMiniAOD.root','pixelEfficiencies_radiation.root','Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root', 'Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root', 'Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root', 'Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root', 'Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root', 'Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root', 'Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root', 'Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root', 'egammaEffi.txt_EGM2D_Tight_UL17.root', 'egammaEffi.txt_Ele_Tight_EGM2D.root', 'egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root', 'egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root','direct_simu_reco2016_cff.py','direct_simu_reco2017_cff.py','direct_simu_reco2018_cff.py','RoccoR2017UL.txt','RoccoR2018UL.txt','MyDataPileupHistogram.root', 'PileupMC.root','RoccoR2016aUL.txt','RoccoR2016bUL.txt' ]
config.JobType.outputFiles = ['out.root']
config.JobType.allowUndistributedCMSSW = True	
config.JobType.maxMemoryMB = 2500
#config.JobType.numCores = 8

config.section_("Data")
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'  #'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/jgomespi' #%s/' % (getUsernameFromSiteDB())
config.Data.publication = False 
#config.Data.outputDatasetTag = 'MC'
	
config.section_("Site")
config.Site.storageSite = 'T2_BR_UERJ'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

tags = [["tag=MC_preVFP","maxEvents=-1"],
        ["tag=MC_postVFP","maxEvents=-1"],
        ["tag=MC_2017","maxEvents=-1"],
        ["tag=MC_2018","maxEvents=-1"],
        ["tag=data_preVFP","maxEvents=-1"],
        ["tag=data_postVFP","maxEvents=-1"],
        ["tag=data_2017","maxEvents=-1"],
        ["tag=data_2018","maxEvents=-1"]]

# Datasets obtained at https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis

MC_preVFP = [['DYJetsToLL_M-50_preVFP','/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['WW_preVFP','/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['WZ_preVFP','/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['ZZ_preVFP','/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GGToMuMu_Elastic_preVFP','/GGToMuMu_Pt-25_Elastic_13TeV-lpair/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GGToMuMu_Inel-El_preVFP','/GGToMuMu_Pt-25_Inel-El_13TeV-lpair/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GGToMuMu_Inel-Inel_preVFP','/GGToMuMu_Pt-25_Inel-Inel_13TeV-lpair/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GGToEE_Elastic_preVFP','/GGToEE_Pt-35_Elastic_13TeV-lpair/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GGToEE_Inel-El_preVFP','/GGToEE_Pt-35_Inel-El_13TeV-lpair/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GGToEE_Inel-Inel_preVFP','/GGToEE_Pt-35_Inel-Inel_13TeV-lpair/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],  
  ['GluGluToWWToENEN_preVFP','/GluGluToWWToENEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GluGluToWWToENMN_preVFP','/GluGluToWWToENMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GluGluToWWToMNEN_preVFP','/GluGluToWWToMNEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['GluGluToWWToMNMN_preVFP','/GluGluToWWToMNMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM'],
  ['WJetsToLNu_preVFP','/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM']]

MC_postVFP = [['DYJetsToLL_M-50_postVFP','/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['WW_postVFP','/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['WZ_postVFP','/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['ZZ_postVFP','/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GGToMuMu_Elastic_postVFP','/GGToMuMu_Pt-25_Elastic_13TeV-lpair/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GGToMuMu_Inel-El_postVFP','/GGToMuMu_Pt-25_Inel-El_13TeV-lpair/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GGToMuMu_Inel-Inel_postVFP','/GGToMuMu_Pt-25_Inel-Inel_13TeV-lpair/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GGToEE_Elastic_postVFP','/GGToEE_Pt-35_Elastic_13TeV-lpair/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GGToEE_Inel-El_postVFP','/GGToEE_Pt-35_Inel-El_13TeV-lpair/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GGToEE_Inel-Inel_postVFP','/GGToEE_Pt-35_Inel-Inel_13TeV-lpair/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GluGluToWWToENEN_postVFP','/GluGluToWWToENEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GluGluToWWToENMN_postVFP','/GluGluToWWToENMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GluGluToWWToMNEN_postVFP','/GluGluToWWToMNEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['GluGluToWWToMNMN_postVFP','/GluGluToWWToMNMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM'],
	   ['WJetsToLNu_postVFP','/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM']]

MC_2017 = [['DYJetsToLL_M-50_2017','/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['WW_2017','/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'],
	   ['WZ_2017','/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'],
	   ['ZZ_2017','/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'],
	   ['GGToMuMu_Elastic_2017','/GGToMuMu_Pt-25_Elastic_13TeV-lpair/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['GGToMuMu_Inel-El_2017','/GGToMuMu_Pt-25_Inel-El_13TeV-lpair/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['GGToMuMu_Inel-Inel_2017','/GGToMuMu_Pt-25_Inel-Inel_13TeV-lpair/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['GGToEE_Elastic_2017','/GGToEE_Pt-35_Elastic_13TeV-lpair/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['GGToEE_Inel-El_2017','/GGToEE_Pt-35_Inel-El_13TeV-lpair/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['GGToEE_Inel-Inel_2017','/GGToEE_Pt-35_Inel-Inel_13TeV-lpair/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'],
	   ['GluGluToWWToENEN_2017','/GluGluToWWToENEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'],
	   ['GluGluToWWToENMN_2017','/GluGluToWWToENMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'],
	   ['GluGluToWWToMNEN_2017','/GluGluToWWToMNEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'],
	   ['GluGluToWWToMNMN_2017','/GluGluToWWToMNMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'], 
	   ['WJetsToLNu_2017','/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM']]

MC_2018 = [['DYJetsToLL_M-50_2018','/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['WW_2018','/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'],
   ['WZ_2018','/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'],
   ['ZZ_2018','/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GGToMuMu_Elastic_2018','/GGToMuMu_Pt-25_Elastic_13TeV-lpair/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GGToMuMu_Inel-El_2018','/GGToMuMu_Pt-25_Inel-El_13TeV-lpair/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GGToMuMu_Inel-Inel_2018','/GGToMuMu_Pt-25_Inel-Inel_13TeV-lpair/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GGToEE_Elastic_2018','/GGToEE_Pt-35_Elastic_13TeV-lpair/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GGToEE_Inel-El_2018','/GGToEE_Pt-35_Inel-El_13TeV-lpair/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GGToEE_Inel-Inel_2018','/GGToEE_Pt-35_Inel-Inel_13TeV-lpair/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GluGluToWWToENEN_2018','/GluGluToWWToENEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GluGluToWWToENMN_2018','/GluGluToWWToENMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GluGluToWWToMNEN_2018','/GluGluToWWToMNEN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['GluGluToWWToMNMN_2018','/GluGluToWWToMNMN_TuneCP5_13TeV_MCFM701_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'],
   ['WJetsToLNu_2018','/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM']]

data_preVFP = [['SingleMuon_B_2','/SingleMuon/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD'],
               ['SingleMuon_C','/SingleMuon/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD'],
               ['SingleElectron_B_2','/SingleElectron/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD'],
               ['SingleElectron_C','/SingleElectron/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD']]

data_postVFP = [['SingleMuon_G','/SingleMuon/Run2016G-UL2016_MiniAODv2-v2/MINIAOD'],
                ['SingleElectron_G','/SingleElectron/Run2016G-UL2016_MiniAODv2-v2/MINIAOD']]

data_2017 = [['SingleMuon_2017B','/SingleMuon/Run2017B-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleMuon_2017C','/SingleMuon/Run2017C-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleMuon_2017D','/SingleMuon/Run2017D-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleMuon_2017E','/SingleMuon/Run2017E-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleMuon_2017F','/SingleMuon/Run2017F-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2017B','/SingleElectron/Run2017B-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2017C','/SingleElectron/Run2017C-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2017D','/SingleElectron/Run2017D-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2017E','/SingleElectron/Run2017E-UL2017_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2017F','/SingleElectron/Run2017F-UL2017_MiniAODv2-v1/MINIAOD']]

data_2018 = [['SingleMuon_2018A','/SingleMuon/Run2018A-UL2018_MiniAODv2-v3/MINIAOD'],
             ['SingleMuon_2018B','/SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD'],
             ['SingleMuon_2018C','/SingleMuon/Run2018C-UL2018_MiniAODv2-v2/MINIAOD'],
             ['SingleMuon_2018D','/SingleMuon/Run2018D-UL2018_MiniAODv2-v3/MINIAOD'],  
             ['SingleElectron_2018A','/EGamma/Run2018A-UL2018_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2018B','/EGamma/Run2018B-UL2018_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2018C','/EGamma/Run2018C-UL2018_MiniAODv2-v1/MINIAOD'],
             ['SingleElectron_2018D','/EGamma/Run2018D-UL2018_MiniAODv2-v2/MINIAOD']]

for tag in tags:
    config.JobType.pyCfgParams=tag
    if (tag[0]=="tag=MC_preVFP") or (tag[0]=="tag=MC_postVFP") or (tag[0]=="tag=MC_2017") or (tag[0]=="tag=MC_2018"):
        config.Data.splitting = 'FileBased'
        config.Data.outputDatasetTag = 'MC'
    else:
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 500
        config.Data.outputDatasetTag = 'data'
        datasets = [0,0]        
    if (tag[0]=="tag=MC_preVFP"):
	    datasets = MC_preVFP
    elif (tag[0]=="tag=MC_postVFP"):
        datasets = MC_postVFP
    elif (tag[0]=="tag=MC_2017"):
        datasets = MC_2017
    elif (tag[0]=="tag=MC_2018"):
        datasets = MC_2018
        
    if (tag[0]=="tag=data_preVFP"):
        config.Data.lumiMask = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Analysis/test/Ntuple_build/CMSSW_10_6_26/src/MakeNtuple/MakeNTuple/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns.txt'
        datasets = data_preVFP
    elif (tag[0]=="tag=data_postVFP"):
        config.Data.lumiMask = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Analysis/test/Ntuple_build/CMSSW_10_6_26/src/MakeNtuple/MakeNTuple/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns.txt' 
        datasets = data_postVFP
    elif (tag[0]=="tag=data_2017"):
        config.Data.lumiMask = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Analysis/test/Ntuple_build/CMSSW_10_6_26/src/MakeNtuple/MakeNTuple/combined_RPIN_CMS.json'
        datasets = data_2017
    elif (tag[0]=="tag=data_2018"):
        config.Data.lumiMask = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Analysis/test/Ntuple_build/CMSSW_10_6_26/src/MakeNtuple/MakeNTuple/CMSgolden_2RPGood_anyarms.json'
        datasets = data_2018

    for dataset in datasets:
        if dataset==0:
           break
        config.General.requestName = dataset[0]
        config.Data.inputDataset = dataset[1]
        p = Process(target=submit, args=(config,)
        p.start()
        p.join()
        submit(config)

