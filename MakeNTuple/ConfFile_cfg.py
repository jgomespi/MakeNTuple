import FWCore.ParameterSet.Config as cms
import copy
import os
import sys
# Define collections
muon             = 'slimmedMuons'
electron        = 'slimmedElectrons'

import FWCore.ParameterSet.VarParsing as VarParsing
# setup 'standard'  options
options = VarParsing.VarParsing ('analysis')
# setup any defaults you want
print("Get the options")
options.tag='MC_signal_preVFP' 
options.maxEvents=100
options.parseArguments()
print("Done")
from Configuration.StandardSequences.Eras import eras
print("Loading the configurations of the given period")
if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP'):
	from Configuration.Eras.Era_Run2_2016_cff import *
	from Validation.CTPPS.simu_config.base_cff import *
	import direct_simu_reco2016_cff as profile
	process = cms.Process('CTPPSTestAcceptance', profile.era)
	profile.LoadConfig(process)
	process.load("CalibPPS.ESProducers.ctppsBeamParametersESSource_cfi")
	process.load("Validation.CTPPS.simu_config.year_2016_cff")
	process.load('direct_simu_reco_cff')
elif (options.tag=='MC_2017') or (options.tag=='MC_signal_2017'):
        #process = cms.Process("CTPPSTestProtonReconstruction", eras.Run2_2017)
	from Configuration.Eras.Era_Run2_2017_cff import *
        from Validation.CTPPS.simu_config.base_cff import *
        import direct_simu_reco2017_cff as profile
        process = cms.Process('CTPPSTestAcceptance', profile.era)
        profile.LoadConfig(process)
        process.load("CalibPPS.ESProducers.ctppsBeamParametersESSource_cfi")
        process.load("Validation.CTPPS.simu_config.year_2017_cff")
        process.load('direct_simu_reco_cff')
elif (options.tag=='MC_2018') or (options.tag=='MC_signal_2018'):
        from Configuration.Eras.Era_Run2_2018_cff import *
        from Validation.CTPPS.simu_config.base_cff import *
        import direct_simu_reco2018_cff as profile
        process = cms.Process('CTPPSTestAcceptance', profile.era)
        profile.LoadConfig(process)
        process.load("CalibPPS.ESProducers.ctppsBeamParametersESSource_cfi")
        process.load("Validation.CTPPS.simu_config.year_2018_cff")
        process.load('direct_simu_reco_cff')
elif (options.tag == "data_preVFP") or (options.tag == "data_postVFP"):
	process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2016)
elif (options.tag == "data_2017"):
        process = cms.Process("CTPPSTestProtonReconstruction", eras.Run2_2017)
elif (options.tag == "data_2018"):
        process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2018)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load('Configuration.StandardSequences.GeometryDB_cff')
# Added in order to get the proton collections (https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsDirectSimulation#Full_step_by_step_recipe_to_setu)
process.load('Configuration.EventContent.EventContent_cff')
print("Done")
print("Processing the maxEvents")
process.maxEvents = cms.untracked.PSet( 
			input = cms.untracked.int32(options.maxEvents) 
)
print("Done")
print("Processing untracked")
process.options=cms.untracked.PSet(
	wantSummary=cms.untracked.bool(True)
	, SkipEvent = cms.untracked.vstring('ProductNotFound')
#	, numberOfThreads = cms.untracked.uint32( 8 )
)
print("Done")
print("Associating the bool vars and the file")
filename = " "
if (options.tag=='MC_signal_2018'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(True)
    DATA=cms.bool(False)
    preVFP=cms.bool(False)
    postVFP=cms.bool(False)
    filename  = "file:signal_2018.root"
if (options.tag=='MC_signal_2017'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(True)
    DATA=cms.bool(False)
    preVFP=cms.bool(False)
    postVFP=cms.bool(False)
    filename  = "file:signal_2017.root"
if (options.tag=='MC_preVFP'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(False)
    preVFP=cms.bool(True)
    postVFP=cms.bool(False)
    filename = "root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL16MiniAODAPVv2/GGToEE_Pt-35_Elastic_13TeV-lpair/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/0A42F294-B631-8A44-9B3D-0D465620AAC5.root"
    #filename = "file:MiniAOD_DY_preVFP.root"
if (options.tag=='MC_postVFP'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(False)
    preVFP=cms.bool(False)
    postVFP=cms.bool(True)
    filename = "root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL16MiniAODv2/WW_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/120000/1E1FCD8B-D656-5642-BFE3-E000EE0D4510.root"
if (options.tag=='MC_2017'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(False)
    preVFP=cms.bool(False)
    postVFP=cms.bool(False)
    #filename =  "root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL16MiniAODv2/WW_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/120000/1E1FCD8B-D656-5642-BFE3-E000EE0D4510.root"
    filename = "root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL17MiniAODv2/GGToMuMu_Pt-25_Elastic_13TeV-lpair/MINIAODSIM/106X_mc2017_realistic_v9-v2/110000/6CF52578-30AF-0140-AF5D-C610145BCEF1.root"
if (options.tag=='MC_2018'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(False)
    preVFP=cms.bool(False)
    postVFP=cms.bool(False)
    #filename =  "root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18MiniAODv2/WW_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/22B9609C-A972-F24B-AB8F-E12BC5000CD1.root"
    filename =  "root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18MiniAODv2/GGToMuMu_Pt-25_Elastic_13TeV-lpair/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/5BCB678B-7C30-4746-BFE2-0B7B167BB56D.root"
if (options.tag=='MC_signal_preVFP'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(True)
    DATA=cms.bool(False)
    preVFP=cms.bool(True)
    postVFP=cms.bool(False)
    filename="file:signal_preVFP.root"
    filename = "file:/eos/home-j/jgomespi/workspace/FullSim_MiniAOD/MuMu_preVFP/merge/MuMu_preVFP_2e-06_0.0_6.root"
    #filename="file:/eos/home-j/jgomespi/workspace/FullSim_MiniAOD/MuMu_preVFP/MiniAOD_MuMu_preVFP_0.0_0.0_833.root"
if (options.tag=='MC_signal_postVFP'):
    MC=cms.bool(True)
    MC_Signal=cms.bool(True)
    DATA=cms.bool(False)
    preVFP=cms.bool(False)
    postVFP=cms.bool(True)
    filename = "file:signal.root"
if (options.tag=='data_preVFP'):
    filename = "root://cmsxrootd.fnal.gov//store/data/Run2016B/SingleMuon/MINIAOD/ver2_HIPM_UL2016_MiniAODv2-v2/120000/0042DCA3-FD73-4641-B984-636AA05DFB55.root"
    MC=cms.bool(False)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(True)
    preVFP=cms.bool(True)
    postVFP=cms.bool(False)
if (options.tag=='data_postVFP'):
    filename = "root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/MINIAOD/UL2016_MiniAODv2-v2/120000/001FDE5F-A989-2F48-A280-D4D0F7766D95.root"
    MC=cms.bool(False)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(True)
    preVFP=cms.bool(False)
    postVFP=cms.bool(True)
if (options.tag=='data_2017'):
    filename = "root://cmsxrootd.fnal.gov//store/data/Run2017B/SingleMuon/MINIAOD/UL2017_MiniAODv2-v1/260000/9032A966-8ED0-B645-97B6-A8EBC1D8D3B9.root"
    MC=cms.bool(False)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(True)
    preVFP=cms.bool(False)
    postVFP=cms.bool(False)
if (options.tag=='data_2018'):
    filename = "root://cmsxrootd.fnal.gov//store/data/Run2018A/SingleMuon/MINIAOD/UL2018_MiniAODv2-v3/2530000/002A113D-FB15-1341-A170-638E53A7261F.root"
    MC=cms.bool(False)
    MC_Signal=cms.bool(False)
    DATA=cms.bool(True)
    preVFP=cms.bool(False)
    postVFP=cms.bool(False)
print("Done")
print("Duplicate check")
process.source = cms.Source("PoolSource"
    , fileNames = cms.untracked.vstring(filename)
    , duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)
print("Done")
print("GlobalTag config")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#Joao Pedro changed:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun2LegacyAnalysis
if (options.tag=='MC_preVFP') or (options.tag=='MC_signal_preVFP'):
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_preVFP_v11','') # MC: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2016preVFPAnalysis
elif (options.tag=='MC_postVFP') or (options.tag=='MC_signal_postVFP'):
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '') # MC: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2016postVFPAnalysis
elif (options.tag=='MC_2017') or (options.tag=='MC_signal_2017'):
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v10', '') # MC: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis
elif (options.tag=='MC_2018') or (options.tag=='MC_signal_2018'):
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '') # MC: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis
else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37') # data
print("Done")
print("Proton configuration")
# local RP reconstruction chain with standard settings

process.load("RecoCTPPS.Configuration.recoCTPPS_cff")


##################      P R O T O N    ####################
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  beamDivergenceVtxGenerator = cms.PSet(initialSeed =cms.untracked.uint32(849))
)

if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_2017') or (options.tag=='MC_2018') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP') or (options.tag=='MC_signal_2017') or (options.tag=='MC_signal_2018'):
	# override LHCInfo source         
	if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP'):
		UseCrossingAngle(185.,process)
	
	if (options.tag=='MC_2017') or (options.tag=='MC_signal_2017'):
		process.load("CalibPPS.ESProducers.ctppsLHCInfoRandomXangleESSource_cfi")
        	process.ctppsLHCInfoRandomXangleESSource.generateEveryNEvents = 1
                process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "xangle.root"
		process.ctppsLHCInfoRandomXangleESSource.xangleHistogramObject = "2017_preTS2"
		process.ctppsLHCInfoRandomXangleESSource.beamEnergy = 6500.
		process.ctppsLHCInfoRandomXangleESSource.betaStar = 0.40
	
	if (options.tag=='MC_2018') or (options.tag=='MC_signal_2018'):
                process.load("CalibPPS.ESProducers.ctppsLHCInfoRandomXangleESSource_cfi")
                process.ctppsLHCInfoRandomXangleESSource.generateEveryNEvents = 1
		process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "xangle.root"
                process.ctppsLHCInfoRandomXangleESSource.xangleHistogramObject = "2018_preTS1"
		process.ctppsLHCInfoRandomXangleESSource.beamEnergy = 6500.
		process.ctppsLHCInfoRandomXangleESSource.betaStar = 0.30
	
	if (options.tag!='MC_preVFP') and (options.tag!='MC_postVFP') and (options.tag!='MC_signal_preVFP') and (options.tag!='MC_signal_postVFP'):
		process.esPreferLHCInfo = cms.ESPrefer("CTPPSLHCInfoRandomXangleESSource", "ctppsLHCInfoRandomXangleESSource")

		# override beam-parameter source
		process.load("CalibPPS.ESProducers.ctppsBeamParametersFromLHCInfoESSource_cfi")

		process.ctppsBeamParametersFromLHCInfoESSource.beamDivX45 = process.ctppsBeamParametersESSource.beamDivX45
		process.ctppsBeamParametersFromLHCInfoESSource.beamDivX56 = process.ctppsBeamParametersESSource.beamDivX56
		process.ctppsBeamParametersFromLHCInfoESSource.beamDivY45 = process.ctppsBeamParametersESSource.beamDivY45
		process.ctppsBeamParametersFromLHCInfoESSource.beamDivY56 = process.ctppsBeamParametersESSource.beamDivY56

		process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX45 = process.ctppsBeamParametersESSource.vtxOffsetX45
		process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX56 = process.ctppsBeamParametersESSource.vtxOffsetX56
		process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY45 = process.ctppsBeamParametersESSource.vtxOffsetY45
		process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY56 = process.ctppsBeamParametersESSource.vtxOffsetY56
		process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ45 = process.ctppsBeamParametersESSource.vtxOffsetZ45
		process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ56 = process.ctppsBeamParametersESSource.vtxOffsetZ56

		process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevX = process.ctppsBeamParametersESSource.vtxStddevX
		process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevY = process.ctppsBeamParametersESSource.vtxStddevY
		process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevZ = process.ctppsBeamParametersESSource.vtxStddevZ
	
	# update settings of beam-smearing module
	process.beamDivergenceVtxGenerator.src = cms.InputTag("")
	process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
  		cms.InputTag("prunedGenParticles")	
		#,cms.InputTag("genPUProtons","genPUProtons")
	)

	# do not apply vertex smearing again                                                                                                                   
	process.ctppsBeamParametersESSource.vtxStddevX = 0
	process.ctppsBeamParametersESSource.vtxStddevY = 0
	process.ctppsBeamParametersESSource.vtxStddevZ = 0

	# undo CMS vertex shift                                                                                                                                
	if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP'):
                process.ctppsBeamParametersESSource.vtxOffsetX45 = -0.09163
                process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.16955
                process.ctppsBeamParametersESSource.vtxOffsetZ45 = -0.9315
		#process.ctppsBeamParametersESSource.vtxOffsetX45 = -0.1077 
		#process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.418
		#process.ctppsBeamParametersESSource.vtxOffsetZ45 = +1.576

        if (options.tag=='MC_2017') or (options.tag=='MC_signal_2017'):
                process.ctppsBeamParametersESSource.vtxOffsetX45 = 0.024755
                process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.069233
                process.ctppsBeamParametersESSource.vtxOffsetZ45 = -0.82054

        if (options.tag=='MC_2018') or (options.tag=='MC_signal_2018'):
                process.ctppsBeamParametersESSource.vtxOffsetX45 = -0.0107682
                process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.041722
                process.ctppsBeamParametersESSource.vtxOffsetZ45 = -0.035748

print("Done")


print("Processing Trigger")
###################    T R I G G E R    ###################
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = copy.deepcopy(hltHighLevel)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
if (options.tag=='MC_2018') or (options.tag=='MC_signal_2018') or (options.tag=='data_2018'):
	process.hltFilter.HLTPaths = ['HLT_Ele28_WPTight_Gsf_v*','HLT_IsoMu24_v*']
else:
	process.hltFilter.HLTPaths = ['HLT_Ele27_WPTight_Gsf_v*','HLT_IsoMu24_v*']

process.hltFilter.throw = cms.bool(False)
process.hltFilter.andOr = cms.bool(True) # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
print("Done")
print("MET and Filter")

## The three lines below should always be included
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

if DATA:
	DATA_=True
else:
	DATA_=False
runMetCorAndUncFromMiniAOD(process,
                           isData=DATA_
                           )
process.load("RecoMET.METFilters.metFilters_cff")

#process.METFilter = process.metFilters

'''
if (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP') or (options.tag=='MC_signal_2017') or (options.tag=='MC_signal_2018'):
	process.metFilters.TriggerResultsTag =  cms.InputTag("TriggerResults","","PAT") # cms.InputTag("TriggerResults","","RECO")
else:
	process.metFilters.TriggerResultsTag =  cms.InputTag("TriggerResults","","RECO")

#process.metFilters.throw = cms.bool(False) # throw exception on unknown path names
#process.metFilters.andOr = cms.bool(False) # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
# Filters from here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
# MET filter recommendations will be updated for UltraLegacy datasets. The filters are currently under review and it is likely that the recommendations for ECAL related filters will change. Updates recommendations to be released as soon as they are available.

if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP'):
	process.metFilters.HLTPaths = ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_BadPFMuonDzFilter','Flag_eeBadScFilter'] #MC
elif (options.tag=='data_preVFP') or (options.tag=='data_postVFP'):
        process.metFilters.HLTPaths = ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_BadPFMuonDzFilter','Flag_eeBadScFilter'] #data
elif (options.tag=='MC_2017') or (options.tag=='MC_2018') or (options.tag=='MC_signal_2017') or (options.tag=='MC_signal_2018'):
        process.metFilters.HLTPaths = ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_BadPFMuonDzFilter','Flag_eeBadScFilter','Flag_ecalBadCalibFilter'] #MC
elif (options.tag=='data_2017') or (options.tag=='data_2018'):
        process.metFilters.HLTPaths = ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_BadPFMuonDzFilter','Flag_eeBadScFilter','Flag_ecalBadCalibFilter'] #data

from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
)
'''
print("Done")
###########################################################

# COMMENTED BY JOAO PEDRO 18/02/21
'''
###################    C L E A N E R    ###################
process.cleanPatMuons = cms.EDProducer("PATMuonCleaner",
    src = cms.InputTag(muon),
    # preselection (any string-based cut for pat::Muon)
    preselection = cms.string("pt() > 20. && abs(eta()) < 2.4"),
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Muon)
    finalCut = cms.string(''),
)

process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    ## pat electron input source
    src = cms.InputTag(electron),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string("pt() > 20. && abs(eta()) < 2.4"), #electronID('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)
'''
print("Prefiring")

#################### P R E F I R I N G ####################
from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
if (options.tag=="MC_preVFP") or (options.tag=="MC_signal_preVFP") or (options.tag=="data_preVFP"):
	pref_string_muons = "2016preVFP"
	pref_string_electrons = "UL2016preVFP"
if (options.tag=="MC_postVFP") or (options.tag=="MC_signal_postVFP") or (options.tag=="data_postVFP"):
        pref_string_muons = "2016postVFP"
        pref_string_electrons = "UL2016postVFP"
if (options.tag=="MC_2017") or (options.tag=="MC_signal_2017") or (options.tag=="data_2017"):
        pref_string_muons = "20172018"
        pref_string_electrons = "UL2017BtoF"
if (options.tag=="MC_2018") or (options.tag=="MC_signal_2018") or (options.tag=="data_2018"):
        pref_string_muons = "20172018"
        pref_string_electrons = "None"
if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_2017') or (options.tag=='MC_2018') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP') or (options.tag=='MC_signal_2017') or (options.tag=='MC_signal_2018'):
	process.prefiringweight = l1PrefiringWeightProducer.clone(
    		TheJets = cms.InputTag("slimmedJets") #this should be the slimmedJets collection with up to date JECs !
    		, DataEraECAL = cms.string(pref_string_electrons) 
    		, DataEraMuon = cms.string(pref_string_muons)
    		, UseJetEMPt = cms.bool(False)
    		, PrefiringRateSystematicUnctyMuon = cms.double(0.2)
    		, PrefiringRateSystematicUnctyECAL = cms.double(0.2)
    		#, SkipWarnings = False
	)
print("Done")
print("EGamma post reco")
############# E / G A M M A   P O S T   R E C O ###########
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if (options.tag=="MC_preVFP") or (options.tag=="MC_signal_preVFP") or (options.tag=="data_preVFP"): 
	era_ = '2016preVFP-UL'
elif (options.tag=="MC_postVFP") or (options.tag=="MC_signal_postVFP") or (options.tag=="data_postVFP"): 
	era_ = '2016postVFP-UL'
elif (options.tag=="MC_2017") or (options.tag=="MC_signal_2017") or (options.tag=="data_2017"):
        era_ = '2017-UL'
elif (options.tag=="MC_2018") or (options.tag=="MC_signal_2018") or (options.tag=="data_2018"):
        era_ = '2018-UL'

setupEgammaPostRecoSeq(
		process
		,era=era_ # change era
) 
print("Done")
# MET
METCollection="slimmedMETsPuppi"
print("Filter")

# Create a EDFilter to choose for 2 oposite charge leptons with |eta|<2.4 and pt>38GeV and loose ID, primary vertex at |z| < 15cm.
process.demo2 = cms.EDFilter('MyUserSelector'
        , muons                    = cms.InputTag(muon)
        , electrons                = cms.InputTag(electron)
        , vertices                 = cms.InputTag('offlineSlimmedPrimaryVertices')
	, Z_vtx			   = cms.double(15.)
	, pT_min		   = cms.double(38.)
	, eta_max		   = cms.double(2.4)	
        , filter = cms.bool(True)
)
tag_ = str(options.tag)
# demo process, which processes MakeNTuple:
print("Done")
print("MakeNTuple")
process.demo = cms.EDAnalyzer('MakeNTuple'
        , muons                    = cms.InputTag(muon)
        , electrons                = cms.InputTag(electron)
        , MET                      = cms.InputTag('slimmedMETs')#cms.InputTag("slimmedMETs","","PAT")
        , PFCand                   = cms.InputTag('packedPFCandidates')
        , PileupSumInfoInputTag = cms.InputTag('slimmedAddPileupInfo')
        , vertices                 = cms.InputTag('offlineSlimmedPrimaryVertices')
	, MC			   = MC
        , MC_Signal                = MC_Signal
        , DATA                     = DATA
        , preVFP_                  = preVFP
        , postVFP_                 = postVFP
	, tag_			   = cms.string(tag_)
        , ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP")
        , ppsRecoProtonMultiRPTag  = cms.InputTag("ctppsProtons", "multiRP")
	, TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
        , genTag  = cms.InputTag('prunedGenParticles')
	, MCEvent = cms.untracked.InputTag('source','generator')
)
print("Done")
print("Output")
# Output
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string("out.root")
)
print("Done")
#process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(False)

# After the modification by JOAO PEDRO 18/02/21:
if (options.tag=='MC_preVFP') or (options.tag=='MC_postVFP') or (options.tag=='MC_2017') or (options.tag=='MC_2018') or (options.tag=='MC_signal_preVFP') or (options.tag=='MC_signal_postVFP') or (options.tag=='MC_signal_2017') or (options.tag=='MC_signal_2018'):
        process.p = cms.Path(
		process.beamDivergenceVtxGenerator
                *process.ctppsDirectProtonSimulation
                *process.reco_local
                *process.ctppsProtons
                *process.hltFilter
                *process.fullPatMetSequence 
#                *process.metFilters
#                *process.cleanPatElectrons
#                *process.cleanPatMuons
                *process.egammaPostRecoSeq
#		*process.BadPFMuonFilterUpdateDz
		*process.prefiringweight
                *process.demo2
                *process.demo
                )
else:
	process.p = cms.Path(
		process.hltFilter
                *process.fullPatMetSequence 
#		*process.metFilters
#		*process.cleanPatElectrons
#		*process.cleanPatMuons
		*process.egammaPostRecoSeq
#		*process.BadPFMuonFilterUpdateDz
		*process.demo2
		*process.demo
		)

