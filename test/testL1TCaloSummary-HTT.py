import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.parseArguments()

def formatLumis(lumistring, run) :
    lumis = (lrange.split('-') for lrange in lumistring.split(','))
    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
    return ['-'.join(l) for l in runlumis]

print 'Getting files for run %d...' % options.runNumber

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('L1Trigger.Configuration.SimL1Emulator_cff')

process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')
process.load('L1Trigger.L1TCaloLayer1.simCaloStage2Layer1Digis_cfi')

process.load('L1Trigger.L1TCaloLayer1.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1BoostedTauStudies_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
				#'file:/eos/user/p/pdas/L1Boosted/ggHbb/MiniAOD/RunIIAutumn18MiniAOD_21Dec_0_5300.root'
                                'file:/eos/user/p/pdas/L1Boosted/ggHtautau/MiniAOD/RunIIAutumn18MiniAOD_21Dec_0_5300.root'
),
                            secondaryFileNames = cms.untracked.vstring(
                                #file:/eos/user/p/pdas/L1Boosted/ggHtautau/DR/RunIIAutumn18DRPremix_21Dec_0_5300.root
                                'file:/eos/user/p/pdas/L1Boosted/ggHtautau/DR/RunIIAutumn18DRPremix_step1_21Dec_0_5300.root'
				#'file:/eos/user/p/pdas/L1Boosted/ggHbb/DR/RunIIAutumn18DRPremix_step1_21Dec_0_5300.root'
                           )
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:734","1:961","1:966","1:982")
#process.source.eventsToProcess = cms.untracked.VEventRange("1:960110","1:965580")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent.root"),
    outputCommands = cms.untracked.vstring('keep *')
    #outputCommands = cms.untracked.vstring('drop *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-Htt.root")
)

process.p = cms.Path(process.RawToDigi*process.simCaloStage2Layer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

process.e = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.p,process.e)
process.schedule = cms.Schedule(process.p)

from L1Trigger.Configuration.customiseReEmul import L1TReEmulMCFromRAWSimHcalTP
process = L1TReEmulMCFromRAWSimHcalTP(process)
from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2018_v1_3
process = L1TSettingsToCaloParams_2018_v1_3(process)

process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")

process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHBQIE11 = 1
process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHEQIE11 = 1
process.simHcalTriggerPrimitiveDigis.weightsQIE11 = {
    "ieta1" :  [-0.47, 1.0],
    "ieta2" :  [-0.47, 1.0],
    "ieta3" :  [-0.47, 1.0],
    "ieta4" :  [-0.47, 1.0],
    "ieta5" :  [-0.47, 1.0],
    "ieta6" :  [-0.47, 1.0],
    "ieta7" :  [-0.47, 1.0],
    "ieta8" :  [-0.47, 1.0],
    "ieta9" :  [-0.47, 1.0],
    "ieta10" : [-0.47, 1.0],
    "ieta11" : [-0.47, 1.0],
    "ieta12" : [-0.47, 1.0],
    "ieta13" : [-0.47, 1.0],
    "ieta14" : [-0.47, 1.0],
    "ieta15" : [-0.47, 1.0],
    "ieta16" : [-0.47, 1.0],
    "ieta17" : [-0.47, 1.0],
    "ieta18" : [-0.47, 1.0],
    "ieta19" : [-0.47, 1.0],
    "ieta20" : [-0.47, 1.0],
    "ieta21" : [-0.43, 1.0],
    "ieta22" : [-0.43, 1.0],
    "ieta23" : [-0.43, 1.0],
    "ieta24" : [-0.43, 1.0],
    "ieta25" : [-0.43, 1.0],
    "ieta26" : [-0.43, 1.0],
    "ieta27" : [-0.43, 1.0],
    "ieta28" : [-0.43, 1.0]
}

process.HcalTPGCoderULUT.contain1TSHB = True
process.HcalTPGCoderULUT.contain1TSHE = True

# Pick one of the pairs of lines below based on the intended scenario for running
process.HcalTPGCoderULUT.containPhaseNSHB = 1.0 # For Run3 MC
process.HcalTPGCoderULUT.containPhaseNSHE = 1.0 # For Run3 MC

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# Multi-threading
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
