import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")
process = cms.Process("Test")

#import EventFilter.L1TXRawToDigi.util as util

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
#if len(options.inputFiles) is 0 and options.inputFileList is '' :
#    inputFiles = util.getFilesForRun(options.runNumber, options.dataStream)
#elif len(options.inputFileList) > 0 :
#    with open(options.inputFileList) as f :
#        inputFiles = list((line.strip() for line in f))
#else :
#    inputFiles = cms.untracked.vstring(options.inputFiles)
#if len(inputFiles) is 0 :
#    raise Exception('No files found for dataset %s run %d' % (options.dataStream, options.runNumber))
#print 'Ok, time to analyze'


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v10', '')

process.hcalDigis.saveQIE10DataNSamples = cms.untracked.vint32( 10)
process.hcalDigis.saveQIE10DataTags = cms.untracked.vstring( "MYDATA" )

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

#process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1BoostedJetStudies_cfi")

#process.l1NtupleProducer.isData = cms.bool(False)
process.l1NtupleProducer.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.l1NtupleProducer.hcalToken = cms.InputTag("l1tCaloLayer1Digis")
process.l1NtupleProducer.useECALLUT = cms.bool(False)
process.l1NtupleProducer.useHCALLUT = cms.bool(False)
process.l1NtupleProducer.useHFLUT   = cms.bool(False)
process.l1NtupleProducer.useLSB     = cms.bool(True)
process.l1NtupleProducer.verbose    = cms.bool(False)

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT   = cms.bool(False)
process.uct2016EmulatorDigis.useLSB     = cms.bool(True)
process.uct2016EmulatorDigis.verbose    = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken  = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken  = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:/hdfs/store/user/ojalvo/ZB-events/mini-1.root'
#                                '/store/data/Run2018E/ZeroBias/MINIAOD/PromptReco-v1/000/325/684/00000/7A23D558-A9D3-A546-95FA-12B56317C0A0.root',
#                                '/store/data/Run2018E/ZeroBias/MINIAOD/PromptReco-v1/000/325/484/00000/C08A230D-9ADC-F14E-AF94-3048FA9BA9F0.root',
#                                '/store/data/Run2018E/ZeroBias/MINIAOD/PromptReco-v1/000/325/449/00000/357A6203-9FE7-2F49-A7D3-04033E1F3DEC.root',
#                                '/store/data/Run2018E/ZeroBias/MINIAOD/PromptReco-v1/000/325/458/00000/E492B615-2F02-944E-A4D2-695290837592.root',
#                                '/store/data/Run2018E/ZeroBias/MINIAOD/PromptReco-v1/000/325/430/00000/76B231DA-491C-E244-B3C6-C89FEE3AD921.root',
#                                '/store/data/Run2018E/ZeroBias/MINIAOD/PromptReco-v1/000/325/283/00000/4A15B34D-ED00-784B-8AE9-CF03D646BF92.root'
                                                          ),
                            secondaryFileNames = cms.untracked.vstring(
                                'file:/hdfs/store/user/ojalvo/ZB-events/raw-1.root'
#                                '/store/data/Run2018E/ZeroBias/RAW/v1/000/325/684/00000/DAC64D8B-241A-2546-9100-65B44FFCCF83.root',
#                                '/store/data/Run2018E/ZeroBias/RAW/v1/000/325/484/00000/AF499180-B00E-6B40-9836-F0C5DC7D8497.root',
#                                '/store/data/Run2018E/ZeroBias/RAW/v1/000/325/449/00000/9D10EB45-D6E9-684E-8489-DDF7BC4A9A3E.root',
#                                '/store/data/Run2018E/ZeroBias/RAW/v1/000/325/458/00000/8CD5C424-E460-9541-94EC-679125C2ECE4.root',
#                                '/store/data/Run2018E/ZeroBias/RAW/v1/000/325/430/00000/7A52F28E-9F80-4845-92E2-1A71F82301F8.root',
#                                '/store/data/Run2018E/ZeroBias/RAW/v1/000/325/283/00000/2C9521F5-8784-FA47-963B-64BE55AB9CA7.root'
                                                                   )

                            )



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
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-Boosted-ZB.root")
)


process.L1TRawToDigi_Stage2 = cms.Task(process.caloLayer1Digis, process.caloStage2Digis)
process.RawToDigi_short = cms.Sequence(process.L1TRawToDigi_Stage2)
process.p = cms.Path(process.RawToDigi_short*process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)
process.e = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.p,process.e)
process.schedule = cms.Schedule(process.p)

#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
