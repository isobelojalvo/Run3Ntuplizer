import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")

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
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1TRegionNtupleProducer_cfi")

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring(inputFiles)#,
                            #secondaryFileNames = cms.untracked.vstring(secondaryMap[options.inputFiles[0]])
                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/MINIAODSIM/106X_upgrade2018_realistic_v4-v2/10000/E6E75F38-FAFC-4242-9C38-52C2A3D451AC.root') ,
                            secondaryFileNames = cms.untracked.vstring(
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/96B24904-342B-394E-B84B-341429E877F1.root',
			    'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/BBA135B4-8449-114C-B6DC-1447340C8FF2.root',
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/E6674DCD-AF07-E24F-A254-4EFEC2DF40E8.root'
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/FB2DCB4B-A9BA-1947-B661-FA03EF8D4ED9.root',
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/BF573032-3292-F045-8F39-5EC7FBEE28BA.root',
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/C587B60B-857E-7C48-90F1-23EC77B4BBF6.root',
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/05BA435A-401F-084B-AA00-4B695B8BAA19.root',
                            'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/3159EC8E-F6A9-BF4B-A6BD-D0B8846A863F.root',
			    'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/CBC10CE6-691E-5C4F-A4BF-3D2AF0F59844.root',
			    'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_6_0/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/106X_upgrade2018_realistic_v4-v2/10000/7A513272-70AD-B149-B9B0-6469897C2417.root'
                                                                   ),
inputCommands = 
cms.untracked.vstring("keep *", "drop patHcalDepthEnergyFractionsedmValueMap_packedPFCandidates_hcalDepthEnergyFractions_RECO")
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("321149:1055","320757:185")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

#outputFile = '/data/' + os.environ['USER'] + '/l1tCaloSummary-' + str(options.runNumber) + '.root'

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent-SPG.root"),
    outputCommands = cms.untracked.vstring('drop *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-SPG.root")
)

process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

#process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.p)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
