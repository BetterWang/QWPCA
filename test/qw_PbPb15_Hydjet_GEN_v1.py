import FWCore.ParameterSet.Config as cms

process = cms.Process("PCA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v13', '')

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
#        fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/q/qwang/work/cleanroomRun2/Ana/CMSSW_7_5_8_patch2/src/QWAna/QWCumuV3/test/HIMinBias_28.root")
	fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/q/qwang/work/cleanroomRun2/Ana/data/hydjet_pixel.root")
)

import HLTrigger.HLTfilters.hltHighLevel_cfi


process.QWPCA = cms.EDAnalyzer('QWPCA'
		, bGen = cms.untracked.bool(True)
		, bSim = cms.untracked.bool(False)
		, bEff = cms.untracked.bool(False)
		, minPt = cms.untracked.double(0.3)
		, maxPt = cms.untracked.double(3.0)
		, centrality = cms.InputTag("centralityBin", "HFtowers")
		, trackTag = cms.untracked.InputTag('genParticles')
		, vertexSrc = cms.untracked.InputTag('hiSelectedVertex', "")
#		, fweight = cms.untracked.InputTag('EffCorrectionsPixel_TT_pt_0_10_v2.root')
		, pterrorpt = cms.untracked.double(0.1)
		, dzdzerror = cms.untracked.double(3.0)
		, d0d0error = cms.untracked.double(3.0)
		, minvz = cms.untracked.double(-1.0)
		, maxvz = cms.untracked.double(15.0)
		, minEta = cms.untracked.double(-2.4)
		, maxEta = cms.untracked.double(2.4)
		, minCent = cms.untracked.int32(-1)
		, maxCent = cms.untracked.int32(500)
		)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('pca.root')
)

process.path= cms.Path(process.QWPCA)

process.schedule = cms.Schedule(
	process.path,
)
