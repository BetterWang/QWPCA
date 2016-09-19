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
        fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/q/qwang/work/cleanroomRun2/Ana/data/ppReco.root")
)

import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltMB = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMB.HLTPaths = [
	"HLT_HIL1MinimumBiasHF2AND_*",
	"HLT_HIL1MinimumBiasHF1AND_*",
]
process.hltMB.andOr = cms.bool(True)
process.hltMB.throw = cms.bool(False)

process.QWPCA = cms.EDAnalyzer('QWPCA'
		, bGen = cms.untracked.bool(False)
		, bSim = cms.untracked.bool(False)
		, bEff = cms.untracked.bool(True)
		, minPt = cms.untracked.double(0.3)
		, maxPt = cms.untracked.double(3.0)
		, centrality = cms.InputTag("centralityBin", "")
		, trackTag = cms.untracked.InputTag('generalTracks')
		, vertexSrc = cms.untracked.InputTag('offlinePrimaryVertices', "")
		, fweight = cms.untracked.InputTag('Hydjet_eff_mult_v1.root')
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

#process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
#process.centralityBin.Centrality = cms.InputTag("hiCentrality")
#process.centralityBin.centralityVariable = cms.string("HFtowers")

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.clusterCompatibilityFilter.clusterPars = cms.vdouble(0.0,0.006)

process.primaryVertexFilter.src = cms.InputTag("offlinePrimaryVertices")

process.eventSelection = cms.Sequence(
        process.hfCoincFilter3
        + process.primaryVertexFilter
#        + process.clusterCompatibilityFilter
)

process.centralityBin = cms.EDProducer('QWPPRecoCentBinProducer')

process.path= cms.Path(process.eventSelection*process.centralityBin*process.QWPCA)

process.schedule = cms.Schedule(
	process.path,
)
