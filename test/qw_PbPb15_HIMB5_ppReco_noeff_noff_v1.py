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


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('pca.root')
)

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.clusterCompatibilityFilter.clusterPars = cms.vdouble(0.0,0.006)

process.primaryVertexFilter.src = cms.InputTag("offlinePrimaryVertices")

process.eventSelection = cms.Sequence(
        process.hfCoincFilter3
        + process.primaryVertexFilter
#        + process.clusterCompatibilityFilter
)


process.QWPCA = cms.EDAnalyzer('QWPCA'
		, trackEta = cms.untracked.InputTag('QWEvent', "eta")
		, trackPhi = cms.untracked.InputTag('QWEvent', "phi")
		, trackWeight = cms.untracked.InputTag('QWEvent', "weight")
		, vertex = cms.untracked.InputTag('offlinePrimaryVertices', "")
		, centrality = cms.untracked.InputTag('Noff', "")
		, minvz = cms.untracked.double(-1.0)
		, maxvz = cms.untracked.double(15.0)
		, nvtx = cms.untracked.int32(100)
		)

process.load('PbPb_HIMB5_ppReco_noeff')

#process.path= cms.Path(process.eventSelection*process.centralityBin*process.QWEvent*process.QWPCA*process.histNoff*process.vectPhi*process.vectEta*process.vectPt*process.vectPhiW*process.vectEtaW*process.vectPtW)
process.path= cms.Path(process.eventSelection*process.makeEvent*process.QWPCA*process.vectMonW)

process.schedule = cms.Schedule(
	process.path,
)
