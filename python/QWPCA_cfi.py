import FWCore.ParameterSet.Config as cms

QWPCA = cms.EDAnalyzer('QWPCA'
		, bGen = cms.untracked.bool(False)
		, bSim = cms.untracked.bool(False)
		, bEff = cms.untracked.bool(False)
		, minPt = cms.untracked.double(0.3)
		, maxPt = cms.untracked.double(3.0)
		, centrality = cms.InputTag("centralityBin", "HFtowers")
		, trackTag = cms.untracked.InputTag('hiGeneralAndPixelTracks')
		, algoParameters = cms.vint32(4,5,6,7)
		, vertexSrc = cms.untracked.InputTag('hiSelectedVertex', "")
#		, fweight_ = cms.untracked.InputTag('PbPb_dijet_TT_5TeV_v2.root')
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

