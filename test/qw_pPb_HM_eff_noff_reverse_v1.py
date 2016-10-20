import FWCore.ParameterSet.Config as cms

process = cms.Process("PCA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.GlobalTag.globaltag = 'GR_P_V43::All'


process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/q/qwang/work/cleanroom2/CMSSW_5_3_20/pPb_HM_1000_1_Bgt.root")
)


import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltHM100 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM100.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
#       "HLT_PAPixelTracks_Multiplicity130_v*",
#       "HLT_PAPixelTracks_Multiplicity160_v*",
#       "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM100.andOr = cms.bool(True)
process.hltHM100.throw = cms.bool(False)

process.hltHM130 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM130.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
#       "HLT_PAPixelTracks_Multiplicity160_v*",
##      "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM130.andOr = cms.bool(True)
process.hltHM130.throw = cms.bool(False)


process.hltHM160 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM160.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
        "HLT_PAPixelTracks_Multiplicity160_v*",
#       "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM160.andOr = cms.bool(True)
process.hltHM160.throw = cms.bool(False)




process.hltHM190 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM190.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
        "HLT_PAPixelTracks_Multiplicity160_v*",
        "HLT_PAPixelTracks_Multiplicity190_v*",
#       "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM190.andOr = cms.bool(True)
process.hltHM190.throw = cms.bool(False)


process.hltHM220 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM220.HLTPaths = [
        "HLT_PAPixelTracks_Multiplicity100_v*",
        "HLT_PAPixelTracks_Multiplicity130_v*",
        "HLT_PAPixelTracks_Multiplicity160_v*",
        "HLT_PAPixelTracks_Multiplicity190_v*",
        "HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM220.andOr = cms.bool(True)
process.hltHM220.throw = cms.bool(False)





process.TFileService = cms.Service("TFileService",
    fileName = cms.string('pca.root')
)


process.QWPCA = cms.EDAnalyzer('QWPCA'
		, trackEta = cms.untracked.InputTag('QWEvent', "eta")
		, trackPhi = cms.untracked.InputTag('QWEvent', "phi")
		, trackWeight = cms.untracked.InputTag('QWEvent', "weight")
		, vertexZ = cms.untracked.InputTag('QWEvent', "vz")
		, centrality = cms.untracked.InputTag('Noff', "")
		, minvz = cms.untracked.double(-1.0)
		, maxvz = cms.untracked.double(15.0)
		, nvtx = cms.untracked.int32(100)
		)

process.load('pPb_HM_eff')

process.QWEvent.bFlip = cms.untracked.bool(True)
process.NoffFilter100 = cms.EDFilter("QWIntFilter",
	selectedBins = cms.vint32(120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149),
	BinLabel = cms.InputTag("Noff")
	)

process.NoffFilter130 = cms.EDFilter("QWIntFilter",
	selectedBins = cms.vint32(
		150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184),
	BinLabel = cms.InputTag("Noff")
	)

process.NoffFilter160 = cms.EDFilter("QWIntFilter",
	selectedBins = cms.vint32(
		185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219),
	BinLabel = cms.InputTag("Noff")
	)

process.NoffFilter190 = cms.EDFilter("QWIntFilter",
	selectedBins = cms.vint32(
		220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259),
	BinLabel = cms.InputTag("Noff")
	)

process.NoffFilter220 = cms.EDFilter("QWIntFilter",
	selectedBins = cms.vint32(
		260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,
		311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349),
	BinLabel = cms.InputTag("Noff")
	)
#process.path= cms.Path(process.eventSelection*process.centralityBin*process.QWEvent*process.QWPCA*process.histNoff*process.vectPhi*process.vectEta*process.vectPt*process.vectPhiW*process.vectEtaW*process.vectPtW)

process.path100 = cms.Path(process.hltHM100*process.makeEvent*process.NoffFilter100*process.QWPCA*process.vectMonW)
process.path130 = cms.Path(process.hltHM130*process.makeEvent*process.NoffFilter130*process.QWPCA*process.vectMonW)
process.path160 = cms.Path(process.hltHM160*process.makeEvent*process.NoffFilter160*process.QWPCA*process.vectMonW)
process.path190 = cms.Path(process.hltHM190*process.makeEvent*process.NoffFilter190*process.QWPCA*process.vectMonW)
process.path220 = cms.Path(process.hltHM220*process.makeEvent*process.NoffFilter220*process.QWPCA*process.vectMonW)

process.schedule = cms.Schedule(
	process.path100,
	process.path130,
	process.path160,
	process.path190,
	process.path220,
)
