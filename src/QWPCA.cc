#include <algorithm>
#include <math.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "QWAna/QWPCA/interface/QWPCA.h"

using namespace std;

QWPCA::QWPCA(const edm::ParameterSet& iConfig):
	trackEta_( iConfig.getUntrackedParameter<edm::InputTag>("trackEta") ),
	trackPt_( iConfig.getUntrackedParameter<edm::InputTag>("trackPt") ),
	trackPhi_( iConfig.getUntrackedParameter<edm::InputTag>("trackPhi") ),
	trackWeight_( iConfig.getUntrackedParameter<edm::InputTag>("trackWeight") ),
	vertexTag_( iConfig.getUntrackedParameter<edm::InputTag>("vertex") ),
	centralityTag_( iConfig.getUntrackedParameter<edm::InputTag>("centrality") )
{
	minvz_ = iConfig.getUntrackedParameter<double>("minvz", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz", 15.);

	mineta_ = iConfig.getUntrackedParameter<double>("mineta", -2.4);
	maxeta_ = iConfig.getUntrackedParameter<double>("maxeta", 2.4);

	minpt_ = iConfig.getUntrackedParameter<double>("minpt", 0.3);
	maxpt_ = iConfig.getUntrackedParameter<double>("maxpt", 0.3);

	bCent_ = iConfig.getUntrackedParameter<bool>("bCent", false);

	cmode_ = iConfig.getUntrackedParameter<int>("cmode", 1);
	nvtx_ = iConfig.getUntrackedParameter<int>("nvtx", 100);

	consumes<reco::VertexCollection>(vertexTag_);
	consumes<std::vector<double> >(trackEta_);
	consumes<std::vector<double> >(trackPt_);
	consumes<std::vector<double> >(trackPhi_);
	consumes<std::vector<double> >(trackWeight_);

	edm::Service<TFileService> fs;
	trV = fs->make<TTree>("trV", "trV");

	trV->Branch("cent", &(Cent), "cent/I");
	trV->Branch("mult", &(Mult), "mult/I");
//	trV->Branch("noff", &(t.Noff), "mult/I");

	trV->Branch("pQeta2",  "std::vector< std::complex<double> >", &pQeta2);
	trV->Branch("pQeta3",  "std::vector< std::complex<double> >", &pQeta3);
	trV->Branch("pQeta4",  "std::vector< std::complex<double> >", &pQeta4);
	trV->Branch("pQetaW", "std::vector<double>", &pQetaW);
}


//////////////////
void QWPCA::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	Handle<std::vector<double> >	hEta;
	Handle<std::vector<double> >	hPt;
	Handle<std::vector<double> >	hPhi;
	Handle<std::vector<double> >	hWeight;

	iEvent.getByLabel(trackEta_,	hEta);
	iEvent.getByLabel(trackPt_,	hPt);
	iEvent.getByLabel(trackPhi_,	hPhi);
	iEvent.getByLabel(trackWeight_, hWeight);

	unsigned int sz = hEta->size();
	if ( sz == 0 ) {
		cout << " --> sz = 0 empty event" << endl;
		return;
	}

	if ( sz != hPt->size() or sz != hPhi->size() or sz != hWeight->size() ) {
		cout << " --> inconsistency" << endl;
		return;
	}

	Handle<VertexCollection> vertexCollection;
	iEvent.getByToken(vertexToken_, vertexCollection);
	VertexCollection recoVertices = *vertexCollection;
	if ( recoVertices.size() < 1 ) return;
	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
			return a.tracksSize() > b.tracksSize() ? true:false;
			});
	int primaryvtx = 0;
	double vz = recoVertices[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		return;
	}

	// centrality

	edm::Handle<int> ch;
	iEvent.getByLabel(centralityTag_,ch);
	Cent = *(ch.product());
	Mult = sz;

	std::vector<Complex> Qeta2(NETA, Complex(0.,0.));
	std::vector<Complex> Qeta3(NETA, Complex(0.,0.));
	std::vector<Complex> Qeta4(NETA, Complex(0.,0.));
	std::vector<double>  QetaW(NETA, 0.);
	double etabinW = 4.8 / NETA;
	for ( int i = 0; i < sz; i++ ) {
		int ieta = ((*hEta)[i] + 2.4)/etabinW;
		if ( ieta < 0 or ieta >= NETA ) continue;
		Qeta2[ieta] += (*hWeight)[i] * Complex(cos(2*(*Phi)[i]), sin(2*(*Phi)[i]));
		Qeta3[ieta] += (*hWeight)[i] * Complex(cos(3*(*Phi)[i]), sin(3*(*Phi)[i]));
		Qeta4[ieta] += (*hWeight)[i] * Complex(cos(4*(*Phi)[i]), sin(4*(*Phi)[i]));

		QetaW[ieta] += (*hWeight)[i];
	}

	pQeta2 = &Qeta2;
	pQeta3 = &Qeta3;
	pQeta4 = &Qeta4;
	pQetaW = &QetaW;
	trV->Fill();

	return;
}

//////////////////
QWPCA::~QWPCA()
{
	return;
}


//////////////////
void QWPCA::beginJob()
{
	return;
}

//////////////////
void QWPCA::endJob()
{
	return;
}

void QWPCA::beginRun(edm::Run const&, edm::EventSetup const&)
{
	return;
}

void QWPCA::endRun(edm::Run const&, edm::EventSetup const&)
{
	return;
}

void QWPCA::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	return;
}

void QWPCA::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(QWPCA);
