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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "QWAna/QWPCA/interface/QWPCA.h"


QWPCA::QWPCA(const edm::ParameterSet& iConfig):
	  bGen_(iConfig.getUntrackedParameter<bool>("bGen", false))
	, bSim_(iConfig.getUntrackedParameter<bool>("bSim", false))
	, bEff_(iConfig.getUntrackedParameter<bool>("bEff", false))
	, minPt_(iConfig.getUntrackedParameter<double>("minPt", 1.0))
	, maxPt_(iConfig.getUntrackedParameter<double>("maxPt", 3.0))
	, centralityToken_( consumes<int>(iConfig.getParameter<edm::InputTag>("centrality")) )
	, trackToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackTag")))
	, vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")) )
	, fweight_( iConfig.getUntrackedParameter<edm::InputTag>("fweight", std::string("NA")) )
{

	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror", 3.);
	d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error", 3.);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt", 0.1);
	minvz_ = iConfig.getUntrackedParameter<double>("minvz", -1.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz", 15.);
	minEta_ = iConfig.getUntrackedParameter<double>("minEta", -2.4);
	maxEta_ = iConfig.getUntrackedParameter<double>("maxEta", 2.4);

	minCent_ = iConfig.getUntrackedParameter<int>("minCent", -1);
	maxCent_ = iConfig.getUntrackedParameter<int>("maxCent", 500);

	std::string streff = fweight_.label();
	if ( streff == std::string("NA") ) {
		std::cout << "!!! eff NA" << std::endl;
		bEff_ = false;
	} else {
		TFile * fEffFak = new TFile(streff.c_str());
		std::cout << "!!! Using particle weight " << streff << std::endl;
		if ( bEff_ ) {
			std::cout << "!!! Apply Eff correction" << std::endl;
			for ( int i = 0; i < 20; i++ ) {
				if ( streff == std::string("Hydjet_eff_mult_v1.root") ) {
					TH2D * h = (TH2D*) fEffFak->Get("rTotalEff3D_1");
					for ( int c = 0; c < 200; c++ ) {
						hEff_cbin[c] = h;
					}
				}
			}
			std::cout << "!!! eff histo done" << std::endl;
		}
	}

	edm::Service<TFileService> fs;
	trV = fs->make<TTree>("trV", "trV");

	trV->Branch("cent", &(t.Cent), "cent/I");
	trV->Branch("mult", &(t.Mult), "mult/I");
//	trV->Branch("noff", &(t.Noff), "mult/I");

	trV->Branch("pQeta2",  "std::vector< std::complex<double> >", &pQeta2);
	trV->Branch("pQeta3",  "std::vector< std::complex<double> >", &pQeta3);
	trV->Branch("pQeta4",  "std::vector< std::complex<double> >", &pQeta4);
	trV->Branch("pQetaW", "std::vector<double>", &pQetaW);
}


//////////////////
void QWPCA::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//std::cout << __LINE__ << "\tstart ana" << std::endl;
	if ( bGen_ ) analyzeMC(iEvent, iSetup);
	else analyzeData(iEvent, iSetup);

	if ( bSim_ ) overRide();

	//std::cout << __LINE__ << "\t" << t.Mult << std::endl;
	if ( t.Mult == 0 ) return;

	std::vector<Complex> Qeta2(NETA, Complex(0.,0.));
	std::vector<Complex> Qeta3(NETA, Complex(0.,0.));
	std::vector<Complex> Qeta4(NETA, Complex(0.,0.));
	std::vector<double>  QetaW(NETA, 0.);

	double etabinW = 4.8 / NETA;
	for ( int i = 0; i < t.Mult; i++ ) {
		int ieta = (t.Eta[i] + 2.4)/etabinW;
		if ( ieta < 0 or ieta >= NETA ) continue;
		Qeta2[ieta] += t.weight[i] * Complex(cos(2*t.Phi[i]), sin(2*t.Phi[i]));
		Qeta3[ieta] += t.weight[i] * Complex(cos(3*t.Phi[i]), sin(3*t.Phi[i]));
		Qeta4[ieta] += t.weight[i] * Complex(cos(4*t.Phi[i]), sin(4*t.Phi[i]));

		QetaW[ieta] += t.weight[i];
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
void QWPCA::overRide()
{
	t.Mult = 7;
	t.Cent = 150;
	t.Charge[0] = 1;
	t.Charge[1] = 1;
	t.Charge[2] = 1;
	t.Charge[3] = -1;
	t.Charge[4] = -1;
	t.Charge[5] = -1;
	t.Charge[6] = -1;

	t.Eta[0] = 0;
	t.Eta[1] = 0;
	t.Eta[2] = 0;
	t.Eta[3] = 0;
	t.Eta[4] = 0;
	t.Eta[5] = 0;
	t.Eta[6] = 0;
	t.Eta[7] = 0;

	t.Phi[0] = 0;
	t.Phi[1] = 3.1;
	t.Phi[2] = 0.1;
	t.Phi[3] = 0.;
	t.Phi[4] = 0.1;
	t.Phi[5] = -3.1;
	t.Phi[6] = 3.1;

	t.weight[0] = 1;
	t.weight[1] = 1;
	t.weight[2] = 1;
	t.weight[3] = 1;
	t.weight[4] = 1;
	t.weight[5] = 1;
	t.weight[6] = 1;
}

//////////////////
//int QWPCA::getNoffCent(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& Noff)
//{
//	// very hard coded Noff track centrality cut
//	using namespace edm;
//	using namespace reco;
//	//      int Noff = 0;
//
//	Handle<VertexCollection> vertexCollection;
//	iEvent.getByToken(vertexToken_, vertexCollection);
//	const VertexCollection * recoVertices = vertexCollection.product();
//
//	if ( recoVertices.size() < 1 ) return;
//	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
//			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
//			return a.tracksSize() > b.tracksSize() ? true:false;
//			});
//
//	int primaryvtx = 0;
//	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
//	double vxError = (*recoVertices)[primaryvtx].xError();
//	double vyError = (*recoVertices)[primaryvtx].yError();
//	double vzError = (*recoVertices)[primaryvtx].zError();
//
//
//	Handle<TrackCollection> tracks;
//	iEvent.getByToken(trackToken_,tracks);
//	for(TrackCollection::const_iterator itTrack = tracks->begin();
//			itTrack != tracks->end();
//			++itTrack) {
//
//		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
//		if ( itTrack->charge() == 0 ) continue;
//		if ( itTrack->pt() < 0.4 ) continue;
//
//		double d0 = -1.* itTrack->dxy(v1);
//		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
//		double dz=itTrack->dz(v1);
//		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
//		if ( fabs(itTrack->eta()) > 2.4 ) continue;
//		if ( fabs( dz/dzerror ) > 3. ) continue;
//		if ( fabs( d0/derror ) > 3. ) continue;
//		if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
//
//		Noff++;
//	}
//
//	int cent = ;
//	while ( [cent] <= Noff ) cent--;
//	return cent;
//}
//////////////////
void QWPCA::analyzeData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	t.Mult = 0;

	Handle<VertexCollection> vertexCollection;
	iEvent.getByToken(vertexToken_, vertexCollection);
	VertexCollection recoVertices = *vertexCollection;
	if ( recoVertices.size() < 1 ) return;
	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
			return a.tracksSize() > b.tracksSize() ? true:false;
			});

	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();

	double vz = recoVertices[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		//std::cout << __LINE__ << std::endl;
		return;
	}
	t.vz = vz;

	// centrality

	edm::Handle<int> ch;
	iEvent.getByToken(centralityToken_,ch);
	t.Cent = *(ch.product());
	if ( t.Cent < 0 or t.Cent >= 200 ) {
		return;
	}

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByToken(trackToken_,tracks);

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {
		if ( itTrack->charge() == 0 ) continue;
		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
		if ( itTrack->pt() > maxPt_ or itTrack->pt() < minPt_ ) continue;
		if ( itTrack->eta() > maxEta_ or itTrack->eta() < minEta_ ) continue;
		int nHits = itTrack->numberOfValidHits();
		if ( itTrack->pt() > 2.4 and nHits < 11 ) continue;
		if ( itTrack->pt() < 2.4 and nHits!=3 and nHits!=4 and nHits!=5 and nHits!=6) continue;

		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
//		if ( itTrack->hitPattern().pixelLayersWithMeasurement() == 0 ) continue;

		if ( itTrack->originalAlgo() != 4 and
			itTrack->originalAlgo() != 5 and
			itTrack->originalAlgo() != 6 and
			itTrack->originalAlgo() != 7
		) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		if ( fabs( d0/derror ) > d0d0error_ ) continue;

		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;

		t.Charge[t.Mult] = itTrack->charge();
		t.Pt[t.Mult] = itTrack->pt();
		t.Eta[t.Mult] = itTrack->eta();
		t.Phi[t.Mult] = itTrack->phi();

		if ( bEff_ ) {
			double eff = hEff_cbin[t.Cent]->GetBinContent( hEff_cbin[t.Cent]->FindBin(itTrack->eta(), itTrack->pt()) );
			t.weight[t.Mult] = 1.0 / eff;
		} else {
			t.weight[t.Mult] = 1.0;
		}

		t.Mult++;
	}

	return;
}

//////////////////
void QWPCA::analyzeMC(const edm::Event&, const edm::EventSetup&)
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
