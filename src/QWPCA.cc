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
	  bGen_(iConfig.getUntrackedParameter<bool>("bGen", false))
	, bSim_(iConfig.getUntrackedParameter<bool>("bSim", false))
	, bEff_(iConfig.getUntrackedParameter<bool>("bEff", false))
	, minPt_(iConfig.getUntrackedParameter<double>("minPt", 1.0))
	, maxPt_(iConfig.getUntrackedParameter<double>("maxPt", 3.0))
	, centralityToken_( consumes<int>(iConfig.getParameter<edm::InputTag>("centrality")) )
	, trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackTag"))
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

	minCent_ = iConfig.getUntrackedParameter<int>("minCent", 0);
	maxCent_ = iConfig.getUntrackedParameter<int>("maxCent", 200);

	if ( trackTag_.label() == "hiGeneralTracks" ) {
		sTrackQuality = HIReco;
	} else if ( trackTag_.label() == "generalTracks" ) {
		sTrackQuality = ppReco;
	} else if ( trackTag_.label() == "hiGeneralAndPixelTracks" ) {
		sTrackQuality = Pixel;
	} else {
		sTrackQuality = trackUndefine;
	}
	if ( trackTag_.label() == "genParticles" ) {
		trackGenToken_ = consumes<reco::GenParticleCollection>(trackTag_);
	} else {
		trackToken_ = consumes<reco::TrackCollection>(trackTag_);
	}

	cout << "!!! using Track cuts " << sTrackQuality << endl;
	std::string streff = fweight_.label();
	if ( streff == std::string("NA") ) {
		std::cout << "!!! eff NA" << std::endl;
		bEff_ = false;
	} else {
		TFile * fEffFak = new TFile(streff.c_str());
		std::cout << "!!! Using particle weight " << streff << std::endl;
		if ( bEff_ ) {
			std::cout << "!!! Apply Eff correction" << std::endl;
			if ( streff == string("PbPb_MB_TT_5TeV_v2.root") or streff == string("PbPb_dijet_TT_5TeV_v2.root") ) {
				TH2D * h = (TH2D*) fEffFak->Get("rTotalEff3D_0_5");
				for ( int c = 0; c < 10; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_5_10");
				for ( int c = 10; c < 20; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_10_30");
				for ( int c = 20; c < 60; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_30_50");
				for ( int c = 60; c < 100; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_50_100");
				for ( int c = 100; c < 200; c++ ) {
					hEff_cbin[c] = h;
				}
			} else if ( streff == std::string("Hydjet_eff_mult_v1.root") ) {
				TH2D * h = (TH2D*) fEffFak->Get("rTotalEff3D_1");
				for ( int c = 0; c < 120; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_2");
				for ( int c = 120; c < 260; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_3");
				for ( int c = 260; c < 400; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_4");
				for ( int c = 400; c < 800; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("rTotalEff3D_5");
				for ( int c = 800; c < 2000; c++ ) {
					hEff_cbin[c] = h;
				}
			} else if ( streff == std::string("EffCorrectionsPixel_TT_pt_0_10_v2.root") ) {
				TH2D * h = (TH2D*) fEffFak->Get("Eff_0_5");
				for ( int c = 0; c < 10; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("Eff_5_10");
				for ( int c = 10; c < 20; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("Eff_10_30");
				for ( int c = 20; c < 60; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("Eff_30_50");
				for ( int c = 60; c < 100; c++ ) {
					hEff_cbin[c] = h;
				}
				h = (TH2D*) fEffFak->Get("Eff_50_100");
				for ( int c = 100; c < 200; c++ ) {
					hEff_cbin[c] = h;
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
///
bool
QWPCA::TrackQuality_ppReco(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
        if ( itTrack->charge() == 0 ) {
                return false;
        }
	if ( fabs(itTrack->eta()) > 2.4 ) return false;
	if ( itTrack->pt() < minPt_ or itTrack->pt() > maxPt_ ) return false;
        if ( !itTrack->quality(reco::TrackBase::highPurity) ) {
                return false;
        }
        if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) {
                return false;
        }
	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();
        double d0 = -1.* itTrack->dxy(v1);
        double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
        if ( fabs( d0/derror ) > d0d0error_ ) {
                return false;
        }
        double dz=itTrack->dz(v1);
        double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
        if ( fabs( dz/dzerror ) > dzdzerror_ ) {
                return false;
        }
        return true;
}

///
bool
QWPCA::TrackQuality_HIReco(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
	if ( itTrack->charge() == 0 ) return false;
	if ( itTrack->pt() < minPt_ or itTrack->pt() > maxPt_ ) return false;
	if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
	if ( fabs(itTrack->eta()) > 2.4 ) return false;
	if ( itTrack->numberOfValidHits() < 11 ) return false;
	if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
		return false;
	}
	if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) {
		return false;
	}
	if (
		itTrack->originalAlgo() != 4 and
		itTrack->originalAlgo() != 5 and
		itTrack->originalAlgo() != 6 and
		itTrack->originalAlgo() != 7
	) {
		return false;
	}

	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();
	double d0 = -1.* itTrack->dxy(v1);
	double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
	if ( fabs( d0/derror ) > d0d0error_ ) {
		return false;
	}

	double dz=itTrack->dz(v1);
	double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
	if ( fabs( dz/dzerror ) > dzdzerror_ ) {
		return false;
	}
	return true;
}

///
bool
QWPCA::TrackQuality_Pixel(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
	if ( itTrack->charge() == 0 ) return false;
	if ( itTrack->pt() < minPt_ or itTrack->pt() > maxPt_ ) return false;
	if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
	if ( fabs(itTrack->eta()) > 2.4 ) return false;
	bool bPix = false;
	int nHits = itTrack->numberOfValidHits();
//	std::cout << __LINE__ << "\tnHits = " << nHits << std::endl;
	if ( itTrack->pt() < 2.4 and (nHits==3 or nHits==4 or nHits==5 or nHits==6) ) bPix = true;
	if ( not bPix ) {
		if ( nHits < 11 ) return false;
		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
			return false;
		}
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) {
			return false;
		}
		if (
			itTrack->pt() > 2.4 and
			itTrack->originalAlgo() != 4 and
			itTrack->originalAlgo() != 5 and
			itTrack->originalAlgo() != 6 and
			itTrack->originalAlgo() != 7
		) {
			return false;
		}

		int primaryvtx = 0;
		math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
		double vxError = recoVertices[primaryvtx].xError();
		double vyError = recoVertices[primaryvtx].yError();
		double vzError = recoVertices[primaryvtx].zError();
		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		if ( fabs( d0/derror ) > d0d0error_ ) {
			return false;
		}

		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs( dz/dzerror ) > dzdzerror_ ) {
			return false;
		}
	}
	return true;
}
///
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

	double vz = recoVertices[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		return;
	}
	t.vz = vz;

	// centrality

	edm::Handle<int> ch;
	iEvent.getByToken(centralityToken_,ch);
	int Cent = *(ch.product());
	if ( minCent_ >= 0 and Cent < minCent_ ) {
		return;
	} else if ( minCent_ < 0 and Cent < -minCent_ ) {
		Cent = -minCent_;
	}
	if ( maxCent_ >= 0 and Cent > maxCent_ ) {
		return;
	} else if ( maxCent_ < 0 and Cent > -maxCent_ ) {
		Cent = -maxCent_;
	}

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByToken(trackToken_,tracks);

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {
		if ( sTrackQuality == HIReco and not TrackQuality_HIReco(itTrack, recoVertices) ) continue;
		else if ( sTrackQuality == ppReco and not TrackQuality_ppReco(itTrack, recoVertices) ) continue;
		else if ( sTrackQuality == Pixel  and not TrackQuality_Pixel (itTrack, recoVertices) ) continue;

		t.Charge[t.Mult] = itTrack->charge();
		t.Pt[t.Mult] = itTrack->pt();
		t.Eta[t.Mult] = itTrack->eta();
		t.Phi[t.Mult] = itTrack->phi();

		if ( bEff_ ) {
			double eff = hEff_cbin[Cent]->GetBinContent( hEff_cbin[Cent]->FindBin(itTrack->eta(), itTrack->pt()) );
			t.weight[t.Mult] = 1.0 / eff;
		} else {
			t.weight[t.Mult] = 1.0;
		}
//		cout << itTrack->eta() << "\t" << t.weight[t.Mult] << endl;

		t.Mult++;
	}
	t.Cent = Cent;
//	std::cout << " t.Mult = " << t.Mult << std::endl;

	return;
}

//////////////////
void QWPCA::analyzeMC(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

	double vz = recoVertices[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		return;
	}
	t.vz = vz;

	// centrality

	edm::Handle<int> ch;
	iEvent.getByToken(centralityToken_,ch);
	int Cent = *(ch.product());
	if ( minCent_ >= 0 and Cent < minCent_ ) {
		return;
	} else if ( minCent_ < 0 and Cent < -minCent_ ) {
		Cent = -minCent_;
	}
	if ( maxCent_ >= 0 and Cent > maxCent_ ) {
		return;
	} else if ( maxCent_ < 0 and Cent > -maxCent_ ) {
		Cent = -maxCent_;
	}

	// track
	Handle< std::vector<GenParticle> > tracks;
	iEvent.getByToken(trackGenToken_,tracks);


	for(std::vector<GenParticle>::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {

		if ( itTrack->status()!=1 ) continue;
		if ( itTrack->charge() == 0 ) continue;
		if ( itTrack->pt() < minPt_ or itTrack->pt() > maxPt_ ) continue;
		if ( fabs(itTrack->eta()) > 2.4 ) continue;

		t.Charge[t.Mult] = itTrack->charge();
		t.Pt[t.Mult] = itTrack->pt();
		t.Eta[t.Mult] = itTrack->eta();
		t.Phi[t.Mult] = itTrack->phi();

		t.weight[t.Mult] = 1.0;

		t.Mult++;
	}
	t.Cent = Cent;

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
