#include <correlations/Types.hh>
#include <correlations/Result.hh>
#include <correlations/QVector.hh>
#include <correlations/closed/FromQVector.hh>
// for the fuxking pp rereco centrality
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include <complex>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TH2D.h"

#pragma link C++ class std::vector< std::complex<double> >+;

typedef std::complex<double> Complex;
// event structure
const int NMAX_TRK = 10000;
typedef struct QWEvent_ {
	int     Cent;
	int     Mult;
	double  vz;
	int     Noff;
	double  Pt[NMAX_TRK];
	double  Eta[NMAX_TRK];
	double  Phi[NMAX_TRK];
	int     Charge[NMAX_TRK];
	double  weight[NMAX_TRK];
	int     RunId;
	int     EventId;
} QWEvent;


class QWQVector : public edm::EDAnalyzer {
public:
	explicit QWQVector(const edm::ParameterSet&);
	~QWQVector();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	////////////////////////////////
	void analyzeData(const edm::Event&, const edm::EventSetup&);
	void analyzeMC(const edm::Event&, const edm::EventSetup&);
//	int  getNoffCent(const edm::Event&, const edm::EventSetup&, int&);

	void overRide();

	bool bGen_;
	bool bSim_;
	bool bEff_;

	double minPt_;
	double maxPt_;

	edm::EDGetTokenT<int>                           centralityToken_;
	edm::EDGetTokenT<reco::TrackCollection>		trackToken_;
	std::vector<int>				algoParameters_;
	edm::EDGetTokenT<reco::VertexCollection>	vertexToken_;
	edm::InputTag					fweight_;

	double	dzdzerror_;
	double	d0d0error_;
	double	pterrorpt_;
	double  minvz_;
	double  maxvz_;
	double  minEta_;
	double  maxEta_;

	int	minCent_;
	int	maxCent_;

	TH2D *	hEff_cbin[200];

	QWEvent t;

	TTree * trV;
	static const int NETA = 48;
	std::vector<std::complex<double>>	* pQeta2;
	std::vector<std::complex<double>>	* pQeta3;
	std::vector<std::complex<double>>	* pQeta4;
	std::vector<double>			* pQetaW;
};


