// for the fuxking pp rereco centrality

#include <complex>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TH2D.h"


typedef std::complex<double> Complex;


class QWPCA : public edm::EDAnalyzer {
public:
	explicit QWPCA(const edm::ParameterSet&);
	~QWPCA();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	////////////////////////////////

		edm::InputTag					trackEta_;
		edm::InputTag					trackPhi_;
		edm::InputTag					trackWeight_;

		edm::InputTag					vertexTag_;
		edm::InputTag					centralityTag_;

		double	minvz_, maxvz_;

		unsigned int	nvtx_;

	TTree * trV;
	static const int NETA = 48;
	std::vector<std::complex<double>>	* pQeta2;
	std::vector<std::complex<double>>	* pQeta3;
	std::vector<std::complex<double>>	* pQeta4;
	std::vector<double>			* pQetaW;
	int	Cent;
	int	Mult;
};


