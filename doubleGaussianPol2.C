#include <cmath>

#include "TMath.h"


struct doubleGaussianPol2 {

	enum ModelComponent {total = 0, signal = 1, gaussian1 = 2, gaussian2 = 3, background = 4};

	doubleGaussianPol2(
		const ModelComponent component               = total,
		const bool           forceCommonGaussianMean = false)
		: _component(component),
			_forceCommonGaussianMean(forceCommonGaussianMean)
	{ }

	double posDefPol2 (
		const double x,
		const double p0,
		const double p1,
		const double p2)
	{
		const double linTerm = p1 + p2 * x;
		return p0 * p0 + linTerm * linTerm;
	}

	double operator() (double* vars, double* pars)
	{
		const double x = vars[0];

		const double p0     = pars[0];
		const double p1     = pars[1];
		const double p2     = pars[2];
		const double A      = pars[3];
		const double r      = sin(pars[4]) * sin(pars[4]);  // ensure [0, 1] range
		const double sigma1 = pars[5];
		const double sigma2 = pars[6];
		const double mu1    = pars[7];
		const double mu2    = (_forceCommonGaussianMean) ? pars[7] : pars[8];

		switch (_component) {
			case total:
				return A * ((1 - r) * TMath::Gaus(x, mu1, sigma1, true) + r * TMath::Gaus(x, mu2, sigma2, true)) + posDefPol2(x, p0, p1, p2);
			case signal:
				return A * ((1 - r) * TMath::Gaus(x, mu1, sigma1, true) + r * TMath::Gaus(x, mu2, sigma2, true));
			case gaussian1:
				return A * (1 - r) * TMath::Gaus(x, mu1, sigma1, true);
			case gaussian2:
				return A * r * TMath::Gaus(x, mu2, sigma2, true);
			case background:
				return posDefPol2(x, p0, p1, p2);
			default:
				return 0;
		}
	}

	ModelComponent _component;
	bool           _forceCommonGaussianMean;

};
