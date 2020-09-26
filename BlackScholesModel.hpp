#ifndef BlackScholesModel.HPP
#define BlackScholesModel.HPP

#include "matrix.hpp"

class BlackScholesModel{
//holds parameters relevant to the Black Scholes Model for a particular stock 

public:	
	double stockPrice;
	double date; //taken to be a real number
	double drift;
	double riskFreeRate; //continuously compounded rate
	double volatility;

	//Methods for generating paths for Monte Carlo Pricing

	Matrix generatePricePath(double toDate, int nSteps) const;
	Matrix generateRiskPricePath(double toDate, int nSteps) const;
	Matrix generateMultiplePricePath(double toDate, int nSteps, const Matrix& initialPrices) const;
	Matrix generateMultipleRiskPricePath(double toDate, int nSteps, const Matrix& initialPrices) const;

	//Methods for getenerating lattice for Binomial Tree Pricing

	Matrix generateBinomialLattice(double toDate, int nSteps) const;

private:
	Matrix generatePricePath(double toDate, int nSteps, double mu) const;
	Matrix generateMultiplePricePath(double toDate, int nSteps, double mu, const Matrix& initialPrices) const;

};

#endif
