#include <cmath>

#include "MonteCarloPricer.hpp"
#include "matrix.hpp"

double MonteCarloPricer::pricer(const European& Eur, const BlackScholesModel& bsm) const{
	//prices a European Option by using Monte Carlo Simulations

	double sum=0; //holds the sum of the payoffs for each path
	for(int i=0;i<nScenerios;i++){
		Matrix v = bsm.generateRiskPricePath(Eur.maturity, 1);
		//since the price of a European option only depends on the terminal price . . .
		//we only generate a path with one time step

		sum+=Eur.payoff(v(0,0));
	}
	double mean = sum/nScenerios; //expected value of the payoffs
	double T = Eur.maturity - bsm.date;
	return exp(-bsm.riskFreeRate*T)*mean; //discounted expected value of the payoffs

}

double price(const American& Am, const BlackScholesModel& bsm, int N) {
	//simulates on stock price path and returns the payoff for the option based on that path
	//the parameter N is the number of time steps in the stock price path

	Matrix v = bsm.generateRiskPricePath(Am.maturity,N);
	double dt = (Am.maturity - bsm.date)/N;
	double discountFactor = exp(-bsm.riskFreeRate*dt);
	
	double currentPrice = Am.payoff(v(N-1,0));
	//this is the price of the option at the terminal time step
	//this variable will be modified as we move back in time through the stock price path

	for (int i = N-2; i>=0; i--){
		double discountedPrice = discountFactor * currentPrice;
		double immediatePayoff = Am.payoff(v(i,0));
		currentPrice = discountedPrice > immediatePayoff ? discountedPrice : immediatePayoff;
		//the current price of the option is the maximum the discounted price and the intrinsic price
	}

	return currentPrice;

}

double MonteCarloPricer::pricer(const American& Am, const BlackScholesModel& bsm, int N) const{
	//prices an American Option by using Monte Carlo Simulations

	//the parameter N is the number of time steps in the stock price path

	double sum=0;//holds the sum of the payoffs for each path
	for (int i=0;i<nScenerios;i++){
		sum+=price(Am,bsm,N);
	}

	double mean = sum/nScenerios; 
	//dicounted expected value
	//note that the discounting happens in the price method of the Monte Carlo class . . .
	//so we do not need to do it here

	return mean;

}

double MonteCarloPricer::pricer(const PathDependent& Ex, const BlackScholesModel& bsm, int N) const{
	//prices a Path dependent Option by using Monte Carlo Simulations

	double sum=0;//holds the sum of the payoffs for each path
	for (int i=0;i<nScenerios;i++){
		Matrix v = bsm.generateRiskPricePath(Ex.maturity,N);
		sum+=Ex.payoff(v);
	}

	double dt = (Am.maturity - bsm.date)/N;
	double discountFactor = exp(-bsm.riskFreeRate*dt);

	double mean = sum/nScenerios; //expected value of the payoffs
	return discountFactor*mean; //discounted expected value of the payoffs
}











