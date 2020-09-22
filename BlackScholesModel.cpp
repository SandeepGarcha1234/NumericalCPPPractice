#include <cmath>

#include "matrix.hpp"
#include "BlackScholesModel.hpp"

Matrix BlackScholesModel::generateMultiplePricePath(double toDate, int nSteps, double mu, const Matrix& initialPrices) const{
	//generates a random path for each initial stock price in the variable initialPrices
	
	//the variable initialPrices should be in column vector form (ie a matrix of size Mx1)
	//where M is the number of initialPrices

	//the drift of the paths is set to mu

	//the paths are returned in a matrix of size M x nSteps
	//each row corresponds to a price path, specific to an initial price

	//Please note that the price paths are correlated as the same vector of normals is used for each path generation

	assert(initialPrices.getCols()==1);

	int M = initialPrices.getRows();

	Matrix paths(M,nSteps);
	double delta_time=(toDate - date)/nSteps;
	double a = (mu - 0.5*volatility*volatility)*delta_time;
	double b = volatility * sqrt(delta_time);
	Matrix epsilon = boxMuller(nSteps);
	Matrix currentLogPrices(M,1);
	for (int i =0;i<M;i++){
		currentLogPrices(i,0)=log(initialPrices(i,0));
	}
	for (int i =0;i<M;i++){
		for (int j=0;j<nSteps;j++){
			currentLogPrices(i,0)+=a+b*epsilon(j,0);
			paths(i,j)=exp(currentLogPrices(i,0));
		}
	}

	return paths;
}

Matrix BlackScholesModel::generateMultiplePricePath(double toDate, int nSteps, const Matrix& initialPrices) const{
	//generates a random path of stock prices with the drift equal to the drift of the stock for each initial stock price

	//the variable initialPrices should be in column vector (ie of size Mx1)
	//where M is the number of initialPrices

	//the paths are returned in a matrix of size M x nSteps
	//each row corresponds to a price path specific to an initial price

	//Please note that the price paths are correlated as the same vector of normals is used for each path generation

	return generateMultiplePricePath(toDate, nSteps, drift, initialPrices);
}

Matrix BlackScholesModel::generateMultipleRiskPricePath(double toDate, int nSteps, const Matrix& initialPrices) const{
	//generates a random path for each initial stock price in the variable initialPrices
	
	//the variable initialPrices should be in column vector form (ie a matrix of size Mx1)
	//where M is the number of initialPrices

	//the drift of the paths is set to the riskFreeRate
	//ie this is the path under the risk neutral measure

	//the paths are returned in a matrix of size M x nSteps
	//each row corresponds to a price path specific to an initial price

	//Please note that the price paths are correlated as the same vector of normals is used for each path generation



	return generateMultiplePricePath(toDate, nSteps, riskFreeRate, initialPrices);
}

Matrix BlackScholesModel::generatePricePath(double toDate, int nSteps, double mu) const{
	//generates a random path of stock prices with drift mu
	//the path is returned in a matrix ("column vector form") of size nStepsx1

	Matrix path(nSteps,1); //this holds the stock price path
	double delta_time=(toDate - date)/nSteps;
	double a = (mu-0.5*volatility*volatility)*delta_time;
	double b = volatility * sqrt(delta_time);
	Matrix epsilon = boxMuller(nSteps); //generate nStep standard normals
	double currentLogPrice=log(stockPrice);
	for (int i=0;i<nSteps;i++){
		currentLogPrice+=(a+b*epsilon(i,0));
		path(i,0) = exp(currentLogPrice);
	}
	return path;
}

Matrix BlackScholesModel::generatePricePath(double toDate, int nSteps) const{
	//generates a random path of stock prices with the drift equal to the drift of the stock
	
	//the path is returned in a matrix ("column vector form") of size nStepsx1
	return generatePricePath(toDate, nSteps, drift);
}

Matrix BlackScholesModel::generateRiskPricePath(double toDate, int nSteps) const{
	//generates a random path of stock prices with the drift equal to the riskFreeRate
	//ie. this is path generation under the risk-neutral measure

	//the path is returned in a matrix ("column vector form") of size nStepsx1
	return generatePricePath(toDate, nSteps, riskFreeRate);
}

Matrix BlackScholesModel::generateBinomialLattice(double toDate, int nSteps) const{
	//returns a matrix (column vector form) of size 2*nSteps+1 x 1
	//this matrix contains the stock prices for the last two time steps
	//note that the number of stock prices (nodes) at time step n equals n+1
	//thus the total number of nodes at the last two prices is 2*nSteps+1

	//the reason for this structure is because we are using the CRR model
	//in this model u=1/d and so the stock prices repeat in the tree
	//hence we do not need to store the entire tree but only the last two layers

	//thus when we go from the ith entry to the i+1th entry, we are either jumping forward . . . 
	//in time or we are jumping back in time

	//the zero index entry in the matrix corresponds the stock price in the last time step . . .
	//where we take the bottom branch each time (ie. the stock price is (d^nSteps)*initialStockPrice)
	//the entry corresponding to index 1 if the stock price in the second to last time step . . .
	//where we take the bottom branch each time (ie. the stock price is d^(nSteps-1)*initialStockPrice)

	//More generally speaking, the even indices correspond to the last time step and . . .
	//the odd indices correspond to the second to last time step

	//we can get the i+1 entry from the ith entry by multiplying by u

	//Below we fill in the middle entry first, which is the initial stock price
	//then we work our way down by multiplying by d and work our way up by multiplying by u


	int N = 2*nSteps+1; //lattice size

	Matrix lattice(N,1);
	double delta_time = (toDate - date)/nSteps;
	double u = exp(volatility*sqrt(delta_time));
	double d = 1/u;

	lattice((N-1)/2,0) = stockPrice;//the stock price in the middle of the tree is equal to the initial stock price

	for (int i=(N-3)/2;i>=0;i--){
		lattice(i,0) = lattice(i+1,0)*d;//these are the stock prices below the middle part of the tree
	}
	for (int i=(N+1)/2;i<N;i++){
		lattice(i,0) = lattice(i-1,0)*u;//these are the stock prices above the middle part of the tree
	}

	return lattice;


}





