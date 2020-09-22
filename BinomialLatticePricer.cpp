#include <cassert>
#include <cmath>

#include "BinomialLatticePricer.hpp"

Tree::Tree(int n) : size(n) {
	start = new double* [n+1];
	for (int i=0;i<n+1;i++){
		start[i] = new double[i+1];
		//start[i] corresponds to the ith time step
		//the ith time step must contain i+1 roots

	}
};

Tree::~Tree() {
	for (int i=0;i<size+1;i++){
		delete[] start[i];
	}

	delete[] start;
}

void Tree::assign(const Tree& T) {
	//used for the copy constructor and the assignment operator

	size = T.getsize();
	start = new double* [size+1];
	for (int i=0;i<size+1;i++){
		start[i] = new double [i+1];
	}

	for (int i =0; i<size+1;i++){
		for (int j=0;j<i+1;j++){
			start[i][j] = T(i,j);
		}
	}
}

Tree::Tree(const Tree& T){
	assign(T);
}

Tree& Tree::operator=(const Tree& T){
	//deletes the previous tree held in the variable . . .
	//before assigning it a new tree

	for (int i=0;i<size+1;i++){
		delete[] start[i];
	}
	delete[] start;
	assign(T);
	return *this;
}

double& Tree::operator()(int i, int j) const{
	//returns the value held in the ith time step in the jth root . . .
	//counting from the bottom, which correspond to j=0
	//ie this is zero-based indexing

	assert(i>=0 && i<size+1 && j<=i);
	//the condition i<size+1 assures i is less than the total number of time steps
	//the condition j<=i is to ensure that we are not accessing a node that is not in the ith time step
	//since there i+1 nodes in the ith time step, and we are using zero-based indexing . . . 
	//the maximum index of the node is i

	return start[i][j];
}

double& Tree::operator()(int i, int j) {
	//returns the value held in the ith time step in the jth root . . .
	//counting from the bottom, which correspond to j=0
	//ie this is zero-based indexing

	assert(i>=0 && i<size+1 && j<=i);
	//the condition i<size+1 assures i is less than the total number of time steps
	//the condition j<=i is to ensure that we are not accessing a node that is not in the ith time step
	//since there i+1 nodes in the ith time step, and we are using zero-based indexing . . . 
	//the maximum index of the node is i

	return start[i][j];
}

double factorial(int n) {
	//return n! where n is a non-negatice integer


	if (n==0){return 1;}
	else {return n*factorial(n-1);}
}

double choose(int n, int k){
	//return n choose k  . . .
	//where n and k are non-negative integers

	double product =1;
	for(int i=0;i<k;i++){
		product*=n-i;
	}
	return product/factorial(k);
}


double BinomialLatticePricer::pricer(const European& Eur, const BlackScholesModel& bsm) const{
	Matrix stockPriceLattice = bsm.generateBinomialLattice(Eur.maturity,steps);
	//generates the stock prices for the last two layers

	double dt = (Eur.maturity - bsm.date)/steps;
	double u = exp(bsm.volatility*sqrt(dt));
	double d = 1/u;
	double p = (exp(bsm.riskFreeRate*dt) - d)/(u-d);
	int n = steps;

	double sum = 0;

	for (int i=0;i<steps+1;i++){
		sum+=choose(n,i)*pow(p,i)*pow(1-p,n-i)*Eur.payoff(stockPriceLattice(2*i,0));
		//calculates the expected value of the option under the risk neutral measure
		//note that this loop only takes even indexed values from the stockPriceLattice . . . 
		//as these are the values corresponding to the last time step
	}

	return sum*exp(-bsm.riskFreeRate * (Eur.maturity - bsm.date)); //discounts the expected value and returns the result
}

double BinomialLatticePricer::pricer(const American& Am, const BlackScholesModel& bsm) const{
	Matrix stockPriceLattice = bsm.generateBinomialLattice(Am.maturity,steps);
	//generates the stock prices for the last two layers

	double dt = (Am.maturity - bsm.date)/steps;
	double u = exp(bsm.volatility*sqrt(dt));
	double d = 1/u;
	double p = (exp(bsm.riskFreeRate*dt) - d)/(u-d);

	double discountFactor = exp(-bsm.riskFreeRate*dt);

	Tree T(steps); 
	//this tree will hold the value of the American option for each time step and node
	//T(i,j) will correspond to the option value at time step i at node j

	for(int i=0;i<steps+1;i++){
		T(steps,i) = Am.payoff(stockPriceLattice(2*i,0));
		//the value of the American option in the last time step is just its intrinsic value
		//note we only take the even indexed values from the stockPriceLattice . . .
		//as these are the values corresponding to the last time step

	}

	for (int i=steps-1;i>=0;i--){
		for (int j=0;j<i+1;j++){
			double discountedPrice = discountFactor*((1-p)*T(i+1,j)+p*T(i+1,j+1));
			//this is the discounted expected value from the two nodes in the nex time step
			//T(i+1,j) is the option value at time step i+1 if the stock goes down(by factor d) when  . . .
			//going forward in time from time step i and node j
			//And T(i+1,j) is the option value at time step i+1 if the stock goes up(by factor u) when  . . .
			//going forward in time from time step i and node j

			int k = steps - i;
			//as we go back in the time layers, the bottom branch and the top branch of the stockPriceLattice
			//are not longer achievable values for the stock


			double immediatePayoff = Am.payoff(stockPriceLattice(k+2*j,0));
			//this calculates the intrinsic value of the option
			//we multiply the j by two since the indices j and j+1 correspond to diffent time layers



			T(i,j) = discountedPrice > immediatePayoff ? discountedPrice : immediatePayoff;
			//compares the dicounted price with the intrinsic price
			//The price of the option at that node is the maximum of the two
		}
		
	}

	return T(0,0);

}


