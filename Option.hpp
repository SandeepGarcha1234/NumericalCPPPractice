#ifndef OPTION.HPP
#define OPTION.HPP

#include "matrix.hpp"

class Option{
public:
	virtual ~Option(){};
	double maturity; //date when the option expires

};

class PathIndependent : public Option {
public:
	double strikePrice;
	virtual ~PathIndependent(){};
	virtual double payoff(double stockPrice) const = 0;

};

class European : public PathIndependent {
public:
	virtual ~European() {};

};

class American : public PathIndependent {
public:
	virtual ~American() {};

};

class PathDependent : public Option {
public:
	double strikePrice;
	virtual double payoff(const Matrix& stockPricePath) const = 0;
	//since the payoff is dependent on the path, this payoff method takes a matrix . . .
	//of size Nx1 (column vector form) and then returns the payoff

};

class EuropeanCall : public European{
public:
	double payoff(double stockPrice) const;
	double price(const BlackScholesModel& bsm) const;

};

class EuropeanPut : public European{
public:
	double payoff(double stockPrice) const;
	double price(const BlackScholesModel& bsm) const;
};

class AmericanPut : public American{
	double payoff(double stockPrice) const;
};




#endif