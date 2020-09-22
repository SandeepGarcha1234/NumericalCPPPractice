#include <cmath>
#include "stats.hpp"

double EuropeanCall::payoff(double stockPrice) const{
    //returns the intrinsic payoff for a European Call Option given the stockPrice

	if (stockPrice > strikePrice) {return stockPrice - strikePrice;}
	else{return 0;}
}

double EuropeanCall::price(const BlackScholesModel& bsm) const{
    //returns the price of a European Call Option

    //takes as a parameter, a variable of the type BlackScholesModel, which holds the stock Price, risk Free rate, . . .
    //and other relevant parameters

    double T=maturity-bsm.date;
    double numerator=log(bsm.stockPrice/strikePrice)+(bsm.riskFreeRate+bsm.volatility*bsm.volatility*0.5)*T;
    double denominator = bsm.volatility*sqrt(T);
    double d1=numerator/denominator;
    double d2=d1 - denominator;
    return bsm.stockPrice*normcdf(d1)-exp(-bsm.riskFreeRate*T)*strikePrice*normcdf(d2);
}

double EuropeanPut::payoff(double stockPrice) const{
    //returns the intrinsic payoff for a European Put Option given the stockPrice

	if (stockPrice<strikePrice){return strikePrice - stockPrice;}
	else{return 0;}
}

double EuropeanPut::price(const BlackScholesModel& bsm) const{
    //returns the price of a European Put Option
    
    //takes as a parameter, a variable of the type BlackScholesModel, which holds the stock Price, risk Free rate, . . .
    //and other relevant parameters
    double T=maturity-bsm.date;
    double numerator=log(bsm.stockPrice/strikePrice)+(bsm.riskFreeRate+bsm.volatility*bsm.volatility*0.5)*T;
    double denominator = bsm.volatility*sqrt(T);
    double d1=numerator/denominator;
    double d2=d1 - denominator;
    return bsm.stockPrice*(normcdf(d1)-1)+exp(-bsm.riskFreeRate*T)*strikePrice*(1-normcdf(d2));
}

double AmericanPut::payoff(double stockPrice) const{
    //returns the intrinsic payoff for an American Put Option given the stockPrice

	if (stockPrice<strikePrice){return strikePrice - stockPrice;}
	else{return 0;}
}