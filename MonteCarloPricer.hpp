#ifndef MonteCarloPricer.hpp
#define MonteCarloPricer.hpp

#include "BlackScholesModel.hpp"
#include "Option.hpp"


class MonteCarloPricer{
public:
	MonteCarloPricer(int N) : nScenerios(N) {};
	
	int nScenerios;
	//this is the number of Scenarios, that we wish to run for the Monte Carlo Simulation

	double pricer(const European& Eur, const BlackScholesModel& bsm) const;
	double pricer(const American& Am, const BlackScholesModel& bsm, int N) const;
	//the parameter N is the number of time steps for each sample stock price path in the simulation



	double pricer(const PathDependent& Ex, const BlackScholesModel& bsm, int N) const;
	//the parameter N is the number of time steps for each sample stock price path in the simulation
};


#endif