#ifndef FIRSTORDERODE
#define FIRSTORDERODE

#include <vector>
#include "matrix.hpp"


class FirstOrderODE{
//holds a first order ODE (vector valued) of the form dv/dt = F(t,v(t))
private:
	std::vector<double (*) (double t, Matrix v)> F;
	//F is function from R x Rn -> Rn
	//ie. it is of the form F(t,v) where t is the time parameter and v is a vector in Rn

	double initialTime; //corresponds to t0
	Matrix initialValue; 
	//corresponds to v0
	//should be in (column vector form)
	//ie it should be nx1 matrix

	//v0 = v(t0)


	int order;
	//this corresponds to the dimension v


public:
	FirstOrderODE(const std::vector<double (*) (double t, Matrix v)>& FF, double t0, const Matrix& y0);
	
	std::vector<double (*) (double t, Matrix v)> getRHS() const {return F;}
	double getinitialTime() const {return initialTime;}
	Matrix getinitialValue() const{return initialValue;}
	int getorder() const {return order;}

	void setRHS(std::vector<double (*) (double t, Matrix v)> FF) {F=FF;}
	void setinitialTime(double t) {initialTime =t;}
	void setinitialValue(std::vector<double> y0) {initialValue = y0;}

};

class ODEsolve{
private:
	double finalTime;
	double numberOfSteps;

	//solves a First Order ODE (that is of type FirstOrderODE)
	//the solver take uniform steps from the inital Time to the final time

public:
	ODEsolve(double tn, double N) : finalTime(tn), numberOfSteps(N) {};

	double getnumberOfSteps() {return numberOfSteps;}
	double getfinalTime() {return finalTime;}

	void setnumberOfSteps(double N) {numberOfSteps = N;}
	void setfinalTime(double tn) {finalTime = tn;}

	virtual Matrix solve(const FirstOrderODE& A) {};
	//method to the first order ODE
	//the result is returned in a form of a matrix of size (dimension of F x numberOfSteps)
	//the jth column in the result corresponds to the solution v evaluated at the jth time step

};

class EulerMethod : public ODEsolve {
public:
	EulerMethod(double tn, double N) : ODEsolve(tn,N) {}; //tn is the final time and N is the number of steps
	Matrix solve(const FirstOrderODE& A);
};

class ImprovedEulerMethod : public ODEsolve{
public:
	ImprovedEulerMethod(double tn, double N) : ODEsolve(tn,N) {}; //tn is the final time and N is the number of steps
	Matrix solve(const FirstOrderODE& A);

};

class MidpointMethod : public ODEsolve{
public:
	MidpointMethod(double tn, double N) : ODEsolve(tn,N) {}; //tn is the final time and N is the number of steps
	Matrix solve(const FirstOrderODE& A);

};

class RK4 : public ODEsolve{
public:
	RK4(double tn, double N) : ODEsolve(tn,N) {}; //tn is the final time and N is the number of steps
	Matrix solve(const FirstOrderODE& A);
};



















#endif