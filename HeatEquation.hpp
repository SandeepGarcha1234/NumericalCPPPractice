#ifndef HEATEQUATION.HPP
#define HEATEQUATION.HPP

#include "matrix.hpp"


typedef double (*fptr) (double); 
//typedefines the type:pointer to a function that takes a single double as an argument . . .
//and gives it the name fptr

class HeatEquation{
//stores the heat equation du/dt= diffusivity* d2u/dx2 . . . where u(t,x) is a function from R2 to R
//t is the time variable
//note that the space variable is in R
//with the x variable being bounded at the left by the variable leftX and at the right by rightX

//the boundary condition at the left boundary(leftX) is given by the function (*leftBoundary)(t)
//note that the variable leftBoundary is a function pointer

//the boundary condition at the right boundary(rightX) is given by the function (*rightBoundary)(t)
//note that the variable rightBoundary is a function pointer

//the initial condition (ie. u(0,x)) is given by the function (*initialCondition)(x)
//note that the variable initialCondition is a function pointer

private:
	double leftX;
	double rightX;

	double diffusivity;

    fptr leftBoundary;
    fptr rightBoundary;
    fptr initialCondition;

public:
	double getleftX() const{return leftX;}
	double getrightX() const{return rightX;}
	double getdiffusivity() const{return diffusivity;}
	fptr getleftBoundary() const {return leftBoundary;}
	fptr getrightBoundary() const{return rightBoundary;}
	fptr getinitialCondition() const{return initialCondition;}

	void setleftX(double x) {leftX=x;}
	void setrightX(double x) {rightX=x;}
	void setdiffusivity(double r) {
		assert(r>0);
		diffusivity=r;
	}
	void setleftBoundary(fptr f) {leftBoundary=f;}
	void setrightBoundary(fptr f) {rightBoundary=f;}
	void setinitialCondition(fptr f) {initialCondition=f;}

};

class HeatSolver{
private:
	double delta_time;
	double delta_x;
	double finalTime;
public:
	HeatSolver(double dt, double dx, double Tn) : delta_time(dt), delta_x(dx), finalTime(Tn) {};

	double getdt() const {return delta_time;}
	double getdx() const {return delta_x;}
	double getfinalTime() const{return finalTime;}

	void setdt(double dt) {delta_time=dt;}
	void setdx(double dx) {delta_x=dx;}
	void setfinalTime(double Tn) {finalTime=Tn;}


	virtual Matrix solve(const HeatEquation& A){};

};

class ExplicitHeatSolver : public HeatSolver{
public:
	ExplicitHeatSolver(double dt, double dx, double Tn) : HeatSolver(dt,dx,Tn) {};
	Matrix solve(const HeatEquation& A);
};

class ImplicitHeatSolver : public HeatSolver{
public:
	ImplicitHeatSolver(double dt, double dx, double Tn) : HeatSolver(dt,dx,Tn) {};
	Matrix solve(const HeatEquation& A);

};

class CrankNicolsonHeatSolver : public HeatSolver{
public:
	CrankNicolsonHeatSolver(double dt, double dx, double Tn) : HeatSolver(dt,dx,Tn) {};
	Matrix solve(const HeatEquation& A);
};














#endif