#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#include "stats.hpp"
#include "matrix.hpp"

const double PI = 3.141592653589793;

double normcdf(double x){
	//returns the cumulative distribution function evaluated at x
    if (x>=0){
        double k=1.0/(1.0+0.2316419*x);
        return 1-(1.0/sqrt(2*PI))*exp(-x*x/2.0)*
        k*(0.319381530+k*(-0.356563782+k*(1.781477937+k*(-1.821255978+1.330274429*k))));
    }
    else{return 1-normcdf(-x);}
}


double horner0(double x, double a0) {return a0;}
double horner1(double x, double a0, double a1) {return a0+x*horner0(x,a1);}
double horner2(double x, double a0, double a1,double a2) {return a0+x*horner1(x,a1,a2);}
double horner3(double x, double a0, double a1,double a2,double a3) {return a0+x*horner2(x,a1,a2,a3);}
double horner4(double x, double a0, double a1,double a2,double a3,double a4) {return a0+x*horner3(x,a1,a2,a3,a4);}
double horner5(double x, double a0, double a1,double a2,double a3,double a4,double a5) {return a0+x*horner4(x,a1,a2,a3,a4,a5);}
double horner6(double x, double a0, double a1,double a2,double a3,double a4,double a5,double a6) {return a0+x*horner5(x,a1,a2,a3,a4,a5,a6);}
double horner7(double x, double a0, double a1,double a2,double a3,double a4,double a5,double a6, double a7) {return a0+x*horner6(x,a1,a2,a3,a4,a5,a6,a7);}
double horner8(double x, double a0, double a1,double a2,double a3,double a4,double a5,double a6, double a7,double a8) {return a0+x*horner7(x,a1,a2,a3,a4,a5,a6,a7,a8);}




double norminv(double x){
	//returns the xth percentile from the normal distribution

	assert(x>=0 && x<=1);


	double y=x-0.5;
	if (abs(y)<=0.42){
		double r=y*y;

		double a0=2.50662823884;
		double a1=-18.61500062529;
		double a2=41.39119773534;
		double a3=-25.44106049637;
		double b1=-8.47351093090;
		double b2=23.08336743743;
		double b3=-21.06224101826;
		double b4=3.13082909833;


		return y*horner3(r,a0,a1,a2,a3)/horner4(r,1.0,b1,b2,b3,b4);
	}

	else if(y<0){
		double r=x;
		double s=log(-log(r));

		double c0=0.3374754822726147;
		double c1=0.9761690190917186;
		double c2=0.1607979714918209;
		double c3=0.0276438810333863;
		double c4=0.0038405729373609;
		double c5=0.0003951896511919;
		double c6=0.0000321767881768;
		double c7=0.0000002888167364;
		double c8=0.0000003960315187;



		double t=horner8(s,c0,c1,c2,c3,c4,c5,c6,c7,c8);
		if (x>0.5) {return t;}
		else {return -t;}
	}

	else {
		double r=1-x;
		double s=log(-log(r));

		double c0=0.3374754822726147;
		double c1=0.9761690190917186;
		double c2=0.1607979714918209;
		double c3=0.0276438810333863;
		double c4=0.0038405729373609;
		double c5=0.0003951896511919;
		double c6=0.0000321767881768;
		double c7=0.0000002888167364;
		double c8=0.0000003960315187;



		double t=horner8(s,c0,c1,c2,c3,c4,c5,c6,c7,c8);
		if (x>0.5) {return t;}
		else {return -t;}
	}
}


Matrix randuniform(int n){
	//returns a nx1 matrix, with entries being random samples
	//from the uniform distribution on [0,1]


    Matrix v(n,1);
    for (int i=0;i<n;i++){
        double x=double(rand())/RAND_MAX;
        v(i,0)=x;
    }
    
    return v;
}

double boxMuller(){
	//returns a random sample from the normal distribution N(0,1)

	//since the Box Muller method produces two normals, one is stored and returned
	//when the function is called for a second time


    double z1;
    static double z2;
    static bool has_stored=false; //this variable determines whether a normal is stored
    
    if (has_stored){
        has_stored=false;
        return z2;
    }
    
    Matrix u(2,1);
    double r2;
    
    while (true){
        u=randuniform(2);
        for(int i=0;i<2;i++){
            u(i,0)=2*u(i,0)-1; 
            //now u has two samples from the uniform distribution between [-1,1]
        }
        r2=u(0,0)*u(0,0)+u(1,0)*u(1,0);
        if (r2!=0 && r2<1){
            break;
            //only breaks out of the loop, if the sample is from the unit disk
        }
    }
    double s=sqrt(-2*log(r2)/r2);
    z1=s*u(0,0);
    z2=s*u(1,0);
    
    has_stored=true;
    return z1;
    
}

Matrix boxMuller(int N){
	//returns a matrix of size Nx1 with the entries being
	//random samples from the normal distribution N(0,1)

    Matrix v(N,1);
    for(int i=0;i<N;i++){
        v(i,0)=boxMuller();
    }
    return v;
}










