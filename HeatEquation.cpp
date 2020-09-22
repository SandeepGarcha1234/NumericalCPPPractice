#include "HeatEquation.hpp"
#include "matrix.hpp"

Matrix ExplicitHeatSolver::solve(const HeatEquation& A){
	double dt = getdt();
	double dx = getdx();

	double mu = A.getdiffusivity() * dt / (dx*dx);
	assert (mu<= 0.5);

	double nu = 1-2*mu;
	
	double leftPoint = A.getleftX();
	double rightPoint = A.getrightX();

	int M = int((rightPoint - leftPoint)/dx) - 1;
	int N = int(getfinalTime()/dt);
	
	Matrix result(M,N);

	double l = leftPoint + dx;
	double c = l+dx;
	double r =c+dx;

	for (int i=0; i<M; i++){
		
		if (i==0){
			double leftValue = (*(A.getleftBoundary()))(0);
			result(i,0) = mu*(leftValue + (*(A.getinitialCondition()))(leftPoint + 2*dx)) + nu*(*(A.getinitialCondition()))(leftPoint+dx);
		}
		else if (i==M-1){
			double rightValue = (*(A.getrightBoundary()))(0);
			result(i,0) = mu *(rightValue + (*(A.getinitialCondition()))(rightPoint - 2*dx)) + nu*(*(A.getinitialCondition()))(rightPoint-dx);
		}

		else {
			result(i,0) = mu*(*(A.getinitialCondition()))(l) + (*(A.getinitialCondition()))(r) + nu*(*(A.getinitialCondition()))(c);
			l+=dx;
			c+=dx;
			r+=dx;
		}

	}

	int timeStep =1;

	while (timeStep<N){
		double left = leftPoint + dx;
		double center = l+dx;
		double right =c+dx;

		for (int i=0; i<M;i++){
			if (i==0){
				double leftValue = (*(A.getleftBoundary()))(timeStep * dt);
				result(i,timeStep) = mu*(leftValue + result(i+1,timeStep-1)) + nu*result(i,timeStep-1);
			}

			else if(i==M-1){
				double rightValue = (*(A.getrightBoundary()))(timeStep * dt);
				result(i,timeStep) = mu *(rightValue + result(i-1,timeStep-1)) + nu*result(i,timeStep-1);
			}

			else{
				result(i,timeStep) = mu*(result(i-1,timeStep-1) + result(i+1,timeStep-1)) + nu*(result(i,timeStep-1));
				left+=dx;
				center+=dx;
				right+=dx;
			}
		}
	
		timeStep++;
	
	}

	return result;

}

Matrix HeatTriagonalSolver(double diag, double offdiag, const Matrix& vec){
	assert(vec.getCols()==1);
	int n = vec.getRows();

	Matrix c(n,1);
	Matrix d(n,1);

	Matrix result(n,1);

	c(0,0) = vec(0,0)/offdiag;
	d(0,0) = -diag/offdiag;

	c(1,0) = vec(1,0)/offdiag - diag*c(0,0)/offdiag;
	d(1,0) = -diag*d(0,0)/offdiag - 1;

	for (int i=2;i<n-1;i++){
		c(i,0) = (vec(i,0) - offdiag*c(i-2,0) - diag*c(i-1,0))/offdiag;
		d(i,0) = -(d(i-2,0) + diag* d(i-1,0)/offdiag);
	}

	c(n-1,0) = (vec(n-1,0) - offdiag * c(n-3,0))/diag;
	d(n-1,0) = -offdiag*d(n-3,0)/diag;

	result(0,0) = (c(n-2,0)-c(n-1,0))/(d(n-1,0)-d(n-2,0));
	for (int i=1; i<n; i++){
		result(i,0) = c(i-1,0) + d(i-1,0)*result(0,0);
	}

	return result;

}

Matrix ImplicitHeatSolver::solve(const HeatEquation& A){
	double dt = getdt();
	double dx = getdx();

	double mu = A.getdiffusivity() * dt / (dx*dx);

	double nu = 1+2*mu;
	
	double leftPoint = A.getleftX();
	double rightPoint = A.getrightX();

	int M = int((rightPoint - leftPoint)/dx) - 1;
	int N = int(getfinalTime()/dt);
	
	Matrix result(M,N);

	Matrix u0(M,0);

	double l = leftPoint + dx;

	for (int i=0;i<M;i++){
		u0(i,0) = (*(A.getinitialCondition()))(l);
		l+=dx;
	}

	u0(0,0) = mu * (*(A.getleftBoundary()))(0);
	u0(M-1,0) = mu * (*(A.getrightBoundary()))(0);

	Matrix u1 = HeatTriagonalSolver(nu,-mu,u0);

	result.insertColumn(u1,0);

	int timeStep=1;
	while (timeStep<N){
		Matrix v = result.column(timeStep-1);
		v(0,0) += (*(A.getleftBoundary()))(timeStep*dt);
		v(M-1,0) += (*(A.getrightBoundary()))(timeStep*dt);

		Matrix w = HeatTriagonalSolver(nu,-mu,v);
		result.insertColumn(w,timeStep);

		timeStep++;
	}

	return result;

}

Matrix HeatTriagonalMultiplier(double diag, double offdiag, const Matrix& vec){
	assert(vec.getCols()==1);
	int n = vec.getRows();

	Matrix result(n,1);
	result(0,0) = diag * vec(0,0) + offdiag * vec(1,0);
	for (int i=1; i<n-1; i++){
		result(i,0) = offdiag*(vec(i-1,0) + vec(i+1,0)) + diag*vec(i,0);
	}
	result(n-1,0) = offdiag*(vec(n-2,0)) + diag*vec(n-1,0);

	return result;
}

Matrix CrankNicolsonHeatSolver::solve(const HeatEquation& A){
	double dt = getdt();
	double dx = getdx();

	double mu = A.getdiffusivity() * dt / (dx*dx);
	
	double leftPoint = A.getleftX();
	double rightPoint = A.getrightX();

	int M = int((rightPoint - leftPoint)/dx) - 1;
	int N = int(getfinalTime()/dt);
	
	Matrix result(M,N);

	Matrix u0(M,0);

	double l = leftPoint + dx;

	for (int i=0;i<M;i++){
		u0(i,0) = (*(A.getinitialCondition()))(l);
		l+=dx;
	}

	u0 = HeatTriagonalMultiplier(1-mu,0.5*mu,u0);

	u0(0,0) = 0.5 * mu * ((*(A.getleftBoundary()))(0) + (*(A.getleftBoundary()))(dt));
	u0(M-1,0) = 0.5 * mu * ((*(A.getrightBoundary()))(0) + (*(A.getrightBoundary()))(dt));

	Matrix u1 = HeatTriagonalSolver(1+mu,-0.5*mu,u0);
	result.insertColumn(u1,0);

	int timeStep = 1;
	while (timeStep < N){
		Matrix v0 = result.column(timeStep-1);
		v0 = HeatTriagonalMultiplier(1-mu,0.5*mu,v0);
		v0(0,0) = 0.5 * mu * ((*(A.getleftBoundary()))((timeStep-1)*dt) + (*(A.getleftBoundary()))(timeStep*dt));
		v0(M-1,0) = 0.5 * mu * ((*(A.getrightBoundary()))((timeStep-1)*dt) + (*(A.getrightBoundary()))(timeStep*dt));
		v0 = HeatTriagonalSolver(1+mu,-0.5*mu,v0);

		result.insertColumn(v0,timeStep);
		timeStep++;
	}

	return result;
}


