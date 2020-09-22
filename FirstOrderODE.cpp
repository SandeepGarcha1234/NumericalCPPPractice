#include <cassert>
#include "FirstOrderODE"

FirstOrderODE::FirstOrderODE(const std::vector<double (*) (double t, Matrix v)>& FF, double t0, const Matrix& y0) :
F(FF), initialTime(t0), initialValue(y0){
	assert(FF.size() == y0.getRows());
	order = y0.getRows();
}

Matrix evaluate(const std::vector<double (*) (double t, Matrix u)>& F, double t, const Matrix& v){
	//evaluates F(t,v) and returns the result in a matrix of size nx1 (column vector form)

	int n = v.getRows();
	assert(F.size() == n);

	Matrix result(n,1);

	for (int i=0; i<n; i++){
		double y = (*(F[i]))(t,v);
		result(i,0) = y;
	}

	return result;
}


Matrix EulerMethod::solve(const FirstOrderODE& A){
	//Uses Euler Method to solve the ODE
	//the value at time step k+1 ie. (v(k+1)) is as follows

	//v(k+1)=v(k)+h*F(tk,v(k))

	int m = A.getorder();
	int n = getnumberOfSteps();

	Matrix result(m,n);
	
	double s = A.getinitialTime(); //this variable represent the current time step
	double h = (getfinalTime() - s)/n; //this is the size of the time step

	Matrix currentVector = A.getinitialValue();


	for (int j=0; j<n; j++){
		

		currentVector = currentVector + (h* evaluate(A.getRHS(),s,currentVector));
		s+=h;//increments to the next time step
		result.insertColumn(currentVector,j); //insert v(j) into the jth index column
	}

	return result;

}


Matrix ImprovedEulerMethod::solve(const FirstOrderODE& A){
	//Uses the Improved Euler Method to solve the ODE
	//the value at time step k+1 ie. (v(k+1)) is as follows

	//v(k+1)=v(k)+0.5*h*[F(tk,v(k)) + F(t+h, v(k)+h*F(tk,v(k)))]
	int m = A.getorder();
	int n = getnumberOfSteps();

	Matrix result(m,n);
	
	double s = A.getinitialTime();
	double h = (getfinalTime() - s)/n; //this is the size of the time step

	Matrix currentVector = A.getinitialValue();

	for (int j=0; j<n; j++){

		Matrix evaluationVector = evaluate(A.getRHS(),s,currentVector); //h*F(tk,v(k))
		Matrix eulerVector = currentVector + (h* evaluationVector); //v(k) + F(tk,v(k))
		currentVector = currentVector + 0.5*h * (evaluationVector + evaluate(A.getRHS(),s+h,eulerVector)); //v(k+1)

		s+=h;
		result.insertColumn(currentVector,j); 

	}

	return result;
}

Matrix MidpointMethod::solve(const FirstOrderODE& A){
	//Uses the Midpoint Method to solve the ODE
	//the value at time step k+1 ie. (v(k+1)) is as follows

	//v(k+1)=v(k)+h*F(tk+0.5*h, v(k)+0.5*h*F(tk,v(k))) 
	int m = A.getorder();
	int n = getnumberOfSteps();

	Matrix result(m,n);
	
	double s = A.getinitialTime();
	double h = (getfinalTime() - s)/n;

	Matrix currentVector = A.getinitialValue();

	for (int j=0; j<n; j++){
		Matrix eulerVector = currentVector + 0.5*h*evaluate(A.getRHS(),s,currentVector);//v(k)+0.5*h*F(tk,v(k))
		currentVector = currentVector + h*evaluate(A.getRHS(),s+0.5*h,eulerVector);

		s+=h;
		result.insertColumn(currentVector,j);
	}

	return result;
}

Matrix RK4::solve(const FirstOrderODE& A){
	//Uses the Runge-Kutta Method of order 4 to solve the ODE
	//the value at time step k+1 ie. (v(k+1)) is as follows

	//v(k+1)=v(k)+(h/6.0)*[F(tk,v(k))+2F(t2k,v2k)+2F(t3k,v3k)+F(t4k,v4k)] where
	//t2k = tk+0.5h		v2k = u(k) + 0.5h*F(tk,vk)
	//t3k = t2k			v3k = u(k) + 0.5h*F(t2k,v2k)
	//t4k = tk+h		v4k = u(k) + h*F(t3k,v3k)

	int m = A.getorder();
	int n = getnumberOfSteps();

	Matrix result(m,n);
	
	double s = A.getinitialTime();
	double h = (getfinalTime() - s)/n;

	Matrix currentVector = A.getinitialValue();

	for (int j=0;j<n;j++){

		double s2 = s+ 0.5*h; //t2k
		double s4 = s+ h;	//t4k

		Matrix e1 = evaluate(A.getRHS(),s,currentVector); //F(tk,vk)
		Matrix u2 = currentVector + 0.5*h*e1; //v2k

		Matrix e2 = evaluate(A.getRHS(),s2,u2); //F(t2k,v2k)
		Matrix u3 = currentVector + 0.5*h*e2; //v3k

		Matrix e3 = evaluate(A.getRHS(),s2,u3); //F(t3k,v2k)
		Matrix u4 = currentVector + h*e3; //v4k

		currentVector = currentVector + (h/6.0) * (e1+2*e2+2*e3+evaluate(A.getRHS(),s4,u4)); //v(k)+(h/6.0)*[F(tk,v(k))+2F(t2k,v2k)+2F(t3k,v3k)+F(t4k,v4k)]
		s=s4; //the next time step should be t+h which is t4k

		result.insertColumn(currentVector,j);

	}

	return result;


}



