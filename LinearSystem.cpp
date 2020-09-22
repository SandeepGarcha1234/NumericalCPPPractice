#include <cassert>
#include "LinearSystem.hpp"
#include "matrix.hpp"


LinearSystem::LinearSystem(const Matrix& A, const Matrix& b):
mat(A), columnVector(b){
	assert(A.getRows()==b.getRows() && b.getCols()==1);
	//ensures that b is a column vector and
	//A and b are of compatible size

}

Matrix forwardSubstitution(const Matrix& L, const Matrix& b){
	//L should be a lower triangular matrix

	//solves the equation Lx = b by forward substitution

	int n = L.getCols();

	Matrix x(n,1); //declaration of "vector" that will store our result
	x(0,0) = b(0,0)/L(0,0);
	for (int i=1;i<n;i++){
		double sum =b(i,0);
		for (int j=0;j<i;j++){
			sum-= L(i,j)*x(j,0);
		}
		x(i,0)= sum/L(i,i);
	}

	return x;
}


Matrix backSubstitution(const Matrix &U, const Matrix& b){
	//U should be an upper traingular matrix

	//solves the equation Ux = b by back substitution

	int n = U.getCols();

	Matrix x(n,1);
	x(n-1,0) = b(n-1,0)/U(n-1,n-1);
	for (int i = n-2; i>=0; i--){
		double sum = b(i,0);
		for (int j=n-1;j>i;j--){
			sum-=U(i,j)*x(j,0);
		}
		x(i,0) = sum/U(i,i);
	}

	return x;
}

Matrix GaussianElimination::solve() const{
	//solves the equation Ax=b for x where A and b are members of the class LinearSystem

	//the equation is solved by first getting the PLU decomposition for the matrix A
	
	//then Ax=b is equivalent to Pb=PAx=LUx
	//below c=Pb and Ux=y
	//first we solve c=Ly for y by forward substitution
	//then we solve Ux=y for x by back subsitution



	int m = getMatrix().getRows();
	Matrix P(m,m);
	Matrix L(m,m);
	Matrix U(m,m);

	getMatrix().PLUdecomposition(P,L,U);
	Matrix c = P*getcolumnVector();

	Matrix y = forwardSubstitution(L,c);
	Matrix x = backSubstitution(U,y);

	return x;

}


Matrix Jacobi::nextVector(const Matrix& v0) const{
	//computes the next vector for the Jacobi iteration method for solving linear equations . . .
	//given the initial vector v0

	//it does this by decomposing A=L+D+U
	//where L is all the entries below the diagonal(excluding the diagonal)
	//D is the matrix of just diagonals of A
	//U is the matrix with the entries above the diagonal(excluding the diagonal)

	//if v1 is the next vector, then
	//D*v1 = b - (L+U)*v0


	int n = getMatrix().getCols();
	assert(v0.getRows() == n);

	Matrix T = getMatrix() - getMatrix().diagonal(); //this is the L+U matrix

	Matrix v1 = getcolumnVector() - (T*v0); //this calculates b-(L+U)*v0
	for (int i =0; i<n;i++){
		v1(i,0)/= getMatrix()(i,i); //solves for v1 by dividing by the D
	}

	return v1;
}

double norm(const Matrix& v){
	//calculates the l2 norm for a "column vector"

	assert (v.getCols()==1);
	int n = v.getRows();
	double sum = 0;
	for (int i =0 ;i<n;i++){
		sum+= v(i,0)*v(i,0);
	}

	return sqrt(sum);
}

Matrix Jacobi::solve(const Matrix& v0){
	//solves the linear equation by using Jacobi iteration
	//v0 is the initial guess


	Matrix previousVector = v0;
	int i =0;
	int N = getIterations();
	
	while (i<N-1){
		//the reason that this loop is ran only N-1 times is 
		//because for the Nth iteration we need to store both
		//of the last vectors in the sequence
		//in order to calculate the normOfError variable

		previousVector = nextVector(previousVector);
		i++;
	}

	Matrix finalVector = nextVector(previousVector);
	//by declaring the finalVector, we can store the result in the finalVector
	//and the 2nd last vector of the sequence in the variable previousVector

	normOfError = norm(finalVector - previousVector);

	return finalVector;

}

Matrix GaussSeidel::nextVector(const Matrix& v0) const{
	//computes the next vector for the Gauss Seidel iteration method for solving linear equations . . .
	//given the initial vector v0

	//it does this by decomposing A=L+D+U
	//where L is all the entries below the diagonal(excluding the diagonal)
	//D is the matrix of just diagonals of A
	//U is the matrix with the entries above the diagonal(excluding the diagonal)

	//if v1 is the next vector, then
	//(L+D)*v1 = b - U*v0
	//solves for v1 by using forward substitution since L+D is a lower diagonal matrix



	int n = getMatrix().getCols();
	assert (v0.getRows() == n);

	Matrix L = getMatrix().lower(); 
	// note that this matrix is L+D from the comment above
	// note that that lower returns the lower diagonal matrix (including the diagonal)

	Matrix U = getMatrix() - L;
	Matrix v1 = getcolumnVector() - U*v0; //this is the b-U*v0
	v1 = forwardSubstitution(L,v1);
	return v1;
}

Matrix GaussSeidel::solve(const Matrix& v0){
	//solves the linear equation by using Gauss Seidel iteration
	//v0 is the initial guess

	Matrix previousVector = v0;
	int i =0;
	int N = getIterations();
	while (i<N-1){
		//the reason that this loop is ran only N-1 times is 
		//because for the Nth iteration we need to store both
		//of the last vectors in the sequence
		//in order to calculate the normOfError variable

		previousVector = nextVector(previousVector);
		i++;
	}

	Matrix finalVector = nextVector(previousVector);
	//by declaring the finalVector, we can store the result in the finalVector
	//and the 2nd last vector of the sequence in the variable previousVector

	normOfError = norm(finalVector - previousVector);

	return finalVector;
}

Matrix SOR::nextVector(const Matrix& v0, double alpha) const{
	//computes the next vector for the SOR iteration method for solving linear equations . . .
	//given the initial vector v0

	//it does this by decomposing A=L+alpha*D - (alpha-1)*D+U
	//where alpha = 1/relaxtion_parameter
	//where L is all the entries below the diagonal(excluding the diagonal)
	//D is the matrix of just diagonals of A
	//U is the matrix with the entries above the diagonal(excluding the diagonal)

	//if v1 is the next vector, then
	//(L+alpha*D)*v1 = b + ((alpha-1)D-U)*v0
	//solves for v1 by using forward substitution since L+alpha*D is a lower diagonal matrix



	int n = getMatrix().getCols();
	assert (v0.getRows() == n);
	
	Matrix D = getMatrix().diagonal();
	Matrix L = getMatrix().lower()-D;
	Matrix U = getMatrix() - L - D;
	L+=alpha*D; //sets up the matrix L+alpha*D
	U = ((alpha-1)*D)-U;	//sets up the matrix (alpha-1)*D-U

	Matrix v1 = getcolumnVector() + U*v0; //this is b+((alpha-1)D-U)*v0
	v1 = forwardSubstitution(L,v1);
	return v1;
}

Matrix SOR::solve(const Matrix& v0, double relaxationParameter){
	//solves the linear equation by using SOR iteration
	//v0 is the initial guess

	Matrix previousVector = v0;
	int i =0;
	int N = getIterations();
	double alpha = 1.0/relaxationParameter;
	while (i<N-1){
		//the reason that this loop is ran only N-1 times is 
		//because for the Nth iteration we need to store both
		//of the last vectors in the sequence
		//in order to calculate the normOfError variable

		previousVector = nextVector(previousVector, alpha);
		i++;
	}

	Matrix finalVector = nextVector(previousVector,alpha);
	//by declaring the finalVector, we can store the result in the finalVector
	//and the 2nd last vector of the sequence in the variable previousVector


	normOfError = norm(finalVector - previousVector);

	return finalVector;
}


