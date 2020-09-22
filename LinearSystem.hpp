#ifndef LINEARSYSTEM.HPP
#define LINEARSYSTEM.HPP

#include <cassert>
#include "matrix.hpp"


class LinearSystem{
//sets up the equation Ax=b

private:
	Matrix mat; //the matrix A
	Matrix columnVector; //the "vector" b

public:
	LinearSystem(const Matrix& A, const Matrix& b); //b must a Matrix of size nx1
	virtual ~LinearSystem(){};
	Matrix getMatrix() const {return mat;}
	Matrix getcolumnVector() const {return columnVector;}
	virtual Matrix solve(){};

};

class GaussianElimination : public LinearSystem {
public:
	GaussianElimination(const Matrix& A, const Matrix& b) : LinearSystem(A,b) {}
	Matrix solve() const;

};

class Iterative : public LinearSystem{
private:
	int numberOfIterations;
public:
	Iterative(const Matrix& A, const Matrix& b, int N = 40) : LinearSystem(A,b), numberOfIterations(N){}
	virtual ~Iterative(){};

	void setIterations(int N) {numberOfIterations = N;}
	int getIterations() const {return numberOfIterations;}
	virtual Matrix solve(){};
};

class Jacobi : public Iterative{
//solves a linear equation by using the Jacobi iteartion method
//Although the class does not check for convergence of the method,
//it does provide the l2 norm of the last two vectors in the sequence of the iteration
//This result is stored in the normOfError variable

private:
	double normOfError;
	Matrix nextVector(const Matrix& v0) const;
public:
	Jacobi(const Matrix& A, const Matrix& b, int N = 40) : Iterative(A,b,N) {normOfError= -1;}
	//assigning normOfError to -1 is to assure that this variable is only accessed
	//after the iteration method has been ran.
	//See the getnormOfError() method

	double getnormOfError() const {
		assert(normOfError!= -1); //ensures that the iteration has been run
		return normOfError;
	}

	Matrix solve(const Matrix& v0); //v0 is the initial guess

};
	
class GaussSeidel : public Iterative{
//solves a linear equation by using the Gauss Seidel iteartion method
//Although the class does not check for convergence of the method,
//it does provide the l2 norm of the last two vectors in the sequence of the iteration
//This result is stored in the normOfError variable

private:
	double normOfError;
	Matrix nextVector(const Matrix& v0) const;
public:
	GaussSeidel(const Matrix& A, const Matrix& b, int N = 40) : Iterative(A,b,N) {normOfError= -1;}
	//assigning normOfError to -1 is to assure that this variable is only accessed
	//after the iteration method has been ran.
	//See the getnormOfError() method

	double getnormOfError() const {
		assert(normOfError!= -1); //ensures that the iteration has been run
		return normOfError;
	}

	Matrix solve(const Matrix& v0); //v0 is the initial guess
};

class SOR : public Iterative{
//solves a linear equation by using the SOR iteartion method
//Although the class does not check for convergence of the method,
//it does provide the l2 norm of the last two vectors in the sequence of the iteration
//This result is stored in the normOfError variable

private:
	double normOfError;
	Matrix nextVector(const Matrix& v0, double relationParameter) const;
public:
	SOR(const Matrix& A, const Matrix& b, int N = 40) : Iterative(A,b,N) {normOfError= -1;}
	//assigning normOfError to -1 is to assure that this variable is only accessed
	//after the iteration method has been ran.
	//See the getnormOfError() method

	double getnormOfError() const {
		assert(normOfError!= -1); //ensures that the iteration has been run
		return normOfError;
	}

	Matrix solve(const Matrix& v0, double relationParameter); //v0 is the initial guess
};



#endif