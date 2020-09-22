#ifndef MATRIX.HPP
#define MATRIX.HPP

#include <iostream>
#include <string>
#include <vector>

class Matrix{

private:
	int nRows;
	int nCols;
	double** start;
	void assign(const Matrix& mat);

public:
	Matrix(int rows, int cols); //declares a rows x cols matrix
	virtual ~Matrix();
	Matrix(const Matrix& mat);
	Matrix(const std::string& s);
	Matrix(const std::vector<double>& vec, bool column=true);

	//Methods to access private members
	int getRows() const {return nRows;}
	int getCols() const {return nCols;}
	double** getStart() const {return start;}

	//Operator overloading

	Matrix& operator=(const Matrix& mat);
	Matrix& operator+=(const Matrix& mat);
	Matrix& operator-=(const Matrix& mat);
	Matrix& operator*=(const Matrix& mat);
	Matrix& operator*=(double a);

	double& operator()(int i, int j) const;
	double& operator()(int i, int j);

	Matrix operator+(const Matrix& mat) const;
	Matrix operator-(const Matrix& mat) const;
	Matrix operator*(const Matrix& mat) const;
	Matrix operator*(double a) const;
	Matrix operator*(std::vector<double>& colVec) const;

	int pivotIndex(int i) const; //finds the row index of the greatest absolute value below entry (i,i) for a square matrix
	void PLUdecomposition(Matrix& P, Matrix& L, Matrix& U) const;
	double determinant() const;
	double trace() const;
	Matrix transpose() const;
	Matrix diagonal() const;
	Matrix lower() const;

	Matrix row(int i) const;
	Matrix column(int j) const;

	void insertRow(const std::vector<double>& v, int i);
	void insertColumn(const std::vector<double>& v, int j);

	void insertRow(const Matrix& v, int i);
	void insertColumn(const Matrix& v, int j);

};

Matrix operator*(double a, const Matrix& mat);
std::ostream& operator<<(std::ostream& out, const Matrix& mat);
Matrix identity(int n);




#endif
