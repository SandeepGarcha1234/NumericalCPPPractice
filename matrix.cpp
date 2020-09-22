#include <cassert>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include "matrix.hpp"



Matrix::Matrix(int rows, int cols):
nRows(rows),nCols(cols){
	start = new double* [rows];
	for (int i=0; i<rows; i++){
		start[i] = new double [cols];
	}
}

Matrix::~Matrix(){
	for(int i=0;i<nRows;i++){
		delete[] start[i];
	}
	delete[] start;
}

void Matrix::assign(const Matrix& mat){
	// This method is to help with the implementation of thecopy constructor and the assignment operator


	nRows = mat.getRows();
	nCols = mat.getCols();
	start = new double* [nRows];
	for (int i=0;i<nRows;i++){
		start[i] = new double [nCols];
	}
	for (int i =0; i<nRows; i++){
		for (int j=0; j<nCols; j++){
			start[i][j] = mat(i,j);
		}
	}
}

Matrix::Matrix(const std::string& s){
	//Function takes a string and creates a matrix
	//The string should begin with [
	//Entries should be entered row by row with each row entry being separated by ","
	//The last entry in the row should be followed by a ";" unless it is the last entry in the matrix
	// in which case it should be followed by ]

	//For example to create the matrix
	/*
	[1,2,3
	 4,5,6
	 7,8,9]
	*/

	// and give the variable name Mat, use the following command
	// Matrix Mat("[1,2,3;4,5,6;7,8,9]");

	int m = std::count(s.begin(),s.end(),';')+1;
	//Note that each row ends with ";" except for the last row and hence the plus one
	

	int n = std::count(s.begin(),s.end(),',')/m+1; 
	//Each entry in a row is going to have a "," after it except for the last entry
	//Hence to get the number of columns take total number of "," and divide by the number of columns and then add 1
	
	

	nRows = m;
	nCols = n;
	start = new double* [m];
	for (int i =0;i<m;i++){
		start[i]= new double [n];
	}
	
	auto pchar = s.begin()+1;
	// Beginning from the second character since the first character is going to be "["

	std::string t;
	double a;
	for (int i=0;i<m;i++){
		for (int j=0;j<n;j++){
			t=""; //this variable hold the digits of the next entry as a string
			while((*pchar)!=',' && (*pchar)!=';' && (*pchar)!=']'){
				t.push_back(*pchar);
				pchar++;
			}
			a = std::stod(t); //converting the string variable to double to insert into matrix 
			start[i][j] = a;
			pchar++; //This is to skip past the delimiter (, ; ])
		}
	}

}

Matrix::Matrix(const Matrix& mat){
	assign(mat);
}

Matrix::Matrix(const std::vector<double>& vec, bool column){
	//Creates nx1 or 1xn matrix which can be treated as a column or row vector respectively
	//if the boolean variable column is true (default), then it creates a nx1 matrix
	//if the boolean variable column is false, then it creates a 1xn matrix




	int n = int(vec.size());
	if (column){
		nRows = n;
		nCols = 1;
		start = new double* [n];
		for (int i =0; i<n; i++){
			start [i] = new double [1];
			start[i][0] = vec[i];
		}
	}

	else{
		nRows = 1;
		nCols = n; 
		start = new double* [1];
		start[0] = new double [n];
		for (int i =0; i<n; i++){
			start[0][i] = vec[i];
		}
	}
}

Matrix& Matrix::operator=(const Matrix& mat){
	
	for(int i =0;i<nRows;i++){
		delete[] start[i];
	}
	delete[] start;
	//The above steps are to make sure that Matrix variable being assigned to, is emptied before new values are put into it
	

	assign(mat);
	return *this;
	
}

double& Matrix::operator()(int i, int j) const{

	//Allows for access to matrix entries
	//Note the zero-based indexing is being used
	//To access the entry in the first row and first column of the matrix variable Mat,
	//use the following command
	//Mat(0,0);

	assert(0<=i && i<nRows && 0<=j && j<nCols); // to ensure that a valid entry is being accessed
	return start[i][j];
}

double& Matrix::operator()(int i, int j) {

	//Allows for access to matrix entries
	//Note the zero-based indexing is being used
	//To access the entry in the first row and first column of the matrix variable Mat,
	//use the following command
	//Mat(0,0);

	assert(0<=i && i<nRows && 0<=j && j<nCols);
	return start[i][j];
}

Matrix Matrix::operator+(const Matrix& mat) const{
	//overloads the + operator to allow for the addition of two matrices
	//to add two matrices named mat1 and mat2
	//use the following command
	//mat1 + mat2;


	assert(nRows==mat.getRows() && nCols==mat.getCols()); // ensure that the matrices are of the same size
	Matrix result(nRows,nCols);
	double** resultStart = result.getStart();
	for (int i=0; i<nRows; i++){
		for (int j=0; j<nCols; j++){
			resultStart[i][j] = start[i][j] + mat(i,j);
		}
	}

	return result;
}

Matrix Matrix::operator-(const Matrix& mat) const{

	//overloads the - operator to allow for the subtraction of two matrices
	//to subtract two matrices named mat1 and mat2
	//use the following command
	//mat1 - mat2;

	assert(nRows==mat.getRows() && nCols==mat.getCols()); // ensure that the matrices are of the same size
	Matrix result(nRows,nCols);
	double** resultStart = result.getStart();
	for (int i=0; i<nRows; i++){
		for (int j=0; j<nCols; j++){
			resultStart[i][j] = start[i][j] - mat(i,j);
		}
	}

	return result;
}

Matrix Matrix::operator*(const Matrix& mat) const{

	//overloads the * operator to allow for the multiplication of two matrices in the linear algebra sense
	//please note that it is not entry-wise multiplication
	//to multiply two matrices named mat1 and mat2
	//use the following command
	//mat1 * mat2;


	assert(nCols==mat.getRows()); 
	//ensure the number of columns of the left matrix equals the number of rows for the right matrix
	

	int n = mat.getCols();
	Matrix result(nRows,n);

	for (int i=0; i<nRows; i++){
		for (int j=0; j<n; j++){
			result(i,j)=0;
			for (int k =0; k<nCols;k++){
				result(i,j)+=((*this)(i,k)*mat(k,j));
			}
		}
	}

	return result;
}

Matrix Matrix::operator*(double a) const{

	//allows to muliply a matrix by a scalar

	Matrix result(nRows,nCols);
	for(int i =0;i<nRows;i++){
		for(int j=0;j<nCols;j++){
			result(i,j) = (*this)(i,j)*a;
		}
	}

	return result;
}

Matrix Matrix::operator*(std::vector<double>& colVec) const{
	//allows right multiplication of a vector to matrix

	Matrix colMat(colVec, true); 
	//converts vector variable to variable of type matrix of size nx1 since it is right multiplication
	//then uses the operloaded * operator to return result
	//no need to check size compatability since this was done in the previous overloading of *

	return (*this)*colMat;
}

Matrix operator*(const std::vector<double>& vec, const Matrix& mat){
	//allows left multiplication of a vector to matrix

	Matrix rowMat(vec, false);
	//converts vector variable to variable of type matrix of size 1xn since it is left multiplication
	//then uses the operloaded * operator to return result
	//no need to check size compatability since this was done in the previous overloading of *

	return rowMat*mat;
}

Matrix& Matrix::operator+=(const Matrix& mat){

	//allows for increment of matrix variable by another matrix variable as long as they are of the same size
	//does this by using the overloaded + operator

	Matrix result = (*this).operator+(mat);
	*this = result;
	return *this;
}

Matrix& Matrix::operator-=(const Matrix& mat){

	//allows for decrement of matrix variable by another matrix variable as long as they are of the same size
	//does this by using the overloaded - operator

	Matrix result = (*this).operator-(mat);
	*this = result;
	return *this;
}

Matrix& Matrix::operator*=(const Matrix& mat){

	//allows to update variable by multiplying it by another matrix as long as the size is compatible
	//for instance, if A,B,C are three variables of type Matrix, then the following
	/*
	C=A*B;
	A=C;
	*/

	//would be equivalent to 

	/*
	A*=B;
	*/

	//in terms of the values held in the Matrix A

	Matrix result = (*this)*mat;
	*this = result;
	return *this;
}

Matrix& Matrix::operator*=(double a){

	//updates matrix by multiplying all of the entries by the value held in variable a

	Matrix result = (*this).operator*(a);
	*this = result;
	return *this;
}


int Matrix::pivotIndex(int i) const{

	//This function finds the row index of the entry with the largest absolute at or below entry (i,i) in column i
	//This function is for the PLUdecomposition function for partial pivoting
	//Note that all indicies are zero-based
	//For example, suppose the Matrix variable Mat, stores the following matrix
	/*
	[1,5,4
	2,-9,2
	-7,8,1]
	*/

	//given the following code

	/*
	
	double x1 = Mat.pivotIndex(0);
	double x2 = Mat.pivotIndex(1);
	double x3 = Mat.pivotIndex(2);

	*/

	//the variable x1 would be 2
	//the variable x2 would be 1
	//the variable x3 would be 2


	int n=getRows();
	assert(i>=0 && i<n);
	int k=i;
	double max = fabs((*this)(i,i));
	for (int j=i+1; j<n; j++){
		double entry = (*this)(j,i);
		if (fabs(entry) > max){
			max = fabs(entry);
			k=j;
		}
	}

	return k;
}



void Matrix::PLUdecomposition(Matrix& P, Matrix& L, Matrix& U) const{
	
	//provide the results for PLU decomposition for square matrices ONLY 
	//result is stored in variable passed by reference
	//at this moment, this function assumes the matrix is invertible
	//this will be corrected in the future

	//matrices passed by reference to store the results, do not have to be of the correct size
	//the size is corrected in the function
	
	//the function uses partial pivoting (ensuring that the pivot entry is the largest from that row)
	//the function used is pivotIndex(int)  ----   defined above

	//to increase numerical accuracy


	int m = getRows();
	int n = getCols();
	assert(m==n); // ensure that the matrix is square
	U=*this;
	L=identity(n);
	P=identity(n);


	double** pStart = P.getStart();
	double** uStart = U.getStart();
	double* ptemp; // pointer to allow for row exchanges for the U and P matrix during pivot
	double temp;
	for (int i=0; i<n; i++){
		int k = U.pivotIndex(i); // find the largest pivot possible
		if (i!=k){
			// if the largest pivot is not on the diagonal
			
			ptemp = pStart[i];
			pStart[i] = pStart[k];
			pStart[k] = ptemp;

			// This switches rows i and k for the matrix P


			ptemp = uStart[i];
			uStart[i] = uStart[k];
			uStart[k] = ptemp;

			//This switches row i and k for the matrix U

			for (int j=0; j<i;j++){
				temp = L(i,j);
				L(i,j) = L(k,j);
				L(k,j) = temp;
			}

			//Note for matrix L, we do not switch entire rows but only the entries of rows i and k below the diagonal

		}

		//below modifies row j by subtracting the coeff * row i from it
		//the coeff is entered as an entry to the L matrix

		for (int j = i+1; j<n; j++){
			double coeff = U(j,i)/U(i,i);
			L(j,i) = coeff;
			U(j,i) =0;
			for(int l=i+1; l<n;l++){
				U(j,l)-= coeff*U(i,l);
			}
		}

	}

	
}



double parityOfPermutationMatrix(const Matrix& mat){

	//calculates the determinant of a permutation matrix
	//it does this by decomposing the permutation into a product of cycles
	//then for each cycle it calculates it parity
	//if the cycle is of even length, then its parity is -1
	//if the cycle is of odd length, then its parity is 1
	//this parity is multiplied to the overall parity



	int n = mat.getRows();
	int* indices = new int [n];
	//the indices variable that will keep track of the cycles
	//it starts with the 0 index, and then adds the index to which 0 is mapped to
	//then that index is mapped to another index which is then added to this array as well
	//as soon as the cycle is complete, we find another index which has not been included in any previous cycle
	//then we calculate the cycle for this index


	indices[0]=0; //we start by getting the cycle involving the 0 index
	
	int counter = 1;
	//in the counter variable, we count the number of indices that have been accounted for in any cycle
	//once the counter is n, that means that each index has been included in one of the cycles

	int parity = 1;
	//the parity variable keeps track of the overall determinant
	//the parity of each cycle will be multiplied to this variable

	int beginningIndex = 0;
	//the beginning variable keeps track of the position where the beginning index of 
	//the current cycle resides in the indicies array



	while (counter<n){
		//note that once all the indices have been
		//ie counter=n
		//then this while loop will not execute and
		//hence the parity of the final cycle will not be accounted for
		//we account for this final cycle below

		for(int j=0;j<n;j++){
			if(mat(indices[counter-1],j)==1){
				//here we are locating the next index of the cycle
				if(j!=indices[beginningIndex]){
					//if the next index is not currently in the cycle 
					//ie(is not equal to the index at the front of the cycle)
					//we add this index to indices array
					indices[counter] = j;
					
					counter++;
					//since we have added another index, we have one less index to account for
					//thus we increment the counter variable by 1
					break;
				}

				else{
					//this is the case where the cycle has been completed
					//ie this index is mapped to the index at the front of the cycle
					if((counter - beginningIndex)%2==0){
						//here the cycle is of odd length

						//to see this, suppose that beginning index is one and the end index is three
						//then the cycle has three element (at index 1, 2 and 3)
						//and 3-1=2 which is even
						parity*=-1;

						//below we find the index that has not yet been included in any of the cycles
						//by comparing them to the indices in the indices array

						int k=1;
						//we start with index 1 since we know that the first cycle started with index 0
						//if this index has been in a previous cycle, then we increment k
						//until we find an index that has not yet been included in a cycle



						while(k<n){
							bool match = false;
							//the match variable tells us whether we have found our index k
							//is already in one of the cycles

							for(int l=0;l<counter;l++){
								if(k==indices[l]){
									//this is the case where the index was already in a previos cycle
									match = true;
									k++; //by updating k, we are going to check if the next index has been accounted for
									break;
								}

							}

							if(!match){
								//this is the case where the index is not in any previous cycle

								indices[counter]=k;
								beginningIndex = counter;
								//we change the beginning index as k will be the beginning of our new cycle

								counter++;
								break;
							}
						}
					}

					else{
						//this is the case the cycle is of even length
						//since the parity of the cycle is 1 
						//we do not need update the overall parity variable

						//below we find the index that has not yet been included in any of the cycles
						//we use the same method as above

						int k=1;
						while(k<n){
							bool match = false;
							for(int l=0;l<counter;l++){
								if(k==indices[l]){
									match = true;
									k++;
									break;
								}

							}

							if(!match){
								indices[counter]=k;
								beginningIndex = counter;
								counter++;
								break;
							}
						}

					}
				}
			}
		}
	}

	if (counter==n){
		//here we calculate the parity of the final cycle
        if((counter- beginningIndex)%2==0){
            parity*=-1;
        }
        
    }

	delete[] indices;
	return parity;

}

double Matrix::determinant() const{
	//calculates the determinant of a matrix by first calculating the PLU decomposition
	//the result is given by multiplying the determinant of U by the determinant of P
	//Note that the determinant of L is 1


	int n = getRows();

	Matrix P(n,n);
	Matrix L(n,n);
	Matrix U(n,n);

	(*this).PLUdecomposition(P,L,U);

	//the code below calculates the determinant of U
	//by multplying the diagonals

	double result=1;
	for (int i =0; i<n;i++){
		result*= U(i,i);
	}

	result *= parityOfPermutationMatrix(P);
	return result;
}

double Matrix::trace() const{
	//calculates the trace of a matrix


	int m = getRows();
	int n = getCols();
	assert(m==n);

	double result = 0;

	for (int i=0; i<n; i++){
		result+=(*this)(i,i);
	}

	return result;
}

Matrix Matrix::transpose() const{
	//returns the transpose of a matrix

	int m = getRows();
	int n = getCols();

	Matrix result(n,m);

	for (int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			result(j,i)=(*this)(i,j);
		}
	}

	return result;
}

Matrix operator*(double a, const Matrix& mat){
	return mat*a;
}

std::ostream& operator<<(std::ostream& out, const Matrix& mat){
	//inserts the matrix into the output stream
	//for example, if we had Matrix mat holding the following matrix
	/*
	[1,2,3
	 4,5,6
	 7,8,9]
	*/

	//then the output to the stream would be
	//[1,2,3;4,5,6;7,8,9]

	int m = mat.getRows();
	int n = mat.getCols();
	out<<'[';
	for(int i =0;i<m;i++){
		for (int j=0;j<n-1;j++){
			out<<mat(i,j)<<',';
			//we are outputting the first n-1 entries for the row 
			//followed by the a comma (,)
		}

		out<<mat(i,n-1);
		if(i==m-1){out<<']';}//this is when we have entered the last entry of the matrix
		else{out<<';';}//the case when it is the last entry of the row but not the matrix
	}
	return out;
}

Matrix identity(int n){
	//returns the nxn identity matrix

	Matrix I(n,n);
	double** start = I.getStart();
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			if(i==j){
				start[i][j]=1;
			}
			else {start[i][j]=0;}
		}
	}

	return I;
}

Matrix Matrix::diagonal() const{
	//returns a matrix of the same size with the same diagonals and
	//zeros on the off-diagonals

	int m = getRows();
	int n = getCols();

	assert(m==n);
	Matrix D(m,n);
	for (int i =0; i<m; i++){
		for (int j =0; j<n;j++){
			if(i!=j){D(i,j)=0;}
			else {D(i,j) = (*this)(i,j);}
		}
	}

	return D;
}

Matrix Matrix::lower() const{
	//returns a matrix with the same entry on the lower diagonal
	//including the diagonal
	//has zeros on the upper diagonal


	int m = getRows();
	int n = getCols();

	assert(m==n);

	Matrix L(m,n);
	for (int i =0;i<m;i++){
		for (int j =0; j<n; j++){
			if(i>=j){L(i,j) = (*this)(i,j);}
			else{L(i,j)=0;}
		}
	}

	return L;
}

Matrix Matrix::row(int i) const{
	//returns the ith row from the matrix in the form of a 1xn matrix
	//note that because of zero-based indexing
	//when using the ith index,
	//the i+1 row is returned

	int m = getRows();
	int n = getCols();

	Matrix rowVector(1,n);
	for (int j=0; j<n; j++) {
		rowVector(0,j) = (*this)(i,j);
	}

	return rowVector;
}

Matrix Matrix::column(int j) const{
	//returns the jth column from the matrix in the form of a nx1 matrix
	//note that because of zero-based indexing
	//when using the jth index,
	//the j+1 column is returned

	int m = getRows();
	int n = getCols();

	Matrix columnVector(m,1);
	for (int i =0; i<m; i++){
		columnVector(i,0) = (*this)(i,j);
	}

	return columnVector;
}

void Matrix::insertRow(const std::vector<double>& v, int i){
	//inserts the vector into the ith indexed row
	//because of zero-based indexing
	//it is inserted into the i+1 row

	int m = getRows();
	int n = getCols();

	assert(i>=0 && i<m && v.size()==n);
	for (int j=0;j<n;j++){
		(*this)(i,j) = v[j];
	}

}

void Matrix::insertColumn(const std::vector<double>& v, int j){
	//inserts the vector into the jth indexed column
	//because of zero-based indexing
	//it is inserted into the j+1 column

	int m = getRows();
	int n = getCols();

	assert(j>=0 && j<n && v.size()==m);
	for (int i =0;i<m;i++){
		(*this)(i,j) = v[i];
	}
}

void Matrix::insertRow(const Matrix& v, int i){
	//inserts the matrix v into the ith indexed row
	//v can be of dimension 1xn or nx1
	//because of zero-based indexing
	//it is inserted into the i+1 row


	int m = getRows();
	int n = getCols();

	assert(i>=0 && i<m);

	int p = v.getRows();
	int q = v.getCols();

	assert(p==1 || q==1);
	//ensures that either the row size or column size is 1

	if (p==1){
		//this is the case where v is a row vector
		assert(q==n);
		for (int j=0;j<n;j++){
			(*this)(i,j) = v(0,j);
		}

	}

	else{
		//this is the case where v is a column vector
		assert(p==n);
		for(int j=0;j<n;j++){
			(*this)(i,j) = v(j,0);
		}
	}
}

void Matrix::insertColumn(const Matrix& v, int j){
	//inserts the matrix v into the jth indexed column
	//v can be of dimension 1xn or nx1
	//because of zero-based indexing
	//it is inserted into the j+1 column

	int m = getRows();
	int n = getCols();

	assert(j>=0 && j<n);
	//ensures that either the row size or column size is 1

	int p = v.getRows();
	int q = v.getCols();

	assert(p==1 || q==1);

	if(p==1){
		//this is the case where v is a row vector
		assert(q==m);
		for (int i=0;i<m;i++){
			(*this)(i,j) = v(0,i);
		}
	}

	else{
		//this is the case where v is a column vector
		assert(p==m);
		for (int i=0;i<m;i++){
			(*this)(i,j) = v(i,0);
		}
	}

}





