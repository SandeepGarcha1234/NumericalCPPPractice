# Numerical CPP Practice
I have recently started to learn C++ and the purpose of this code is for me to gain practice while working on some projects relevant to quantitative finance. This project is by no means complete and I will continue to add to this project. 

## The Matrix Class

This class stores two dimensional arrays. There are some linear algebra functionalities, however more need to be added. It should be noted that this class uses zero-based indexing. In addition, this project does not have a vector class and any computations involving vectors are done with the Matrix class, using matrices of size Nx1 or 1XN depending on the context.

### Constructing A Matrix

In order to create, an empty matrix name M of size mxn, just use the command Matrix M(m,n). For example, if the M should have 3 rows and 2 columns, then use the following command:

```
Matrix M(3,2);
```

If the entries of the matrix are already known, then the matrix can be created from a string. For example, to create a matrix M with the following entries

[ 1, 2

   3, 4]
  
use the following command:

```
Matrix M("[1,2;3,4]";

```

Please note that each entry in a specific row is separated by a `,` and the end of a row is followed a `;` (expect the last row. The begininning of the matrix begins with a `[` and ends with a `]`. There are no spaces between any of the characters.

A matrix can also be created from a standard C++ vector type. The following code shows how to do this.

```
std::vector<double> v({1,2});
Matrix M(v,true); //will be a 2x1 matrix
Matrix N(v,false); //will be a 1x2 matrix

```
Please note that the second argument specifies if the matrix is in column vector form(true) or if the matrix is in row vector form(false). This argument is true by default.

The identity matrix of size n can be created by using the identity function. For instance, the following code creates a 2x2 identity matrix and stores it in I.

`Matrix I = identity(2);`

Please note the identity function exists outside of the class.

### Methods on Matrices

The methods getRows() and getCols() return the number of rows and columns, respectively. For instance, consider the following code:

```
Matrix M("[1,2;3,4;5,6]");
int m = M.getRows();
int n = M.getCols();
```

Here m will take the value 3 and n will take the value 2.

The PLUdecomposition returns the PLU decomposition for square invertible matrices(this will be expanded in the future). The result is returned in matrices that are passed by as arguments and they must be passed in the correct order (permutation matrix first, then the lower diagonal matrix, and finally the upper diagonal matrix. For example, suppose we have a matrix under the variable A as follows

[1&nbsp;&nbsp;&nbsp;&nbsp;2&nbsp;&nbsp;&nbsp;&nbsp;-1&nbsp;&nbsp;&nbsp;&nbsp;0

 2&nbsp;&nbsp;&nbsp;&nbsp;    4&nbsp;&nbsp;&nbsp;&nbsp;     -2&nbsp;&nbsp;&nbsp;&nbsp;    -1

-3&nbsp;&nbsp;&nbsp;&nbsp;   -5&nbsp;&nbsp;&nbsp;&nbsp;     -6&nbsp;&nbsp;&nbsp;&nbsp;    1

-1&nbsp;&nbsp;&nbsp;&nbsp;    2&nbsp;&nbsp;&nbsp;&nbsp;      8&nbsp;&nbsp;&nbsp;&nbsp;    -2]

Here, the following code can be used to get the PLU decomposition.

```
Matrix A("[1,2,-1,0;2,4,-2,-1;-3,-5,6,1;-1,2,8,-2]");
Matrix P(4,4);
Matrix L(4,4);
Matrix U(4,4);
    
A.PLUdecomposition(P,L,U);

```

In this example, the permuation matrix would be stored in P; the lower diagonal matrix would be stored in L; And the upper diagonal matrix would be stored in U.

In order to calculate the determinant and the trace of a matrix, the determinant() and trace() methods can be used. For the above example, the code would be as follows,

```
Matrix A("[1,2,-1,0;2,4,-2,-1;-3,-5,6,1;-1,2,8,-2]");
double det = A.determinant();
double tr = A.trace();

```
Since the determinant is calculated using the PLU decomposition, it currently only works invertible matrices(this will be corrected in the future).

To get the transpose of A, use the code `A.transpose();`

The method diagonal() returns a matrix of the same size with the same entries on the diagonals and zeros on the off-diagonals. For example, consider the following code:

```
Matrix A("[1,2,-1,0;2,4,-2,-1;-3,-5,6,1;-1,2,8,-2]");
Matrix D = A.diagonal();
```
Then the matrix D is 4x4 and has (1, 4, 6, -2) on the diagonals and zeros on the off-diagonals.

The lower method returns a matrix of the same size, with the same entries on the lower diagonal(including the diagonal) and with zeros on the off-diagonal. The call for the method would be as follows:

```
Matrix A("[1,2,-1,0;2,4,-2,-1;-3,-5,6,1;-1,2,8,-2]");
Matrix L = A.lower();
```
This class has two insertRow methods and two insertColumn methods. In one version of the methods (ie. one of insertRow methods and one of the insertColumn methods), the methods take two arguments. The first argument is a C++ standard vector<double>, and the second in an integer. the vector contains the values we want to insert into the matrix (either as a row or column, depending on whether we are calling insertRow or insertColumn) and the second argument(the integer) is the row/column in which we wish to insert these values. It should be noted that this class uses zero-based indexing and so the code `A.insertRow(v,2);` would insert the vector v into the 3rd column of A. Finally, the vector has to be of the correct size. For instance, when using the insertRow method, the size of the vector has to equal the number of columns. As an example, the following code:

```
std::vector<double> v({0,0});
Matrix A("[1,2;3,4]");
A.insertColumn(v,0);
```
inserts 0s into the first column of A. The other vesrion of the insertRows and insertColumn methods work in the same except that the first argument they take is a Matrix. This matrix must be of size NX1 or 1XN and N must be of appropriate size.

### Overloaded Operators

The operators `+` , `-`, and `*` all have been overloaded and are compatible with linear algebra. This means that `*` corresponds to the linear algebra version of multiplication and is not entry-wise multiplication. However, `*` can also be used to multiply a matrix by scalar. Please also not that, `*` is compatible with multiplication the standard C++ vector<double> type. Please consider the following code:
  
```
Matrix A("[1,2;3,4]");
std::vector<double> v({1,2});
    
Matrix C = A*v;
Matrix D = v*A;
```
Here, C is 2x1 matrix containing the values 5 and 11. And D is a 1x2 matrix containing the values 7 and 10.

Please note that `<<` also has been overloaded. Thus matrices can be sent to an output stream. For instance, the following code:

```
Matrix I= identity(2);
std::cout<<I<<'\n';
```

would display [1,0;0,1]. Note that this has the same format as the Matrix constructor that takes string as an input.

The entries of a matrix can be accessed using round brackets. For instance to access the entry in the first row and in the first column of a matrix A, the following code would be used `double a = A(0,0);`. Please note the use of zero-based indexing. This round brackets can also be used to change the values in the matrix. For instance, `A(0,0)=5` would change the value in the first row and first column of A to 5.

Finally, matrices can be updated using the overloaded operators `+=`, `-=`, and `*=`. For instance, the following code,
```
Matrix I= identity(2);
Matrix A("[1;2]");
I*=A;
```

would change the value of I from a 2x2 identity matrix to 2x1 matrix containing the values 1 and 2.

## Linear Systems

In this project, there are four ways to solve a linear equation Ax=b for x where A is a matrix and b is column matrix(matrix of size nx1). Please note that at this stage, the linear solvers require square invertible matrices(this limitation will be eliminated in the future).

The first of these methods uses the class GaussianElimination. As the name suggests, this class uses Gaussian Elimination to solve the linear equation. Let us consider the following example, by declaring and defining the Matrices A and b as follows:

```
Matrix A("[1,2,-1,0;2,4,-2,-1;-3,-5,6,1;-1,2,8,-2]");
Matrix b("[1;-1;3;0]");
```
Please note that size b is 4x1 and not 1x4. Then we can set up the linear system to be solved by Gaussian Elimination by declaring G to be of type GaussianElimination as follows:

`GaussianElimination G(A,b);`

Then we can solve the system and store the result in the variable x as follows:

`Matrix x=G.solve();`

In this example, x end up being 4x1 matrix with entries (2,0,1,3). The GaussianElimination Class has two other methods, getMatrix() and getColumnVector(), both of which take no arguments and return Left hand side(A) and right hand side(b) of the linear system.

The second of these methods is the Jacobi Iterative method. To use this method to solve a linear equation, we must use the Jacobi class. The idea is very similar to that of when we are using the Gaussian Elimination Class. For example, we start defining A and B as follows:

```
Matrix A("[4,1,0,1,0;1,4,1,0,1;0,1,4,1,0;1,0,1,4,1;0,1,0,1,4]");
Matrix b("[1;2;-1;2;1]");
```
Then we need to declare a variable J of type Jacobi as follows:

`Jacobi J(A,b,45);`

Unlike the GaussianElimination class, the Jacobi class constructor can take an extra argument. This argument is an integer that determines the number of iterations for the Jacobi method. The default is set to 40. It should be noted that this class does not directly check for convergence of this method. However, it does provide an indirect method to check convergence. This class has getnormOfError() method that returns the l2 norm of the difference of the last two vectors in the iteration. Please note that this method can only be ran after the iteration has been run (ie. we have solved for x using the Jacobi iteration). From this code, one can determine the convergence of the method.

Then in order to solve this equation, we must provide J with an initial guess. This is done as follows:

```
Matrix v0("[0;0;0;0;0]");
Matrix x = J.solve(v0);
```
Here our initial guess is the zero "vector" (5x1 Matrix containing all zeros). Please note that the solve method takes an initial guess as an argument in the form of column matrix. In this example, we get that x = (-0.1,0.7,-0.6,0.7,-0.1). 

Only after we have run the solve method, can we call the getnormOfError() method. Entering the code `double err = J.getnormOfError();` result in err holding the value 3.06929e-10.

The other two methods, I have implemented so far, for solving linear equations are the Gauss Seidel method and SOR method. These methods can be accessed using the GaussSeidel class and the SOR class. Both of these classes function, in exactly the same way as the Jacobi class and have the same methods.

## Statistics

The statistical functions in this project are very brief. I have only included the functions that I would need.

The normcdf(x) function takes one argument(x) - a double  and returns the normal cumulative distribution function evaluated at x. For instance `double y = normcdf(1.95)` sets the value y to be 0.974412. The norminv(x) function takes one argument(x) - a double and returns the double corresponding to the xth percentile in the normal distribution. Please note that x must be between 0 and 1. For instance, `double y = norminv(0.975)` sets the value y to be 1.9591.

The randuniform function takes 1 argument (n) - an integer, and return a nx1 Matrix with the entries being "random" samples from [0,1]. The boxMuller() function without any arguments returns a random sample from the standard normal distribution. The boxMuller(int n) function that takes one argument (n) - an integer, returns a nx1 Matrix with the entries being random samples from the standard normal.

## Black Scholes Model


