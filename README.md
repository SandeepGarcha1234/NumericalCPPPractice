# NumericalCPPPractice
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


