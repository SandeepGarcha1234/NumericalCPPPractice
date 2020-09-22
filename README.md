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

The PLUdecomposition returns the PLU decomposition for square invertible matrices(this will be expanded in the future). The result is returned in matrices that are passed by reference and they must be passed in the correct order (permutation matrix first, then the lower diagonal matrix, and finally the upper diagonal matrix.







