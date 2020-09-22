# NumericalCPPPractice
I have recently started to learn C++ and the purpose of this code is for me to gain practice while working on some projects relevant to quantitative finance 

## The Matrix Class

This class stores two dimensional arrays. There are some linear algebra functionalities, however more need to be added. It should be noted that this class uses zero-based indexing.

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
