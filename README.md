# Benny Sparse

BennySparse is a C++ implementation of the sparse linear equation solver 
presented by Timothy A. Davis in 
[_Direct Methods for Sparse Linear Systems_](https://epubs.siam.org/doi/book/10.1137/1.9780898718881).

[YouTube playlist](https://youtube.com/playlist?list=PL5EvFKC69QIyRLFuxWRnH6hIw6e1-bBXB) of lectures from the book.

## Compiling at Command Line

```console
g++ main.cpp benny_sparse.cpp benny_symbolic.cpp benny_numeric.cpp
```

## Minimal Examples

Store a 3x3 matrix in triplet form, then compress.

```cpp
#include "benny_sparse.hpp"

int main (int argc, char **argv)
{
  benny_sparse kb(3,3,5);
  kb.benny_entry(0,0,1.0);
  kb.benny_entry(1,1,4.0);
  kb.benny_entry(1,2,2.0);
  kb.benny_entry(2,1,2.0);
  kb.benny_entry(2,2,4.0);

  kb.benny_print(std::cout);

  benny_sparse *kbsparse = kb.benny_compress();

  kbsparse->benny_print(std::cout);

  delete kbsparse;

  return 0;
}
```

Load a dense matrix from file, then compress and solve.

```cpp
#include "benny_sparse.hpp"

#include <fstream>

int main (int argc, char **argv)
{
  ifstream input("matrix_1767.txt");
  benny_sparse T(input,false);
  input.close();
  
  T.benny_print(std::cout);

  benny_sparse *Tc = T.benny_compress();

  Tc->benny_print(std::cout);

  double *x = new double[Tc->ncols()];

  for (int i = 0; i < Tc->ncols(); i++) x[i] = 0.0;
  x[0] = 1;
  int ok = Tc->benny_cholsol(0,x);
  if (ok < 0) {
    std::cout << "Cholesky decomposition failed" << std::endl;
  } 
  else {
    std::cout << x[0] << std::endl;
  }

  for (int i = 0; i < Tc->ncols(); i++) x[i] = 0.0;
  x[0] = 1;
  ok = Tc->benny_lusol(0,x,0.001);  
  if (ok < 0) {
    std::cout << "LU decomposition failed" << std::endl;
  }
  else {
    std::cout << x[0] << std::endl;
  }

  delete [] x;
  delete Tc;

  return 0;
}
```
