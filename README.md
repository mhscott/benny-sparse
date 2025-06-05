# Benny Sparse

BennySparse is a C++ implementation of the sparse linear equation solver 
presented in 
[_Direct Methods for Sparse Linear Systems_](https://epubs.siam.org/doi/book/10.1137/1.9780898718881) 
by Timothy A. Davis.

## Minimal Examples

Store a 3x3 matrix in triplet form, then compress.

```cpp
#include "benny_sparse.hpp"
#include "benny_utility.hpp"

int main (int argc, char **argv)
{
  benny_sparse kb(3,3,5);
  kb.benny_entry(0,0,1.0);
  kb.benny_entry(1,1,4.0);
  kb.benny_entry(1,2,2.0);
  kb.benny_entry(2,1,2.0);
  kb.benny_entry(2,2,4.0);

  kb.benny_print();

  benny_sparse *kbsparse = kb.benny_compress();

  kbsparse->benny_print();

  delete kbsparse;

  return 0;
}
```

Load a dense matrix from file, then compress.

```cpp
#include "benny_sparse.hpp"
#include "benny_utility.hpp"

#include <fstream>

int main (int argc, char **argv)
{
  ifstream input("matrix_1767.txt");
  benny_sparse T(input,false);
  input.close();
  
  T.benny_print();

  benny_sparse *Tc = T.benny_compress();

  Tc->benny_print();

  delete Tc;

  return 0;
}
```

## Sparse Matrix Class Definition

```cpp
class benny_sparse {
 private:
  int nzmax; // maximum number of entries
  int m; // number of rows
  int n; // number of columns
  int *p; // column pointers
  int *i; // row indices
  double *x; // values
  int nz; // number of entries in triplet matrix  (=-1 if compressed)

  bool is_triplet(void) const {return nz >= 0;}
  bool is_csc(void) const {return nz == -1;}

 public:
  // Construct m x n matrix with nzmax entries
  benny_sparse(int m, int n, int nzmax, bool values = true, bool triplet = true);
  // Construct from input stream, e.g., file of triplets or file of full matrix
  benny_sparse(std::istream &input, bool triplets = false, double tol = 1e-16, bool minus1 = false):
  // Copy constructor
  benny_sparse(const benny_sparse& A);
  // Default constructor, uses nzmax = 1
  benny_sparse();
  // Destructor
  ~benny_sparse()
};
```
