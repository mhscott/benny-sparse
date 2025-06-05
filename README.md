# Benny Sparse

BennySparse is a C++ implementation of the sparse linear equation solver 
presented in 
[_Direct Methods for Sparse Linear Systems_](https://epubs.siam.org/doi/book/10.1137/1.9780898718881) 
by Timothy A. Davis.

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
