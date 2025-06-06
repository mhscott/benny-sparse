#ifndef BENNY_SPARSE_HPP
#define BENNY_SPARSE_HPP

#include "benny_utility.hpp"

#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

class benny_sparse {
private:
 int nzmax; // maximum number of entries
 int m; // number of rows
 int n; // number of columns
 int *p; // column pointers
 int *i; // row indices
 double *x; // values
 int nz; // number of entries in triplet matrix  (=-1 if compressed)

 int benny_sprealloc(int newnzmax)
 {
  if (newnzmax <= 0) newnzmax = is_csc() ? p[n] : nz;

  // Resize the row index vector
  int *newi = new int[newnzmax];
  for (int i = 0; i < nzmax; i++) newi[i] = this->i[i];
  delete [] this->i;
  this->i = newi;

  // Resize the column index vector
  if (is_triplet()) {
    int *newp = new int[newnzmax];
    for (int i = 0; i < nzmax; i++) newp[i] = this->p[i];
    delete [] this->p;
    this->p = newp;
  }

  // Resize the values vector
  if (this->x != 0) {
    double *newx = new double[newnzmax];
    for (int i = 0; i < nzmax; i++) newx[i] = this->x[i];
    delete [] this->x;
    this->x = newx;
  }

  nzmax = newnzmax;

  return 0;
 }

public:
 bool is_triplet(void) const {return nz >= 0;}
 bool is_csc(void) const {return nz == -1;}
 int nrows(void) const {return m;}
 int ncols(void) const {return n;}
 int ntriplets(void) const {return nz;}

 /// @brief Print the matrix in triplet or compressed format, C++ implementation of cs_print on page 23
 /// @param brief - limits the amount of printing for large matrices
 /// @return 0 success, < 0 failure
 int benny_print(bool brief = true) const
 {
  std::cout << "BennySparse Version 1.0" << std::endl;
  if (nz < 0) {
    std::cout << m << "-by-" << n
    << ", nzmax: " << nzmax
    << ", nnz: " << this->p[n] << std::endl;
    for (int j = 0; j < n; j++) {
      std::cout << "   col " << j << " : locations "
      << p[j] << " to " << p[j+1]-1 << std::endl;
      for (int p = this->p[j]; p < this->p[j+1]; p++) {
        std::cout << "      " << i[p] << " : " << (this->x ? this->x[p] : 1) << std::endl;
        if (brief && p > 20) {
          std::cout << std::endl;
          return 0;
        }
      }
    }
  }
  else {
    std::cout << "triplet: " << m << "-by-" << n
    << ", nzmax: " << nzmax
    << ", nnz: " << nz << std::endl;
    for (int p = 0; p < nz; p++) {
        std::cout << "      " << i[p] << ' ' << this->p[p] << " : " << (this->x ? this->x[p] : 1) << std::endl;
        if (brief && p > 20) {
          std::cout << std::endl;
          return 0;
        }
      }
  }
  std::cout << std::endl;
  return 0;
 }
 
 /// @brief C++ implementation of cs_spalloc on page 11
 /// @param m
 /// @param n
 /// @param nzmax
 /// @param values 
 /// @param triplet 
 benny_sparse(int m, int n, int nzmax, bool values = true, bool triplet = true):
 nzmax(nzmax),m(m),n(n),p(0),i(0),x(0),nz(0)
 {
  if (nzmax < 1) nzmax = 1;
  nz = triplet ? 0 : -1;

  p = new int[triplet ? nzmax : n+1];
  i = new int[nzmax];
  x = values ? new double[nzmax] : 0;
 }

 benny_sparse():nzmax(1),m(0),n(0),p(0),i(0),x(0),nz(0)
 {
  // Create object with storage for 1 value
  this->p = new int[1];
  this->i = new int[1];
  this->x = new double[1];
 }

 benny_sparse(std::istream &input, bool triplets = false, double tol=1e-16, bool minus1 = false):
 nzmax(1),m(0),n(0),p(0),i(0),x(0),nz(0)
 {
   this->p = new int[1];
   this->i = new int[1];
   this->x = new double[1];

   int i, j;
   double val;
   if (triplets) {
    while (input >> i) {
      input >> j >> val;
      if (fabs(val) < tol) continue;
      if (minus1) {
        i--;
        j--;
      }
      this->benny_entry(i,j,val);
    }
  } else {
   i = 0;
   std::string str;
   while (getline(input,str)) {
    istringstream ss(str);
    j = 0;
    while (ss >> val) {
      if (fabs(val) >= tol) this->benny_entry(i,j,val);
      j++;
    }
    i++;
   }
  }
}

 benny_sparse(const benny_sparse& A):x(0)
 {
  nzmax = A.nzmax;
  m = A.m;
  n = A.n;
  nz = A.nz;

  p = new int[nzmax];
  for (int i = 0; i < nzmax; i++) p[i] = A.p[i];

  i = new int[nzmax];
  for (int i = 0; i < nzmax; i++) this->i[i] = A.i[i];

  if (A.x != 0) {
    x = new double[nzmax];
    for (int i = 0; i < nzmax; i++) x[i] = A.x[i];
  }
 }

 /// @brief C++ implementation of cs_spfree on page 11
 ~benny_sparse()
 {
   if (p != 0) delete [] p;
   if (i != 0) delete [] i;
   if (x != 0) delete [] x;
 }

 /// @brief Compute y = A*x + y where A is this matrix
 int benny_gaxpy(const double *x, double *y) const
{
if (!is_csc() || x == 0 || y == 0) return -1;

for (int j = 0; j < n; j++)
  for (int p = this->p[j]; p < this->p[j+1]; p++)
  y[this->i[p]] += this->x[p]*x[j];

return 0;
}

 /// @brief Add an entry to this matrix
 /// @param r 
 /// @param c 
 /// @param val 
 /// @return 
 int benny_entry(int r, int c, double val)
 {
   // Check if matrix has been compressed or (r,c) indices out of bounds
   if (!is_triplet() || r < 0 || c < 0) return -1;

   // Need to check size
   if (nz >= nzmax) this->benny_sprealloc(2*nzmax);

   // Assign triplet
   if (x != 0) x[nz] = val;
   if (i != 0) i[nz] = r;
   if (p != 0) p[nz++] = c;

   // Dimensions
   m = (m >= r+1) ? m : r+1;
   n = (n >= c+1) ? n : c+1;

   return 0;
 }

 /// @brief Permute a sparse matrix, C=P*A*Q
 /// C++ implementation of cs_permute on page 21
 /// @param pinv Inverse row permutation vector
 /// @param q Column permutation vector
 /// @param values Boolean on whether or not to copy values (default = true)
 /// @return new benny_sparse object -- calling object should delete
 benny_sparse* benny_permute(const int *pinv, const int *q, bool values = true)
 {
  if (!is_csc()) return 0; // Not compressed

  benny_sparse *C = new benny_sparse(this->m, this->n, this->p[n], values && (this->x != 0), false);

  int nz = 0;
  for (int k = 0; k < n; k++) {
    C->p[k] = nz;
    int j = q ? q[k] : k;
    for (int t = this->p[j]; t < this->p[j+1]; t++) {
      if (C->x != 0) C->x[nz] = this->x[t];
      C->i[nz++] = pinv ? pinv[this->i[t]] : this->i[t];
    }
  }
  C->p[n] = nz;

  return C;
 }

 /// @brief Transforms triplets into compressed matrix format
 /// C++ implementation of cs_compress on page 13
 /// @param void
 /// @return new benny_sparse object -- calling object should delete
 benny_sparse* benny_compress(void) {
   if (!is_triplet()) return 0; // Already compressed

   benny_sparse *C = new benny_sparse(*this);

   int *w = new int[n];
   for (int i = 0; i < n; i++) w[i] = 0;
   for (int k = 0; k < nz; k++) w[p[k]]++;

   int *Cp = C->p;
   benny_cumsum(Cp, w, n);

   int *Ci = C->i;
   double *Cx = C->x;
   int p;
   for (int k = 0; k < nz; k++) {
     Ci[p = w[this->p[k]]++] = this->i[k];
     if (Cx != 0) Cx[p] = this->x[k];
   }

   C->nz = -1;

   delete [] w;

   return C;
 }

 /// @brief Combines duplicate entires in a compressed matrix
 /// C++ implementation of cs_dupl on page 15
 /// @param void
 /// @return 0 if success, < 0 if failure
 int benny_dupl(void)
 {
  if (!is_csc()) return 0; // Not compressed

  int *w = new int[m];
  for (int i = 0; i < m; i++) w[i] = -1;
  int nz = 0;
  for (int j = 0; j < n; j++) {
    int q = nz;
    for (int p = this->p[j]; p < this->p[j+1]; p++) {
      int i = this->i[p];
      if (w[i] >= q) this->x[w[i]] += this->x[p];
      else {
        w[i] = nz;
        this->i[nz] = i;
        this->x[nz++] = this->x[p];
      }
    }
    this->p[j] = q;
  }
  this->p[n] = nz;
  delete [] w;

  return 0;
 }

 /// @brief Solve Lx=b where L is lower triangular
 /// @param x On input, x contains b
 /// @return 
 int benny_lsolve(double *x) const
 {
  if (!is_csc() || x == 0) return -1;

  int *Lp = this->p;
  int *Li = this->i;
  double *Lx = this->x;
  for (int j = 0; j < n; j++) {
    x[j] /= Lx[Lp[j]];
    for (int p = Lp[j]+1; p < Lp[j+1]; p++) x[Li[p]] -= Lx[p]*x[j];
  }

  return 0;
 }

 /// @brief Solve L'x=b where L is lower triangular
 /// @param x On input, x contains b
 /// @return 
 int benny_ltsolve(double *x) const
 {
  if (!is_csc() || x == 0) return -1;

  int *Lp = this->p;
  int *Li = this->i;
  double *Lx = this->x;
  for (int j = n-1; j >= 0; j--) {
    for (int p = Lp[j]+1; p < Lp[j+1]; p++) x[j] -= Lx[p]*x[Li[p]];
    x[j] /= Lx[Lp[j]];
  }

  return 0;
 } 

 /// @brief Solve Ux=b where U is upper triangular
 /// @param x On input, x contains b
 /// @return 
 int benny_usolve(double *x) const
 {
  if (!is_csc() || x == 0) return -1;

  int *Up = this->p;
  int *Ui = this->i;
  double *Ux = this->x;
  for (int j = 0; j < n; j++) {
    x[j] /= Ux[Up[j+1]-1];
    for (int p = Up[j]; p < Up[j+1]-1; p++) x[Ui[p]] -= Ux[p]*x[j];
  }

  return 0;
 }

 /// @brief Solve U'x=b where L is uppser triangular
 /// @param x On input, x contains b
 /// @return 
 int benny_utsolve(double *x) const
 {
  if (!is_csc() || x == 0) return -1;

  int *Up = this->p;
  int *Ui = this->i;
  double *Ux = this->x;
  for (int j = 0; j < n; j++) {
    for (int p = Up[j]; p < Up[j+1]-1; p++) x[j] -= Ux[p]*x[Ui[p]];
    x[j] /= Ux[Up[j+1]-1];
  }

  return 0;
 } 

};

#endif
