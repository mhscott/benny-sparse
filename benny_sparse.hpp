#ifndef BENNY_SPARSE_HPP
#define BENNY_SPARSE_HPP

#include <iostream>
#include <sstream>
#include <cmath>

#define BENNY_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define BENNY_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define BENNY_FLIP(i) (-(i)-2)
#define BENNY_UNFLIP(i) (((i) < 0) ? BENNY_FLIP(i) : (i))
#define BENNY_MARKED(w,j) (w[j] < 0)
#define BENNY_MARK(w,j) {w[j] = BENNY_FLIP(w[j]); }

using namespace std;

class benny_sparse;

class benny_numeric {
 private:
  benny_sparse *L;
  benny_sparse *U;
  int *pinv;
  double *B;
 public:
  benny_numeric():L(0), U(0), pinv(0), B(0) {}
  ~benny_numeric();
  friend class benny_sparse;
};

class benny_symbolic {
 private:
  int *pinv;
  int *q;
  int *parent;
  int *cp;
  int *leftmost;
  int n2;
  double lnz;
  double unz;
 public:
  benny_symbolic():pinv(0), q(0), parent(0), cp(0), leftmost(0), n2(0), lnz(0.0), unz(0.0) {}
  ~benny_symbolic() {
    if (pinv != 0) delete [] pinv;
    if (q != 0) delete [] q;
    if (parent != 0) delete [] parent;
    if (cp != 0) delete [] cp;
    if (leftmost != 0) delete [] leftmost;
  }
  friend class benny_sparse;    
};

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
  for (int i = 0; i < BENNY_MIN(nzmax,newnzmax); i++) newi[i] = this->i[i];
  delete [] this->i;
  this->i = newi;

  // Resize the column index vector
  if (is_triplet()) {
    int *newp = new int[newnzmax];
    for (int i = 0; i < BENNY_MIN(nzmax,newnzmax); i++) newp[i] = this->p[i];
    delete [] this->p;
    this->p = newp;
  }

  // Resize the values vector
  if (this->x != 0) {
    double *newx = new double[newnzmax];
    for (int i = 0; i < BENNY_MIN(nzmax,newnzmax); i++) newx[i] = this->x[i];
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

 /// @brief Query a value from compressed matrix (mainly for testing)
 /// @param r Row index, 0-based
 /// @param c Column index, 0-based
 /// @return A(r,c) value if present, otherwise 0.0
 /// Also returns zero if matrix not compressed or r or c out of bounds
 double value(int r, int c) const
 {
  if (!is_csc() || r+1 > m || c+1 > n) return 0.0;

  for (int j = p[c]; j < p[c+1]; j++) {
    if (i[j] == r) return x[j];
  }

  return 0.0;
}

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
 benny_sparse(int m, int n, int nzmax = 1, bool values = true, bool triplet = true):
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

 double benny_norm()
 {
   if (!is_csc() || x == 0) return -1.0;
   double norm = 0.0;
   for (int j = 0; j < n; j++) {
    double s = 0.0;
    for (int p = this->p[j]; p < this->p[j+1]; p++) s += fabs(this->x[p]);
    norm = BENNY_MAX(norm, s);
   }
   return norm;
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

 int benny_scatter(int j, double beta, int *w, double *x, int mark, benny_sparse *C, int nz) const
 {
  if (!is_csc() || w == 0 || !C->is_csc()) return -1;

  int *Ap = p; int *Ai = i; double *Ax = this->x; int *Ci = C->i;
  for (int p = Ap[j]; p < Ap[j+1]; p++) {
    int i = Ai[p];
    if (w[i] < mark) {
      w[i] = mark;
      Ci[nz++] = i;
      if (x != 0) x[i] = beta*Ax[p];
    }
    else {
      if (x != 0) x[i] += beta*Ax[p];
    }
  }

  return nz;
 }

 benny_sparse *benny_add(const benny_sparse *B, double alpha, double beta)
 {
  if (!is_csc() | !B->is_csc()) return 0;

  int anz = p[n];
  int *Bp = B->p; int *Bi = B->i; double *Bx = B->x; int bnz = Bp[B->n];
  int *w = new int[m];
  for (int i = 0; i < m; i++) w[i] = 0;
  bool values = (this->x != 0) && (Bx != 0);
  double *x = values ? new double[m] : 0;
  benny_sparse *C = new benny_sparse(m, B->n, anz + bnz, values, false);
  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  int nz = 0;
  for (int j = 0; j < B->n; j++) {
    Cp[j] = nz;
    nz = this->benny_scatter(j, alpha, w, x, j+1, C, nz);
    nz = B->benny_scatter(j, beta, w, x, j+1, C, nz);
    if (values) {
      for (int p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
    }
  }
  Cp[B->n] = nz;
  C->benny_sprealloc(0);

  delete [] w;
  delete [] x;

  return C;
 }

 benny_sparse *benny_multiply(const benny_sparse *B)
 {
  if (!is_csc() || !B->is_csc()) return 0;

  int anz = p[n];
  int *Bp = B->p; int *Bi = B->i; double *Bx = B->x; int bnz = Bp[B->n];
  int *w = new int[m];
  for (int i = 0; i < m; i++) w[i] = 0;
  bool values = (this->x != 0) && (Bx != 0);
  double *x = values ? new double[m] : 0;
  benny_sparse *C = new benny_sparse(m, B->n, anz + bnz, values, false);
  int *Cp = C->p;
  int nz = 0;
  for (int j = 0; j < B->n; j++) {
    if (nz + m > C->nzmax && !C->benny_sprealloc(2*(C->nzmax)+m)) {
      delete [] w;
      delete [] x;
      delete C;
      return 0;
    }
    int *Ci = C->i; double *Cx = C->x;
    Cp[j] = nz;
    for (int p = Bp[j]; p < Bp[j+1]; p++) {
      nz = this->benny_scatter(Bi[p], Bx ? Bx[p] : 1, w, x, j+1, C, nz);
    }
    if (values) {
      for (int p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
    }
  }
  Cp[B->n] = nz;
  C->benny_sprealloc(0);

  delete [] w;
  delete [] x;

  return C;
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

 static double benny_cumsum(int *p, int *c, int n)
{
  if (p == 0 || c == 0) return -1;

  int nz = 0;
  double nz2 = 0.0;
  for (int i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i];
    c[i] = p[i];
  }
  p[n] = nz;
  
  return nz2;
}
 
static int benny_pvec(const int *p, const double *b, double *x, int n)
{
  if (x == 0 || b == 0) return -1;
  for (int k = 0; k < n; k++) x[k] = b[p ? p[k] : k];
  return 0;
}

static int benny_ipvec(const int *p, const double *b, double *x, int n)
{
  if (x == 0 || b == 0) return -1;
  for (int k = 0; k < n; k++) x[p ? p[k] : k] = b[k];
  return 0;
}

static int* benny_pinv(const int *p, int n)
{
  int *pinv = new int[n];
  for (int k = 0; k < n; k++) pinv[p[k]] = k;
  return pinv;
}

 /// @brief Permute a sparse matrix, C=P*A*Q
 /// C++ implementation of cs_permute on page 21
 /// @param pinv Inverse row permutation vector
 /// @param q Column permutation vector
 /// @param values Boolean on whether or not to copy values (default = true)
 /// @return new benny_sparse object -- calling object should delete
 benny_sparse *benny_permute(const int *pinv, const int *q, bool values = true)
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

 benny_sparse *benny_symperm(const int *pinv, bool values = true)
 {
  if (!is_csc()) return 0;

  int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
  benny_sparse *C = new benny_sparse(n, n, Ap[n], values && Ax != 0, false);
  if (C == 0) {
    std::cerr << "benny_symperm - ran out of memory creating C" << std::endl;
    return 0;
  }
  int *w = new int[n];
  if (w == 0) {
    std::cerr << "benny_symperm - ran out of memory creating w" << std::endl;
    return 0;
  }
  for (int i = 0; i < n; i++) w[i] = 0;
  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  for (int j = 0; j < n; j++) {
    int j2 = pinv ? pinv[j] : j;
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      int i = Ai[p];
      if (i > j) continue;
      int i2 = pinv ? pinv[i] : i;
      w[BENNY_MAX(i2,j2)]++;
    }
  }
  benny_cumsum(Cp, w, n);
  for (int j = 0; j < n; j++) {
    int j2 = pinv ? pinv[j] : j;
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      int i = Ai[p];
      if (i > j) continue;
      int i2 = pinv ? pinv[i] : i;
      int q;
      Ci[q = w[BENNY_MAX(i2,j2)]++] = BENNY_MIN(i2,j2);
      if (Cx != 0) Cx[q] = Ax[p];
    }
  }

  delete [] w;

  return C;
 }

 /// @brief Transforms triplets into compressed matrix format
 /// C++ implementation of cs_compress on page 13
 /// @param void
 /// @return new benny_sparse object -- calling object should delete
 benny_sparse *benny_compress(void)
 {
   if (!is_triplet()) return 0; // Already compressed

   int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
   benny_sparse *C = new benny_sparse(m, n, nz, Ax != 0, false);

   int *w = new int[n];
   for (int i = 0; i < n; i++) w[i] = 0;
   for (int k = 0; k < nz; k++) w[p[k]]++;

   int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
   benny_cumsum(Cp, w, n);

   int p;
   for (int k = 0; k < nz; k++) {
     Ci[p = w[Ap[k]]++] = Ai[k];
     if (Cx != 0) Cx[p] = Ax[k];
   }

   C->nz = -1;

   delete [] w;

   return C;
 }

 benny_sparse *benny_transpose(bool values = true)
 {
  if (!is_csc()) return 0;

  int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
  benny_sparse *C = new benny_sparse(this->m, this->n, Ap[n], values && Ax != 0, false);

  int *w = new int[n];
  for (int i = 0; i < n; i++) w[i] = 0;

  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  for (int p = 0; p < Ap[n]; p++) w[Ai[p]]++;
  benny_cumsum(Cp, w, m);

  int q;
  for (int j = 0; j < n; j++) {
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      Ci[q = w[Ai[p]]++] = j;
      if (Cx != 0) Cx[q] = Ax[p];
    }
  }

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
  for (int j = n-1; j >= 0; j--) {
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

private:

 public:
 /// @brief C++ implementation of cs_reach function on page 33
 /// @param G 
 /// @param B 
 /// @param k 
 /// @param xi 
 /// @param pinv 
 /// @return 
 int benny_reach(benny_sparse *G, const benny_sparse *B, int k, int *xi, const int *pinv)
 {
   if (!G->is_csc() || !B->is_csc() || xi == 0) return -1;

   int n = G->n; int *Bp = B->p; int *Bi = B->i; int *Gp = G->p;
   int top = G->n;
   for (int p = Bp[k]; p < Bp[k+1]; p++)
     if (!BENNY_MARKED(Gp,Bi[p])) top = benny_dfs(Bi[p],G,top,xi,xi+n,pinv);
   for (int p = top; p < n; p++) BENNY_MARK(Gp, xi[p]);

   return top;
 }

 /// @brief C++ implementation of cs_dfs function on page 33
 /// @param j 
 /// @param G 
 /// @param top 
 /// @param xi 
 /// @param pstack 
 /// @param pinv 
 /// @return 
 int benny_dfs(int j, benny_sparse *G, int top, int *xi, int *pstack, const int *pinv)
 {
   if (!G->is_csc() || xi == 0 || pstack == 0) return -1;

   int *Gp = G->p; int *Gi = G->i;
   int head = 0;
   xi[0] = j;
   while (head >= 0) {
    j = xi[head];
    int jnew = (pinv != 0) ? pinv[j] : j;
    if (!BENNY_MARKED(Gp,j)) {
      BENNY_MARK(Gp,j);
      pstack[head] = (jnew < 0) ? 0 : BENNY_UNFLIP(Gp[jnew]);
    } 
    bool done = true;
    int p2 = (jnew < 0) ? 0 : BENNY_UNFLIP(Gp[jnew+1]);
    for (int p = pstack[head]; p < p2; p++) {
      int i = Gi[p];
      if (BENNY_MARKED(Gp,i)) continue;
      pstack[head] = p;
      xi[++head] = i;
      done = false;
      break;
    }
    if (done) {
      head--;
      xi[--top] = j;
    }
   }

   return top;
 }

 int benny_spsolve(benny_sparse *G, const benny_sparse *B, int k, int *xi, double *x, const int *pinv, bool lo = true)
 {
   if (!G->is_csc() || !B->is_csc() || xi == 0 || x == 0) return -1;

   int *Gp = G->p; int *Gi = G->i; double *Gx = G->x; int n = G->n;
   int *Bp = B->p; int *Bi = B->i; double *Bx = B->x; 
   int top = benny_reach(G, B, k, xi, pinv);
   for (int p = top; p < n; p++) x[xi[p]] = 0;
   for (int p = Bp[k]; p < Bp[k+1]; p++) x[Bi[p]] = Bx[p];
   for (int px = top; px < n; px++) {
    int j = xi[px];
    int J = (pinv != 0) ? pinv[j] : j;
    if (J < 0) continue;
    x[j] /= Gx[lo ? Gp[J] : Gp[J+1]-1];
    int p = lo ? Gp[J]+1 : Gp[J];
    int q = lo ? Gp[J+1] : Gp[J+1]-1;
    for ( ; p < q; p++) x[Gi[p]] -= Gx[p] * x[j];
   }
   
   return top;
 }

public:
 /// @brief C++ implementation of cs_etree on page 42
 /// @param ata true - elimination tree of Cholesky factorization of A'*A
 /// false - elimination tree of Cholesky factorization of A
 /// @return new elimination tree array -- calling object should delete
 int *benny_etree(bool ata = false)
 {
  if (!is_csc()) return 0;

  int *Ap = this->p;
  int *Ai = this->i;

  int *parent = new int[n];
  int *w = new int[n + (ata ? m : 0)];

  int *ancestor = w;
  int *prev = w + n;

  if (ata) for (int i = 0; i < m; i++) prev[i] = -1;
  for (int k = 0; k < n; k++) {
    parent[k] = -1;
    ancestor[k] = -1;
    for (int p = Ap[k]; p < Ap[k+1]; p++) {
      int i = ata ? prev[Ai[p]] : Ai[p];
      int inext;
      for ( ; i != -1 && i < k; i = inext) {
        inext = ancestor[i];
        ancestor[i] = k;
        if (inext == -1) parent[i] = k;
      }
      if (ata) prev[Ai[p]] = k;
    }
  }

  delete [] w;

  return parent;
 }

 /// @brief C++ implementation of cs_ereach on page 43
 /// @param k 
 /// @param parent 
 /// @param s 
 /// @param w 
 /// @return 
 int benny_ereach(int k, const int *parent, int *s, int *w)
 {
  if (!is_csc() || parent == 0 || s == 0 || w == 0) return -1;

  int top = n;
  int *Ap = this->p; int *Ai = this->i;
  BENNY_MARK(w,k);
  for (int p = Ap[k]; p < Ap[k+1]; p++) {
    int i = Ai[p];
    if (i > k) continue;
    int len;
    for (len = 0; !BENNY_MARKED(w,i); i = parent[i]) {
      s[len++] = i;
      BENNY_MARK(w,i);
    }
    while (len > 0) s[--top] = s[--len];
  }
  for (int p = top; p < n; p++) BENNY_MARK(w,s[p]);
  BENNY_MARK(w,k);
  return top;
 }

private:
 /// @brief Depth-first search of a tree, C++ implementation of cs_tdfs on page 45
 /// @param j 
 /// @param k 
 /// @param head 
 /// @param next 
 /// @param post 
 /// @param stack 
 /// @return 
 int benny_tdfs(int j, int k, int *head, const int *next, int *post, int *stack)
 {
  if (head == 0 || next == 0 || post == 0 || stack == 0) return -1;

  int top = 0;
  stack[0] = j;
  while (top >= 0) {
    int p = stack[top];
    int i = head[p];
    if (i == -1) {
      top--;
      post[k++] = p;
    }
    else {
      head[p] = next[i];
      stack[++top] = i;
    }
  }
  return k;
 }

public:
 /// @brief Postorder a tree, C++ implementation of cs_post on page 45
 /// @param parent 
 /// @param n 
 /// @return new postorder array -- calling object should delete
 int *benny_post(const int *parent, int n)
 {
   if (parent == 0) return 0;
   
   int *post = new int[n];
   int *w = new int[3*n];
   int *head = w; int *next = w+n; int *stack = w+2*n;
   for (int j = 0; j < n; j++) head[j] = -1;
   for (int j = n-1; j >= 0; j--) {
     if (parent[j] == -1) continue;
     next[j] = head[parent[j]];
     head[parent[j]] = j;
   }
   int k = 0;
   for (int j = 0; j < n; j++) {
     if (parent[j] != -1) continue;
     k = benny_tdfs(j, k, head, next, post, stack);
   }

   delete [] w;

   return post;
 }

private:
 /// @brief Find level and first descendant of a node in an elimination tree
 /// C++ implmentation of firstdesc on page 47
 /// @param n 
 /// @param parent 
 /// @param post 
 /// @param first 
 /// @param level 
 void benny_firstdesc(int n, int *parent, int *post, int *first, int *level)
 {
   for (int i = 0; i < n; i++) first[i] = -1;
   for (int k = 0; k < n; k++) {
     int i = post[k];
     int len = 0;
     int r;
     for (r = i; r != -1 && first[r] == -1; r = parent[r], len++) first[r] = k;
     len += (r == -1) ? -1 : level[r];
     for (int s = i; s != r; s = parent[s]) level[s] = len--;
   }  
 }

 /// @brief C++ implementation of cs_leaf on page 51
 /// @param i 
 /// @param j 
 /// @param first 
 /// @param maxfirst 
 /// @param prevleaf 
 /// @param ancestor 
 /// @param jleaf 
 /// @return 
 int benny_leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf)
 {
  if (first == 0 || maxfirst == 0 || prevleaf == 0 || ancestor == 0 || jleaf == 0) return -1;

  *jleaf = 0;
  if (i <= j || first[j] <= maxfirst[i]) return -1;
  maxfirst[i] = first[j];
  int jprev = prevleaf[i];
  prevleaf[i] = j;
  *jleaf = (jprev == -1) ? 1 : 2;
  if (*jleaf == 1) return i;
  int q;
  for (q = jprev; q != ancestor[q]; q = ancestor[q]) {}
  int sparent;
  for (int s = jprev; s != q; s = sparent) {
    sparent = ancestor[s];
    ancestor[s] = q;
  }

  return q;
 }

public:
 int *benny_rowcnt(int *parent, int *post)
 {
  int *Ap = this->p; int *Ai = this->i;
  int *w = new int[5*n];
  int *ancestor = w;
  int *maxfirst = w+n;
  int *prevleaf = w+2*n;
  int *first = w+3*n;
  int *level = w+4*n;

  int *rowcount = new int[n];

  benny_firstdesc(n, parent, post, first, level);

  for (int i = 0; i < n; i++) {
    rowcount[i] = 1;
    prevleaf[i] = -1;
    maxfirst[i] = -1;
    ancestor[i] = i;
  }
  for (int k = 0; k < n; k++) {
    int j = post[k];
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      int i = Ai[p];
      int jleaf;
      int q = benny_leaf(i,j,first,maxfirst,prevleaf,ancestor,&jleaf);
      if (jleaf != 0) rowcount[i] += level[j] - level[q];
    }
    if (parent[j] != -1) ancestor[j] = parent[j];
  }

  delete [] w;

  return rowcount;
 }

private:
 void init_ata(const int *post, int *w, int **head, int **next)
 {
  int *ATp = this->p; int *ATi = this->i;
  *head = w+4*n; *next = w+5*n+1;
  for (int k = 0; k < n; k++) w[post[k]] = k;
  for (int i = 0; i < n; i++) {
    int k = n, p;
    for (int p = ATp[i]; p < ATp[i+1]; p++) k = (k <= w[ATi[p]]) ? k : w[ATi[p]];
    (*next)[i] = (*head)[k];
    (*head)[k] = i;
  }
 }

public:
 int *benny_counts(const int *parent, const int *post, bool ata = false)
 {
  if (!is_csc() || parent == 0 || post == 0) return 0;

  int s = 4*n + (ata ? n+m+1 : 0);
  int *colcount = new int[n];
  int *delta = colcount;
  int *w = new int[s];
  benny_sparse *AT = this->benny_transpose(false);

  int *ancestor = w;
  int *maxfirst = w+n;
  int *prevleaf = w+2*n;
  int *first = w+3*n;
  for (int k = 0; k < s; k++) w[k] = -1;
  for (int k = 0; k < n; k++) {
    int j = post[k];
    delta[j] = (first[j] == -1) ? 1 : 0;
    for ( ; j != -1 && first[j] == -1; j = parent[j]) first[j] = k;
  }

  int *ATp = AT->p; int *ATi = AT->i;

  int *head = 0;
  int *next = 0;
  if (ata) AT->init_ata(post, w, &head, &next);
  for (int i = 0; i < n; i++) ancestor[i] = i;
  for (int k = 0; k < n; k++) {
    int j = post[k];
    if (parent[j] != -1) delta[parent[j]]--;
    for (int J = ata ? head[k] : j; J != -1; J = ata ? next[J] : -1) {
      for (int p = ATp[J]; p < ATp[J+1]; p++) {
        int i = ATi[p];
        int jleaf;
        int q = benny_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
        if (jleaf >= 1) delta[j]++;
        if (jleaf == 2) delta[q]--;
      }
    }
    if (parent[j] != -1) ancestor[j] = parent[j];
  }

  for (int j = 0; j < n; j++) if (parent[j] != -1) colcount[parent[j]] += colcount[j];

  delete AT;
  delete [] w;

  return colcount;
 }

public:

 benny_symbolic *benny_schol(int order)
 {
  if (!is_csc()) return 0;

  benny_symbolic *S = new benny_symbolic();

  // Natural ordering
  int *P = new int[this->n];
  for (int i = 0; i < this->n; i++) P[i] = i;

  S->pinv = benny_pinv(P, this->n);
  delete [] P;
  if (order != 0 && S->pinv == 0) {
    delete S;
    return 0;
  }
  benny_sparse *C = this->benny_symperm(S->pinv, false);
  S->parent = C->benny_etree(false);
  int *post = benny_post(S->parent, this->n);
  int *c = C->benny_counts(S->parent, post, false);
  delete [] post;
  delete C;
  S->cp = new int[n+1];
  S->unz = S->lnz = benny_cumsum(S->cp, c, n);
  delete [] c;

  if (S->lnz >= 0) return S;
  else {
    delete S;
    return 0;
  }
 }

 benny_numeric *benny_chol(const benny_symbolic *S)
 {
   if (!is_csc() || S == 0 || S->cp == 0 || S->parent == 0) return 0;

  benny_numeric *N = new benny_numeric();
  int *c = new int[2*n];
  double *x = new double[n];

  int *cp = S->cp; int *pinv = S->pinv; int *parent = S->parent;
  benny_sparse *C = pinv ? this->benny_symperm(pinv, true) : this;
  benny_sparse *E = pinv ? C : 0;
  int *s = c+n;
  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  
  benny_sparse *L;
  N->L = L = new benny_sparse(n, n, cp[n], true, false);
  int *Lp = L->p; int *Li = L->i; double *Lx = L->x;
  for (int k = 0; k < n; k++) Lp[k] = c[k] = cp[k];
  for (int k = 0; k < n; k++) {
    // Non-zero pattern of L(k,:)
    int top = C->benny_ereach(k, parent, s, c);
    x[k] = 0.0;
    for (int p = Cp[k]; p < Cp[k+1]; p++) {
      if (Ci[p] <= k) x[Ci[p]] = Cx[p];
    }
    double d = x[k];
    x[k] = 0.0;

    // Triangular solve
    for ( ; top < n; top++) {
      int i = s[top];
      double lki = x[i] / Lx[Lp[i]];
      x[i] = 0.0;
      for (int p = Lp[i]+1; p < c[i]; p++) x[Li[p]] -= Lx[p]*lki;
      d -= lki*lki;
      int p = c[i]++;
      Li[p] = k;
      Lx[p] = lki;
    }

    // Compute L(k,k
    if (d <= 0.0) {
      delete E;
      delete [] c;
      delete [] x;
      delete N;
      return 0;
    }
    int p = c[k]++;
    Li[p] = k;
    Lx[p] = sqrt(d);
  }
  Lp[n] = cp[n];

  delete E;
  delete [] c;
  delete [] x;

  return N;
 }

 int benny_cholsol(int order, double *b)
 {
  if (!is_csc() || b == 0) return -1;

  benny_symbolic *S = this->benny_schol(order);
  benny_numeric *N = this->benny_chol(S);
  double *x = new double[n];
  bool ok = S != 0 && N != 0 && x != 0;
  if (ok) {
    benny_ipvec(S->pinv, b, x, n); // x = P*b
    N->L->benny_lsolve(x); // x = L\x
    N->L->benny_ltsolve(x); // x= L'\x
    benny_pvec(S->pinv, x, b, n); // b = P'*x
  }

  delete [] x;
  delete S;
  delete N;

  return ok ? 0 : -1;
 }
};

// Down here to get rid of a compiler warning
benny_numeric::~benny_numeric() {
    if (L != 0) delete L;
    if (U != 0) delete U;
    if (pinv != 0) delete [] pinv;
    if (B != 0) delete [] B;
}
#endif
