#ifndef BENNY_UTILITY_HPP
#define BENNY_UTILITY_HPP

double benny_cumsum(int *p, int *c, int n)
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
 
int benny_pvec(const int *p, const double *b, double *x, int n)
{
  if (x == 0 || b == 0) return -1;
  for (int k = 0; k < n; k++) x[k] = b[p ? p[k] : k];
  return 0;
}

int benny_ipvec(const int *p, const double *b, double *x, int n)
{
  if (x == 0 || b == 0) return -1;
  for (int k = 0; k < n; k++) x[p ? p[k] : k] = b[k];
  return 0;
}

int* benny_pinv(const int *p, int n)
{
  int *pinv = new int[n];
  for (int k = 0; k < n; k++) pinv[p[k]] = k;
  return pinv;
}

#endif