#include "benny_sparse.hpp"

int main (int argc, char **argv)
{
    benny_sparse T(3,3,1);
    T.benny_entry(0,0,1.0);
    T.benny_entry(1,1,4.0);
    T.benny_entry(1,2,2.0);
    T.benny_entry(2,1,2.0);
    T.benny_entry(2,2,4.0);

    T.benny_print(cout);

    benny_sparse *Tc = T.benny_compress();

    double x[3] = {0,0,0};
    x[1] = 1;
    
    if (Tc->benny_cholsol(1,x)) {
        std::cerr << "failed to solve Ax=b by cholesky decomposition" << std::endl;
        return -1;
    }
    for (int i = 0; i < 3; i++) std::cout << x[i] << std::endl;

    delete Tc;

    return 0;
}