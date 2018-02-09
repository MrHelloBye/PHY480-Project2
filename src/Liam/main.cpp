#include <cstdio>
#include <cmath>
#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xeval.hpp"
#include <chrono>
#include <algorithm>

#include <time.h>
double get_cpu_time()
{
    return (double)clock() / CLOCKS_PER_SEC;
}

using namespace std;
using namespace xt;

double error(xtensor<double,1>&, xtensor<double,1>&);
void tridiag_sym(xtensor<double,1>, xtensor<double,1>, xtensor<double,1>&);

int main()
{
    int n = 10000;
    
    xtensor<double, 1> alpha = ones<double>({n});
    xtensor<double, 1> beta = ones<double>({n-1});
    xtensor<double, 1> b = zeros<double>({n});
    xtensor<double, 1> analyt = ones<double>({n});
    
    alpha *= 2;
    beta *= -1;
    
    
    double N = static_cast<double>(n);
    double j = 0;
    for (unsigned long i = 0; i<n; i++)
    {
        j = static_cast<double>(i);
        b.at(i) = 100*exp(-10*j/N);
    }
    
    double t1 = get_cpu_time();
    
    tridiag_sym(alpha, beta, b);
    b /= (n+2);
    b /= (n+2);
    
    double t2 = get_cpu_time();
    
    printf("It took me %f seconds.\n",t2-t1);
    
    /* ----- Data output ----- */
    
    FILE *fp = fopen("solution.tsv","w");
    int mod = 1000;
    
    for(unsigned long i = 0; i<n; i++)
    {
        if (!(i%mod))
            fprintf(fp, "%f\n", b.at(i)); 
    }
    
    fclose(fp);
    return 0;
}


double error(xtensor<double,1> &approx, xtensor<double,1>&analyt)
{
    const unsigned long n = approx.shape()[0];
    if (n != analyt.shape()[0])
        throw "Vector sizes are different.";
    
    xtensor<double, 1> errors = zeros<double>({n});
    
    for (unsigned long i = 0; i<n; i++)
    {
        errors.at(i) = fabs(approx.at(i)-analyt.at(i))/analyt.at(i);
    }
    
    return *max_element(errors.begin(),errors.end());
}

void tridiag_sym(xtensor<double, 1> alpha,
    xtensor<double, 1> beta, xtensor<double, 1> &b)
{
    unsigned long n = alpha.shape()[0];
    
    auto temp = beta.at(0);
    
    for (unsigned int k = 1; k<n; k++)
    {
        temp = beta.at(k-1);
        beta.at(k-1) = temp/alpha.at(k-1);
        alpha.at(k) -= temp*beta.at(k-1);
    }
    
    for (unsigned long k = 1; k<n; k++)
    {
        b.at(k) -= beta.at(k-1)*b.at(k-1);
    }
    
    b.at(n-1) /= alpha.at(n-1);
    
    for (unsigned long k = n-2; k>0; k--)
    {
        b.at(k) = b.at(k)/alpha.at(k) - beta.at(k)*b.at(k+1);
    }
    
    
}