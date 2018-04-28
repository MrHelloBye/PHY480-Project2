#include <cstddef>

#ifndef EIGEN_H
#define EIGEN_H

using namespace std;

void symSchur(double**, const size_t&, const size_t&, double &, double&);

void classicalJacobi(double**, double**, const size_t &, const double &);

void cyclicJacobi(double**, double**, const size_t &, const double &);

double maxOff(double **A, const size_t &dim, size_t &p, size_t &q);

void jacobiSimilarityRotation(double**, double**, const size_t &,
    const size_t &, const size_t &, const double &, const double &);
    
void jacobiRightRotation(double **V, double **temp, const size_t &dim,
    const size_t &p, const size_t &q, const double &c, const double &s);

void jacobiLeftRotation(double **V, double **temp, const size_t &dim,
    const size_t &p, const size_t &q, const double &c, const double &s);

void householderVector(double *x, size_t const &dim, double *v, double &b);

void householderTridiag(double **A, size_t const &dim);

void householderOrthoMat(double **T, double **Q, size_t const &dim);
    
#endif //EIGEN_H