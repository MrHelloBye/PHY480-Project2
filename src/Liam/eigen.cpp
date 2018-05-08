#include <cstddef>
#include <cmath>
#include "eigen.h"
#include <stdio.h>

using namespace std;

void symSchur(double **A, const size_t &p, const size_t &q,
    double &c, double &s)
{
    if (A[p][q] != 0)
    {
        double tau = (A[q][q]-A[p][p])/(2*A[p][q]);
        double t;
        
        if (tau >= 0)
            t = 1/(tau+sqrt(1+tau*tau));
        else
            t = 1/(tau-sqrt(1+tau*tau));
        
        c = 1/sqrt(1+t*t);
        s = t*c;
    }
    
    else
    {
        c = 1;
        s = 0;
    }
}

void classicalJacobi(double **A, double **V,
    const size_t &dim, const double &tol)
{
    //Initialize orthogonal matrix V to Identity
    for (int i = 0; i<dim; i++)
    {
        for (int j = 0; j<dim; j++)
        {
            if (i==j)
                V[i][j] = 1;
            else
                V[i][j] = 0;
        }
    }
    
    size_t p; size_t q;
    double c; double s;
    
    double **temp = new double*[dim];
    for(size_t i = 0; i<dim; i++)
        temp[i] = new double[dim];
    
    int counter = 0;
    double max = maxOff(A, dim, p, q);
    while ( max >= tol)
    {
        symSchur(A, p, q, c, s);
        jacobiSimilarityRotation(A, temp, dim, p, q, c, s);
        jacobiRightRotation(V, temp, dim, p, q, c, s);
        
        max = maxOff(A, dim, p, q);
        printf("%d %f\n", counter, max);
        counter++;
    }
}

void cyclicJacobi(double **A, double **V, const size_t &dim, const double &tol)
{
    //Initialize orthogonal matrix V to Identity
    for (int i = 0; i<dim; i++)
    {
        for (int j = 0; j<dim; j++)
        {
            if (i==j)
                V[i][j] = 1;
            else
                V[i][j] = 0;
        }
    }
    
    size_t p; size_t q;
    double c; double s;
    
    double **temp = new double*[dim];
    for(size_t i = 0; i<dim; i++)
        temp[i] = new double[dim];
    
    int counter = 0;
    double max = maxOff(A, dim, p, q);
    while ( fabs(max) >= tol)
    {
        max = 0.0;
        for(p = 1; p<dim; p++)
        {
            for(q = p+1; q<dim; q++)
            {
                if(fabs(A[p][q]) > max)
                    max = A[p][q];
                symSchur(A, p, q, c, s);
                jacobiSimilarityRotation(A, temp, dim, p, q, c, s);
                jacobiRightRotation(V, temp, dim, p, q, c, s);
            }
            printf("p: %d q: %d\n", p, q);
        }
        
        
        printf("Pass: %d Max: %f\n", counter, max);
        counter++;
    }
}

double maxOff(double **A, const size_t &dim, size_t &p, size_t &q)
{
    double max = 0;
    
    for (int i = 0; i<dim; i++)
    {
        for (int j = i+1; j<dim; j++)
        {
            if (fabs(A[i][j]) > max)
            {
                max = fabs(A[i][j]);
                p = i;
                q = j;
            }
        }
    }
    
    return max;
}

void jacobiSimilarityRotation(double **A, double **temp, const size_t &dim,
    const size_t &p, const size_t &q, const double &c, const double &s)
{
    jacobiLeftRotation(A, temp, dim, p, q, c, s);
    jacobiRightRotation(A, temp, dim, p, q, c, s);
}

void jacobiRightRotation(double **V, double **temp, const size_t &dim,
    const size_t &p, const size_t &q, const double &c, const double &s)
{
    //Calculate new values and store in temp array
    for (int i = 0; i<dim; i++)
    {
        temp[i][p] = c*V[i][p] - s*V[i][q];
        temp[i][q] = s*V[i][p] + c*V[i][q];
    }
    
    //Update old array
    for (int i = 0; i<dim; i++)
    {
        V[i][p] = temp[i][p];
        V[i][q] = temp[i][q];
    }
}

void jacobiLeftRotation(double **V, double **temp, const size_t &dim,
    const size_t &p, const size_t &q, const double &c, const double &s)
{
    //Calculate new values and store in temp array
    for (int j = 0; j<dim; j++)
    {
        temp[p][j] = c*V[p][j] - s*V[q][j];
        temp[q][j] = s*V[p][j] + c*V[q][j];
    }
    
    //Update old array
    for (int j = 0; j<dim; j++)
    {
        V[p][j] = temp[p][j];
        V[q][j] = temp[q][j];
    }
}

void householderVector(double *x, size_t const &dim, double *v, double &b)
{
    double s = 0;
    v[0] = 1;
    for(size_t i = 1; i<dim; i++)
    {
        s += x[i]*x[i];
        v[i] = x[i];
    }
    
    // b = 2/mag(v)^2
    if (s == 0 && x[0] >= 0)
        b = 0;
    else if (s == 0 && x[0] <0)
        b = -2;
    else
    {
        double m = sqrt(x[0]*x[0]+s);
        if (x[0] <= 0)
            v[0] = x[0] - m;
        else
            v[0] = -s/(x[0]+m);
        b = 2*v[0]*v[0]/(s+v[0]*v[0]);
        
        for(size_t i = dim-1; i<dim; i--)
            v[i] = v[i]/v[0];
    }
}

void householderTridiag(double **A, size_t const &dim)
{
    double *v = new double[dim];
    double b;
    double *p = new double[dim];
    size_t sub_dim; // the length of the householder vector
    double pTv;
    double *w = new double[dim];
    double sub_norm;
    
    /*
    printf("\nInitial Matrix\n");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", A[i][j]);
        printf("\n");
    }
    printf("\n");
    */

    for (size_t k = 0; k<dim-2; k++)
    {
        //sub_dim is the length of the householder vector, it is length
        // dim-(k+1) because we are zeroing all elements below the subdiagonal
        sub_dim = dim-(k+1);
        householderVector(A[k]+(k+1), sub_dim, v, b);
        
        /*
        printf("dim: %zu, sub_dim: %zu\n",dim,sub_dim);
        printf("Beta: %f\n", b);
        printf("Vector: ");
        for(size_t i = 0; i<sub_dim; i++)
            printf("%f\t",v[i]);
        printf("\n");
        */
        
        
        //Calculate p vector
        for(size_t i = 0; i<sub_dim; i++)
        {
            p[i] = 0;
            for(size_t j = 0; j<sub_dim; j++)
            {
                p[i] += b*A[k+1+i][k+1+j]*v[j];
                //printf("i: %zu, k: %zu, k+1+j: %lu, j: %zu",i,k,k+1+j,j);
                //printf("\tA[k][k+1+j]: %f, v[j]: %f\n", A[k+1][k+1+j], v[j]);
            }
        }
        
        //Calculate pTv
        pTv = 0;
        for(size_t i = 0; i<sub_dim; i++)
            pTv += p[i]*v[i];
        
        //Calculate w vector
        for(size_t i = 0; i<sub_dim; i++)
            w[i] = p[i] - (b*pTv/2)*v[i];
        
        //Calculate superdiagonal and zero other off diagonals
        sub_norm = A[k][k+1]*A[k][k+1];
        for(size_t i = k+2; i<dim; i++)
        {
            //This could be reused from householderVector()
            sub_norm += A[k][i]*A[k][i];
            A[k][i] = 0;
            A[i][k] = v[i-k-1];
        }
            
        sub_norm = sqrt(sub_norm);
        A[k][k+1] = sub_norm;
        A[k+1][k] = b;
        
        for(size_t i = 0; i<sub_dim; i++)
        {
            for(size_t j = 0; j<sub_dim; j++)
                A[k+1+i][k+1+j] -= v[i]*w[j] + v[j]*w[i];
        }
        
        /*
        for(size_t i = 0; i<dim; i++)
        {
            for(size_t j = 0; j<dim; j++)
                printf("%f\t", A[i][j]);
            printf("\n");
        }
        printf("\n");
        */
    }
}

//Should change to c++11 style array template instead of passing dim
void householderOrthoMat(double **T, double **Q, size_t const &dim)
{
    double *QTempColumn = new double[dim];
    double beta = 0;
    double *v = new double[dim];
    
    //Initialize beta and v from factored form
    beta = T[dim-2][dim-3];
    for(int i = 0; i<dim; i++)
    {
        if(i < dim-2)
            v[i] = 0;
        else if(i == dim-2)
            v[i] = 1;
        else
            v[i] = T[i][dim-3];
    }
    
    //Start Q matrix calculation with
    // smallest Householder vector
    
    for(int i = dim-2; i<dim; i++)
    {
        for(int j = dim-2; j<dim; j++)
        {
            Q[i][j] = -beta*v[i]*v[j];
        }
        Q[i][i] += 1;
    }
    /*
    for(int i = 0; i<dim; i++)
    {
        for(int j = 0; j<dim; j++)
            Q[i][j] = 0;
        Q[i][i] = 1;
    }
    */
        
    //Do matrix multiplication to calculate Q
    printf("Initial Q matrix\n");
    for(int i = 0; i<dim; i++)
    {
        for(int j = 0; j<dim; j++)
        {
            printf("%f\t",Q[i][j]);
        }
        printf("\n");
    }
    
    for(int k = dim-3; k>0; k--)
    {
        printf("k: %d\n",k);
        beta = T[k][k-1];
        v[k] = 1;
        double temp;
        for(int i = k+1; i<dim-(k-1); i++)
        {
            printf("i: %d k: %d dim: %d\n", i,k,dim);
            v[i] = T[k-1+i][k-1];
            //printf("v[%d]: %f\n",i,v[i]);
        }
        
        for(int j = k; j<dim; j++)
        {
            if(j == k)
            {
                for(int i = k; i<dim; i++)
                    Q[i][j] = -beta*v[i]*v[j];
                Q[k][k] += 1;
            }
            else
            {
                for(int i = k; i<dim; i++)
                {
                    QTempColumn[i] = Q[i][j];
                    //printf("QTempColumn[%d]: %f\n",i,QTempColumn[i]);
                    Q[i][j] = 0;
                }
            
                for(int i = k; i<dim; i++)
                {
                    for(int l = k+1; l<dim; l++)
                    {
                        //printf("i: %d l: %d j: %d QTempColumn[l]: %f\n",
                        //i, l, j, QTempColumn[l]);
                        if(i==l)
                            Q[i][j] += (1-beta*v[i]*v[l])*QTempColumn[l];
                        else
                            Q[i][j] += -beta*v[i]*v[l]*QTempColumn[l];
                    }
                }
            }
        }
        
        printf("Q matrix\n");
        for(int i = 0; i<dim; i++)
        {
            for(int j = 0; j<dim; j++)
                printf("%f\t",Q[i][j]);
            printf("\n");
        }
    }
    
    // Symmetrize T matrix
    for(int i = 0; i<dim; i++)
    {
        for(int j = 0; j<i-1; j++)
        {
            T[i][j] = 0;
        }
    }
    for(int i = 0; i<dim-1; i++)
        T[i+1][i] = T[i][i+1];
}