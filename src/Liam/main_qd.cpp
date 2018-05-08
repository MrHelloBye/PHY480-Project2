#include "eigen.h"
#include <stdio.h>
#include <iostream>
#include <chrono>

int main()
{
    const size_t dim = 500;
    //Initialize matrices to zero and then change
    //diagonals to allow vectorization

    //Create an identity matrix
    double **V = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        V[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            V[i][j] = 0;
        V[i][i] = 1; //Done this way to allow vectorization
    }
    
    //Initialize the matrix to diagonalize
    double **A = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
        {
            A[i][j] = 0;
        }
    }
    
    double rmax = 15;
    double dr = rmax/dim;
    for(size_t i = 0; i<dim-1; i++)
    {
        A[i][i+1] = -1/(dr*dr);
        A[i+1][i] = -1/(dr*dr);
    }
    for(size_t i = 1; i<=dim; i++)
    {
        A[i-1][i-1] = 2/(dr*dr) + 0.5*0.5*(i*dr)*(i*dr) + 1/(i*dr);
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    cyclicJacobi(A, V, dim, 0.01);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << "ns" << std::endl;

    FILE* fp = fopen("QD/A.csv","w");

    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim-1; j++)
            fprintf(fp, "%f,", A[i][j]);
        fprintf(fp, "%f\n", A[i][dim-1]);
    }
    fclose(fp);

    fp = fopen("QD/V.csv","w");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim-1; j++)
            fprintf(fp, "%f,", V[i][j]);
        fprintf(fp, "%f\n", V[i][dim-1]);
    }
    fclose(fp);
    
    return 0;
}
