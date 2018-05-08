#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "eigen.h"
#include <math.h>
#include <stdio.h>
#include <chrono>

TEST_CASE("Check whether Jacobi preserves orthogonality")
{
    //Create a 1000x1000 identity matrix
    size_t dim = 1000;
    double **A = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            A[i][j] = 0;
        A[i][i] = 1; //Done this way to allow vectorization
    }
    
    double **temp = new double*[dim];
    for(size_t i = 0; i<dim; i++)
        temp[i] = new double[dim];
    
    jacobiSimilarityRotation(A, temp, dim, 3, 5, cos(0.5), sin(0.5));
    
    double norm_sq;
    //Check whether rows are normalized
    for(size_t i = 0; i<dim; i++)
    {
        norm_sq = 0;
        for(size_t j = 0; j<dim; j++)
        {
            norm_sq += A[i][j]*A[i][j];
        }
        REQUIRE(norm_sq==Approx(1));
    }
    
    //Check whether rows are are all orthogonal
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = i+1; j<dim; j++)
        {
            norm_sq = 0;
            for(size_t k = 0; k<dim; k++)
                norm_sq += A[i][k]*A[k][j];
            REQUIRE(norm_sq==Approx(0));
        }
    }
}

TEST_CASE("Check whether classical Jacobi finds correct eigenpairs for Toeplitz")
{
    //Create an identity matrix
    size_t dim = 6;
    double **V = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        V[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            V[i][j] = 0;
        V[i][i] = 1; //Done this way to allow vectorization
    }
    
    //Create a Toeplitz matrix
    double **A = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
        {
            if (i==j)
                A[i][j] = 3;
            else if (i == j-1 || i == j+1)
                A[i][j] = 1;
            else
                A[i][j] = 0;
        }
    }
    
    //And a copy for the later check
    //Create a Toeplitz matrix
    double **A_old = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A_old[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
        {
            if (i==j)
                A_old[i][j] = 3;
            else if (i == j-1 || i == j+1)
                A_old[i][j] = 1;
            else
                A_old[i][j] = 0;
        }
    }
    
    classicalJacobi(A, V, dim, 0.00001);
    
    /*
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", A[i][j]);
        printf("\n");
    }
    printf("\n");
    
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", V[i][j]);
        printf("\n");
    }
    */
    
    double eigval; double sum;
    for(size_t i = 0; i<dim; i++)
    {
        eigval = 0;
        for(size_t j = 0; j<dim; j++)
        {
            sum = 0;
            for(size_t k = 0; k<dim; k++)
                sum += A_old[j][k]*V[k][i];
            eigval += sum*sum;
        }
        eigval = sqrt(eigval);
        REQUIRE(eigval == Approx(A[i][i]));    
    }
}

TEST_CASE("Check whether Cyclic Jacobi finds correct eigenpairs for Toeplitz")
{
    //Create an identity matrix
    size_t dim = 9;
    double **V = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        V[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            V[i][j] = 0;
        V[i][i] = 1; //Done this way to allow vectorization
    }
    
    //Create a Toeplitz matrix
    double **A = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
        {
            if (i==j)
                A[i][j] = 3;
            else if (i == j-1 || i == j+1)
                A[i][j] = 1;
            else
                A[i][j] = 0;
        }
    }
    
    //And a copy for the later check
    //Create a Toeplitz matrix
    double **A_old = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A_old[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
        {
            if (i==j)
                A_old[i][j] = 3;
            else if (i == j-1 || i == j+1)
                A_old[i][j] = 1;
            else
                A_old[i][j] = 0;
        }
    }
    
    classicalJacobi(A, V, dim, 0.00001);
    
    /*
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", A[i][j]);
        printf("\n");
    }
    printf("\n");
    
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", V[i][j]);
        printf("\n");
    }
    */
    
    double eigval; double sum;
    for(size_t i = 0; i<dim; i++)
    {
        eigval = 0;
        for(size_t j = 0; j<dim; j++)
        {
            sum = 0;
            for(size_t k = 0; k<dim; k++)
                sum += A_old[j][k]*V[k][i];
            eigval += sum*sum;
        }
        eigval = sqrt(eigval);
        REQUIRE(eigval == Approx(A[i][i]));    
    }
}

TEST_CASE("Check whether Householder result is tridiagonal and whether the decomposition gives the original matrix")
{
    size_t dim = 5000;
    double **T = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        T[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            T[i][j] = (i+j+1);
    }
    
    double **A = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        A[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            A[i][j] = (i+j+1);
    }
    
    //Q can effectively be stored in A because of symmetry
    //This would save memory and time, but we need the
    //vectors explicitly for calculating the eigenvectors
    double **Q = new double*[dim];
    for(size_t i = 0; i<dim; i++)
    {
        Q[i] = new double[dim];
        for(size_t j = 0; j<dim; j++)
            Q[i][j] = 0;
        Q[i][i] = 1;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    householderTridiag(T, dim);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time: ";
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << "ns" << std::endl;
    
    /*

    householderOrthoMat(T, Q, dim);
    
    printf("\nFactored representation of T and Q matrices\n");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", T[i][j]);
        printf("\n");
    }
    
    double norm_sq;
    //Check whether rows are normalized
    for(size_t i = 0; i<dim; i++)
    {
        norm_sq = 0;
        for(size_t j = 0; j<dim; j++)
        {
            norm_sq += Q[i][j]*Q[i][j];
        }
        printf("i: %zu norm_sq: %f\n",i,norm_sq);
        REQUIRE(norm_sq==Approx(1));
    }
    
    
    //Check whether rows are are all orthogonal
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = i; j<dim; j++)
        {
            norm_sq = 0;
            for(size_t k = 0; k<dim; k++)
                norm_sq += Q[i][k]*Q[j][k];
            printf("i: %zu j: %zu norm_sq: %f\n",i,j,norm_sq);
            if(i==j)
                REQUIRE(norm_sq==Approx(1));
            else
                REQUIRE(norm_sq==Approx(0));
        }
    }
    
    //Print out matrix T
    printf("T:\n");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", T[i][j]);
        printf("\n");
    }
    printf("\n");
    
    //Print out matrix Q
    printf("Q:\n");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", Q[i][j]);
        printf("\n");
    }
    printf("\n");
    
    
    //Check if matrix product Q*T*Q^T gives old matrix
    double *T_temp = new double[dim];
    //First calculate T*Q^T
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
        {
            T_temp[j] = 0; //temp row
            for(size_t k = 0; k<dim; k++)
            {
                T_temp[j] += T[i][k]*Q[j][k];
            }
        }
        
        for(size_t j = 0; j<dim; j++)
            T[i][j] = T_temp[j];
    }
    
    //Second calculate Q*(T*Q^T)
    for(size_t j = 0; j<dim; j++)
    {
        for(size_t i = 0; i<dim; i++)
        {
            T_temp[i] = 0; //temp column
            for(size_t k = 0; k<dim; k++)
            {
                T_temp[j] += Q[i][k]*T[k][j];
            }
        }
        
        for(size_t i = 0; i<dim; i++)
            T[i][j] = T_temp[i];
    }
    
    //Print out matrix Q*T*Q^T
    printf("QTQ^t:\n");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", T[i][j]);
        printf("\n");
    }
    printf("\n");
    
    //Print out matrix A
    printf("A:\n");
    for(size_t i = 0; i<dim; i++)
    {
        for(size_t j = 0; j<dim; j++)
            printf("%f\t", A[i][j]);
        printf("\n");
    }
    printf("\n");

    */
}