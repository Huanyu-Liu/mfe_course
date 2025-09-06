//
//  matrix_cal.cpp
//  Project5
//
//  Created by Huanyu Liu on 2/18/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include "matrix_cal.hpp"

double * matrix_cal::addition(double * a, double * b, int k){
    double * result = new double[k];
    for (int i = 0; i < k; ++i){
        result[i] = a[i] + b[i];
    }
    return result;
}

double * matrix_cal::matrix_multi(double **a, double * b, int k){
    double * result = new double[k];
    for (int i = 0; i != k; ++i){
        result[i] = 0;
        for (int j = 0; j != k; ++j){
            result[i] += a[i][j] * b[j];
        }
    }
    return result;
}


double **matrix_cal::LUPDecompose(double **A, int * P, int N, double Tol) {
    
    int i, j, k, imax;
    double maxA, absA, *ptr;
    //double * temp;
    double **updated_A = new double *[N];
    for (int i = 0; i != N; ++i){
        updated_A[i] = new double[N];
        for (int j = 0; j != N; ++j){
            updated_A[i][j] = A[i][j];
        }
    }
    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    
    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;
        
        for (k = i; k < N; k++)
            if ((absA = std::abs(updated_A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }
        
        //if (maxA < Tol) return; //failure, matrix is degenerate
        
        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            //pivoting rows of A
            ptr = updated_A[i];
            //temp = A[i];
            updated_A[i] = updated_A[imax];
            updated_A[imax] = ptr;
            
            //counting pivots starting from N (for determinant)
            P[N]++;
        }
        
        for (j = i + 1; j < N; j++) {
            //std::cout << updated_A[i][i] << std::endl;
            updated_A[j][i] = updated_A[j][i] / updated_A[i][i];
            
            for (k = i + 1; k < N; k++)
                updated_A[j][k] = updated_A[j][k] - updated_A[j][i] * updated_A[i][k];
        }
    }
    return updated_A;
    //decomposition done
}
///* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
// * OUTPUT: x - solution vector of A*x=b
// */
double * matrix_cal::LUPSolve(double **A, int *P, double * b, int N) {
    double *result = new double[N];
    for (int i = 0; i < N; i++) {
        result[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            result[i] -= A[i][k] * result[k];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            result[i] -= A[i][k] * result[k];

        result[i] /= A[i][i];
    }
    return result;
}
//
double * matrix_cal::linear_equation(double **A, double * b, int N){
    int * p = new int[N];
    double **updated_A = LUPDecompose(A, p, N, SMALL);
    double * result = LUPSolve(updated_A, p, b, N);
    delete [] p;
    for (int i = 0; i != N; ++i){
        delete [] updated_A[i];
    }
    delete [] updated_A;
    return result;
}
