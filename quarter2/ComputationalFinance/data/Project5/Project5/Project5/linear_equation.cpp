//
//  linear_equation.cpp
//  Project5
//
//  Created by Huanyu Liu on 2/21/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include "linear_equation.hpp"

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */
int * LUPDecompose(double **A, int N, double Tol) {
    
    int * P = new int[N];
    int i, j, k, imax;
    double maxA, *ptr, absA;
    
    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    
    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;
        
        for (k = i; k < N; k++)
            if ((absA = std::abs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }
        
        if (maxA < Tol) return 0; //failure, matrix is degenerate
        
        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;
            
            //counting pivots starting from N (for determinant)
            P[N]++;
        }
        
        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];
            
            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    
    return P;  //decomposition done
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
double * LUPSolve(double **A, int *P, double *b, int N) {
    double * x = new double[N];
    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];
        
        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }
    
    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];
        
        x[i] = x[i] / A[i][i];
    }
    return x;
}
