//
//  matrix_cal.cpp
//  Project5
//
//  Created by Huanyu Liu on 2/18/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include "matrix_cal.hpp"

vector<double> matrix_cal::addition(vector<double> a, vector<double> b, int k){
    vector<double> result;
    for (int i = 0; i < k; ++i){
        result.push_back(a[i] + b[i]);
    }
    return result;
}

vector<double> matrix_cal::matrix_multi(vector<vector<double>> a, vector<double> b, int k){
    vector<double> result(k);
    for (int i = 0; i != k; ++i){
        result[i] = 0;
        for (int j = 0; j != k; ++j){
            result[i] += a[i][j] * b[j];
        }
    }
    return result;
}


vector<vector<double>> matrix_cal::LUPDecompose(vector<vector<double>> A, int * P, int N, double Tol) {
    
    int i, j, k, imax;
    double maxA, absA;
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
        
        //if (maxA < Tol) return; //failure, matrix is degenerate
        
        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            //pivoting rows of A
            //ptr = A[i];
            vector<double> temp = A[i];
            A[i] = A[imax];
            A[imax] = temp;
            
            //counting pivots starting from N (for determinant)
            P[N]++;
        }
        
        for (j = i + 1; j < N; j++) {
            A[j][i] = A[j][i] / A[i][i];
            
            for (k = i + 1; k < N; k++)
                A[j][k] = A[j][k] - A[j][i] * A[i][k];
        }
    }
    return A;
    //decomposition done
}
///* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
// * OUTPUT: x - solution vector of A*x=b
// */
vector<double> matrix_cal::LUPSolve(vector<vector<double>> A, int *P, vector<double> b, int N) {
    vector<double> result;
    for (int i = 0; i < N; i++) {
        result.push_back(b[P[i]]);

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
vector<double> matrix_cal::linear_equation(vector<vector<double>> A, vector<double> b, int N){
    int * p = new int[N];
    vector<vector<double>> updated_A = LUPDecompose(A, p, N, SMALL);
    vector<double> result = LUPSolve(updated_A, p, b, N);
    delete [] p;
    return result;
}
