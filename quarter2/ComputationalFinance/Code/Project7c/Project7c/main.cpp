//
//  main.cpp
//  Project7c
//
//  Created by Huanyu Liu on 3/20/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include <iostream>
#include <vector>
#include "Option.hpp"
#include "AmericanOption.hpp"

using std::vector;
using std::cout;
using std::endl;
int main(int argc, const char * argv[]) {
    double s0 = 10;
    double r = 0.04;
    double sigma = 0.2;
    double k = 10;
    double t = 0.5;
    double delta = 0.002;
    double delta_x = sigma * sqrt(4 * delta);
    double delta_s = 0.25;
    int N = 200;
    //int dimension = 2 * N + 1;
    //double x0 = log(s0);
    Option op(s0,r,sigma,k,t);
    vector<double> F = op.efd(N, delta, delta_x);
//    for (int i = 4; i != 17; ++i){
//        double temp = int(round((log(i) - x0) / delta_x));
//
//        cout << F[N - temp] << endl;
//    }
//    F = op.ifd(N, delta, delta_x);
//    //cout << F[N] << endl;
//
//    for (int i = 4; i != 17; ++i){
//        double temp = int(round((log(i) - x0) / delta_x));
//
//        cout << F[N - temp] << endl;
//    }
//    F = op.cnfd(N, delta, delta_x);
//    for (int i = 4; i != 17; ++i){
//        double temp = int(round((log(i) - x0) / delta_x));
//
//        cout << F[N - temp] << endl;
//    }
    
    N = (int)(s0 / delta_s);
    AmericanOption american_op(s0,r,sigma,k,t);
    F = american_op.american_efd(delta, delta_s, false);
    for (int i = 4; i != 17; ++i){
        double temp = int(round((i - s0) / delta_s));
        cout << F[N - temp] << endl;
    }
    F = american_op.american_ifd(delta, delta_s, false);
    //cout << F[N] << endl;
    for (int i = 4; i != 17; ++i){
        double temp = int(round((i - s0) / delta_s));
        cout << F[N - temp] << endl;
    }
    F = american_op.american_cnfd(delta, delta_s, false);
    for (int i = 4; i != 17; ++i){
        double temp = int(round((i - s0) / delta_s));
        cout << F[N - temp] << endl;
    }
    //cout << F[N] << endl;
//    double A[4][4] = {{1,1,1,1},{2,3,0,-1},{-3,4,1,2},{1,2,-1,1}};
//    vector<vector<double>> B = {{1,1,1,1},{2,3,0,-1},{-3,4,1,2},{1,2,-1,1}};
//    vector<double> b = {13,-1,10,1};
    //double a[4];
//    double ** B = new double * [4];
//    for (int i = 0; i < 4; ++i){
//        B[i] = A[i];
//    }
//    vector<double> a = matrix_cal::linear_equation(B, b, 4);
//    for (int i = 0; i < 4; ++i){
//        for (int j = 0; j < 4; ++j){
//            cout << B[i][j] << ", ";
//        }
//        cout << endl;
//    }

//    for (int i = 0; i < 4; ++i){
//        for (int j = 0; j < 4; ++j){
//            cout << A[i][j] << ", ";
//        }
//        cout << endl;
//    }


//    for (int i = 0; i < 4; ++i){
//        cout << a[i] << ", ";
//    }


    return 0;
}
