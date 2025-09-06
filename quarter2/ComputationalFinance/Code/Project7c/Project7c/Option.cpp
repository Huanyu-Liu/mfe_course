//
//  Option.cpp
//  Project7
//
//  Created by Huanyu Liu on 2/28/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include "Option.hpp"

Option::Option(double s0, double r, double sigma, double k, double t){
    this->s0 = s0;
    this->r = r;
    this->sigma = sigma;
    this->k = k;
    this->t = t;
}

vector<double> Option::log_stock(int N, double delta_x){
    double x0 = std::log(s0);
    vector<double> stock_price;
    for (int i = 0; i < 2 * N + 1; ++i){
        stock_price.push_back(x0 + (N - i) * delta_x);
    }
    return stock_price;
}

void Option::stock(double * stock_price, int N, double delta_s){
    for (int i = 0; i < 2 * N + 1; ++i){
        stock_price[i] = s0 + (N - i) * delta_s;
    }
}

vector<double> Option::efd(int N, double delta, double delta_x){
    vector<double> F;
    double temp1 = sigma * sigma / (delta_x * delta_x);
    double temp2 = (r - sigma * sigma / 2) / (2 * delta_x);
    int m = (int)(t / delta);
    double pu = delta * (temp1 / 2 + temp2);
    double pm = 1 - delta * temp1 - r * delta;
    double pd = delta * (temp1 / 2 - temp2);
    
    int dimension = 2 * N + 1;
    vector<double> rows(dimension);
    vector<vector<double>> A(dimension, rows);
    A[0][0] = pu;
    A[0][1] = pm;
    A[0][2] = pd;
    A[dimension - 1][dimension - 3] = pu;
    A[dimension - 1][dimension - 2] = pm;
    A[dimension - 1][dimension - 1] = pd;
    vector<double> B(dimension);
    vector<double> stock_price(dimension);
    //double * stock_price = new double[dimension];
    stock_price = log_stock(N,delta_x);
    B[dimension - 1] = exp(stock_price[dimension - 2]) - exp(stock_price[dimension - 1]);
    for (int i = 1; i != dimension - 1; ++i){
        A[i][i - 1] = pu;
        A[i][i] = pm;
        A[i][i+1] = pd;
    }
    for (int i = 0; i < dimension; ++i){
        double payoff = k - exp(stock_price[i]);
        F.push_back(payoff > 0 ? payoff : 0);
    }
    vector<double> result;
    for (int i = 0; i < m; ++i){
        result = matrix_cal::matrix_multi(A, F, dimension);
        //matrix_cal::matrix_multi(A, F, result, dimension);
        F = matrix_cal::addition(B, result, dimension);
    }
    return F;
}



vector<double> Option::ifd(int N, double delta, double delta_x){
    vector<double> F;
    double temp1 = sigma * sigma / delta_x / delta_x;
    double temp2 = (r - sigma * sigma / 2 ) / delta_x;
    double pu = -0.5 * delta * (temp1 + temp2);
    double pm = 1 + delta * temp1 + r * delta;
    double pd = -0.5 * delta * (temp1 - temp2);
    int m = (int)(t / delta);
    int dimension = 2 * N + 1;
    vector<double> temp(dimension);
    vector<vector<double>> A(dimension, temp);
//    double ** A = new double * [dimension];
//    for (int i = 0; i < dimension; ++i){
//        A[i] = new double[dimension]();
//    }
    A[0][0] = 1;
    A[0][1] = -1;
    A[dimension - 1][dimension - 2] = 1;
    A[dimension - 1][dimension - 1] = -1;
    vector<double> B(dimension);
    vector<double> stock_price;
//    double * B = new double[dimension]();
//    double * stock_price = new double[dimension];
    stock_price = log_stock(N,delta_x);
    B[dimension - 1] = exp(stock_price[dimension - 2]) - exp(stock_price[dimension - 1]);
    for (int i = 1; i < dimension - 1; ++i){
        A[i][i - 1] = pu;
        A[i][i] = pm;
        A[i][i+1] = pd;
    }
    for (int i = 0; i < dimension; ++i){
        double payoff = k - exp(stock_price[i]);
        F.push_back(payoff > 0 ? payoff : 0);
    }
    //int * p = new int[dimension];
    //matrix_cal::LUPDecompose(A, p, dimension, SMALL);
    for (int i = 0; i < m; ++i){
        for (int j = 1; j < dimension - 1; ++j){
            B[j] = F[j];
        }
        F = matrix_cal::linear_equation(A, B, dimension);
    }
    return F;
}

vector<double> Option::cnfd(int N, double delta, double delta_x){
    vector<double> F;
    double temp1 = sigma * sigma / delta_x / delta_x;
    double temp2 = (r - sigma * sigma / 2) / delta_x;
    double pu = -0.25 * delta * (temp1 + temp2);
    double pm = 1 + delta * temp1 / 2 + r * delta / 2;
    double pd = -0.25 * delta * (temp1 - temp2);
    int m = (int)(t / delta);
    int dimension = 2 * N + 1;
    vector<double> temp(dimension);
    vector<vector<double>> A(dimension, temp);
    A[0][0] = 1;
    A[0][1] = -1;
    A[dimension - 1][dimension - 2] = 1;
    A[dimension - 1][dimension - 1] = -1;
    //double * B = new double[dimension]();
    vector<double> B(dimension);
    //double * stock_price = new double[dimension];
    vector<double> stock_price = log_stock(N,delta_x);
    B[dimension - 1] = exp(stock_price[dimension - 2]) - exp(stock_price[dimension - 1]);
    for (int i = 1; i < dimension - 1; ++i){
        A[i][i - 1] = pu;
        A[i][i] = pm;
        A[i][i+1] = pd;
    }
    for (int i = 0; i < dimension; ++i){
        double payoff = k - exp(stock_price[i]);
        F.push_back(payoff > 0 ? payoff : 0);
    }
    for (int i = 0; i < m; ++i){
        for (int j = 1; j < dimension - 1; ++j){
            B[j] = F[j - 1] * -pu - (pm - 2) * F[j] - pd * F[j + 1];
        }
        F = matrix_cal::linear_equation(A, B, dimension);
    }
    return F;
}
