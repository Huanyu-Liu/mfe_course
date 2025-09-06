//
//  AmericanOption.cpp
//  Project7c
//
//  Created by Huanyu Liu on 3/20/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include "AmericanOption.hpp"

AmericanOption::AmericanOption(double s0, double r, double sigma, double k, double t){
    this->s0 = s0;
    this->r = r;
    this->sigma = sigma;
    this->k = k;
    this->t = t;
}

vector<double> AmericanOption::stock(int N, double delta_s){
    vector<double> stock_price;
    for (int i = 0; i < 2 * N + 1; ++i){
        stock_price.push_back(s0 + (N - i) * delta_s);
    }
    return stock_price;
}

vector<double> AmericanOption::american_efd(double delta, double delta_s, bool is_call){
    vector<double> F;
    int N = (int)(s0 / delta_s);
    int dimension = 2 * N + 1;
    int m = (int)(t / delta);
    //double * stock_price = new double[dimension];
    vector<double> stock_price = stock(N, delta_s);
    vector<double> B(dimension);
    vector<double> exercise;
    if (is_call){
        B[0] = stock_price[0] - stock_price[1];
        for (int i = 0; i < dimension; ++i){
            exercise.push_back(stock_price[i] - k > 0 ? stock_price[i] - k : 0);
            F.push_back(exercise[i]);
        }
    }
    else{
        B[dimension - 1] = stock_price[dimension - 2] - stock_price[dimension - 1];
        for (int i = 0; i < dimension; ++i){
            exercise.push_back(k - stock_price[i] > 0 ? k - stock_price[i] : 0);
            F.push_back(exercise[i]);
            //std::cout << F[i] << std::endl;
        }
    }
    vector<double> pu;
    vector<double> pm;
    vector<double> pd;
    for (int i = 0; i < dimension - 2; ++i){
        double temp1 = r * (i + 1) / 2;
        double temp2 = sigma * sigma * (i + 1) * (i + 1);
        pu.push_back(delta * (temp1 + temp2 / 2));
        pm.push_back(1 - delta * (temp2 + r));
        pd.push_back(delta * (-temp1 + temp2 / 2));
    }
    vector<double> rows(dimension);
    vector<vector<double>> A(dimension, rows);
    A[0][0] = pu[dimension - 3];
    A[0][1] = pm[dimension - 3];
    A[0][2] = pd[dimension - 3];
    A[dimension - 1][dimension - 3] = pu[0];
    A[dimension - 1][dimension - 2] = pm[0];
    A[dimension - 1][dimension - 1] = pd[0];
    for (int i = 1; i < dimension - 1; ++i){
        A[i][i - 1] = pu[2 * N - 1 - i];
        A[i][i] = pm[2 * N - 1 - i];
        A[i][i + 1] = pd[2 * N - 1 - i];
    }
    //double * ecv = new double[dimension];
    vector<double> ecv;
    for (int i = 0; i < m; ++i){
        ecv = matrix_cal::matrix_multi(A, F, dimension);
        for (int j = 0; j < dimension; ++j){
            //std::cout << B[j] << std::endl;
            ecv[j] += B[j];
            F[j] = exercise[j] > ecv[j] ? exercise[j] : ecv[j];

        }
    }
    return F;
}

vector<double> AmericanOption::american_ifd(double delta, double delta_s, bool is_call){
    vector<double> F;
    int N = (int)(s0 / delta_s);
    int dimension = 2 * N + 1;
    int m = (int)(t / delta);
    vector<double> stock_price = stock(N, delta_s);
//    double * B = new double[dimension]();
//    double * exercise = new double[dimension];
    vector<double> B(dimension);
    vector<double> exercise;
    if (is_call){
        B[0] = stock_price[0] - stock_price[1];
        for (int i = 0; i < dimension; ++i){
            exercise.push_back(stock_price[i] - k > 0 ? stock_price[i] - k : 0);
            F.push_back(exercise[i]);
        }
    }
    else{
        B[dimension - 1] = stock_price[dimension - 2] - stock_price[dimension - 1];
        for (int i = 0; i < dimension; ++i){
            exercise.push_back(k - stock_price[i] > 0 ? k - stock_price[i] : 0);
            F.push_back(exercise[i]);
            //std::cout << F[i] << std::endl;
        }
    }

    vector<double> a1, a2, a3;
    for (int i = 0; i < dimension - 2; ++i){
        double temp1 = r * (i + 1) / 2;
        double temp2 = sigma * sigma * (i + 1) * (i + 1);
        a1.push_back(-delta * (temp1 + temp2 / 2));
        a2.push_back(1 + delta * (temp2 + r));
        a3.push_back(-delta * (-temp1 + temp2 / 2));
    }
    vector<double> rows(dimension);
    vector<vector<double>> A(dimension, rows);
    A[0][0] = 1;
    A[0][1] = -1;
    A[dimension - 1][dimension - 2] = 1;
    A[dimension - 1][dimension - 1] = -1;
    for (int i = 1; i < dimension - 1; ++i){
        A[i][i - 1] = a1[2 * N - 1 - i];
        A[i][i] = a2[2 * N - 1 - i];
        A[i][i + 1] = a3[2 * N - 1 - i];
    }
    vector<double> ecv;
    for (int i = 0; i < m; ++i){
        for (int j = 1; j < dimension - 1; ++j){
            B[j] = F[j];
        }
        ecv = matrix_cal::linear_equation(A, B, dimension);

        for (int j = 0; j < dimension; ++j){
            //std::cout << B[j] << std::endl;
            F[j] = exercise[j] > ecv[j] ? exercise[j] : ecv[j];

        }
    }
    return F;
}

vector<double> AmericanOption::american_cnfd(double delta, double delta_s, bool is_call){
    vector<double> F;
    int N = (int)(s0 / delta_s);
    int dimension = 2 * N + 1;
    int m = (int)(t / delta);
    vector<double> stock_price;
    stock_price = stock(N, delta_s);
    vector<double> B(dimension);
    vector<double> exercise;
    if (is_call){
        B[0] = stock_price[0] - stock_price[1];
        for (int i = 0; i < dimension; ++i){
            exercise.push_back(stock_price[i] - k > 0 ? stock_price[i] - k : 0);
            F.push_back(exercise[i]);
        }
    }
    else{
        B[dimension - 1] = stock_price[dimension - 2] - stock_price[dimension - 1];
        for (int i = 0; i < dimension; ++i){
            exercise.push_back(k - stock_price[i] > 0 ? k - stock_price[i] : 0);
            F.push_back(exercise[i]);
            //std::cout << F[i] << std::endl;
        }
    }
    //delete [] stock_price;
    vector<double> a1, a2, a3, b1, b2, b3;
    for (int i = 0; i < dimension - 2; ++i){
        double temp1 = r * (i + 1);
        double temp2 = sigma * sigma * (i + 1) * (i + 1);
        a1.push_back(-delta * (temp1 + temp2) / 4);
        a2.push_back(1 + delta * (temp2 + r) / 2);
        a3.push_back(-delta * (-temp1 + temp2) / 4);
        b1.push_back(-a1[i]);
        b2.push_back(2 - a2[i]);
        b3.push_back(-a3[i]);
    }

    vector<double> rows(dimension);
    vector<vector<double>> A(dimension, rows);
    A[0][0] = 1;
    A[0][1] = -1;
    A[dimension - 1][dimension - 2] = 1;
    A[dimension - 1][dimension - 1] = -1;
    for (int i = 1; i < dimension - 1; ++i){
        A[i][i - 1] = a1[2 * N - 1 - i];
        A[i][i] = a2[2 * N - 1 - i];
        A[i][i + 1] = a3[2 * N - 1 - i];
    }
    vector<double> ecv;
    for (int i = 0; i < m; ++i){
        for (int j = 1; j < dimension - 1; ++j){
            B[j] = F[j - 1] * b1[j - 1] + F[j] * b2[j - 1] + F[j + 1] * b3[j - 1];
        }
        ecv = matrix_cal::linear_equation(A, B, dimension);

        for (int j = 0; j < dimension; ++j){
            //std::cout << B[j] << std::endl;
            F[j] = exercise[j] > ecv[j] ? exercise[j] : ecv[j];

        }
    }
    return F;
}
