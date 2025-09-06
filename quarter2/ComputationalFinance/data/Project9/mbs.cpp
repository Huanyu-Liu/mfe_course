//
//  mbs.cpp
//  Project9
//
//  Created by Huanyu Liu on 3/3/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include "mbs.hpp"

MBS::MBS(double r0, double k, double r_bar, double sigma, double T, double wac, double loan){
    this->r0 = r0;
    this->k = k;
    this->r_bar = r_bar;
    this->sigma = sigma;
    this->T = T;
    this->wac = wac;
    this->loan = loan;
}

double MBS::sum(vector<double> v, int start, int end){
    double sum = 0;
    for (int i = start; i != end; ++i){
        sum += v[i];
    }
    return sum;
}

vector<double> MBS::rf_rate_path(int size, double delta){
    random_generator rg;
    vector<double> bm = rg.brownian_motion(size, delta);
    vector<double> r;
    r.push_back(r0);
    for (int i = 0; i != size - 1; ++i){
        r.push_back(r[i] + k * (r_bar - r[i]) * delta + sigma * sqrt(r[i]) * bm[i]);
    }
    return r;
    
}

double MBS::P_theory(double T, double FV, double r0, double rbar, double sigma, double k){
    double h1 = sqrt(k * k + 2 * sigma * sigma);
    double h2 = (k + h1) / 2;
    double h3 = 2 * k * rbar / (sigma * sigma);
    double B = (exp(h1 * T) - 1) / (h2 * (exp(h1 * T) - 1) + h1);
    double A = pow((h1 * exp(h2 * T)/(h2 * (exp(h1 * T) - 1) + h1)), h3);
    return FV * A * exp(-B * r0);
}


double MBS::CPR(double t, int month, double pvt, double r10){
    double RI = 0.28 + 0.14 * atan(-8.57 + 430 * (wac - r10));
    double BU = 0.3 + 0.7 * pvt / loan;
    double SG = t / 30 < 1 ? t / 30 : 1;
    double SY[12] = {0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.1,1.18,1.22,1.23,0.98};
    return RI * BU * SG * SY[month - 1];
}

double MBS::price(int month, double oas, int simulation_count){
    double delta = 1.0 / 12;
    double pvt;
    double TPP, ct;
    double temp1, temp2;
    double r = wac / 12;
    int mortgage_step = (int)(T * 12);
    double result = 0;
    vector<double> r_path, r10_path;
    random_generator rg;
    for (int simulate = 0; simulate != simulation_count; ++simulate){
        int start_month = month;
        pvt = loan;
        r_path = rg.cir_stock_path(mortgage_step, r0, k, r_bar, delta, sigma);
        double r10;
        double R = 0;
        double sum = 0;
        for (int i = 0; i != mortgage_step; ++i){
            r10 = P_theory(10, 1, r_path[i], r_bar, sigma, k);
            r10 = - log(r10) / 10;
            temp1 = 1 / (1 - pow(1 + r, i - mortgage_step)) - 1;
            temp2 = 1 - pow(1 - CPR(i + 1, start_month, pvt, r10), delta);
            TPP = pvt * r * temp1 + (pvt - pvt * r * temp1) * temp2;
            ct = TPP + pvt * r;
            pvt -= TPP;
            R += (r_path[i] + oas);
            sum += ct * exp(-R * delta);
            start_month = start_month % 12 + 1;
        }
        //std::cout << sum << std::endl;
        r_path.clear();
        result += sum;
    }
    return result / simulation_count;
}
tuple<double, double, double> MBS::oas(int month, double market_price, int simulation_count){
    double duration, convexity;
    //int counter = 0;
    double x = -0.01;
    double y = 0.0005;
    double left = -1.0;
    double right = 0.0;
    double mbs = price(month, x, simulation_count);
    vector<double> r_path;
    double price_plus, price_minus;
    while (abs(mbs - market_price) > SMALL) {
        //++counter;
        //std::cout << counter << std::endl;
        if (mbs < market_price){
            right = x;
        }
        else{
            left = x;
        }
        x = (left + right) / 2;
        mbs = price(month, x, simulation_count);
    }
    double p0 = price(month, x, simulation_count);
    price_plus = price(month, x + y, simulation_count);
    price_minus = price(month, x - y, simulation_count);
    duration = (price_minus - price_plus) / (2 * y * p0);
    convexity = (price_plus + price_minus - 2 * market_price) / (2 * p0 * y * y);
    return std::make_tuple(x, duration, convexity);

}
