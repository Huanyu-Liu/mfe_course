//
//  AmericanOption.hpp
//  Project7c
//
//  Created by Huanyu Liu on 3/20/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#ifndef AmericanOption_hpp
#define AmericanOption_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include "matrix_cal.hpp"
#endif /* AmericanOption_hpp */

using std::vector;

class AmericanOption{
    double s0, r, sigma, k, t;
public:
    AmericanOption(double s0, double r, double sigma, double k, double t);
    vector<double> efd(int N, double delta, double delta_x);
    vector<double> ifd(int N, double delta, double delta_x);
    vector<double> cnfd(int N, double delta, double delta_x);
    vector<double> log_stock(int N, double delta_x);
    vector<double> stock(int N, double delta_s);
    vector<double> american_efd(double delta, double delta_s, bool is_call);
    vector<double> american_ifd(double delta, double delta_s, bool is_call);
    vector<double> american_cnfd(double delta, double delta_s, bool is_call);
};
