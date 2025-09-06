//
//  matrix_cal.hpp
//  Project5
//
//  Created by Huanyu Liu on 2/18/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#ifndef matrix_cal_hpp
#define matrix_cal_hpp
#define SMALL 0.000000001

#pragma once

#include <stdio.h>
#include <cmath>
#include <vector>
#endif /* matrix_cal_hpp */

using std::vector;
namespace matrix_cal {
    vector<double> matrix_multi(vector<vector<double>> a, vector<double> b, int k);
    vector<double> addition(vector<double> a, vector<double> b, int k);
    vector<vector<double>> LUPDecompose(vector<vector<double>> A, int * P, int N, double Tol);
    vector<double> LUPSolve(vector<vector<double>> A, int *P, vector<double> b, int N);
    vector<double> linear_equation(vector<vector<double>> A, vector<double> b, int N);
}
