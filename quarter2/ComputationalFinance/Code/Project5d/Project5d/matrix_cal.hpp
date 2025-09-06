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
#include <iostream>
#endif /* matrix_cal_hpp */

using std::vector;
namespace matrix_cal {
    double * matrix_multi(double **a, double * b, int k);
    double * addition(double * a, double * b, int k);
    double **LUPDecompose(double **A, int * P, int N, double Tol);
    double * LUPSolve(double **A, int *P, double * b, int N);
    double * linear_equation(double **A, double * b, int N);
}
