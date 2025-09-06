//
//  main.cpp
//  Project5c
//
//  Created by Huanyu Liu on 3/20/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "option.hpp"

using std::cout;
using std::endl;
using std::cin;

int main(int argc, const char * argv[]) {
    // insert code here...
    int size = 100;
    int path_count = 100000;
    double r = 0.06;
    double s0 = 36;
    double sigma = 0.2;
    double t = 0.5;
    int k;
    int method = 0;
    double x = 40;
    
    // Problem 1
    // (a)
    
    cout << "Problem1" << endl << endl;
    
    option op1(s0,r,sigma,x,t);
    cout << "Laguerre polynomial:" << endl;
//    k = 2;
//    cout << "s0 = 36, t = 0.5, k = 2: " << op1.price(k, method, path_count, size) << endl;
    k = 3;
    cout << "s0 = 36, t = 0.5, k = 3: " << op1.price(k, method, path_count, size) << endl;
    k = 4;
    cout << "s0 = 36, t = 0.5, k = 4: " << op1.price(k, method, path_count, size) << endl;
    return 0;
}
