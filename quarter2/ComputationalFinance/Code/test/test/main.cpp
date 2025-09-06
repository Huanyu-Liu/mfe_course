//
//  main.cpp
//  test
//
//  Created by Huanyu Liu on 3/20/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include <iostream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;

vector<double> test(vector<double> v){
    v[2] = 5;
    return v;
}

int main(int argc, const char * argv[]) {
    vector<double> b(5);
    vector<vector<double>> a(5, b);
//    for (int i = 0; i != 5; ++i){
//        for (int j = 0; j != 5; ++j){
//            cout << a[i][j] << endl;
//        }
//    }
    b[3] = 4;
    
    double A[4][4] = {{1,1,1,1},{2,3,0,-1},{-3,4,1,2},{1,2,-1,1}};
    //    vector<vector<double>> B = {{1,1,1,1},{2,3,0,-1},{-3,4,1,2},{1,2,-1,1}};
    //    vector<double> b = {13,-1,10,1};
    //double a[4];
    double ** B = new double * [4];
    for (int i = 0; i < 4; ++i){
        B[i] = A[i];
    }
    A[0][0] = 5;
    cout << B[0][0] << endl;
    
//    vector<double> result = test(b);
//    cout << b[2] << endl;
//    cout << result[2] << endl;
    return 0;
}
