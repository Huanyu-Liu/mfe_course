//
//  main.cpp
//  Final_1
//
//  Created by Huanyu Liu on 3/21/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include <iostream>
using std::cout;
using std::endl;

class random_generator{
    unsigned long seed;
    unsigned long m;
    unsigned long a;
public:
    random_generator(unsigned long x = (1 << 10) + 1);
    double uniform_generator();
};

int main(int argc, const char * argv[]) {
    // insert code here...
    int simulation_count = 10000;
    double x = 1.1;
    double result = 0;
    double temp;
    random_generator rg;
    for (int i = 0; i!= simulation_count; ++i){
        double s = 0;
        int k = 0;
        while (s < x) {
            s += rg.uniform_generator();
            k++;
        }
        temp = 4.54 > k ? 4.54 - k : 0;
        result += temp;
    }
    result /= simulation_count;
    cout << result << endl;
    return 0;
}

random_generator::random_generator(unsigned long x){
    seed = x;
    m = (1 << 31) - (unsigned long)1;
    a = 168019;
}

double random_generator::uniform_generator(){
    seed = (a * seed) % m;
    return (double)seed / m;
}
