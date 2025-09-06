//
//  main.cpp
//  Final_3
//
//  Created by Huanyu Liu on 3/21/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>
using std::vector;
using std::cout;
using std::endl;
class random_generator{
    unsigned long seed;
    unsigned long m;
    unsigned long a;
public:
    random_generator(unsigned long x = (1 << 10) + 1);
    ~random_generator();
    unsigned long get_seed(){return seed;}
    double uniform_generator();
    vector<double> normal_generator(int size);
    vector<double> brownian_motion(int size, double t);
    vector<double> cir_stock_path(int size, double r0, double kappa, double r_bar, double delta, double sigma);
    vector<double> stock_path(int size, double t, double r, double sigma, double s0);
    void set_seed(unsigned long seed);
    void delfault();
};

int main(int argc, const char * argv[]) {
    // insert code here...
    double s0 = 100, r = 0.05, sigma = 0.35, T = 5, k = 100;
    int time_interval = 100, simulation_count = 10000;
    double delta = T / time_interval;
    double tau;
    random_generator rg;
    double result = 0;
    int exercise_count = 0, exercise_L = 0;
    for (int i = 0; i != simulation_count; ++i){
        vector<double> stock_path = rg.stock_path(time_interval, T, r, sigma, s0);
        //vector<double> bm = rg.brownian_motion(time_interval, delta);
        //double s = s0;
        tau = T;
        double payoff = 0;
        double Lt, Ut;
        for (int j = 0; j != time_interval; ++j){
            Lt = 50 * exp(0.138629 * delta * j);
            Ut = 200 - Lt;
            if (stock_path[j] < Lt){
                tau = j * delta;
                //cout << stock_path[j] << endl;
                payoff = k - stock_path[j];
                exercise_count++;
                exercise_L++;
                break;
            }
            if (stock_path[j] > Ut){
                tau = j * delta;
                payoff = stock_path[j] - k;
                exercise_count++;
                break;
            }
            //s += s * r * delta + sigma * s * bm[j];
            payoff = 0;
        }
        result += payoff * exp(-r * tau);
    }
    result /= simulation_count;
    double prob = 1.0 * exercise_L / exercise_count;
    cout << result << endl;
    //cout << exercise_count << endl;
    cout << prob << endl;
    return 0;
}


random_generator::random_generator(unsigned long x){
    seed = x;
    m = (1 << 31) - (unsigned long)1;
    a = 168019;
}
random_generator::~random_generator(){
    
}
double random_generator::uniform_generator(){
    seed = (a * seed) % m;
    return (double)seed / m;
}

vector<double> random_generator::normal_generator(int size){
    vector<double> box_muller_norm;
    for (int i = 0; i < size; i+=2){
        double u1 = uniform_generator();
        double u2 = uniform_generator();
        box_muller_norm.push_back(sqrt(-2 * log(u1)) * cos(2 * M_PI * u2));
        box_muller_norm.push_back(sqrt(-2 * log(u1)) * sin(2 * M_PI * u2));
    }
    return box_muller_norm;
}

vector<double> random_generator::brownian_motion(int size, double t){
    vector<double> w;
    w = normal_generator(size);
    for (int i = 0; i < size; i++){
        w[i] *= sqrt(t);
    }
    return w;
}

vector<double> random_generator::cir_stock_path(int size, double r0, double kappa, double r_bar, double delta, double sigma){
    vector<double> path;
    path.push_back(r0);
    vector<double> bm = brownian_motion(size, delta);
    for (int i = 0; i != size - 1; ++i){
        path.push_back(path[i] + kappa * (r_bar - path[i]) * delta + sigma * sqrt(abs(path[i])) * bm[i]);
    }
    return path;
}

vector<double> random_generator::stock_path(int size, double t, double r, double sigma, double s0){
    vector<double> path;
    double delta = t/size;
    //double* bm = brownian_motion(size - 1, delta);
    vector<double> bm = brownian_motion(size, delta);
    path.push_back(s0);
    for (int i = 1; i < size; i++){
        if (path[i - 1] > 0){
            path.push_back(path[i-1] + r * delta * path[i-1] + sigma * path[i - 1] * bm[i - 1]);
        }
        else{
            path.push_back(0);
        }
    }
    return path;
}
