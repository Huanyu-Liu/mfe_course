//
//  main.cpp
//  Final_2
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

class Option{
    double s0, r, t, v0, alpha, beta, gamma;
public:
    Option(double s0, double r, double t, double v0, double alpha, double beta, double gamma);
    double price(double rho, int time_interval, int simulation_count);
};

double max(double x, double y = 0){
    return x > y ? x : y;
}

int main(int argc, const char * argv[]) {
    // insert code here...
    double v0 = 0.06, alpha = 0.45, beta = -5.105, gamma = 0.25, s0 = 20, r = 0.05, T = 2;
    double rho = -0.75;
    int time_interval = 100, simulation_count = 10000;
    Option op(s0, r, T, v0, alpha, beta, gamma);
    cout << op.price(rho, time_interval, simulation_count) << endl;
    rho = 0;
    cout << op.price(rho, time_interval, simulation_count) << endl;
    rho = 0.75;
    cout << op.price(rho, time_interval, simulation_count) << endl;
    return 0;
}

Option::Option(double s0, double r, double t, double v0, double alpha, double beta, double gamma){
    this->s0 = s0;
    this->r = r;
    this->v0 = v0;
    this->alpha = alpha;
    this->t = t;
    this->beta = beta;
    this->gamma = gamma;
}
double Option::price(double rho, int time_interval, int simulation_count){
    double delta = t / time_interval;
    random_generator rg;
    double result = 0;
    for (int i = 0; i != simulation_count; ++i){
        double average = 0;
        double v = v0;
        double s = s0;
        vector<double> bm = rg.brownian_motion(time_interval * 2, delta);
        for (int j = 0; j != time_interval; ++j){
            double temp = (alpha + beta * max(v)) * delta + gamma * sqrt(max(v)) * bm[j];
            double temp2 = r * s * delta + s * sqrt(max(v)) * (rho * bm[j] + sqrt(1 - rho * rho) * bm[time_interval + j]);
            v += temp;
            s += temp2;
            average += s;
        }
        average /= time_interval;
        result += max(s - average);
    }
    return result * exp(-r * t) / simulation_count;
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
