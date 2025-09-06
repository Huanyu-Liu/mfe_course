//
//  main.cpp
//  Final_4
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

class Bond{
    double sigma, alpha, beta, gamma;
public:
    Bond(double sigma, double alpha, double beta, double gamma);
    double discount_bond(double r0, double T, double face, int time_interval, int path_count = 500);
    double option(double r0, double S, double face, double T, double x, int path_count, int time_interval);
};

double max(double x, double y = 0){
    return x > y ? x : y;
}


int main(int argc, const char * argv[]) {
    double r0 = 0.05, alpha = 0.36, beta = -5.86, sigma = 0.36, gamma = 2, T = 0.5, S = 1, face = 10000, strike = 9800;
    int path_count = 10000, time_interval = 100;
    Bond bond(sigma, alpha, beta, gamma);
    double option_price = bond.option(r0, S, face, T, strike, path_count, time_interval);
    cout << option_price << endl;
    return 0;
}

Bond::Bond(double sigma, double alpha, double beta, double gamma){
    //this->r0 = r0;
    this->sigma = sigma;
    this->alpha = alpha;
    this->beta = beta;
    this->gamma = gamma;
}

double Bond::discount_bond(double r0, double T, double face, int time_interval, int path_count){
    double r;
    double delta = T / time_interval;
    random_generator rg;
    vector<double> bm;
    double sum = 0;
    double R;
    for (int i = 0; i != path_count; ++i){
        bm = rg.brownian_motion(time_interval, delta);
        r = r0;
        R = 0;
        for (int j = 0; j != time_interval; ++j){
            r += (alpha + beta * r) * delta + sigma * pow(r, gamma) * bm[j];
            R += r;
        }
        sum += exp(- delta * R);
    }
    return face * sum / path_count;
}

double Bond::option(double r0, double S, double face, double T, double x, int path_count, int time_interval){
    double r;
    //double original_r = r0;
    double delta = T / time_interval;
    random_generator rg;
    vector<double> bm;
    double sum = 0;
    double R;
    double bond, option;
    for (int i = 0; i != path_count; ++i){
        bm = rg.brownian_motion(time_interval, delta);
        r = r0;
        R = 0;
        for (int j = 0; j != time_interval; ++j){
            r += (alpha + beta * r) * delta + sigma * pow(r, gamma) * bm[j];
            R += r;
        }
        bond = discount_bond(r, S - T, face, time_interval);
        option = max(x - bond);
        sum += option * exp(-delta * R);
    }
    //this->set_r0(original_r);
    return sum / path_count;
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
