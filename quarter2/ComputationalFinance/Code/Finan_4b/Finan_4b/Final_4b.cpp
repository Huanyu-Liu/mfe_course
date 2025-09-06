//
//  main.cpp
//  Finan_4b
//
//  Created by Huanyu Liu on 3/21/19.
//  Copyright © 2019 Huanyu Liu. All rights reserved.
//

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

class random_generator{
    unsigned long seed;
    unsigned long m;
    unsigned long a;
public:
    random_generator(unsigned long x = (1 << 10) + 1);
    unsigned long get_seed();
    double uniform_generator();
    double* normal_generator(int size);
    double** bivariate_normal(int size,double x_sigma, double y_sigma, double rho);
    double* brownian_motion(int size, double t);
    double* geometric_brownian_motion(int size, double t, double r, double sigma, double s);
    double* stock_path(int size, double t, double r, double sigma, double s0);
    double * exponential(int size, double lambda);
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
    // insert code here...
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
    //double *bm;
    double sum = 0;
    double R;
    for (int i = 0; i != path_count; ++i){
        double *bm = rg.brownian_motion(time_interval, delta);
        r = r0;
        R = 0;
        for (int j = 0; j != time_interval; ++j){
            r += (alpha + beta * r) * delta + sigma * pow(r, gamma) * bm[j];
            R += r;
        }
        sum += exp(- delta * R);
        delete [] bm;
    }
    return face * sum / path_count;
}

double Bond::option(double r0, double S, double face, double T, double x, int path_count, int time_interval){
    double r;
    //double original_r = r0;
    double delta = T / time_interval;
    random_generator rg;
    //double *bm;
    double sum = 0;
    double R;
    double bond, option;
    for (int i = 0; i != path_count; ++i){
        double *bm = rg.brownian_motion(time_interval, delta);
        r = r0;
        R = 0;
        for (int j = 0; j != time_interval; ++j){
            r += (alpha + beta * r) * delta + sigma * pow(r, gamma) * bm[j];
            R += r;
        }
        bond = discount_bond(r, S - T, face, time_interval);
        option = max(x - bond);
        sum += option * exp(-delta * R);
        delete [] bm;
    }
    //this->set_r0(original_r);
    return sum / path_count;
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

double* random_generator::normal_generator(int size){
    double* box_muller_norm = new double[size];
    for (int i = 0; i < size; i+=2){
        double u1 = uniform_generator();
        double u2 = uniform_generator();
        box_muller_norm[i] = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        box_muller_norm[i+1] = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    }
    return box_muller_norm;
}

double** random_generator::bivariate_normal(int size, double x_sigma, double y_sigma, double rho){
    double** xy = new double*[2];;
    xy[0] = new double[size];
    xy[1] = new double[size];
    double* normal;
    normal = normal_generator(size * 2);
    for (int i = 0; i < size; i++){
        xy[0][i] = x_sigma * normal[i];
        xy[1][i] = y_sigma * rho * normal[i] + y_sigma * sqrt(1 - rho * rho) * normal[i+size];
    }
    return xy;
}

double* random_generator::brownian_motion(int size, double t){
    double* w = new double[size];
    w = normal_generator(size);
    for (int i = 0; i < size; i++){
        w[i] *= sqrt(t);
    }
    return w;
}

double* random_generator::geometric_brownian_motion(int size, double t, double r, double sigma, double s0){
    double* w = brownian_motion(size, t);
    double* st = new double[size];
    for (int i = 0; i < size; i++){
        st[i] = s0 * exp(sigma * w[i] + (r - sigma * sigma / 2) * t);
    }
    return st;
}
double* random_generator::stock_path(int size, double t, double r, double sigma, double s0){
    double* path = new double[size];
    double delta = t/size;
    double* bm = brownian_motion(size - 1, delta);
    path[0] = s0;
    for (int i = 1; i < size; i++){
        if (path[i - 1] > 0){
            path[i] = path[i-1] + r * delta * path[i-1] + sigma * path[i - 1] * bm[i - 1];
        }
        else{
            path[i] = 0;
        }
    }
    return path;
}

double * random_generator::exponential(int size, double lambda){
    double * result = new double[size];
    double uniform;
    for (int i = 0; i < size; ++i){
        uniform = uniform_generator();
        result[i] = -lambda * log(uniform);
    }
    return result;
}

unsigned long random_generator::get_seed(){
    return seed;
}
void random_generator::set_seed(unsigned long seed){
    this->seed = seed;
}
void random_generator::delfault(){
    set_seed(1025);
}
