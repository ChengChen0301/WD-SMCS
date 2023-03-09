#ifndef SIMULATION_H_INCLUDED
#define SIMULATION_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <tuple>
// #include <omp.h>
#include <limits>
#include <iomanip>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include "variables.h"

using namespace std;

void initialize();
void savedata(int arg);
double generateGaussianNoise(double mu, double sigma);
void computeForce(double paraA, double paraB, double k_bump, double K_friction, double v);
void update();
double vary_func(int t);
// double *get_observation(double test_point[mcNum][2]);
// void storestate(int arg);
// double generateUniform(double a, double b);
// void clear_grid();
// void assignstate(int arg);
// void assignstate2(int arg);
void loadData(int k, char* infile);
// double compute_kl_distance(double q[mcNum], double p[mcNum]);
// double compute_eu_distance(double initData[Np][2]);
// double disturbPara(double para, double a, double b, double level);
// double *normalize_weight(double p[mcNum]);
// double compute_weights(double d);
// int selection(double cumFitness[parNum]);
// int *resample(double cumFitness[parNum]);
// void initialize_from_file(char* infile);
// float dist(feature_t *F1, feature_t *F2);
// float compute_wass_distance(int *count1, int *count2, int ds1, int ds2, int sum1, int sum2, double test_point[mcNum][2]);
// int *count_grid1(double gridLength);
// std::tuple<int*, int, int>count_grid(double gridLength);


#endif // SIMULATION_H_INCLUDED
