#ifndef SIMULATION_H_INCLUDED
#define SIMULATION_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <tuple>
#include <limits>
#include <iomanip>
#include <cstring>
#include <sys/stat.h>
#include <opencv2/opencv.hpp> 
#include <sys/types.h>
#include "variables.h"

using namespace std;
using namespace cv;

void initialize();
double generateGaussianNoise(double mu, double sigma);
void computeForce(double paraA, double paraB, double k_bump, double K_friction, double v);
void update();
double vary_func(int t);
void storestate(int arg);
double generateUniform(double a, double b);
void clear_grid();
void assignstate(int arg);
void assignstate2(int arg);
void loadData(int k, char* infile);
double disturbPara(double para, double a, double b, double level);
double *normalize_weight(double p[mcNum]);
double compute_weights(double d);
int selection(double cumFitness[parNum]);
int *resample(double cumFitness[parNum]);
void initialize_from_file(char* infile);
float compute_wass_distance(int *count1, int *count2, int ds1, int ds2, int sum1, int sum2, double test_point[mcNum][2]);
int *count_grid1(double gridLength);
std::tuple<int*, int, int>count_grid(double gridLength);


#endif // SIMULATION_H_INCLUDED
