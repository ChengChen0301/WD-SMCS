#ifndef SIMULATION_H_INCLUDED
#define SIMULATION_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include<tuple>
// #include <omp.h>
#include <limits>
#include <iomanip>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include "variables.h"
#include "emd.h"

using namespace std;

void initialize();
double generateGaussianNoise(double mu, double sigma);
void computeForce(double paraA, double paraB, double k_bump, double K_friction, double v);
void update();
double generateUniform(double a, double b);
void clear_grid();
void assignstate2(int arg);
void loadData(int k, char* infile);
double disturbPara(double para, double a, double b, double level);
double *normalize_weight(double p[mcNum]);
double compute_weights(double d);
int selection(double cumFitness[parNum]);
int *resample(double cumFitness[parNum]);
void initialize_from_file(char* infile);
float dist(feature_t *F1, feature_t *F2);


#endif // SIMULATION_H_INCLUDED
