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
void loadData(int k, char* infile);


#endif // SIMULATION_H_INCLUDED
