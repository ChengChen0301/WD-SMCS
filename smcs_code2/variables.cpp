#include "variables.h"

Point pt[Np];
Line wall[Nw];
vec door;
double c; // noise ratio
double len = 2*_r + 5*paraB;
int number;
list<int> grids[Ng_x][Ng_y];
double initData[dataScale][Np][2];
double paraA, paraB, k_bump, K_friction, v;
double store_state[parNum][4*Np];

double _left = -5.0;
double _right = 5.0;
double _up = 5.0;
double _down = -5.0;
double width = 1.0;
double theta = 90;


