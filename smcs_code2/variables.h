#ifndef VARIABLES_H
#define VARIABLES_H

#include <list>
#include <string>

#include "Point.h"

using namespace std;

// Constants & Input
#define Np  (100)      // number of particles
#define Nw  (5)        // number of walls, 5 for a square room with a door
#define _r   (0.3)      // radius of the particle
// #define m   (80.0)     // mass of a particle
#define tau (0.5)      // acceleration time
#define deltaT (0.001)   // record interval
#define Ng_x (70) // number of grids per column, more than actual size
#define Ng_y (50) // number of grids per row
#define parNum (500)  // number of ensembles
#define mcNum (100)
#define staDim (3)
#define dataScale (30)

extern Point pt[Np];
extern Line wall[Nw];
extern vec door;
extern double c;       // noise ratio
extern double len;     // side length of grids
extern int number;     // number of points in the room
extern list<int> grids[Ng_x][Ng_y]; // to store the labels of points in each grid
extern double initData[dataScale][Np][2];
extern double paraA, paraB, k_bump, K_friction, v;
extern double store_state[parNum][4*Np];

extern double _left;   // "left" conflicts with the namespace
extern double _right;
extern double _up;
extern double _down;
extern double width;   // width of the door
extern double theta;

#endif // VARIABLES_H
