#include<cstring>

#include "simulation.h"
#include "variables.h"

using namespace std;

#include<iostream>  
#include<fstream>
#include<sstream>
#include<cmath>  
#include<cstdlib>
#include<vector>
#include<cstring>
#include<tuple>
using namespace std;
 

void loadData(int k, char* infile){
    wall[0].set_X(_right, _up, _left, _up);
    wall[1].set_X(_left, _up, _left, _down);
    wall[2].set_X(_left, _down, _right, _down);
    wall[3].set_X(_right, _down, _right, -width);
    wall[4].set_X(_right, width, _right, _up);
    door.set_X(_right,0.0);

    char index[7]={0};
    char infolder[100];
    vec temp;
    double y,z;
    sprintf(index,"%d",(int)(k*nstep));
    strcpy(infolder, infile);
    strcat(infolder, "/obsv");
    strcat(infolder, index);
    strcat(infolder, ".txt");
    fstream myfile;
    myfile.open(infolder);
    int i=0;
    if(myfile) {
        while (myfile >> y >> z) {
            temp.set_X(y, z);
            pt[i].set_C(temp);
            initData[k][i][0] = y;
            initData[k][i][1] = z;
            i = i + 1;
        }
    }
    myfile.close();
    myfile.clear();
}

// generate non-overlap points
void initialize(){
    /*------------------set walls and door-----------------*/
    wall[0].set_X(_right, _up, _left, _up);
    wall[1].set_X(_left, _up, _left, _down);
    wall[2].set_X(_left, _down, _right, _down);
    wall[3].set_X(_right, _down, _right, -width);
    wall[4].set_X(_right, width, _right, _up);
    door.set_X(_right,0.0);
    /*------------------set walls and door-----------------*/


    /*------------------set point position-----------------*/
    number = Np; // initial number of points
    vec temp;
    int idx, idy;
    double _left_ = -5.0;
    double _down_ = -5.0;
    double _right_ = 5.0;
    double _up_ = 5.0;

    for(int i=0; i<Np; i++){ // rest particles
        temp.set_X(initData[0][i][0], initData[0][i][1]);
        pt[i].set_C(temp);
        idx = (int)((initData[0][i][0] - _left)/len);
        idy = (int)((initData[0][i][1] - _left)/len);
        pt[i].set_G(idx, idy);
        grids[idx][idy].push_back(i); 
        pt[i].set_flag(1);
    }
    /*------------------set point position-----------------*/
};


double generateGaussianNoise(double mu, double sigma){
    const double epsilon = numeric_limits<double>::min();
    const double two_pi = 2.0*3.14159265358979323846;

    static double z0, z1;
    static bool generate;
    generate = !generate;

    if(!generate)
        return z1 * sigma + mu;

    double u1, u2;
    do{
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    }
    while ( u1 <= epsilon );
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
};


void computeForce(double paraA, double paraB, double k_bump, double K_friction, double v){
    vec targetV, force, dir, noise, temp;
    double dist, touch, u1, u2;
    int iloc, jloc;

    //calculate the force for each ball
    // #pragma omp parallel num_threads(4)
    // #pragma omp for
    for(int i=0; i<Np; i++){
        if(pt[i].get_flag() == 1){ // still in the room
            force.set_X(0.0,0.0);

            iloc = pt[i].get_G()[0];
            jloc = pt[i].get_G()[1]; // find its grid
            // go through the adjacent 9 grids
            for(int j=-1; j<2; j++){
                for(int k=-1; k<2; k++){
                    // compute points whose labels are in the grids
                    if(iloc+j>=0 && jloc+k>=0 && (!grids[iloc+j][jloc+k].empty())){
                        for(list<int>::iterator iter=grids[iloc+j][jloc+k].begin(); iter!=grids[iloc+j][jloc+k].end();iter++){
                            int label = *iter;
                            if(label!=i && pt[label].get_flag()==1){
                                dist = pt[i].get_C().distance_wrt_point(pt[label].get_C());
                                touch = pt[i].get_C().touch_point(pt[label].get_C());
                                dir = pt[i].get_C().direction_from_point(pt[label].get_C());
                                force = force + dir*(paraA*exp((2.0*_r - dist)/paraB)) + dir*(k_bump*touch)\
                                        + dir.normalvector()*(K_friction*touch*(pt[label].get_V()-pt[i].get_V()).inner_with(dir.normalvector()));
                            }
                        }
                    }
                }
            }
            pt[i].set_F3(force);

            for(int j=0; j<Nw; j++){ //force between ball and walls
                dist = pt[i].get_C().distance_wrt_line(wall[j]);
                touch = pt[i].get_C().touch_line(wall[j]);
                dir = pt[i].get_C().direction_from_line(wall[j]);
                force = force + dir*(paraA*exp((_r - dist)/paraB)) + dir*(k_bump*touch)\
                        - dir.normalvector()*(K_friction*touch*pt[i].get_V().inner_with(dir.normalvector()));
            }
            pt[i].set_F1(force - pt[i].get_F3());

            targetV = door.direction_from_point(pt[i].get_C())*v;
            pt[i].set_F2((targetV - pt[i].get_V())/tau*pt[i].get_weight());
            pt[i].set_F(force + (targetV - pt[i].get_V())/tau*pt[i].get_weight());
        }
        else{
            force.set_X(0.0,0.0);
            pt[i].set_F(force);
            pt[i].set_F1(force);
            pt[i].set_F2(force);
            pt[i].set_F3(force);
        }
    }
};


void update(){
    int idx, idy, iloc, jloc;
    for(int i=0; i<Np; i++){
        if(pt[i].get_flag() == 1){ // still in the room
            iloc = pt[i].get_G()[0];
            jloc = pt[i].get_G()[1];

            pt[i].set_A(pt[i].get_F() / pt[i].get_weight()); //update acceleration
            pt[i].set_V(pt[i].get_V() + pt[i].get_A()*delta); // update speed
            pt[i].set_C(pt[i].get_C() + pt[i].get_V()*delta); // update position

            if(pt[i].get_C().get_x() >= door.get_x()){ // update flag
                pt[i].set_flag(0);
                number--;
            }

            idx = (int)((pt[i].get_C().get_x() - _left)/len);
            idy = (int)((pt[i].get_C().get_y() - _down)/len);

            if(idx != iloc || idy != jloc){ // if grid information changes

                for(list<int>::iterator iter=grids[iloc][jloc].begin(); iter!=grids[iloc][jloc].end();){
                    if(*iter == i){
                        iter = grids[iloc][jloc].erase(iter); //remove it from the old grid
                        continue;
                    }
                    iter++;
                } // delete the point from the old grid
                if(pt[i].get_flag() == 1){
                    grids[idx][idy].push_back(i); // add the point into the new grid
                    pt[i].set_G(idx,idy); // update grid
                }
            }
        }
        else{
            pt[i].set_C(pt[i].get_C() + pt[i].get_V()*delta);
        }
    }
};


void clear_grid(){
    for(int p=0; p<Ng_x; p++){
        for(int q=0; q<Ng_y; q++){
            grids[p][q].clear();
        }
    }
};


void assignstate2(int arg){
    vec temp;
    int idx, idy;
    int sum = 0;
    clear_grid();
    for(int j=0; j<Np; j++){
        temp.set_X(initData[arg][j][0], initData[arg][j][1]);
        pt[j].set_C(temp);
        temp.set_X(0.0, 0.0);
        pt[j].set_V(temp);
        idx = (int)((pt[j].get_C().get_x() - _left)/len);
        idy = (int)((pt[j].get_C().get_y() - _down)/len);
        if(pt[j].get_C().get_x() < door.get_x()){
            pt[j].set_flag(1);
            pt[j].set_G(idx, idy);
            grids[idx][idy].push_back(j); 
            sum += 1;
        }
        else
            pt[j].set_flag(0);
    }
    number = sum;
};


float dist(feature_t *F1, feature_t *F2){
    float dX = F1->X - F2->X, dY = F1->Y - F2->Y;
    return sqrt(dX*dX + dY*dY);   
};


double compute_weights(double d){
    const double two_pi = 2.0*3.14159265358979323846;
    return 1.0 / sqrt(two_pi) * exp(-(d*d)/2.0);
};


int selection(double cumFitness[parNum]){ // Roulette Wheel Selection
    double roll = rand()*(1.0 / RAND_MAX);
    int pick;
    for(int i=0; i<parNum; i++){
        if(cumFitness[i] > roll){
            pick = i;
            break;
        }
    }
    return pick;
};


int *resample(double cumFitness[parNum]){
    int* pickindex = new int[parNum];
    for(int i=0; i<parNum; i++)
        pickindex[i] = selection(cumFitness);
    return pickindex;
}


double cutOffGaussianNoise(double mu, double sigma, int k){
    int flag = 1;
    double cgn;
    while(flag == 1){
        cgn = generateGaussianNoise(mu, sigma);
        if(cgn < mu+k*sigma && cgn > mu-k*sigma)
            flag = 0;
    }
    return cgn;
};


double disturbPara(double para, double a, double b, double level){
    int flag = 1;
    double u, temp;
    while(flag == 1){
        // u = cutOffGaussianNoise(0, 1, 1);
        u = generateGaussianNoise(0, 1);
        temp = para+level*u;
        if(temp>=a && temp<=b){
            flag = 0;
            para = temp;
        }
    }
    return para;
}


double generateUniform(double a, double b){
    double r;
    r = (double)rand() / RAND_MAX *(b - a) + a;
    return r;
};


double *normalize_weight(double p[mcNum]){
    double sum = 0.0;
    for(int i=0; i<mcNum; i++){
        sum += p[i];
    }
    for(int i=0; i<mcNum; i++){
        p[i] /= sum;
    }
    return p;
};
