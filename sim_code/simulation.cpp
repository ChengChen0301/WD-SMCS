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
 
// void loadData(double obsv[2000][81]) {
//     ifstream in("round3/obsv.txt");
//     string line;
//     if (in.fail()){
//         cout << "no file!" << endl;
//         getchar();
//         exit(0);
//     }
//     cout << "there is  file!" << endl;
//     int count = 0;
//     int i = 0, j = 0;
//     while (getline(in, line)) {
//         stringstream ss(line);
//         float x;
//         while (ss >> x) {
//             count++;
//             if (count < 81){
//                 obsv[i][j] = x;
//                 j++;
//             } 
//         }
//         i++;
//         j = 0;
//         count = 0;   
//     }
// }

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
    sprintf(index,"%d",(int)(k*100));
    strcpy(infolder, infile);
    strcat(infolder, "/truth");
    strcat(infolder, index);
    strcat(infolder, ".txt");
    fstream myfile;
    myfile.open(infolder);
    int i=0;
    if(myfile) {
        while (myfile >> y >> z) {
            temp.set_X(y, z);
            pt[i].set_C(temp);
            // if(k == 0){
            initData[k][i][0] = y;
            initData[k][i][1] = z;
            // }
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

    // double x = (double)rand() / RAND_MAX *(_right_ - _left_ - 2.0*_r) + _left_ + _r; // the first particle
    // double y = (double)rand() / RAND_MAX *(_up_ - _down_ - 2.0*_r) + _down_ + _r;
    // idx = (int)((x - _left)/len);
    // idy = (int)((y - _down)/len);
    // temp.set_X(x,y); pt[0].set_C(temp);
    // pt[0].set_G(idx, idy);
    // grids[idx][idy].push_back(0); // put the label "0" into the corresponding grid
    // pt[0].set_flag(1); // in the room

    // for(int i=1; i<Np; i++){ // rest particles
    //     int flag = 1;
    //     while(flag == 1){
    //         x = (double)rand() / RAND_MAX *(_right_ - _left_ - 2.0*_r) + _left_ + _r;
    //         y = (double)rand() / RAND_MAX *(_up_ - _down_ - 2.0*_r) + _down_ + _r;
    //         temp.set_X(x,y);
    //         int index = 0;
    //         for (int j=0; j<i; j++){
    //             if(temp.distance_wrt_point(pt[j].get_C()) <= 2.0*_r){
    //                 index = j;
    //                 break;
    //             }
    //             index = j+1;
    //         }
    //         if(index >= i && temp.distance_wrt_point(pt[i-1].get_C()) > 2.0*_r){
    //             flag = 0;
    //             pt[i].set_C(temp);
    //             idx = (int)((x - _left)/len);
    //             idy = (int)((y - _down)/len);
    //             pt[i].set_G(idx, idy);
    //             grids[idx][idy].push_back(i); 
    //             pt[i].set_flag(1);
    //         }
    //     }
    // }

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


double vary_func(int t){
    const double two_pi = 2.0*3.14159265358979323846;
    // return abs(cos(two_pi / 120.0 * double(t) * 0.001)) + 0.01;
    return cos(two_pi / 100.0 * 0.001);
}


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
        // else{
        //     pt[i].set_C(pt[i].get_C() + pt[i].get_V()*delta);
        // }
      
    }
};

void savedata(int arg){
    // arg = arg/100;
    for(int i=0; i<Np; i++){
        save_position[i][0][arg] = pt[i].get_C().get_x();
        save_position[i][1][arg] = pt[i].get_C().get_y();
        save_position[i][2][arg] = pt[i].get_flag();
    }
};
