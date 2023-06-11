#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "sys/timeb.h"

#include "vec.h"
#include "Point.h"
#include "variables.h"
#include "Line.h"
#include "simulation.h"

using namespace std;

int main(int argc, char *argv[]){

    char *k = argv[1];
    string dir = "../observation";
    mkdir(dir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    srand(time(NULL));
    struct timeb timeSeed;
    ftime(&timeSeed);
    srand(timeSeed.time * 1000 + timeSeed.millitm); 
    double obsv[2000][81];
    char infolder[100]={0};
    double noiseX, noiseY;
    double sigma = 0.04;
    int noiseFlag = 1;
    //-----------------------create the folders------------------------
    string folder;
    folder += dir;
    folder += "/round";
    folder += k;
    mkdir(folder.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);

    string PosName = folder + "/coord.txt";
    string VelName = folder + "/velocity.txt";
    string ForceName = folder + "/force.txt";

    string ETName = folder + "/evacuation_time.txt";
    string OPName = folder + "/outside_point.txt";
    string obsvName = folder + "/obsv.txt";
    //-----------------------create the folders------------------------

    initialize();

    int t=0;
    savedata(t);

    paraA = 2000.0;
    paraB = 0.08;
    k_bump = 120000.0;
    K_friction = 240000.0;
    v = 1.0;
    computeForce(paraA, paraB, k_bump, K_friction, v);
    
    while(number != 0){
        t++;
        update();
        computeForce(paraA, paraB, k_bump, K_friction, v);
        if(t%100 == 0 && number>2){ // record every 0.1 seconds
            savedata(t/100);
        }

    }

    ofstream write;
    ifstream read;

   // //-----------------------to save the positions every other 0.1s in a .txt file------------------------
    for(int j=0; j<=t/100; j++) {
        string index = to_string(j*100);
        string PosName = folder + "/obsv" + index + ".txt";
        write.open(PosName.c_str());
        for (int i = 0; i < Np; i++) {
            if(j == 0){
                noiseX = 0;
                noiseY = 0;
            }
            else{
                if(save_position[i][2][j] == 1){
                    noiseFlag = 1;
                    while(noiseFlag == 1){
                        noiseX = generateGaussianNoise(0, sigma);
                        noiseY = generateGaussianNoise(0, sigma);
                        if((save_position[i][0][j]+noiseX<_right) && (save_position[i][0][j]+noiseX>_left) && (save_position[i][1][j]+noiseY<_up) && (save_position[i][1][j]+noiseY>_down)){
                            noiseFlag = 0;
                        }
                    }
                }
                else{
                    noiseX = 0;
                    noiseY = 0;
                }
            }
            write << left << setw(13) << save_position[i][0][j]+noiseX << '\t' << left << setw(13) << save_position[i][1][j]+noiseY <<"\n";
        }
        write.close();
    }

    for(int j=0; j<=t/100; j++){
        string index = to_string(j*100);
        string PosName = folder + "/truth" + index + ".txt";
        write.open(PosName.c_str());
        for (int i = 0; i < Np; i++) {
             write << left << setw(13) << save_position[i][0][j]+noiseX << '\t' << left << setw(13) << save_position[i][1][j]+noiseY <<"\n";
        }
        write.close();
    }

    read.close();

    cout<<"The "<<k<<"th running ends"<<"\n";
    return 0;
}


