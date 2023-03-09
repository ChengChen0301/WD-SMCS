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
    string dir = "../posterior";
    mkdir(dir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    c = 0.0;
    srand(time(NULL));
    struct timeb timeSeed;
    ftime(&timeSeed);
    srand(timeSeed.time * 1000 + timeSeed.millitm); 
    double obsv[2000][81];
    char infolder[100]={0};
    // strcpy(infolder, "/samples/round2/output0.txt");
    double noiseX, noiseY;
    double sigma = 0.0;
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

    char infile[100] = "../observation/round6";
    loadData(0, infile);
    initialize();
    // initialize_from_file(infolder);
    int t=0;
    // get_observation(t);
    savedata(t);

    // paraA = 500.0*vary_func(t) + 1700.0;
    // paraB = 0.09*vary_func(t) + 0.11;
    // k_bump = 120000.0*vary_func(t);
    // K_friction = 240000.0*vary_func(t);
    // v = 1.5*vary_func(t) + 1.5;
    paraA = 1981.173;
    paraB = 0.082;
    k_bump = 120000.0;
    K_friction = 240000.0;
    v = 0.995;
    computeForce(paraA, paraB, k_bump, K_friction, v);
    
    while(number != 0){
        t++;
        update();
        // paraA = 500.0*vary_func(t) + 1700.0;
        // paraB = 0.09*vary_func(t) + 0.11;
        // k_bump = 120000.0;
        // K_friction = 240000.0;
        // v = 1.5*vary_func(t) + 1.5;
        computeForce(paraA, paraB, k_bump, K_friction, v);
        if(t%100 == 0 && number>2){ // record every 0.1 seconds
            // get_observation(t/100);
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

   //  //-----------------------to save all the positions in a .txt file------------------------
   //  int nFrame = t/100;

   //  write.open(PosName.c_str(), ios::app);
   //  for(int j=0; j<nFrame; j++){
   //      for(int i=0; i<Np; i++)
   //          write<<left<<setw(13)<<save_position[i][0][j]<<'\t'<<left<<setw(13)<<save_position[i][1][j]<<'\n';
   //  }
   //  write.close();

   //  write.open(VelName.c_str(), ios::app);
   //  for(int j=0; j<nFrame; j++){
   //      for(int i=0; i<Np; i++)
   //          write<<left<<setw(13)<<save_velocity[i][0][j]<<'\t'<<left<<setw(13)<<save_velocity[i][1][j]<<'\n';
   //  }
   //  write.close();

   //  write.open(OPName.c_str());
   //  for(int j=0; j<nFrame; j++)
   //      write<<left<<setw(7)<<j*0.1<<'\t'<<left<<setw(7)<<Out_number[j]<<"\n";
   //  write.close();

    // write.open(obsvName.c_str());
    // for(int j=0; j<nFrame; j++){
    //     for(int i=0; i<80; i++){
    //         write<<left<<setw(7)<<save_kde[j][i]<<'\t';
    //     }
    //     write<<left<<setw(7)<<save_kde[j][80]<<'\n';
    // }
    // write.close();
    read.close();


    // std::array<double, 2> p1 = {0, 2}; // a point of data (x, y)
    // std::array<double, 2> p2 = {0.15, 2.3};
    // std::array<double, 2> p3 = {-0.1, 2.5};

    // std::vector<std::array<double, 2>> data = {p1, p2, p3};

    // kdepp::Kde2d<std::array<double,2>> kernel(data);

    // std::array<double, 2> test_point = {0.1, 2.5};
    // double result = kernel.eval(test_point);
    // std::cout << result << std::endl;

    // std::vector<float> data = {-0.1, 0, 0.1, 0.1, 0.2};
    // kdepp::Kde1d<float> kernel(data);

    // auto result = kernel.eval(0.05);
    // std::cout << result << std::endl;

    cout<<"The "<<k<<"th running ends"<<"\n";
    return 0;
}

