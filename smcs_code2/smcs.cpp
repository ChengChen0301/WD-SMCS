#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "sys/timeb.h"
#include <time.h>
#include <thread>

#include "vec.h"
#include "Point.h"
#include "variables.h"
#include "Line.h"
#include "simulation.h"

using namespace std;

int main(int argc, char *argv[]){

    char *k = argv[1];
    string dir = "../result";
    char infile[100] = "../observation/round5";
    mkdir(dir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    srand(time(NULL));
    struct timeb timeSeed;
    ftime(&timeSeed);
    srand(timeSeed.time * 1000 + timeSeed.millitm); 
    clock_t start, end, start0, end0;
    unsigned int numThreads = std::thread::hardware_concurrency();

    //-----------------------create the folders------------------------
    string folder;
    folder += dir;
    folder += "/round5_";
    folder += k;
    mkdir(folder.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);

    double a1=1200.0, b1=2200.0, a2=0.02, b2=0.2, a3=100000, b3=200000, a4=200000, b4=300000, a5=0, b5=1.5;
    double c1=2000.0, c2=0.08, c3=120000, c4=240000, c5=1.0;
    double weights[parNum], cumFitness[parNum], predSamp[parNum][staDim], test_point[mcNum][2], tempSamp[parNum][staDim], nextSamp[parNum][staDim];
    float tempWass1, tempWass2, acceptance;
    double save_ensemble[dataScale][parNum][9], save_eff[dataScale];
    double total_weight = 0.0, cum_weight = 0.0, total_weight_square = 0.0, Neff = 0.0, Nth = 1.0/2.0*parNum;
    int index_i, index_j, *count1, *count2, *count3; 
    int ds1, ds2, ds3, sum1, sum2, sum3, ttt, flag0[parNum];
    double gridLength = (_right - _left) / sqrt(mcNum);
    double u, p1, p2;
    float save_each_wass[dataScale][parNum][2*dataScale], minDist[dataScale], accp[parNum], temp_each_wass[dataScale][parNum][2*dataScale];
    start0 = clock();
    

    ofstream write;
    string pointName = folder + "/test_point.txt";
    write.open(pointName.c_str());
    for(int i=0; i<mcNum; i++){
        index_i = i / sqrt(mcNum);
        index_j = i % (int)(sqrt(mcNum));
        test_point[i][0] = _left + 0.5*(2*index_i+1)*(_right - _left) / sqrt(mcNum);
        test_point[i][1] = _down + 0.5*(2*index_j+1)*(_up - _down) / sqrt(mcNum);
        write << left << setw(13) << test_point[i][0] << '\t' << left << setw(13) << test_point[i][1] <<"\n";
    }
    write.close();

    for(int t=0; t<dataScale; t++) // read all the observation files
        loadData(t, infile);

    for(int i=0; i<parNum; i++){ // initialize the particles
        predSamp[i][0] = generateUniform(a1, b1);
        predSamp[i][1] = generateUniform(a2, b2);
        predSamp[i][2] = generateUniform(a5, b5);
        weights[i] = 1.0 / parNum;
        save_ensemble[0][i][0] = predSamp[i][0];
        save_ensemble[0][i][1] = predSamp[i][1];
        save_ensemble[0][i][2] = predSamp[i][2];
        save_ensemble[0][i][3] = weights[i];
    }
    save_eff[0] = parNum;

    for(int t=1; t<dataScale; t++){
        start = clock();
        cout << t << "\t";
        total_weight = 0.0;
        total_weight_square = 0.0;

        for(int tt=0; tt<dataScale; tt++)
            minDist[tt] = 100.0;

        for(int i=0; i<parNum; i++){
            // generate x* ~ k( | x_{t-1})
            nextSamp[i][0] = disturbPara(predSamp[i][0], a1, b1, 10);
            nextSamp[i][1] = disturbPara(predSamp[i][1], a2, b2, 0.001);
            nextSamp[i][2] = disturbPara(predSamp[i][2], a5, b5, 0.1);

            for(int tt=1; tt<=t; tt++){
                // the observation density y_tt
                assignstate2(tt);
                std::tie(count1, ds1, sum1) = count_grid(gridLength);

                if(tt < t){
                    if(flag0[i] == 0)
                        tempWass1 = save_each_wass[t-1][i][tt];
                    else
                        tempWass1 = save_each_wass[t-1][i][tt+dataScale];
                }
                else{
                    // generate y_tt from x_{t-1}
                    assignstate2(0);
                    ttt = 0;
                    while(ttt < tt){
                        for(int l=0; l<100; l++){
                            computeForce(predSamp[i][0], predSamp[i][1], c3, c4, predSamp[i][2]);
                            update();
                        }
                        ttt++;
                    }
                    // compute p(y_tt | x_{t-1})
                    std::tie(count2, ds2, sum2) = count_grid(gridLength);
                    tempWass1 = compute_wass_distance(count1, count2, ds1, ds2, sum1, sum2, test_point);
                }
                save_each_wass[t][i][tt] = tempWass1;
                if(tempWass1 < minDist[tt]){
                    if(tempWass1 > 0)
                        minDist[tt] = tempWass1;
                }

                // generate y_tt from x*
                assignstate2(0);
                ttt = 0;
                while(ttt < tt){
                    for(int l=0; l<100; l++){
                        computeForce(nextSamp[i][0], nextSamp[i][1], c3, c4, nextSamp[i][2]);
                        update();
                    }
                    ttt++;
                }
                // compute p(y_tt | x*)
                std::tie(count3, ds3, sum3) = count_grid(gridLength);
                tempWass2 = compute_wass_distance(count1, count3, ds1, ds3, sum1, sum3, test_point);
                save_each_wass[t][i][tt+dataScale] = tempWass2;
                if(tempWass2 < minDist[tt]){
                    if(tempWass2 > 0)
                        minDist[tt] = tempWass2;
                }
            }
        }

        total_weight = 0.0;
        for(int i=0; i<parNum; i++){
            p1 = 1.0;
            p2 = 1.0;
            for(int tt=1; tt<=t; tt++){
                p1 *= compute_weights(save_each_wass[t][i][tt]/minDist[tt]);
                p2 *= compute_weights(save_each_wass[t][i][tt+dataScale]/minDist[tt]);
            }
            if(p1 > p2)
                acceptance = p2 / p1;
            else
                acceptance = 1.0;
            u = generateUniform(0, 1);
            if(u < acceptance){
                predSamp[i][0] = nextSamp[i][0];
                predSamp[i][1] = nextSamp[i][1];
                predSamp[i][2] = nextSamp[i][2];
                flag0[i] = 1;
            }
            else
                flag0[i] = 0;
            weights[i] *= compute_weights(save_each_wass[t][i][t]/minDist[t]);
            accp[i] = acceptance;
            total_weight += weights[i];
        }
        // compute ess
        cum_weight = 0.0;
        for(int i=0; i<parNum; i++){
            weights[i] /= total_weight;
            save_ensemble[t][i][0] = predSamp[i][0];
            save_ensemble[t][i][1] = predSamp[i][1];
            save_ensemble[t][i][2] = predSamp[i][2];
            save_ensemble[t][i][3] = weights[i];
            save_ensemble[t][i][4] = accp[i];
            save_ensemble[t][i][5] = flag0[i];
            save_ensemble[t][i][6] = nextSamp[i][0];
            save_ensemble[t][i][7] = nextSamp[i][1];
            save_ensemble[t][i][8] = nextSamp[i][2];
            cum_weight += weights[i];
            cumFitness[i] = cum_weight;
            total_weight_square += weights[i]*weights[i];
        }
        Neff = 1.0 / total_weight_square;
        cout << Neff << "\t";
        save_eff[t] = Neff;
        if(Neff < Nth){
            int *pickindex = resample(cumFitness);
            for(int i=0; i<parNum; i++){
                weights[i] = 1.0 / parNum;
                for(int j=0; j<staDim; j++)
                    tempSamp[i][j] = predSamp[pickindex[i]][j];
                for(int tt=1; tt<=t; tt++){
                    temp_each_wass[t][i][tt] = save_each_wass[t][pickindex[i]][tt];
                    temp_each_wass[t][i][tt+dataScale] = save_each_wass[t][pickindex[i]][tt+dataScale];
                }
            }
            for(int i=0; i<parNum; i++){
                for(int j=0; j<staDim; j++)
                    predSamp[i][j] = tempSamp[i][j];
                for(int tt=1; tt<=t; tt++){
                    save_each_wass[t][i][tt] = temp_each_wass[t][i][tt];
                    save_each_wass[t][i][tt+dataScale] = temp_each_wass[t][i][tt+dataScale];
                }
            }
        }
        end = clock();
        cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s, total time = "<<double(end-start0)/CLOCKS_PER_SEC<<"s"<<endl;

        string wassName = folder + "/wass" + to_string(t) + ".txt";
        write.open(wassName.c_str());
        for(int i=0; i<parNum; i++){
            for(int j=0; j<2*dataScale-1; j++){
                write << left << setw(13) << save_each_wass[t][i][j] << '\t';
            }
            write << left << setw(13) << save_each_wass[t][i][2*dataScale-1] << '\n';
        }
        write.close();
    }

    string ensembleName = folder + "/ensembles.txt";
    write.open(ensembleName.c_str());
    for(int t=0; t<dataScale; t++){
        for(int i=0; i<parNum; i++){
            for(int j=0; j<8; j++){
                write << left << setw(13) << save_ensemble[t][i][j] << '\t';
            }
            write << left << setw(13) << save_ensemble[t][i][8] << '\n';
        }
    }
    write.close();

    string effName = folder + "/eff.txt";
    write.open(effName.c_str());
    for(int t=0; t<dataScale; t++){
        write << save_eff[t] << '\n';
    }
    write.close();

    cout<<"The "<<k<<"th running ends"<<"\n";
    return 0;
}

