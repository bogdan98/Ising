#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <time.h>
#include <ctime>
#include <math.h>
#include <random>
#include <chrono>
using namespace std;

double JT = 0.5; // J/T ratio
int **conf;
bool **inc;
int Lx = 100; //size of the lattice in x-direction
int Ly = 100; //size of the lattice in y-direction
int N1 = 100000; //number of steps for Metropolis
int N2 = 1000; //number of steps for Wolff
int Nb = 2000; //number of blocks
int Zb = 1000; //number of steps in the block
//sum of elements in the array size Lx by Ly
double sum(int **arr, int x, int y){
  double result = 0;
  for (int i = 0; i < x; i++){
    for (int j = 0; j < y; j++){
      result = result + arr[i][j];
    }
  }
  return result;
}

//sum of elements in 1D array from i0 to i1 positions inclusively
double sum2(double *arr, int i0, int i1){
  double res = 0.0;
  for (int i = i0; i <= i1; i++){
    res = res + arr[i];
  }
  return res;
}

void initialize(){
  conf = new int*[Lx];
  for (int i = 0; i < Lx; i++){
    conf[i] = new int[Ly];
    for (int j = 0; j < Ly; j++){
      conf[i][j] = 1;
    }
  }
  inc = new bool*[Lx];
  for (int i = 0; i < Lx; i++){
    inc[i] = new bool[Ly];
  }
  for (int i = 0; i < Lx; i++){
    for (int j = 0; j < Ly; j++){
      inc[i][j] = false;
    }
  }
}

//function that returns the value of M
double net_magn(){
  return sum(conf, Lx, Ly);
}

//function that returns the value of M^2
//both can be used in Metropolis and Wolff (for values accummulation)
double net_magn_sq(){
  return pow(sum(conf, Lx, Ly), 2.0);
}

//one step of Metropolis algorithm
void MetropolisStep(){
  int x; int y; //coordinates of the spin to be flipped
  int neighbor;
  double dE; double R;
  x = rand() % Lx;
  y = rand() % Ly;
  dE = 0;
  neighbor = x - 1;
  if (neighbor > -1){
    dE = dE + conf[neighbor][y];
  }
  neighbor = x + 1;
  if (neighbor < Lx){
    dE = dE + conf[neighbor][y];
  }
  neighbor = y - 1;
  if (neighbor > -1){
    dE = dE + conf[x][neighbor];
  }
  neighbor = y + 1;
  if (neighbor < Ly){
    dE = dE + conf[x][neighbor];
  }
  dE = dE*2*JT*conf[x][y];
  R = exp(-dE);
  if ((R >= 1) || (((double) rand())/((double) RAND_MAX) < R)){
    conf[x][y] = -conf[x][y];
  }
}

//Metropolis algorithm
double Metropolis(){
  double A = 0;
  int Z = 0;
  srand(time(NULL));
  //to give new random values every time the algorithm
  //is called
  while(Z<=N1){
    A = A + abs(net_magn());
    Z++;
    MetropolisStep();
  }
  return A/((double) Z*Lx*Ly);
}

//recursive function that adds new spins to the cluster to be flipped
//here, using random number generation as in Metropolis (i.e. using seed)
//did not work well
void grow_cluster(int state, int xi, int yi){
  inc[xi][yi] = true;
  conf[xi][yi] = - conf[xi][yi];
  if (xi < Lx - 1){
    if ((inc[xi + 1][yi] == false) && (conf[xi + 1][yi] == state) && (exp(-2.0*JT) < ((double) rand())/((double) RAND_MAX))){
      grow_cluster(state, xi + 1, yi);
    }
  }
  if (xi > 0){
    if ((inc[xi - 1][yi] == false) && (conf[xi - 1][yi] == state) && (exp(-2.0*JT) < ((double) rand())/((double) RAND_MAX))){
      grow_cluster(state, xi - 1, yi);
    }
  }
  if (yi < Ly - 1){
    if ((inc[xi][yi + 1] == false) && (conf[xi][yi + 1] == state) && (exp(-2.0*JT) < ((double) rand())/((double) RAND_MAX))){
      grow_cluster(state, xi, yi + 1);
    }
  }
  if (yi > 0){
    if ((inc[xi][yi - 1] == false) && (conf[xi][yi - 1] == state) && (exp(-2.0*JT) < ((double) rand())/((double) RAND_MAX))){
      grow_cluster(state, xi, yi - 1);
    }
  }
}

//function that updates the cluster
void update_cluster(){
  int state; int startx; int starty;
  startx = rand() % Lx;
  starty = rand() % Ly;
  state = conf[startx][starty];
  //int k = rand(); //to make sure the numbers in grow_cluster are more random
  grow_cluster(state, startx, starty);

}

//Wolff algorithm. Returns the value of magnetization per spin
double Wolff(){
  double A = 0;
  int Z = 0;
  srand(time(NULL));
  while (Z<=N2){
    Z++;
    A = A + abs(net_magn());
    update_cluster();
    for (int i = 0; i < Lx; i++){
      for (int j = 0; j < Ly; j++){
        inc[i][j] = false;
      }
    }
  }
  return A/((double) Z*Lx*Ly);
}

//function that takes in the values of J/T ratios and gives the corresponding magnetization per spin values
//can take both Wolff() and Metropolis() functions, i.e. work using either algorithm
double *MagnetizationCurve(double *JTs, int l, int size, int N, double(*f)()){
  double *res = new double[l];
  for (int i = 0; i < l; i++){
    res[i] = f();
  }
  return res;
}

//autocorrelation function using blocking method with Nb blocks and Zb steps in each block
//can be measured for both Metropolis and Wolff algorithms
double* Autocorrelation(){
  double M = sum(conf, Lx, Ly);
  double* blockAv = new double[Nb];
  double *res = new double[Nb];
  for (int i = 0; i < Nb; i++){
    blockAv[i] = 0.0;
    for (int j = 0; j < Zb; j++){
      blockAv[i] = blockAv[i] + abs(M);
      MetropolisStep();
      M = sum(conf, Lx, Ly);
    }
    blockAv[i] = blockAv[i]/(Lx*Ly*Zb);
    printf("%d - %f\n", i, blockAv[i]);
  }
  double avg = sum2(blockAv, 0, Nb - 1)/Nb;
  double num = 0.0; double den = 0.0;
  for (int j = 0; j < Nb; j++){
    num = 0.0; den = 0.0;
    res[j] = 0;
    for (int i = 0; i < Nb - j; i++){
      num = num + (blockAv[i] - avg)*(blockAv[i + j] - avg);
      den = den + (blockAv[i] - avg)*(blockAv[i] - avg);
    }
    res[j] = num/den;
  }
  return res;
}

int main(){
  initialize();
  //everything is ready for work after initialization
  double m1;
  m1 = Wolff();
  printf("%f", m1);
}
