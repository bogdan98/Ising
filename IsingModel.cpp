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
double hT = 0.0; //h is B*magnetic moment, T is temperature, hT is h/T
double T = 100; //temperature
int **conf;
bool **inc;
int Lx = 100; //size of the lattice in x-direction
int Ly = 100; //size of the lattice in y-direction
int Nb = 100; //number of blocks
int Ns = 10; //number of superblocks used for error estimate
int Zb = 1000; //number of steps in the block
int Ntotal = Nb*Zb;
double *values; //array containing calculated values after each MC step.
//Used in error and autocorrelation estimates
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
  values = new double[Ntotal];
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

//function that returns the value of the total energy of the lattice
//can be estimated just by putting it in Metropolis/Wolff algorithm
//useful for magnetic susceptibility/heat capacity/entropy estimates
double total_energy(){
  double energy = 0.0;
  for (int i = 0; i < Lx - 1; i++){
    for (int j = 0; j < Ly; j++){
      energy = energy + conf[i][j]*conf[i + 1][j];
    }
  }

  for (int i = 0; i < Ly - 1; i++){
    for (int j = 0; j < Lx; j++){
      energy = energy + conf[j][i]*conf[j][i + 1];
    }
  }

  return -T*JT*energy - T*hT*net_magn();
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
  dE = dE*2*JT*conf[x][y] + 2*hT*conf[x][y];
  R = exp(-dE);
  if ((R >= 1) || (((double) rand())/((double) RAND_MAX) < R)){
    conf[x][y] = -conf[x][y];
  }
}

//Metropolis algorithm
double Metropolis(){
  int Z = 0;
  double A = 0;
  srand(time(NULL));
  //to give new random values every time the algorithm
  //is called
  while(Z < Ntotal){
    MetropolisStep();
    values[Z] = (double) abs(net_magn())/((double) Lx*Ly); //filling the array for later use in error estimates
    A = A + abs(net_magn());
    Z++;
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
  while (Z < Ntotal){
    update_cluster();
    values[Z] = abs(net_magn())/((double) Lx*Ly); //filling in the array for later use in error estimates
    A = A + abs(net_magn());
    for (int i = 0; i < Lx; i++){
      for (int j = 0; j < Ly; j++){
        inc[i][j] = false;
      }
    }
    Z++;
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

//function that gives the error bars using the method of superblocks
double error(){
  double B[Ns];
  double error = 0.0;
  double avg = 0.0;
  for (int i = 0; i < Ns; i++){
    B[i] = 0.0;
    for (int j = 0; j < Nb/Ns; j++){
      for (int k = 0; k < Zb; k++){
        B[i] = B[i] + values[(i + j*Ns)*Zb + k];
      }
    }
    B[i] = B[i]/(Zb*Nb/Ns);
    avg += B[i];
  }
  avg = avg/Ns;
  for (int i = 0; i < Ns; i++){
    error = error + pow((B[i] - avg), 2.0);
  }
  error = sqrt(error)/Ns;
  return error;
}

//autocorrelation function using blocking method with Nb blocks and Zb steps in each block
//can be measured for both Metropolis and Wolff algorithms
double* Autocorrelation(){
  double* blockAv = new double[Nb];
  double *res = new double[Nb];
  for (int i = 0; i < Nb; i++){
    blockAv[i] = 0.0;
    for (int j = 0; j < Zb; j++){
      blockAv[i] = blockAv[i] + values[i*Zb + j];
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
}
