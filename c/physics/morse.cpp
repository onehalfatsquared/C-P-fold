#include <math.h>
#include <algorithm>
#include "bDynamics.h"
#include "../defines.h"
namespace bd {

double morseP(double r, double rho, double E){
  //evaluate the morse potential derivative
  double Y = exp(-(rho) * ((r) - 1));
  return -2 * (rho) * (E) * ( Y*Y -Y);
}

double morseEval(double* particles, int rho, double* E, int N, int* P) {
  //compute total energy of system with pairwise morse potential //not necc?
	double S = 0; int i, j; int rep = 0.0;
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(j!=i){
                double* Z = new double[DIMENSION];
                double r = euDist(particles, i, j, N, Z);
                int minimum = std::min(i,j);
                int maximum = std::max(i,j);
                if(P[N*maximum+minimum]==1){
                    double Eeff = E[N*maximum+minimum];
                    S += morseP(r, rho, Eeff); 
                }
                else if(r<1){
                    S -= rep/r; 
                }
              delete []Z;
            }
        }
    }
    return S;
}

void morseGrad(double* particles, int rho, double* E, int N, int* P, double* g) {
  //compute gradient of energy of system with pairwise morse potential
	int i, j; int rep = 250.0;
	double* S = new double[DIMENSION];
	for(i=0;i<N;i++){
    S[0] = S[1] = 0;
#if (DIMENSION == 3)
    S[2] = 0;
#endif

    for(j=0;j<N;j++){
      if(j!=i){
        double* Z = new double[DIMENSION];
        double r = euDist(particles, i, j, N, Z);
        int minimum = std::min(i,j);
        int maximum = std::max(i,j);
        if(P[N*maximum+minimum]==1){
          double Eeff = E[N*maximum+minimum];
          for(int k=0; k<DIMENSION;k++) {
            if (r < 1 || abs(i-j) == 1) {
              S[k] += 2.0 * rho * rho * Eeff * (r-1.0) * Z[k];
            }
            else {
              S[k] = S[k] + morseP(r, rho, Eeff) * Z[k];
            }
          }
        }
        else if(r<1.1){
          //printf("repelling\n");
          for(int k=0; k<DIMENSION;k++)
            S[k] = S[k] - rep * Z[k] / (r*r);
        }
      	delete []Z;
    	}
		}
    for(int k=0; k<DIMENSION;k++)
      g[DIMENSION*i+k] = S[k];
  }
  delete []S;
}

void morseGradR(double* particles, int rho, double* E, int N, int* P, double* g) {
  //compute gradient of energy of system with pairwise morse potential
  //includes repulsion force to keep particles from binding
  int i, j; 
  int rep = 45.0;
  double* S = new double[DIMENSION];
  for(i=0;i<N;i++){
    S[0] = S[1] = 0;
#if (DIMENSION == 3)
    S[2] = 0;
#endif

    for(j=0;j<N;j++){
      if(j!=i){
        double* Z = new double[DIMENSION];
        double r = euDist(particles, i, j, N, Z);
        int minimum = std::min(i,j);
        int maximum = std::max(i,j);
        if(P[N*maximum+minimum]==1){
          double Eeff = E[N*maximum+minimum];
          for(int k=0; k<DIMENSION;k++) {
            if (r < 1 || abs(i-j) == 1 || (i == 4 && j == 6) || (i==6 && j ==4)) {
              S[k] += 2.0 * rho * rho * Eeff * (r-1.0) * Z[k];
            }
            else {
              double mult = morseP(r, rho, Eeff);
              S[k] = S[k] + mult * Z[k];
              if (r < 1.2) {
                //printf("p1 %d, p2 %d, Distance is %f, mult is %f\n", i,j,r, mult);
                S[k] = S[k] - mult * Z[k];
                S[k] = S[k] - rep * Z[k] / (r*r);
              }
            }
          }
        }
        delete []Z;
      }
    }
    for(int k=0; k<DIMENSION;k++)
      g[DIMENSION*i+k] = S[k];
  }
  delete []S;
}


}