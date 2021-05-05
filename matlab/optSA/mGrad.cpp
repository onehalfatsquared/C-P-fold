/*********************************************************************
 * mGrad.cpp
 *
 * This mex file computes the gradient of a system of N particles in 3D 
 * with a Morse interaction potential with given range and depth. 
 * Includes functionality for up to 3 particle types and variable
 * interaction strength between particles.   
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <math.h>
#include <algorithm>

/*
 * Function:  euDist
 * --------------------
 * computes euclidean distance b/w particles i and j. Also gives unit vector
 *
 *  particles: N by 3 array of positions. Z: stores unit vector. 
 *
 */

double euDist(double* particles, int i, int j, int N, double* Z){
    Z[0] = particles[i]-particles[j];
    Z[1] = particles[N+i]-particles[N+j];
    Z[2] = particles[2*N+i]-particles[2*N+j];

    double R = sqrt(Z[0]*Z[0]+Z[1]*Z[1]+Z[2]*Z[2]);
    Z[0] /= R; Z[1] /= R; Z[2] /= R;
    return R;
}

/*
 * Function:  morseEval
 * --------------------
 * computes and returns value of Morse potential for distance r, range rho
 * and well-depth E. 
 *
 */

double morseEval(double r, double rho, double E){
    double Y = exp(-(rho) * ((r) - 1));
    return -2 * (rho) * (E) * ( Y*Y -Y);
}

// Main function 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *particles_m, *rho_m, *E_m, *N_m, *P_m, *g_m;
    double *particles, *rho, *E, *P, *g;
    const double rep = 200;
    int N;
    int i,j;

//associate inputs
    particles_m = mxDuplicateArray(prhs[0]);
    rho_m = mxDuplicateArray(prhs[1]);
    E_m = mxDuplicateArray(prhs[2]);
    N_m = mxDuplicateArray(prhs[3]);
    P_m = mxDuplicateArray(prhs[4]);

//associate pointers
    particles = mxGetPr(particles_m);
    rho = mxGetPr(rho_m);
    E = mxGetPr(E_m);
    N = round(*(mxGetPr(N_m)));
    P = mxGetPr(P_m);

//associate outputs
    g_m = plhs[0] = mxCreateDoubleMatrix(3*N,1,mxREAL);
    g = mxGetPr(g_m);

//compute the gradient
    for(i=0;i<N;i++){
        double* S = new double[3]; S[0] = S[1] = S[2] = 0;
        for(j=0;j<N;j++){
            if(j!=i){
                double* Z = new double[3];
                double r = euDist(particles, i, j, N, Z);
                int minimum = std::min(i,j);
                int maximum = std::max(i,j);
                if(P[N*maximum+minimum]==1){
                    double Eeff = E[N*maximum+minimum];
                    for(int k=0; k<3;k++)
                        S[k] = S[k] + morseEval(r, *rho, Eeff) * Z[k];
                }
                else if(r<1){
                    for(int k=0; k<3;k++)
                        S[k] = S[k] - rep * Z[k] / (r*r);
                }
              delete []Z;
            }
        }
      for(int k=0; k<3;k++)
            g[3*i+k] = S[k];
      delete []S;
    }
}
