#pragma once


//Simulation Parameters
#define DIMENSION      2       // number of spatial dimensions, 2 or 3
#define BETA           1       // inverse temp
#define POTENTIAL      0       // (-1,0,1) -> (MC, Morse, Lennard Jones)
#define RANGE         40			 // range paremeter for potentials
#define EULER_TS    5e-6       // time step for EM integrator  5e-6 for bd, 50e-6 for sampling
#define RK_TS       1e-6       // time step for RK integrator
#define SAMPLES      1000       // number of samples to obtain (per prcoessor)
#define EQ           500       // number of samples to equilibrate for (change to time?) 

#if (DIMENSION == 2)
	#define BOND_CUTOFF  1.04       // if two particles are less than this distance, bonded
#elif (DIMENSION == 3)
	#define BOND_CUTOFF  1.05      //original was 1.029
#endif

//Newton's Method Parameters
#define N_ITER        20       // Max iterations for NM
#define N_TOL       1e-6       // tolerance for norm(dx)

//mcmc parameters
#define SIG          0.15       // std for isotropic gaussian proposal
#define KAP          2         // sticky parameter for particles
#define EPS          1         // base energy for lattice proteins
#define MAX_TRY     1000       // if this many MCMC samples dont transition, ignore
#define MAX_ITS SAMPLES*300    // max number of mcmc steps before sim is cut off

//visual parameters
#define PTOL        0.1       // exclude nodes from plot with less than this probability