#include <random>
#include <math.h>
#include <chrono>
#include "bDynamics.h"
#include "../defines.h"
namespace bd{


void EM(double* X0, int N, int Nt, double k, 
					int rho, double* E, int* P, double beta, int pot) {
	//apply the EM method to solve the SDE
	//initialize particle and gradient storage
	double* g = new double[DIMENSION*N];
	double* particles = new double[DIMENSION*N];

	//initialize the random number generator - Normal(0,1) //fix seed
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0,1.0);

	//apply the EM scheme
	for (int i = 0; i < Nt; i++) {
		c2p(X0, particles, N);
		if (pot == -1) {
			morseGradR(particles, rho, E, N, P, g);
		}
		if (pot == 0) {//use morse potential
			morseGrad(particles, rho, E, N, P, g);
		}
		else if (pot == 1) {//use lennard jones potential
			ljGrad(particles, rho, E, N, P, g);
		}
		for (int j = 0; j < DIMENSION*N; j++) {
			X0[j] += -g[j]*k + sqrt(2.0*k/beta)*distribution(generator);
		}
	}

	//free the memory
	delete []g; delete []particles;
}

void solveSDE(double* X0, int N, double T, int rho, double beta,
												 double* E, int* P, int method, int pot) {
	if (method == 1) {
		//set time step
		double k = EULER_TS; int Nt = T/k; 
		//solve the sde
		EM(X0,  N, Nt, k, rho, E, P, beta, pot);
	}
}

void setupChain(double* X, int N) {
	//construct a linear chain of particle positions
	for (int i = 0; i < DIMENSION*N; i++) {
		if (i % DIMENSION == 0) {
			X[i] = i/DIMENSION+1;
		}
		else{
			X[i] = 0;
		}
	}
}

void printCluster(double* X, int N) {
	//print the cluster in X
	for (int i = 0; i < DIMENSION*N; i++) printf("%f\n",X[i]);
}


void makeKappaMap(int numTypes, double* kappaVals, 
									std::map<std::pair<int,int>,double>& kappa) {
	//fill the map with values in kappaVals

	int counter = 0; 

	for (int i = 0; i < numTypes; i++) {
		for (int j = i; j < numTypes; j++) {
			kappa[std::make_pair(i,j)] = kappaVals[counter];
			kappa[std::make_pair(j,i)] = kappaVals[counter];
			counter++;
		}
	}
}

void fillP(int N, int* particleTypes, int* P, double* E,
					 std::map<std::pair<int,int>,double>& kappa) {
	//fill interaction matrices with values

	double range = 40;
	for (int i = 0; i < N-1; i++) {
		P[toIndex(i,i+1,N)] = 1; E[toIndex(i,i+1,N)] = 14; 
		P[toIndex(i+1,i,N)] = 1; E[toIndex(i+1,i,N)] = 14; 
		for (int j = i+2; j < N; j++) {
			int p1 = particleTypes[i]; 
			int p2 = particleTypes[j];
			double kval = kappa[{p1,p2}];
			if (kval > 0.01) {
				double energy = stickyNewton(8.0, range, kval, BETA); 
				//printf("Energy is %f\n", energy);
				P[toIndex(i,j,N)] = 1; E[toIndex(i,j,N)] = energy; 
				P[toIndex(j,i,N)] = 1; E[toIndex(j,i,N)] = energy; 
			}
			else{
				P[toIndex(i,j,N)] = 0; E[toIndex(i,j,N)] = 0.1; 
				P[toIndex(j,i,N)] = 0; E[toIndex(j,i,N)] = 0.1; 
			}
		}
	}
}

void readKappaFile(int numInteractions, double* kappa) {
	//read the file to get particle identities

	//file location
	std::string filename = "input/design/kappa.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return;
	}

	//store the entries
	double k; int count = 0; 
	while (in_str >> k) {
		//increment the counter for the particle number (1 based indexing)
		count++; 

		//if the particle number doesnt exceed the total, store it
		if (count <= numInteractions) {
			kappa[count-1] = k;
		}
	}

	//check if the number of entries in file matches the db we recieved
	if (count != numInteractions) {
		fprintf(stderr, "The supplied design file is incompatible with the database file.\n");
		fprintf(stderr, "You supplied %d strengths, but the system needs %d.\n", count, numInteractions);
		abort();
		return;
	}

}

int readDesignFile(int N, int* particleTypes) {
	//read the file to get particle identities

	int max = -1; //use to determine number of types
	int count = 0; //use to determine if input mathces the system

	//file location
	std::string filename = "input/design/particles.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return -1;
	}

	//store the entries
	int id; 
	while (in_str >> id) {
		//increment the counter for the particle number (1 based indexing)
		count++; 

		//if the particle number doesnt exceed the total, store it
		if (count <= N) {
			particleTypes[count-1] = id;
		}

		//get the max of the particle ids
		if (id > max) {
			max = id;
		}
	}

	//check if the number of entries in file matches the db we recieved
	if (count != N) {
		fprintf(stderr, "The supplied design file is incompatible with the database file.\n");
		fprintf(stderr, "You supplied %d types, but the system has %d particles.\n", count, N);
		abort();
		return -1;
	}

	//return number of particle types
	return max+1;
}


}