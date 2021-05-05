#include "database.h"
#include "nauty.h"
#include "hydro.h"
#include "sampling.h"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include "../defines.h"

#define HYDRO_CUT 1.06 //was 1.06 at last working value - 1.04 works for rho=50
namespace bd {

// ********************************************************************** //
// ********************* Class Functions ******************************** //
// ********************************************************************** //

HydroCluster::HydroCluster() {
	//constructor

	cNum = 0;
	
}

HydroCluster::~HydroCluster() {
	//deconstructor

	delete []clusters;
}

HCC::HCC(int num_clusters_, int N_, int maxT_) {
	//constructor

	num_clusters = num_clusters_;
	N = N_;
	maxT = maxT_;

	hc = new HydroCluster[num_clusters];

	for (int i = 0; i < num_clusters; i++) {
		(hc)[i].N = N_; (hc)[i].maxT = maxT_; 
		if (N_ == 6) {
			(hc)[i].cNum = 1;
		}
		if (N == 7) {
			(hc)[i].cNum = 0;
		}
		(hc)[i].clusters = new double[N_*DIMENSION*maxT_];
		for (int j = 0; j < N_* DIMENSION * maxT_; j++) {
			(hc)[i].clusters[j] = 0;
		}
	}
}

HCC::~HCC() {
	//deconstructor

	delete []hc;
}

HCC* extractData(std::string& filename, int N, int maxT) {
	//read data from HD file and put into a HC class

	//the last particle coordinates for n=7 have been deleted. ignore these
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return NULL;
	}

	//read first line, N = number of particles
	int num_clusters; int file_clusters; bool bug_flag = false;
	in_str >> file_clusters; num_clusters = file_clusters;
	if (N == 7) { //addresses the output bug for n=7
		//num_clusters = 3087;
	}
	if (file_clusters != num_clusters) {
		bug_flag = true;
	}
	//std::cout << "First line reads " << num_clusters << "\n";
	num_clusters /= N;
	if (bug_flag) {
		num_clusters--;
	}
	//std::cout << "Number of clusters is " << num_clusters << "\n";

	std::cout << "Number of clusters in file: " << num_clusters << "\n";

	//create an array of HC objects
	HCC* hc = new HCC(num_clusters, N, maxT);

	double x, y;
	double temp;
	int particle = 0;
	int cluster = 0;
	int timestep = 0;


	//fill the hcc with clusters
	while (in_str >> x) {
		//x is stored already, store y
		in_str >> y;

		//determine where in the cluster array these coordinates goes
		int placement = DIMENSION * particle + timestep * N * DIMENSION;
		(*hc)[cluster].clusters[placement] = x;
		(*hc)[cluster].clusters[placement+1] = y;

		//printf("particle %d, x is %f and y is %f\n", particle, x, y);
		

		//cycle through the useless data
		for (int i = 0; i < 5; i++) {
			in_str >> temp;
		}

		//update indices
		particle++;
		if (particle == N) {              //increment cluster by 1 every N particles
			cluster++; particle = 0;
			if (cluster == num_clusters) {  //increment timestep by 1 after all clusters
				timestep++; cluster = 0;
				if (bug_flag) {//skip last particle data - 6 lines
					for (int line = 0; line < 6; line++) {
						for (int i = 0; i < 7; i++) {
							in_str >> temp;
						}
					}
				}
				in_str >> temp;               //extra value when timestep changes
				if (timestep == maxT) {       //break if run out of storage
					break;                      
				}
			}
		}


		//std::cout << "Time step " << timestep << "\n";

	}

	//close the input file
	in_str.close();

	//return the collection
	return hc;
}

// ********************************************************************** //
// ********************* Transition Detection *************************** //
// ********************************************************************** //


void getCoordinates(HydroCluster& hc, double* X, int N, int time) {
	//place the coordinates of the cluster at given time into X

	int first_element = time * N * DIMENSION; 
	for (int i = 0; i < N * DIMENSION; i++) {
		X[i] = hc.clusters[first_element + i];
	}
}

void getAdj(double* X, int N, int* M, double cutoff) {
	//get the adjacnecy matrix from a cluster X - specify cutoff

	//init a particle array to 0, then fill it from X
	double* particles = new double[DIMENSION*N]; double* Z = new double[DIMENSION];
	for (int p = 0; p < DIMENSION*N; p++) {
		particles[p] = 0;
	}
	c2p(X, particles, N);

	/*
	for (int q = 0; q < 2*N; q++) {
		std::cout << particles[q] << "\n";
	}
	*/

	//compute the distance between particles, construct adj matrix
	int i, j; 
	for (int index = 0; index < N*N; index++) {
		index2ij(index, N, i, j);
		if (j > i) {
			double d = euDist(particles, i, j, N, Z);
			if (d < cutoff) {
				M[index] = 1; M[toIndex(j,i,N)] = 1;
			}
		}
	}

	delete []particles; delete []Z;
}

void checkState(int N, double* X, int state, Database* db,
							  bool& reset, int& new_state) {
	//check if a new state has been reached

	//declare reset false, only true if the sample is invalid
	reset = false;

	//get adjacency matrix of the current state
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(X, N, M, HYDRO_CUT);

	//compare old and new state, check if same
	int* old = new int[N*N]; for (int i = 0; i < N*N; i++) old[i]=0;
	int old_bonds = (*db)[state].getBonds();
	extractAM(N, state, old, db);
	int success = checkSame(old, M, N);

	int B = 0;
	for (int i = 0; i < N*N; i++) {
		B += M[i];
	}
	B /= 2;
	//std::cout << B << "\n";

	if (success) {//same matrix
		delete []M; delete []old;
		return;
	}
	else {//not the same, find matrix in database
		success = mcm::findMatrix(M, old, old_bonds, N, db,  reset,  new_state);
		if (!success) {
			//may output an unphysical state. refine with newton, check again
			refine(N, X, M);
			success = mcm::findMatrix(M, old, old_bonds, N, db,  reset,  new_state);
			if (!success) { //state still not found after refine, ignore this sample
				reset = true; 
				printf("State not found in database, B = %d \n\n\n\n\n", B);
			}
		}
		delete []M; delete []old;
	}
}

//maybe parralelize this later?
void determineTransitions(HCC* hc, Database* db, double tps) {
	//takes a collection of cluster trajectories, determines the folding pathway and
	//times, constructs an adjacency matrix and time distribution

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config
	int timer = 0;                          //store time since last bond
	//array of vectors to compute a probability distribution, and time distribution
	std::vector<int>* transition_times = new std::vector<int>[ns*ns];
	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;

	int folded = 0;

	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);
			//printf("Old state %d, new state %d\n", old_state, new_state);
			/*
			for (int x = 0; x < DIMENSION*N; x++) {
				std::cout << X[x] << ' ';
			}
			std::cout << "\n";
			*/
			


			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.
				printf("Trajectory %d is broken at timestep %d\n", i, time);
				printf("Current state is %d\n", new_state);
				for (int coord = 0; coord < N*DIMENSION; coord++) std::cout << X[coord] << "\n";
				break;
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				transition_times[toIndex(old_state, new_state, ns)].push_back(timer);
				old_state = new_state;
				timer = 0;

				//check if ground state has been reached
				bonds = (*db)[new_state].getBonds();
				if (bonds >= 2*N-3) {
					folded++;
					break;
				}
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}
		}
	}

	//transition_times now has all the data, put it into db
	//loop over transitions from every initial -> final state
	std::vector<Pair> PM;
	for (int initial_state = 0; initial_state < ns; initial_state++) {
		//compute the mfpt and number of transitions to each final state
		double mfpt = 0; double Z = 0;
		PM.clear();

		for (int final_state = 0; final_state < ns; final_state++) {
			//get the vector of transition times
			int index = toIndex(initial_state, final_state, ns);
			std::vector<int> transitions = transition_times[index];

			//get the vector size and add to mfpt if non-zero
			int num_transitions = transitions.size();
			if (num_transitions > 0) {
				PM.push_back(Pair(final_state, num_transitions));
				Z += num_transitions;
				for (int entry = 0; entry < num_transitions; entry++) {
					mfpt += transitions[entry];
				}
			}
		}

		//divide by Z to get mean time
		if (Z > 0) {
			mfpt /= Z;
		}

		//place data in the db
		//make a temp vector with same num of elements as P
		std::vector<bd::Pair> temp; 
		for (int i = 0; i < PM.size(); i++) {
			temp.push_back(bd::Pair(PM[i].index, 0));
		}

		//update database
		(*db)[initial_state].mfpt = mfpt*tps;
		(*db)[initial_state].num = 0;
		(*db)[initial_state].denom = 0;
		(*db)[initial_state].num_neighbors = PM.size();
		(*db)[initial_state].P = PM;
		(*db)[initial_state].Z = temp;
		(*db)[initial_state].Zerr = temp;
		(*db)[initial_state].sigma = 0;

	}

	printf("%d of the %d trajectories succesfully folded\n", folded, nc);

	//free the memory
	delete []X; delete []Xold; delete []transition_times;


}

void determineTransitionTimes(HCC* hc, Database* db, double tps, std::ostream& ofile) {
	//takes a collection of cluster trajectories, determines the folding pathway and
	//times, outputs times of transition

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config
	int timer = 0;                          //store time since last bond
	//array of vectors to compute a probability distribution, and time distribution
	std::vector<int>* transition_times = new std::vector<int>[nc];
	bool* broken = new bool[nc];
	for (int i = 0; i < nc; i++) broken[i] = false;
	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;

	int folded = 0;
	int errored = 0;

	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);
			//printf("Old state %d, new state %d\n", old_state, new_state);
			/*
			for (int x = 0; x < DIMENSION*N; x++) {
				std::cout << X[x] << ' ';
			}
			std::cout << "\n";
			*/
			


			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.

				//at first, check if problem gets fixed next time step
				getCoordinates(currentCluster, X, N, time+1);
				reset = false;
				checkState(N, X, old_state, db, reset, new_state);

				if (reset) { //this only gets called when two bonds form at once
					broken[i] = true;
					for (int i = 0; i < N*DIMENSION; i++) {
						std::cout << X[i] << "\n";
					}
					printf("Cluster %d errored. Ignoring trajectory\n", i);
					errored++;
					

					break;
				}
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				transition_times[i].push_back(time);
				old_state = new_state;
				timer = 0;

				//check if ground state has been reached
				bonds = (*db)[new_state].getBonds();
				if (bonds >= 2*N-3) {
					folded++;
					break;
				}
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}
		}
	}

	printf("Fully Folded Trajectories: %d of %d\n", folded, nc);
	printf("Errored Trajectories: %d of %d\n", errored, nc);
	printf("Incomplete Trajectories: %d of %d\n", nc - folded - errored, nc);

	//transition_times now has all the data, put it into db
	//loop over transitions from every initial -> final state
	for (int i = 0; i < nc; i++) {
		if (!broken[i]) {
			for (int t = 0; t < transition_times[i].size(); t++) {
				ofile << transition_times[i][t] << " ";
			}
			ofile << "\n";
		}
	}


	//free the memory
	delete []X; delete []Xold; delete []transition_times; delete []broken;


}

void determineTransitionStates(HCC* hc, Database* db, double tps, std::ostream& ofile) {
	//takes a collection of cluster trajectories, determines the folding pathway and
	//times, outputs times of transition

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config
	int timer = 0;                          //store time since last bond
	//array of vectors to compute a probability distribution, and time distribution
	std::vector<int>* transition_times = new std::vector<int>[nc];
	bool* broken = new bool[nc];
	for (int i = 0; i < nc; i++) broken[i] = false;
	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;

	int folded = 0;
	int errored = 0;

	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);
			//printf("Old state %d, new state %d\n", old_state, new_state);
			/*
			for (int x = 0; x < DIMENSION*N; x++) {
				std::cout << X[x] << ' ';
			}
			std::cout << "\n";
			*/
			


			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.
				//at first, check if problem gets fixed next time step
				getCoordinates(currentCluster, X, N, time+1);
				reset = false;
				checkState(N, X, old_state, db, reset, new_state);

				if (reset) { //this only gets called when two bonds form at once
					broken[i] = true;
					for (int i = 0; i < N*DIMENSION; i++) {
						std::cout << X[i] << "\n";
					}
					printf("Cluster %d errored. Ignoring trajectory\n", i);
					errored++;
					

					break;
				}
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				transition_times[i].push_back(db->lumpMap[new_state]);
				old_state = new_state;
				timer = 0;

				//check if ground state has been reached
				bonds = (*db)[new_state].getBonds();
				if (bonds >= 2*N-3) {
					folded++;
					break;
				}
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}
		}
	}

	printf("Fully Folded Trajectories: %d of %d\n", folded, nc);
	printf("Errored Trajectories: %d of %d\n", errored, nc);
	printf("Incomplete Trajectories: %d of %d\n", nc - folded - errored, nc);

	//transition_times now has all the data, put it into db
	//loop over transitions from every initial -> final state
	for (int i = 0; i < nc; i++) {
		if (!broken[i]) {
			for (int t = 0; t < transition_times[i].size(); t++) {
				ofile << transition_times[i][t] << " ";
			}
			ofile << "\n";
		}
	}


	//free the memory
	delete []X; delete []Xold; delete []transition_times; delete []broken;


}


// ********************************************************************** //
// ********************* Get statistics ********************************* //
// ********************************************************************** //

//this is also defined in sampling.cpp for other purposes
/*
double gyrationRadius(int N, double* X) {
	//compute the radius of gyration

	double* x = new double[N];
	double* y = new double[N];
	double* r2 = new double[N];

	double xCM = 0; double yCM = 0;

	for (int i = 0; i < N; i++) {
		x[i] = X[DIMENSION*i]; xCM += x[i];
		y[i] = X[DIMENSION*i+1]; yCM += y[i];
	}

	double rg = 0;
	for (int i = 0; i < N; i++) {
		r2[i] = (x[i]-xCM/N)*(x[i]-xCM/N) + (y[i]-yCM/N)*(y[i]-yCM/N);
		rg += r2[i];
	}

	rg /= N;

	//free memory
	delete []x; delete []y; delete []r2;

	return sqrt(rg);
}
*/

void timeAverageFHT(HCC* hc, Database* db, std::vector<double>& q, int which) {
	//compute time averages of the quantity q until the first hitting time

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config

	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;


	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);
		int timer = 0;     //store time since last bond
		double avg = 0;       //store running average

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			//get a sample of some quantity
			double quantity;
			if (which == 0) {
				quantity = gyrationRadius(N, X);
				//std::cout << quantity << "\n";
			}
			avg += quantity;

			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);


			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.
				break;
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				//std::cout << "Final " << avg / timer << "\n";
				q.push_back(avg/timer);
				break;
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}
		}
	}


	//free the memory
	delete []X; delete []Xold; 


}

void distributionFHT(HCC* hc, Database* db, std::vector<double>& q, int which) {
	//compute the quantity q at the first hitting time

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config

	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;


	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);
		int timer = 0;      //store time since last bond

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			
			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);


			
			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.
				//at first, check if problem gets fixed next time step
				getCoordinates(currentCluster, X, N, time+1);
				reset = false;
				checkState(N, X, old_state, db, reset, new_state);

				if (reset) { //this only gets called when two bonds form at once
					//broken[i] = true;
					for (int i = 0; i < N*DIMENSION; i++) {
						std::cout << X[i] << "\n";
					}
					printf("Cluster %d errored. Ignoring trajectory\n", i);
					//errored++;
					
					break;
				}
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				double quantity;
				if (which == 0) {
					quantity = gyrationRadius(N, X);
				}
				else if (which == 1) {
					quantity = boop2d(N, X);
				}
				else if (which == 2) {
					quantity = end2end(N, X);
				}
				q.push_back(quantity);
				break;
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}
			old_state = new_state;
		}
	}


	//free the memory
	delete []X; delete []Xold; 
}

void clustersFHT(HCC* hc, Database* db, std::ostream& ofile) {
	//write all clusters at first hitting time to file

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config

	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;


	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);
		int timer = 0;      //store time since last bond

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			
			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);


			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.
				break;
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				
				for (int p = 0; p < DIMENSION*N; p++) {
					ofile << X[p] << ' ';
				}
				ofile << "\n";
				break;
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}
			old_state = new_state;
		}
	}


	//free the memory
	delete []X; delete []Xold; 
}

void distributionFHT2(HCC* hc, Database* db, std::vector<double>& q, int which) {
	//compute the quantity q at the second hitting time

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config

	
	//variables for state change check
	int old_state;
	int new_state;
	bool reset;
	int bonds;


	//loop over each cluster in HHC
	for (int i = 0; i < nc; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d of %d\n", i+1, nc);
		int timer = 0;      //store time since last bond

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - increment timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) { //check if we reached the end of the time series
				//std::cout << time << "\n";
				break;
			}
			timer++;

			
			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);
			int new_bonds = (*db)[new_state].getBonds();


			if (reset) { //either a bond broke, or two bonds formed at once - ignore rest of traj.
				break;
			}

			if (new_state != old_state && new_bonds == 8) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				double quantity;
				if (which == 0) {
					quantity = gyrationRadius(N, X);
				}
				else if (which == 1) {
					quantity = boop2d(N, X);
				}
				else if (which == 2) {
					quantity = end2end(N, X);
				}
				q.push_back(quantity);
				break;
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
				old_state = new_state;
			}
		}
	}


	//free the memory
	delete []X; delete []Xold; 
}

void sampleStats(std::vector<double> X, double& M, double& V) {
	//return the sample mean of the vector X

	//init the mean and variance at 0, get sample size
	M = 0; V = 0;
	int N = X.size();

	//compute the mean
	for (int i = 0; i < N; i++) M += X[i];
	M /= float(N);

	//compute the variance 
	for (int i = 0; i < N; i++) V += (X[i]-M) * (X[i]-M);
	V /= (N-1);
}


// ********************************************************************** //
// ********************* Testing Area  ********************************** //
// ********************************************************************** //


void testExtract(HCC* hc) {
	//test if the extraction is working

	int test_cluster = 0;
	int time = 6;

	for (int i = 0; i < 12; i++) {
		std::cout << (*hc)[test_cluster].clusters[time*12+i] << "\n";
	}
}










}