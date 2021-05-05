#include "design.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <deque>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/Eigen/Dense>
#include <cstdlib>   // rand and srand
#include <ctime>     // For the time function
#include "point.h"
#include "pair.h"
#include "adjacency.h"
#include "database.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "../defines.h"
namespace bd{

/****************************************************/
/*********** General Functions    *******************/
/****************************************************/

void getBondTypes(int N, int* particleTypes, Database* db, std::vector<int> targets) {
	//for each state in targets, check how many of each bond type each state has
	//hard-coded for 2 types

	for (int index = 0; index < targets.size(); index++) {
		int AA = 0; int AB = 0; int BB = 0;
		int state = targets[index];
		for (int i = 0; i < N; i++) {
			for (int j = i+2; j < N; j++) {
				if ( (*db)[state].isInteracting(i,j,N)) {
			    //printf("%d %d ", i,j);
					int p1 = particleTypes[i];
					int p2 = particleTypes[j];
					if (p1 == 0 && p2 == 0) {
						AA++;
					}
					else if (p1 == 1 && p2 == 1) {
						BB++;
					}
					else {
						AB++;
					}
				}
			}
		}
		//printf("\n");
		printf("\n State: %d, AA: %d, AB %d, BB %d\n", state, AA, AB, BB);
		double* X = new double[N*DIMENSION];
		Cluster c = (*db)[state].getRandomIC();
		c.makeArray3d(X,N);
		printCluster(X,N);
		delete []X;
	}
}

void checkPerm(std::string word, std::deque<std::string>& perms) {
	//check if reflected and swapped word is in perms, if not add it

	//check reflections
	std::string rword = word; std::reverse(rword.begin(), rword.end());
	std::deque<std::string>::iterator itR = std::find(perms.begin(), perms.end(), rword);
	if (itR != perms.end()) {
		return;
	}

	//check swaps
	std::string swapWord = word; 
	for (int i = 0; i < word.length(); i++) {
		if (word[i] == '0') swapWord[i] = '1';
		if (word[i] == '1') swapWord[i] = '0';
	}
	std::deque<std::string>::iterator itS = std::find(perms.begin(), perms.end(), swapWord);
	if (itS != perms.end()) {
		return;
	}

	//check reflections of swaps
	std::string rword2 = swapWord; std::reverse(rword2.begin(), rword2.end());
	std::deque<std::string>::iterator itR2 = std::find(perms.begin(), perms.end(), rword2);
	if (itR2 != perms.end()) {
		return;
	}

	perms.push_back(word);
}

void allPerms(int N, std::deque<std::string>& perms) {
	//generate all perms of N 0s and 1s

	std::string p0 = "0"; std::string p1 = "1";
	perms.push_back(p0); perms.push_back(p1);
	int elements = 2;
	int final_size = pow(2,N); 
	while (elements < final_size) {
		std::string f = perms.front(); perms.pop_front();
		perms.push_back(f+p0); perms.push_back(f+p1);
		elements += 1;
	}
}

void distinctPerms(int N, std::deque<std::string>& perms) {
	//generate all permutations of N 0s and 1s w/o reflections

	//first generate all (N-1) length possibilities, fill preliminary deque
	std::string p0 = "0"; std::string p1 = "1";
	std::deque<std::string> prelim;
	prelim.push_back(p0); prelim.push_back(p1);
	int elements = 2;
	int final_size = pow(2,N-1); 
	while (elements < final_size) {
		std::string f = prelim.front(); prelim.pop_front();
		prelim.push_back(f+p0); prelim.push_back(f+p1);
		elements += 1;
	}

	//next, add 0 and 1 to each (N-1) length permuation, only add it to perms
	//if its reflection is not already in it
	for (int i = 0; i < prelim.size(); i++) {
		std::string f = prelim[i];
		checkPerm(f+p0, perms); checkPerm(f+p1, perms);
	}

	printf("Found %ld distinct permutations\n", perms.size());

}

void setTypes(int N, int* particleTypes, std::deque<std::string> perms, int perm) {
	//select a permutation from the perms vector and set particle types

	std::string p = perms[perm];
	for (int i = 0; i < N; i++) {
		particleTypes[i] = p[i]-'0';
		//std::cout << particleTypes[i] << "\n";
	}
}


int setTypes(int N, int* particleTypes, int IC) {
	//set which particles are which type

	int numTypes;

	if (IC == 0) { //all the same type
		for (int i = 0; i < N; i++) {
			particleTypes[i] = 0;
		}
		numTypes = 1;
	}
	else if (IC == 1) { //alternating type a and b
		for (int i = 0; i < N; i++) {
			particleTypes[i/2] = 0;
		}
		for (int i = 0; i < N; i++) {
			particleTypes[i/2+1] = 1;
		}
		numTypes = 2;
	}

	return numTypes;
}

double getStickyProduct(int N, int state,  Database* db, int* particleTypes, 
												std::map<std::pair<int,int>,double> kappa) {
	//return the product of sticky param ^ num bonds for each type

	double stickyProd = 1.0;

	//loop over the adjacency matrix, determine bond types, get factors
	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ((*db)[state].isInteracting(i,j,N)) {
				int p1 = particleTypes[i]; int p2 = particleTypes[j];
				stickyProd *= kappa[{p1,p2}];
			}
		}
	}

	return stickyProd;
}

void getReweightMaps(std::map<std::pair<int,int>,double> kappa,
										 std::map<std::pair<int,int>,double>& a_new,
										 std::map<std::pair<int,int>,double>& gamma_new) {
	//create seperate maps for a values and gamma values
	// a = U''(d) = 2*range^2*E_{ij}
	//gamma = exp(beta*E)

	double beta = BETA;                //inverse temp
	double r = RANGE;                  //range of potential

	//loop over kappa map and fill other two maps appropriately
	std::map<std::pair<int,int>,double>::iterator it;
	for (it = kappa.begin(); it != kappa.end(); it++) {
		std::pair<int,int> key = it->first;
		double k_val           = it->second;
		double E = stickyNewton(8.0, r, k_val, beta);
		double a = 2*r*r*E; double gam = exp(beta*E);
		//std::cout << k_val - gam/(sqrt(a)) << "\n";
		a_new[key] = a; gamma_new[key] = gam;
	}
}

double getPdet(int N, int state, Database* db, int* particleTypes, 
							 double a, std::map<std::pair<int,int>,double> a_new) {
	//get sqrt of ratio of pseudo determinant of dynamical matrix for cluster

	//get the number of bonds, adjacency matrix, and coordinates from db
	int b = (*db)[state].getBonds();
	int* M = new int[N*N]; 
	extractAM(N, state, M, db);
	const Cluster& c = (*db)[state].getRandomIC();
	//cluster structs to arrays
	double* X = new double[DIMENSION*N];
#if (DIMENSION == 2) 
	c.makeArray2d(X, N);
#endif
#if (DIMENSION == 3)
	c.makeArray3d(X, N);
#endif

	//init the jacobian matrix
	Eigen::MatrixXd J(b, 2*N); J.fill(0.0);
	Eigen::MatrixXd J_new(b, 2*N); J_new.fill(0.0);
	double m = sqrt(a);

	//init distances, bond counter
	double XD, YD, ZD;
	int count = 0;

	//loop over bonds and construct rigidity matrices, J and J_new
	for (int i = 0; i < N; i ++) {
		for (int j = i+1; j < N; j++) {
			if (M[bd::toIndex(i,j,N)] == 1) {
				//compute distances
				XD = X[DIMENSION*i] - X[DIMENSION*j];
				YD = X[DIMENSION*i+1] - X[DIMENSION*j+1];
#if (DIMENSION == 3) 
				ZD = X[DIMENSION*i+2] - X[DIMENSION*j+2];
#endif 

				//get particle types and row multiplier
				int p1 = particleTypes[i];
				int p2 = particleTypes[j];
				double m_new = sqrt(a_new[{p1,p2}]);

				//compute jacobian entries
				J(count, DIMENSION*i) = 2*XD*m; J(count, DIMENSION*j) = -2*XD*m;
				J(count, DIMENSION*i+1) = 2*YD*m; J(count, DIMENSION*j+1) = -2*YD*m;
				J_new(count, DIMENSION*i) = 2*XD*m_new; J_new(count, DIMENSION*j) = -2*XD*m_new;
				J_new(count, DIMENSION*i+1) = 2*YD*m_new; J_new(count, DIMENSION*j+1) = -2*YD*m_new;
#if (DIMENSION == 3) 
				J(count, DIMENSION*i+2) = 2*ZD*m; J(count, DIMENSION*j+2) = -2*ZD*m;
				J_new(count, DIMENSION*i+2) = 2*ZD*m_new; J_new(count, DIMENSION*j+2) = -2*ZD*m_new;
#endif

				count +=1;
			}
		}
	}

	//compute the pseudo-determinant of each dynamical matrix
	Eigen::VectorXd s = J.jacobiSvd().singularValues();
	Eigen::VectorXd s_new = J_new.jacobiSvd().singularValues();
	double det = 1; double det_new = 1;
  for (int i = 0; i < s.size(); i++) {
    if (abs(s(i)) > N_TOL) {  
      det *= abs(s(i))*abs(s(i));
    }
  }
  for (int i = 0; i < s_new.size(); i++) {
    if (abs(s_new(i)) > N_TOL) {  
      det_new *= abs(s_new(i))*abs(s_new(i));
    }
  }

	//take sqrt of ratio of determinants
	double D = sqrt(det/det_new);

	//free memory and return D
	delete []M; delete []X;

	return D;
}

void reweight7(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa) {
	//reweight in the case of 7 particles

	double kap0 = KAP;                 //sticky parameter for initial measurement
	double beta = BETA;                //inverse temp
	double r = RANGE;                  //range of potential
	double E = stickyNewton(8.0, r, kap0, beta); //energy corresponding to kap0
	double a = 2*r*r*E;                // parameter multiplying rigidity matrix
	double alpha2 = sqrt(beta*a);      //alpha^2 parameter
	double gamma = exp(beta*E);        //gamma parameter

	double Z = 0;                     //normalizing constant for eq

	//loop over states, get the equilibrium measure entry
	for (int i = 0; i < num_states; i++) {
		//get initial eq prob and number of bonds
		int b = (*db)[i].getBonds();
		int ntBonds = b - (N-1);               //number of non-trivial bonds
		double prob = (*db)[i].getFrequency();
		double new_prob;

		//determine which reweight to use
		if (b != 12) { //do standard re-weight
			//get the denominator of reweight factor - KAP^b_i
			double denom = pow(kap0, ntBonds);

			//get the numerator of the re-weight factor
			double num = getStickyProduct(N, i, db, particleTypes, kappa); 

			//compute the new probability
			new_prob = prob * num / denom;
		}
		else if (b == 12) { //special reweight
			//get the denominator of reweight factor - KAP^b_i
			double denom = pow(gamma, ntBonds);

			//get a and gamma parameters for each new kappa
			std::map<std::pair<int,int>,double> a_new;
			std::map<std::pair<int,int>,double> gamma_new;
			getReweightMaps(kappa, a_new, gamma_new);

			//get the numerator of reweight factor
			double num = getStickyProduct(N, i, db, particleTypes, gamma_new); 

			//get the pseudo-determinant factor
			double det = 0;
			det = getPdet(N, i, db, particleTypes, a, a_new);

			//debug the flower reweight
			//printf("prob = %f\n det = %f\n num = %f\n denom = %f\n\n", prob, det, num, denom);

			//compute new prob
			new_prob = prob * det * num / denom;
		}

		//debug line for corect reweight
		//printf("State %d, Old prob %f, New prob %f\n", i, prob, new_prob);

		//increment Z and add to array
		Z += new_prob;
		eq[i] = new_prob;
	}

	//re-normalize
	for (int i = 0; i < num_states; i++) eq[i] /= Z;

}
void reweight(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa) {
	//perform the re-weighting of the eq measure for new sticky params

	if (N == 6) {
		double kap0 = KAP;                 //sticky parameter for initial measurement
		double beta = BETA;                //inverse temp

		double Z = 0;                     //normalizing constant for eq

		//loop over states, get the equilibrium measure entry
		for (int i = 0; i < num_states; i++) {
			//get initial eq prob and number of bonds
			int b = (*db)[i].getBonds();
			int ntBonds = b - (N-1);               //number of non-trivial bonds
			double prob = (*db)[i].getFrequency();

			//get the denominator of reweight factor - KAP^b_i
			double denom = pow(KAP, ntBonds);

			//get the numerator of the re-weight factor
			double num = getStickyProduct(N, i, db, particleTypes, kappa); 
			//std::cout << num << "\n";

			//compute the new probability
			double new_prob = prob * num / denom;

			//debug line for corect reweight
			//printf("Test: Old %f, New %f\n", prob, new_prob);

			//increment Z and add to array
			Z += new_prob;
			eq[i] = new_prob;
		}

		//re-normalize
		for (int i = 0; i < num_states; i++) eq[i] /= Z;
	}
	else if (N == 7) {
		reweight7(N, num_states, db, particleTypes, eq, kappa);
	}
	
}

void initKappaVals(int numInteractions, double* kappaVals) {
	//initially set the valyes to all 1 for each interaction

	for (int i = 0; i < numInteractions; i++) {
		kappaVals[i] = 3;
	}
}

void checkPositive(int numInteractions, double* kappaVals) {
	//check if each kappa is positive, if not set to 0.01

	for (int i = 0; i < numInteractions; i++) {
		if (kappaVals[i] < 0.05) {
			kappaVals[i] = 0.05;
		}
	}
}

void applyMax(int numInteractions, double* kappaVals, double M) {
	//check if each kappa is below max, M. set to M-0.01 if not

	for (int i = 0; i < numInteractions; i++) {
		if (kappaVals[i] > M) {
			kappaVals[i] = M-0.01;
		}
	}
}

void constructSurfaceTOY(int N, Database* db, int initial, int target, bool useFile) {
	//plot the surface of hitting probabilities for initial to target states
	//toy model where 3 sticky parameters is hard coded - 2 unique

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* P = new double[num_states*num_states]; //probability transition matrix
	double* U = new double[num_states*num_states]; //hitting probability matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//declare outfile
	std::ofstream ofile;
	ofile.open("hittingProbSurface.txt");

	//do hitting probability calculation
	int M = 60; //num points in each dimension
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < M; y++) {
			//set kappa
			kappaVals[0] = 0.5+ (double)x / 4; 
			kappaVals[1] = 0.5+ (double)y / 4; 
			kappaVals[2] = 0.5+ (double)x / 4; 

			//make the map
			makeKappaMap(numTypes, kappaVals, kappa);
			//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
			//do rewieght
			reweight(N, num_states, db, particleTypes, eq, kappa);
			//reset transition matrices
			for (int i = 0; i < num_states*num_states; i++) {
				P[i] = U[i] = 0;
			}
			//copy Tconst into T
			std::copy(Tconst, Tconst+num_states*num_states, T);
			//fill in transposed entries such that T satisfies detailed balance
			satisfyDB(T, num_states, db, eq);
			//fill in diagonal with negative sum of row entries
			fillDiag(T, num_states);
			//use the filled rate matrix to compute probability matrix
			createProbabilityMatrix(T, num_states, P);
			//solve for hitting probabilities to endStates states
			computeHittingProbability(P, num_states, ground, U);
			//get hit prob to targets
			double p = getHitProbability(num_states, initial, targets, U);

			//output to file
			ofile << kappaVals[0] << ' ' << kappaVals[1] << ' ' << p << "\n";

		}
		std::cout << x << "\n";
	}

	//close file
	ofile.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
	delete []T; delete []P;  delete []U; delete []Tconst;
}

void constructScatterTOY(int N, Database* db, int initial, int target, bool useFile) {
	/*Construct a scatter plot of avg rate vs equilibrium probability using 
	  a 3d grid of kappa values. 
	  Keep track of the kappa values that lead to both the highest eq prob 
	  and the highest rate.                                         */


	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqScatter.txt");

	//construct an array of kappa vals to use
	int M;
	//multiples of some base value
	if (N == 6) {
		M = 33; //num points in each dimension
	}
	if (N == 7 || DIMENSION == 3) {
		M = 20;
	}
	double base = 0.05; //smallest kappa to use - 0.05
	double mult = 1.6; //multiplies base to get next value - 1.3
	if (DIMENSION == 3) {
		mult = 2.0;
	}
	if (N == 7) {
		base = 0.1;
	}
	double* Ks = new double[M]; Ks[0] = base;
	for (int i = 1; i < M; i++) {
		Ks[i] = Ks[i-1] * mult;
	}

	//store max values
	double eqM = 0; double rM = 0;
	double eRate;   double rateE;
	double* ek = new double[numInteractions];
	double* rk = new double[numInteractions];

	//do hitting probability calculation
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < M; y++) {
			for (int z = 0; z < M; z++) {
				//set kappa
				kappaVals[0] = Ks[x];
				kappaVals[1] = Ks[y];
				kappaVals[2] = Ks[z];

				//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";

				//make the map
				makeKappaMap(numTypes, kappaVals, kappa);
				//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
				//do rewieght
				reweight(N, num_states, db, particleTypes, eq, kappa);
				//get eq prob
				double eqProb = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
				//copy Tconst into T
				std::copy(Tconst, Tconst+num_states*num_states, T);
				//fill in transposed entries such that T satisfies detailed balance
				satisfyDB(T, num_states, db, eq);
				//fill in diagonal with negative sum of row entries
				fillDiag(T, num_states);
				//get the transition rate
				if (num_states > 400) { //sparse solve
					computeMFPTsSP(num_states, T, targets, m);
				}
				else { //dense solve
					computeMFPTs(num_states, T, targets, m);
				}
				
				double rate = 1/m[initial];

				//update maxima
				if (eqProb > eqM) {
					eqM = eqProb; eRate = rate;
					ek[0] = kappaVals[0]; ek[1] = kappaVals[1]; ek[2] = kappaVals[2]; 
				}
				if (rate > rM) {
					rM = rate; rateE = eqProb;
					rk[0] = kappaVals[0]; rk[1] = kappaVals[1]; rk[2] = kappaVals[2]; 
				}
				if (eqProb > 0.43 && eqProb < 0.51 && rate > 0.64) {
					//std::cout  << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";
				}

				//get the max index from eq
				int maxEqIndex; double maxEq = 0;
				for (int i = 0; i < num_states; i++) {
					if (eq[i] > maxEq) {
						maxEq = eq[i]; maxEqIndex = i;
					}
				}



				//output to file
				ofile << eqProb << ' ' << rate << ' ' << maxEqIndex << "\n";

			}
		}
		std::cout << x << "\n";
	}

	//print out the maxima and argmax
	printf("Max EQ is %f with %f, %f, %f. The rate is %f\n", eqM, ek[0], ek[1], ek[2], eRate);
	printf("Max rate is %f with %f, %f, %f. The EQ is %f\n", rM, rk[0], rk[1], rk[2], rateE);

	//close file
	ofile.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
	delete []T; delete []Tconst; delete []m; delete []Ks;
	delete []ek; delete []rk;
}

void constructScatterTOY1(int N, Database* db, int initial, int target) {
	/*Construct a scatter plot of avg rate vs equilibrium probability using 
	  a 1d grid of kappa values. 
	  Keep track of the kappa values that lead to both the highest eq prob 
	  and the highest rate.                                         */


	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int IC = 0; 
	int numTypes = setTypes(N, particleTypes, IC);
	int numInteractions = 1;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqScatter.txt");

	//construct an array of kappa vals to use
	int M;
	//multiples of some base value
	if (N == 6) {
		M = 33; //num points in each dimension
	}
	if (N == 7) {
		M = 25;
	}
	double base = 0.05; //smallest kappa to use - 0.05
	double mult = 1.6; //multiplies base to get next value - 1.3

	base = 14; mult = 1.1; M = 100;

	double* Ks = new double[M]; Ks[0] = base;
	for (int i = 1; i < M; i++) {
		Ks[i] = Ks[i-1] * mult;
	}

	//store max values
	double eqM = 0; double rM = 0;
	double eRate;   double rateE;
	double* ek = new double[numInteractions];
	double* rk = new double[numInteractions];

	//do hitting probability calculation
	for (int x = 0; x < M; x++) {
		//set kappa
		kappaVals[0] = Ks[x];

		//make the map
		assert(numTypes == 1);
		makeKappaMap(numTypes, kappaVals, kappa);
		//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
		//do rewieght
		reweight(N, num_states, db, particleTypes, eq, kappa);
		//get eq prob
		double eqProb = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
		//copy Tconst into T
		std::copy(Tconst, Tconst+num_states*num_states, T);
		//fill in transposed entries such that T satisfies detailed balance
		satisfyDB(T, num_states, db, eq);
		//fill in diagonal with negative sum of row entries
		fillDiag(T, num_states);
		//get the transition rate
		computeMFPTs(num_states, T, targets, m);
		double rate = 1/m[initial];

		//update maxima
		if (eqProb > eqM) {
			eqM = eqProb; eRate = rate;
			ek[0] = kappaVals[0]; ek[1] = kappaVals[1]; ek[2] = kappaVals[2]; 
		}
		if (rate > rM) {
			rM = rate; rateE = eqProb;
			rk[0] = kappaVals[0]; rk[1] = kappaVals[1]; rk[2] = kappaVals[2]; 
		}
		if (eqProb > 0.43 && eqProb < 0.51 && rate > 0.64) {
			//std::cout  << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";
		}


		//output to file
		ofile << eqProb << ' ' << rate << "\n";
		std::cout << x << "\n";
	}

	//print out the maxima and argmax
	printf("Max EQ is %f with %f, %f, %f. The rate is %f\n", eqM, ek[0], ek[1], ek[2], eRate);
	printf("Max rate is %f with %f, %f, %f. The EQ is %f\n", rM, rk[0], rk[1], rk[2], rateE);

	//close file
	ofile.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
	delete []T; delete []Tconst; delete []m; delete []Ks;
	delete []ek; delete []rk;
}

void perturbDB(Database* db, double freq_perturb_frac, double rate_perturb_frac) {
	//perturb each estimate in the database

	int num_states = db->getNumStates();
	double Z = 0; //new normalizer

	// Get the system time.
	unsigned seed = time(0);
	 
	// Seed the random number generator.
	srand(seed);

	for (int i = 0; i < num_states; i++) {
		//get the original values
		double original_rate = (*db)[i].mfpt;
		double original_eq   = (*db)[i].freq;

		//get some noise terms
		double noise1 = 2.0 * (double(rand()) / RAND_MAX) - 1.0;
		double noise2 = 2.0 * (double(rand()) / RAND_MAX) - 1.0;

		std::cout << noise1 << ' ' << noise2 << "\n";

		//perturb the estimates
		double new_rate = original_rate * (1+noise1*rate_perturb_frac);
		double new_eq   = original_eq   * (1+noise2*freq_perturb_frac);
		double rate_diff = 100.0 * (original_rate-new_rate)/original_rate;
		double eq_diff   = 100.0 * (original_eq-new_eq)/original_eq;

		printf("Old mfpt: %f, new mfpt %f, diff %f %% \n", original_rate, new_rate, rate_diff);
		printf("Old freq: %f, new freq %f, diff %f %% \n", original_eq, new_eq, eq_diff);
		if (new_rate > 0) {
			(*db)[i].mfpt = new_rate;
		}
		else {
			(*db)[i].mfpt = 0.01;
		}
		if (new_eq > 0) {
			(*db)[i].freq = new_eq;
			Z += new_eq;
		}

	}

	for (int i = 0; i < num_states; i++) {
		(*db)[i].freq /= Z;
	}
}

void perturbDBconstP(Database* db) {
	//perturb the database by making all transitions uniformly likely among possibilities

	int num_states = db->getNumStates();

	for (int i = 0; i < num_states; i++) {
		std::vector<Pair> P = (*db)[i].getP();
		for (int j = 0; j < P.size(); j++) {
			P[j].value = 1.0; 
		}
		(*db)[i].P = P;
	}
}

void constructScatterTOYsensitivity(int N, Database* db, int initial, int target, bool useFile) {
	/*Construct a scatter plot of avg rate vs equilibrium probability using 
	  a 3d grid of kappa values. 

	  Keep track of the kappa values that lead to both the highest eq prob 
	  and the highest rate.      

	  This function perturbs the rates or eq probs to determine stability to 
	  sampling error.                            
	*/


	//get database info
	int num_states = db->getNumStates(); 

	//begin by perturbing the data in the database
	double freq_perturb_frac = 0.25; //perturb the eq prob data by up to this fraction
	double rate_perturb_frac = 0.25; //perturb the rate data by up to this fraction
	perturbDB(db, freq_perturb_frac, rate_perturb_frac);
	//perturbDBconstP(db);


	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqScatter.txt");

	//construct an array of kappa vals to use
	int M;
	//multiples of some base value
	if (N == 6) {
		M = 33; //num points in each dimension
	}
	if (N == 7) {
		M = 25;
	}
	double base = 0.05; //smallest kappa to use - 0.05
	double mult = 1.6; //multiplies base to get next value - 1.3
	double* Ks = new double[M]; Ks[0] = base;
	for (int i = 1; i < M; i++) {
		Ks[i] = Ks[i-1] * mult;
	}

	//store max values
	double eqM = 0; double rM = 0;
	double eRate;   double rateE;
	double* ek = new double[numInteractions];
	double* rk = new double[numInteractions];

	//do hitting probability calculation
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < M; y++) {
			for (int z = 0; z < M; z++) {
				//set kappa
				kappaVals[0] = Ks[x];
				kappaVals[1] = Ks[y];
				kappaVals[2] = Ks[z];

				//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";

				//make the map
				makeKappaMap(numTypes, kappaVals, kappa);
				//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
				//do rewieght
				reweight(N, num_states, db, particleTypes, eq, kappa);
				//get eq prob
				double eqProb = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
				//copy Tconst into T
				std::copy(Tconst, Tconst+num_states*num_states, T);
				//fill in transposed entries such that T satisfies detailed balance
				satisfyDB(T, num_states, db, eq);
				//fill in diagonal with negative sum of row entries
				fillDiag(T, num_states);
				//get the transition rate
				computeMFPTs(num_states, T, targets, m);
				double rate = 1/m[initial];

				//update maxima
				if (eqProb > eqM) {
					eqM = eqProb; eRate = rate;
					ek[0] = kappaVals[0]; ek[1] = kappaVals[1]; ek[2] = kappaVals[2]; 
				}
				if (rate > rM) {
					rM = rate; rateE = eqProb;
					rk[0] = kappaVals[0]; rk[1] = kappaVals[1]; rk[2] = kappaVals[2]; 
				}
				if (eqProb > 0.43 && eqProb < 0.51 && rate > 0.64) {
					//std::cout  << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";
				}

				//get the max index from eq
				int maxEqIndex; double maxEq = 0;
				for (int i = 0; i < num_states; i++) {
					if (eq[i] > maxEq) {
						maxEq = eq[i]; maxEqIndex = i;
					}
				}



				//output to file
				ofile << eqProb << ' ' << rate << ' ' << maxEqIndex << "\n";

			}
		}
		std::cout << x << "\n";
	}

	//print out the maxima and argmax
	printf("Max EQ is %f with %f, %f, %f. The rate is %f\n", eqM, ek[0], ek[1], ek[2], eRate);
	printf("Max rate is %f with %f, %f, %f. The EQ is %f\n", rM, rk[0], rk[1], rk[2], rateE);

	//close file
	ofile.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
	delete []T; delete []Tconst; delete []m; delete []Ks;
	delete []ek; delete []rk;
}

/****************************************************/
/******** Hitting Probability Optimization   ********/
/****************************************************/

double getHitProbability(int num_states, int initial, std::vector<int> targets, double* U) {
	//get the hitting probability for the group of states in targets

	double prob = 0; 
	for (int i = 0; i < targets.size(); i++) {
		double p = U[toIndex(initial, targets[i], num_states)];
		prob += p;
	}

	return prob;
}

double solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
													int numTypes, double* Tconst, std::vector<int> ground, 
													std::vector<int> targets) {
	//set up and solve the absorbtion probability problem

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* P = new double[num_states*num_states]; //probability transition matrix
	double* U = new double[num_states*num_states]; //hitting probability matrix
	double* eq = new double[num_states];           //equilibrium measure

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the hitting probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);
	//reset transition matrices
	for (int i = 0; i < num_states*num_states; i++) {
		P[i] = U[i] = 0;
	}
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);
	//use the filled rate matrix to compute probability matrix
	createProbabilityMatrix(T, num_states, P);
	//solve for hitting probabilities to endStates states
	computeHittingProbability(P, num_states, ground, U);
	//get hit prob to targets
	double p = getHitProbability(num_states, initial, targets, U);

	delete []T; delete []P; delete []U; delete []eq;

	return p;
}



double computeGradHP(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, double* Tconst, std::vector<int> ground, 
								 std::vector<int> targets, double* g) {
	//compute the gradient of the hitting probability in kappa
	//returns the hitting probability at current point

	//declare parameters
	double h = 1e-5; //step size for finite differences

	//get the hitting probability at the current kappa vals
	double f0 = solveAbsorbProblem(initial, kappaVals, db, particleTypes, numTypes,
																 Tconst, ground, targets);

	//std::cout << "Hitting Prob = " << f0 << "\n";

	//take a step h in each direction and compute new hitting prob
	for (int k = 0; k < numInteractions; k++) {
		if (kappaVals[k] < 0.011) {
			g[k] = 0;
		}
		else {
			kappaVals[k] += h;
			double f1 = solveAbsorbProblem(initial, kappaVals, db, particleTypes, numTypes,
																		 Tconst, ground, targets);
			kappaVals[k] -= h;

			g[k] = (f1 - f0) / h;
		}
	}

	return f0;
}


void hittingProbMaxTOYperms(int N, Database* db, int initial, int target) {
	//use optimization scheme to get max hitting probability for target in toy model
	//does optimization over all particle perms and outputs max probability

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 2000; double tol = 2e-3; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//loop over all permutations
	std::deque<std::string> perms; 
	//allPerms(N,perms);
	distinctPerms(N, perms);
	int num_perms = perms.size();
	double* permProb = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double hit;

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			hit = computeGradHP(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												Tconst,  ground, targets, g);

			//update kappa - set to 0.1 if it goes negative
			for (int k = 0; k < numInteractions; k++) {
				kappaVals[k] += step*g[k];
				if (kappaVals[k] < 0) {
					kappaVals[k] = 0.01;
				}
			}

			//check for gradient tolerance
			double res = 0;
			for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
			std::cout << "Iteration " << it << ", Res: " << res << "\n";
			if (res < tol) {
				break;
			}
		}

		permProb[p] = hit;
	}

	//print out the hitting probability with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Hitting Prob: %f, ID: %s\n", permProb[p], perms[p].c_str());
	}

	//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []g; delete []permProb;
}

void hittingProbMaxTOY(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max hitting probability for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 2000; double tol = 2e-3; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//get the permutation
	initKappaVals(numInteractions, kappaVals);
	double hit;

	//optimization - steepest ascent
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of kappa
		hit = computeGradHP(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
											Tconst,  ground, targets, g);

		//update kappa - set to 0.1 if it goes negative
		for (int k = 0; k < numInteractions; k++) {
			kappaVals[k] += step*g[k];
			if (kappaVals[k] < 0) {
				kappaVals[k] = 0.01;
			}
		}

		//check for gradient tolerance
		double res = 0;
		for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
		std::cout << "Iteration " << it << ", Res: " << res << "\n";
		if (res < tol) {
			break;
		}
	}

	std::cout << "Maximum Probability: " << hit << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	double R = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
	printf("The average transition rate with this kappa is %f\n", R);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []g; 
}



/****************************************************/
/******** Equilinrium Probability Optimization   ****/
/****************************************************/

double printEqProb(int initial, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, std::vector<int> targets) {
	//do re-weighting and get the eqiulibrium probability sum of all target states

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* eq = new double[num_states];           //equilibrium measure

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the eq probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//get the new equilibrium probability of target state
	double prob = 0;
	for (int i = 0;  i< targets.size(); i++) {
		printf("State %d, Eq prob %f\n", targets[i], eq[targets[i]]);
		prob += eq[targets[i]];
	}

	//free memory
	delete []eq;

	return prob;
}

double getEqProb(int initial, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes,  std::vector<int> targets) {
	//do re-weighting and get the eqiulibrium probability sum of all target states

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* eq = new double[num_states];           //equilibrium measure

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the eq probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//get the new equilibrium probability of target state
	double prob = 0;
	for (int i = 0;  i< targets.size(); i++) {
		prob += eq[targets[i]];
	}

	//free memory
	delete []eq;

	return prob;
}

double computeGradEQ(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, std::vector<int> targets, double* g) {
	//compute the gradient of the eqiluibrium probability in kappa
	//returns the eq probability at current point

	//declare parameters
	double h = 1e-6; //step size for finite differences

	//get the hitting probability at the current kappa vals
	double f0 = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);

	//std::cout << "Hitting Prob = " << f0 << "\n";

	//take a step h in each direction and compute new hitting prob
	for (int k = 0; k < numInteractions; k++) {
		if (kappaVals[k] < 0.011) {
			g[k] = 0;
		}
		else if (kappaVals[k] > 100000) {
			g[k] = 0;
		}
		else {
			kappaVals[k] += h;
			double f1 = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
			kappaVals[k] -= h;

			g[k] = (f1 - f0) / h;
		}
	}

	return f0;
}

void eqProbMaxTOY(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max equilibrium probability for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 120000; double tol = 1e-6; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//get the permutation
	initKappaVals(numInteractions, kappaVals);
	//readKappaFile(numInteractions, kappaVals);
	double eq;

	//optimization - steepest ascent
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of kappa
		eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
											targets, g);

		//update kappa - set to 0.1 if it goes negative
		for (int k = 0; k < numInteractions; k++) {
			kappaVals[k] += step*g[k];
			if (kappaVals[k] < 0) {
				kappaVals[k] = 0.01;
			}
			else if (kappaVals[k] > 100) {
				kappaVals[k] = 99.9;
			}
		}

		//check for gradient tolerance
		double res = 0;
		for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
		std::cout << "Iteration " << it << ", Res: " << res << "\n";
		if (res < tol) {
			break;
		}
	}

	std::cout << "Maximum Eq Probability: " << eq << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; 
}

void eqProbMaxTOYperms(int N, Database* db, int initial, int target) {
	//use optimization scheme to get max hitting probability for target in toy model
	//does optimization over all particle perms and outputs max probability

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 60000; double tol = 1e-6; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//loop over all permutations
	std::deque<std::string> perms; 
	//allPerms(N, perms);
	distinctPerms(N,perms);
	int num_perms = perms.size();
	double* permProb = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double eq;

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												targets, g);

			//update kappa - set to 0.1 if it goes negative
			for (int k = 0; k < numInteractions; k++) {
				kappaVals[k] += step*g[k];
				if (kappaVals[k] < 0) {
					kappaVals[k] = 0.01;
				}
				else if (kappaVals[k] > 100) {
				kappaVals[k] = 99.9;
				}
			}

			//check for gradient tolerance
			double res = 0;
			for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
			//std::cout << "Iteration " << it << ", Res: " << res << "\n";
			if (res < tol) {
				break;
			}
		}

		permProb[p] = eq;
	}

	//print out the hitting probability with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Eq Prob: %f, ID: %s\n", permProb[p], perms[p].c_str());
	}

	//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; delete []permProb;
}


/****************************************************/
/******** Transition Rate  Optimization   ***********/
/****************************************************/

double getRate(int initial, double* kappaVals, Database* db, int* particleTypes, 
							 int numTypes, double* Tconst, std::vector<int> targets) {
	//get transition rate as 1/mfpt to targets

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//initialize the transition rate matrix
	double* T = new double[num_states*num_states];
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);

	//init array for equilibrium distribution and compute it
	double* eq = new double[num_states];
	//declare and fill kappa map
	std::map<std::pair<int,int>,double> kappa;
	makeKappaMap(numTypes, kappaVals, kappa);   
	//do reweight to fill eq
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);

	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);

	//construct array of mfpts and caluclate it
	double* m = new double[num_states];
	computeMFPTsSP(num_states, T, targets, m);
	double mfpt = m[initial];

	//free memory
	delete []T; delete []eq; delete []m;

	return 1.0/mfpt;
}

double getRateTPT(int initial, double* kappaVals, Database* db, int* particleTypes, 
							 		int numTypes, double* Tconst, std::vector<int> targets) {
	//do re-weighting and get the transition rate to the target state - reactive traj

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//initialize the transition rate matrix
	double* T = new double[num_states*num_states];
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);

	//init array for equilibrium distribution and compute it
	double* eq = new double[num_states];
	//declare and fill kappa map
	std::map<std::pair<int,int>,double> kappa;
	makeKappaMap(numTypes, kappaVals, kappa);   
	//do reweight to fill eq
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);

	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);

	//solve for the committor
	//initialize committor in q
	double* q = new double[num_states];
	for (int i = 0; i < num_states; i++) q[i] = 0;
	
	//solve dirichlet problem for committor, q 
	computeCommittor(q, T, num_states, initial, targets);

	//get the average transition rate
	double R = computeTransitionRateTPT(num_states, q, T, eq);

	//free memory 
	delete []T; delete []eq; delete []q;

	return R;
}

double computeGradRate(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, double* Tconst, std::vector<int> targets, 
								 double* g) {
	//compute the gradient of the avg transition rate in kappa
	//returns the avg transition rate at current point

	//declare parameters
	double h = 1e-5; //step size for finite differences

	//get the hitting probability at the current kappa vals
	double f0 = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);

	//std::cout << "Hitting Prob = " << f0 << "\n";

	//take a step h in each direction and compute new hitting prob
	for (int k = 0; k < numInteractions; k++) {
		if (kappaVals[k] < 0.011) {
			g[k] = 0;
		}
		else {
			kappaVals[k] += h;
			double f1 = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
			kappaVals[k] -= h;

			g[k] = (f1 - f0) / h;
		}
	}

	return f0;
}



void rateMaxTOY(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max rate for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 3000; double tol = 1e-5; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);
	double R;

	//optimization - steepest ascent
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of kappa
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
											Tconst, targets, g);

		//update kappa - set to 0.1 if it goes negative
		for (int k = 0; k < numInteractions; k++) {
			kappaVals[k] += step*g[k];
			if (kappaVals[k] > 10) {
				kappaVals[k] = 1000;
			}
		}
		checkPositive(numInteractions, kappaVals);
		//applyMax(numInteractions, kappaVals, 100);



		//check for gradient tolerance
		double res = 0;
		for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
		std::cout << "Iteration " << it << ", Res: " << res << "\n";
		if (res < tol) {
			break;
		}
	}

	std::cout << "Maximum Rate: " << R << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	double eq = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
	double hp = solveAbsorbProblem(initial, kappaVals, db, particleTypes, numTypes, Tconst, 
																 ground, targets);
	printf("The equilibrium probability at this kappa is %f\n", eq);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; delete []Tconst;
}

void rateMaxTOYperms(int N, Database* db, int initial, int target) {
	//use optimization scheme to get max rate for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 5000; double tol = 1e-10; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//loop over all permutations
	std::deque<std::string> perms; 
	//allPerms(N,perms);
	distinctPerms(N, perms);
	int num_perms = perms.size();
	double* permRate = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double R;

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												Tconst, targets, g);

			//update kappa - set to 0.1 if it goes negative
			for (int k = 0; k < numInteractions; k++) {
				kappaVals[k] += step*g[k];
			}
			checkPositive(numInteractions, kappaVals);
			applyMax(numInteractions, kappaVals, 100);

			//check for gradient tolerance
			double res = 0;
			for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
			//std::cout << "Iteration " << it << ", Res: " << res << "\n";
			if (res < tol) {
				break;
			}
		}

		permRate[p] = R;
	}

	//print out the hitting probability with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Average Rate: %f, ID: %s\n", permRate[p], perms[p].c_str());
	}

	//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; delete []permRate; delete []Tconst;
}


/****************************************************/
/*** Constrained Optimization Problems **************/
/****************************************************/

void applyRateConstraint(double R, double c, double r, int numInteractions, 
												 double* gH, double* gR) {
	//determine if the constraint needs to be applied and return appropriate gradient

	if (R < c) {//constraint active
		for (int i = 0; i < numInteractions; i++) gH[i] += r*gR[i];
	}

}


double lineSearchHit(int initial, int numInteractions, double* kappaVals, Database* db, 
										 int* particleTypes, int numTypes, double* Tconst, 
										 std::vector<int> ground, std::vector<int> targets, 
										 double c, double r, double H, double R, double* g, double& step) {
	//perform a line search to get a step size for the objective fn
	//step is passed by reference, returns difference in obj fn

	//define initial step size and reducing factor
	double step0 = 2; double alpha = 0.5;

	//see if the constraint needs to be applied
	int apply = 0;
	if (R < c) {
		apply = 1;
	}

	//evaluate objective fn
	double obj = H + r*apply*(R-c);

	//printf("Original Objective: %f\n", obj);

	//declare array for test kappaVals
	double* testKappa = new double[numInteractions];


	//perform line searches
	double objNew; double delta = -1; 
	int maxReductions = 20; 
	for (int i = 0; i < maxReductions; i++) {
		//get the new kappa values
		step0 *= alpha;
		for (int k = 0; k < numInteractions; k++) testKappa[k] = kappaVals[k] + g[k]*step0;
		checkPositive(numInteractions, testKappa);

		//evaluate the objective
		double Hnew = solveAbsorbProblem(initial, testKappa, db, particleTypes, numTypes, Tconst, 
																	 ground, targets);
		double Rnew = getRate(initial, testKappa, db, particleTypes, numTypes, Tconst, targets);
		int applyNew = 0;
		if (Rnew < c) {
			applyNew = 1;
		}
		objNew = Hnew + r*applyNew*(Rnew-c);
		//printf("New obj: %f\n", objNew);

		delta = objNew - obj;
		if (delta > 0) {
			break;
		}
	}

	step = step0;


	//free memory
	delete []testKappa;

	return delta;
}

double lineSearchEq(int initial, int numInteractions, double* kappaVals, Database* db, 
										int* particleTypes, int numTypes, double* Tconst, 
										std::vector<int> targets, double c, double r, double H, double R, 
										double* g, double& step) {
	//perform a line search to get a step size for the objective fn
	//step is passed by reference, returns difference in obj fn

	//define initial step size and reducing factor
	double step0 = 2; double alpha = 0.5;

	//see if the constraint needs to be applied
	int apply = 0;
	if (R < c) {
		apply = 1;
	}

	//evaluate objective fn
	double obj = H + r*apply*(R-c);

	//printf("Original Objective: %f\n", obj);

	//declare array for test kappaVals
	double* testKappa = new double[numInteractions];


	//perform line searches
	double objNew; double delta = -1; 
	int maxReductions = 20; 
	for (int i = 0; i < maxReductions; i++) {
		//get the new kappa values
		step0 *= alpha;
		for (int k = 0; k < numInteractions; k++) testKappa[k] = kappaVals[k] + g[k]*step0;
		checkPositive(numInteractions, testKappa);

		//evaluate the objective
		double Enew = getEqProb(initial, testKappa, db, particleTypes, numTypes, targets);
		double Rnew = getRate(initial, testKappa, db, particleTypes, numTypes, Tconst, targets);
		int applyNew = 0;
		if (Rnew < c) {
			applyNew = 1;
		}
		objNew = Enew + r*applyNew*(Rnew-c);
		//printf("New obj: %f\n", objNew);

		delta = objNew - obj;
		if (delta > 0) {
			break;
		}
	}

	step = step0;


	//free memory
	delete []testKappa;

	return delta;
}

double lineSearchRate(int initial, int numInteractions, double* kappaVals, Database* db, 
											int* particleTypes, int numTypes, double* Tconst, 
											std::vector<int> targets, double c, double r, double eq, double R, 
											double* g, double& step) {
	//perform a line search to get a step size for the objective fn
	//step is passed by reference, returns difference in obj fn

	//define initial step size and reducing factor
	double step0 = 2; double alpha = 0.5;

	//evaluate objective fn
	double obj = R - r * (eq - c) * (eq - c);

	//printf("Original Objective: %f\n", obj);

	//declare array for test kappaVals
	double* testKappa = new double[numInteractions];


	//perform line searches
	double objNew; double delta = -1; 
	int maxReductions = 20; 
	for (int i = 0; i < maxReductions; i++) {
		//get the new kappa values
		step0 *= alpha;
		for (int k = 0; k < numInteractions; k++) testKappa[k] = kappaVals[k] + g[k]*step0;
		checkPositive(numInteractions, testKappa);

		//evaluate the objective
		double Enew = getEqProb(initial, testKappa, db, particleTypes, numTypes, targets);
		double Rnew = getRate(initial, testKappa, db, particleTypes, numTypes, Tconst, targets);
		
		objNew = Rnew - r * (Enew - c) * (Enew - c);
		//printf("New obj: %f\n", objNew);

		delta = objNew - obj;
		if (delta > 0) {
			break;
		}
	}

	step = step0;


	//free memory
	delete []testKappa;

	return delta;
}

void eqProbMaxTOYc(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max hitting probability for target in toy model
	//has constraint for transition rate
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 5000; double objTol = 1e-6; 
	double c = 0.8; //rate must be > c (target dependent)
	double r = 200;    //cost 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* gH = new double[numInteractions];
	double* gR = new double[numInteractions];

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);
	double eq;
	double R; 

	//optimization - steepest ascent w/ line search
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of hitting prob and rate
		eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												targets, gH);
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												Tconst, targets, gR);

		double obj = eq;
		if (R < c) {
			obj += r*(R-c);
		}


		//apply the rate penalty if needed
		applyRateConstraint(R, c, r, numInteractions, gH, gR);

		//find a step size to give decrease in obj fn
		double step = 0; 
		double delta = lineSearchEq(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
															Tconst, targets, c, r, eq, R, gH, step);

		//take the step
		for (int i = 0; i < numInteractions; i++) {
			kappaVals[i] += step * gH[i];
		}
		checkPositive(numInteractions, kappaVals);

		//output progress
		printf("Iteration: %d, OBJ: %f, Delta: %f\n", it, obj, delta);


		//check for objective function tolerance 
		if (delta < objTol) {
			break;
		}
	}

	std::cout << "Maximum Equilibrium Probability: " << eq << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	double Rf = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
	printf("The average transition rate with this kappa is %f\n", Rf);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []gH; delete []gR; 
}

void hittingProbMaxTOYc(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max hitting probability for target in toy model
	//has constraint for transition rate
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 2000; double objTol = 1e-5; 
	double c = 0.3; //rate must be > c (target dependent)
	double r = 200;    //cost 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* gH = new double[numInteractions];
	double* gR = new double[numInteractions];

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);
	double hit;
	double R; 

	//optimization - steepest ascent w/ line search
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of hitting prob and rate
		hit = computeGradHP(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
											Tconst,  ground, targets, gH);
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
												Tconst, targets, gR);

		double obj = hit;
		if (R < c) {
			obj += r*(R-c);
		}


		//apply the rate penalty if needed
		applyRateConstraint(R, c, r, numInteractions, gH, gR);

		//find a step size to give decrease in obj fn
		double step = 0; 
		double delta = lineSearchHit(initial, numInteractions, kappaVals, db, particleTypes, 
															numTypes, Tconst, ground, targets, c, r, hit, R, gH, step);

		//take the step
		for (int i = 0; i < numInteractions; i++) {
			kappaVals[i] += step * gH[i];
		}
		checkPositive(numInteractions, kappaVals);

		//output progress
		printf("Iteration: %d, OBJ: %f, Delta: %f\n", it, obj, delta);


		//check for objective function tolerance 
		if (delta < objTol) {
			break;
		}
	}

	std::cout << "Maximum Probability: " << hit << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	double Rf = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
	printf("The average transition rate with this kappa is %f\n", Rf);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []gH; delete []gR; 
}

void rateMaxTOYc(int N, Database* db, double c, int initial, int target, bool useFile) {
	//use optimization scheme to get max rate for target in toy model
	//has equality constraint for equilibrium prob
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 5000; double objTol = 1e-8; 
	double r = 200;    //cost 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* gH = new double[numInteractions];
	double* gR = new double[numInteractions];

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);
	double eq;
	double R; 

	//optimization - steepest ascent w/ line search
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of hitting prob and rate
		eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												targets, gH);
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												Tconst, targets, gR);

		//define the objevtive function
		double obj = R - r * (eq - c) * (eq-c);  //maybe use quadratic penalty

		//get the gradient of obj
		for (int i = 0; i < numInteractions; i++) {
			gR[i] -= 2 * r * (eq - c) * gH[i]; 
		}

		//find a step size to give decrease in obj fn
		double step = 0; 
		double delta = lineSearchRate(initial, numInteractions, kappaVals, db, particleTypes, 
															numTypes, Tconst, targets, c, r, eq, R, gR, step);

		//take the step
		for (int i = 0; i < numInteractions; i++) {
			kappaVals[i] += step * gR[i];
		}
		checkPositive(numInteractions, kappaVals);

		//output progress
		printf("Iteration: %d, OBJ: %f, Delta: %f\n", it, obj, delta);


		//check for objective function tolerance 
		if (delta < objTol) {
			break;
		}
	}

	std::cout << "Maximum Transition Rate: " << R << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	printf("The eq prob with this kappa is %f\n", eq);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []gH; delete []gR; 
}

/****************************************************/
/*********** Some sampling stuff *******************/
/****************************************************/

int findMatrix(int N, double* X, Database* db) {
	//find the adjacency matrix for state X in the database

	//construct the AM
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(X,N,M);
	//refine(N, X, M); //refines X and fills M
	int* dbM = new int[N*N];

	//define the return variable to be -1 if not found
	int state = -1;

	//loop over the database to search
	for (int i = 0; i < db->getNumStates(); i++) {
		extractAM(N, i, dbM, db);
		bool S = checkSame(dbM, M, N);
		if (S) { //found matrix
			state = i;
			break;
		}
	}

	//free memory
	delete []M; delete []dbM;

	return state;
}

void estimateHittingProbability(int N, Database* db, int target) {
	//perform brownian dynamics with file specified interactions to estimate
	//previously determined hitting probabilities

	//set parameters
	int rho = RANGE;
	double beta = BETA;
	int method = 1;
	int num_states = db->getNumStates();
	int num_trajectories = 50;
	int T_cut = 20;
	double DT = 0.01;

	//initialize interaction matrices
	int* types = new int[N];
	int numTypes = bd::readDesignFile(N, types);
	int numInteractions = numTypes*(numTypes+1)/2;
	double* kappa = new double[numInteractions];
	readKappaFile(numInteractions, kappa);
	std::map<std::pair<int,int>, double> kmap; 
	makeKappaMap(numTypes, kappa, kmap);
	int* P = new int[N*N];
	double* E = new double[N*N];
	bd::fillP(N, types, P, E, kmap);

	//set the initial and final state storage
	double* X0 = new double[DIMENSION*N];

	//set potential type
	int pot = POTENTIAL; //0 morse, 1 lj

	//get ground states
	std::vector<int> ground; //vector to hold all ground states
	double* Tconst = new double[num_states*num_states];
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//construct a map from state to number of hits
	std::map<int,double> hitSamples;
	for (int i = 0; i < ground.size(); i++) {
		hitSamples.insert(std::pair<int,double>(ground[i],0.0));
	}

	//vector of hitting times
	std::vector<double> avgT;

	//loop over some number of trajectories to get first hit probs
	int traj = 0;
	while (traj < num_trajectories) {
		setupChain(X0,N);

		double T = 0; 
		while (T < T_cut) {
			//solve the sde, determine new state
			solveSDE(X0, N, DT, rho, beta, E, P, method, pot);
			T += DT;
			int new_state = findMatrix(N, X0, db);

			//check if new state is a ground state, if yes, increment and break
			std::map<int, double>::iterator it = hitSamples.find(new_state); 
			if (it != hitSamples.end()) {
			  it->second ++;
			  traj++;
			  printf("Finished trajectory %d. Hit took %f seconds\n", traj, T);
			  printCluster(X0, N);
			  std::vector<int>::iterator itV = std::find(targets.begin(), targets.end(), new_state);
				if (itV != targets.end()) {
					avgT.push_back(T);
				}
			  break;
			}
		}

		if (T > T_cut) {
			std::cout << "Did not find a ground state in the allotted time\n";
		}
	}

	//normalize the empirical probabilities by num_trajectories
	//get total probability for the target
	double targetProb = 0;
	for (auto it = hitSamples.begin(); it != hitSamples.end(); it++) {
		it->second /= num_trajectories;
		double prob = it->second;
		int state = it->first;
		printf("Ground State %d, Probability %f\n", state, prob);

		std::vector<int>::iterator itV = std::find(targets.begin(), targets.end(), state);
		if (itV != targets.end()) {
			targetProb += prob;
		}
	}

	printf("Total Target Probability: %f\n", targetProb);

	//get average time for transitions
	double S = 0; 
	for (int i = 0; i < avgT.size(); i++) {
		S += avgT[i];
	}
	S /= avgT.size();

	printf("Average transition time: %f\n", S);

	//free memory
	delete []types; delete []kappa; delete []P; delete []E; delete []X0;
	delete []Tconst;

}

/****************************************************/
/*********** Extraneous Testing stuff ***************/
/****************************************************/

void printHitProbability(int num_states, int initial, std::vector<int> targets, double* U) {
	//print hitting probability among one group of ground states

	double prob = 0; 
	std::vector<double> probs;
	for (int i = 0; i < targets.size(); i++) {
		double p = U[toIndex(initial, targets[i], num_states)];
		prob += p;
		probs.push_back(p);
	}

	for (int i = 0; i < targets.size(); i++) {
		printf("State %d, hit prob %f\n", targets[i], probs[i]/prob);
	}
}

void solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
								 				int numTypes, double* Tconst, std::vector<int> targets) {
	//set up and solve the absorbtion probability problem - only to targets

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* P = new double[num_states*num_states]; //probability transition matrix
	double* U = new double[num_states*num_states]; //hitting probability matrix
	double* eq = new double[num_states];           //equilibrium measure

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the hitting probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);
	//reset transition matrices
	for (int i = 0; i < num_states*num_states; i++) {
		P[i] = U[i] = 0;
	}
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);
	//use the filled rate matrix to compute probability matrix
	createProbabilityMatrix(T, num_states, P);
	//solve for hitting probabilities to endStates states
	computeHittingProbability(P, num_states, targets, U);
	//get hit prob to targets
	printHitProbability(num_states, initial, targets, U);

	delete []T; delete []P; delete []U; delete []eq;
}

void evalStats(int N, Database* db, int initial, int target, bool useFile) {
	//evaluate the eq prob, hitting prob, rate, and configurations for the
	//properties in the input files

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	std::cout << "All ground states:\n";
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}
	std::cout << "\n";

	//find all target states consistent with input target
	std::vector<int> targets; //targets.push_back(target);
	// /*
	findIsomorphic(N, num_states, target, db, targets);
	std::cout << "All permutations of target state:\n";
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	std::cout << "\n";
	// */

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);

	//init a gradient array
	double* gE = new double[numInteractions];
	double* gR = new double[numInteractions];
	
	//eval eq, hit prob, rate
	double hit = solveAbsorbProblem(initial, kappaVals, db, particleTypes, numTypes, Tconst, 
																	 ground, targets);
	std::cout << "hello\n";
	double eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
												targets, gE);
	std::cout << "hello\n";
	double rate = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												Tconst, targets, gR);
	std::cout << "hello\n";

	double dydx = gR[0]/gE[0] + gR[1]/gE[1] + gR[2]/gE[2];

	//print info about the current run
	printf("Kappa Values: %f, %f, %f\n", kappaVals[0], kappaVals[1], kappaVals[2]);
	printf("Total Equilibrium Probability: %f\n", eq);
	printf("Total Hitting Probability: %f\n", hit);
	printf("Total Average Folding Rate: %f\n", rate);
	printf("Equilibrium Gradient: %f, %f, %f\n", gE[0], gE[1], gE[2]);
	printf("Rate Gradient: %f, %f, %f\n", gR[0], gR[1], gR[2]);
	printf("Rate vs. p derivative: %f\n", dydx);
	std::cout << "\n";

	//print the eq prob of each permutation individually
	std::cout << "Equilibrium Probabilities:\n";
	printEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
	std::cout << "\n";

	//print the hitting prob of each permutation individually
	std::cout << "Hitting Probabilities:\n";
	solveAbsorbProblem(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
	std::cout << "\n";

	//print the bond types
	std::cout << "Bond distributions:\n";
	getBondTypes(N, particleTypes, db, targets);
	std::cout << "\n";

	//free memory
	delete []particleTypes; delete []Tconst; delete []kappaVals;
	delete []gE; delete []gR;


}

void computeParetoFront(int N, Database* db, int initial, int target, bool useFile) {
	/*use gradient of eq and rate map at rate maxima to determine a pareto front
	  using a discrete approximation.
	  Must provide the kappa value of the rate max in the input file */

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);

	//set the grid spacing (in pi space)
	double h = 0.003;

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqPareto.txt");

	//init a gradient array
	double* gE = new double[numInteractions];
	double* gR = new double[numInteractions];

	//move along the discrete pareto front?
	int M = 250;
	for (int i = 0; i < M; i++) {
		double eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												targets, gE);
		double rate = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
													Tconst, targets, gR);

		printf("rate grad: %f, %f, %f\n", gR[0], gR[1], gR[2]);
		printf("eq grad: %f, %f, %f\n", gE[0], gE[1], gE[2]);


		double dydx = 0;
		for (int i = 0; i < 3; i++) {
			if (kappaVals[i] > 0.3 && kappaVals[i] < 100) {
				dydx += gR[i]/gE[i];
			}
		}
		std::cout << dydx << "\n";
		double new_r = rate + dydx*h;
		double new_e = eq + h;
		printf("r = %f, e = %f \n", rate, eq);

		if (dydx > 0.1 || abs(dydx)>1) {
			dydx = 0;
			gR[0] = 0; gR[1] = 0; gR[2] = 0;
		}

		Eigen::MatrixXd Q(2,3); Q.fill(0.0);
		Eigen::VectorXd b(2);   b.fill(0.0);
		Eigen::VectorXd dk(3);  dk.fill(0.0);

		Q(0,0) = gR[0]; Q(0,1) = gR[1]; Q(0,2) = gR[2];
		Q(1,0) = gE[0]; Q(1,1) = gE[1]; Q(1,2) = gE[2];
		b(0) = dydx*h; b(1) = h;

		dk = Q.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

		kappaVals[0] += dk(0); kappaVals[1] += dk(1); kappaVals[2] += dk(2);
		checkPositive(numInteractions, kappaVals);
		printf("Kappa Values: %f, %f, %f\n", kappaVals[0], kappaVals[1], kappaVals[2]);
		ofile << eq << ' ' << rate << "\n";

	}



	//free memory
	delete []particleTypes; delete []Tconst; delete []kappaVals;
	delete []gE; delete []gR;

	ofile.close();
}

double rateMaxPareto(double& c, int initial, Database* db, int num_states, int* particleTypes, 
										 int numTypes, int numInteractions, double* kappaVals, double* Tconst, 
										 std::vector<int> ground, std::vector<int> targets) {
	//use optimization scheme to get max rate for target in toy model
	//has equality constraint for equilibrium prob
	//uses the file for the permutation

	//set iteration settings
	int max_its = 20; double objTol = 1e-8; 
	double r = 200;    //cost 

	//init a gradient array
	double* gH = new double[numInteractions];
	double* gR = new double[numInteractions];

	//set free variables
	double eq;
	double R; 

	//optimization - steepest ascent w/ line search
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of hitting prob and rate
		std::cout << it << "\n";
		eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes,
												targets, gH);
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
												Tconst, targets, gR);

		//define the objevtive function
		double obj = R - r * (eq - c) * (eq-c);  //maybe use quadratic penalty

		//get the gradient of obj
		for (int i = 0; i < numInteractions; i++) {
			gR[i] -= 2 * r * (eq - c) * gH[i]; 
		}

		//find a step size to give decrease in obj fn
		double step = 0; 
		double delta = lineSearchRate(initial, numInteractions, kappaVals, db, particleTypes, 
															numTypes, Tconst, targets, c, r, eq, R, gR, step);

		//take the step
		for (int i = 0; i < numInteractions; i++) {
			kappaVals[i] += step * gR[i];
		}
		checkPositive(numInteractions, kappaVals);


		//check for objective function tolerance 
		if (delta < objTol) {
			break;
		}
	}


	//free memory
	delete []gH; delete []gR; 

	c = eq;

	//return the rate
	return R;
}

void computeParetoFrontGD(int N, Database* db, int initial, int target, bool useFile) {
	/* try to compute a pareto front using gradient descent with constraints to 
		 optimize to the max rate at a given value of eq prob */

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);

	//set the grid spacing (in pi space)
	double h = 0.01;

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqPareto.txt");

	//get the initial rate and eq, output it
	double r = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
	double eq = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
	//ofile << eq << ' ' << r << ' ' << kappaVals[0] << ' ' << kappaVals[1] << 
	//	' ' << kappaVals[2] <<  "\n";
	eq = 0.01; 
	double eq_grid = eq;

	//increment eq, maximize r
	while (eq_grid < 0.3) {
		eq_grid += h;
		eq = eq_grid;
	//while (eq > 0.35) {
		//eq -= h;

		//find the new r
		r = rateMaxPareto(eq, initial, db, num_states, particleTypes, numTypes, numInteractions, kappaVals, 
											Tconst, ground, targets);

		//output to console for progress check and file for analysis
		std::cout << eq << ' ' << r << "\n";
		ofile << eq << ' ' << r << ' ' << kappaVals[0] << ' ' << kappaVals[1] << 
		' ' << kappaVals[2] <<  "\n";
	}


	//free memory
	delete []particleTypes; delete []Tconst; delete []kappaVals;
	ofile.close();
}


void miskinOptimize(int N, Database* db, int initial, int target, bool useFile) {
	//uses the lumped miskin scheme to optimize parameters to form target

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	std::cout << "All ground states:\n";
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}
	std::cout << "\n";

	//find all target states consistent with input target
	std::vector<int> targets; //targets.push_back(target);
	// /*
	findIsomorphic(N, num_states, target, db, targets);
	std::cout << "All permutations of target state:\n";
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	std::cout << "\n";
	// */

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);

	//init a gradient array
	double* gE = new double[numInteractions];
	double* gR = new double[numInteractions];
	double* eqA = new double[num_states];

	//tests the analytical derivative function
	/*
	double eq1 = analyticEqGradient(db, particleTypes, numTypes, numInteractions, kappaVals, eqA, gE, targets);
	printf("Val %f, Gradients: %f, %f, %f\n", eq1, gE[0], gE[1], gE[2]);
	double eq2 = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
												targets, gE);
	printf("Val %f, Gradients: %f, %f, %f\n", eq2, gE[0], gE[1], gE[2]);
	abort();
	*/
	
	
	
	//solve the ode to give the parameters 
	double dt = 0.01;
	int max_its = 5000;
	double m = 0.15;
	double M = 1e6;
	double input_tol = 1e-9;
	double tol = 1e-7;
	double prev = 0;
	int jump = 0; int max_jump = 10;
	double maxEq = 0;
	double limit_tol = 1e-4;

	for (int i = 0; i < max_its; i++) {
		//evaluate the equilibrium probability of state target, as well as the gradient wrt kappa
		//double eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
												//targets, gE);
		double eq = analyticEqGradient(db, particleTypes, numTypes, numInteractions, kappaVals, eqA, gE, targets);
		double R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, numTypes, 
												Tconst, targets, gR);

		//double prefactor = 1.0-eq*eq;
		double prefactor = 1.0-eq*eq;

		printf("Iteration: %d, Eq: %f, k1 %f, k2 %f, k3 %f\n", i, eq, kappaVals[0], kappaVals[1], kappaVals[2]);

		for (int j = 0; j < 3; j++) {
			if (true) {
				if (i % 2 == 0) {
					if (abs(gE[j]) > tol) {
						double update = dt * prefactor * 1.0 / gE[j];
						kappaVals[j] = kappaVals[j] + update;
						//printf("Update %d, val %f   ", j, update);
					}
				}
				else if (i % 2 == 1) {
					//do a rate update
					//kappaVals[j] = kappaVals[j] + 0.8 * gR[j];

					//do another eq update

					if (abs(gE[j]) > tol) {
						double update = dt * prefactor * 1.0 / gE[j];
						kappaVals[j] = kappaVals[j] + update;
						//printf("Update %d, val %f   ", j, update);
					}
				}
			}
			if (kappaVals[j] < m) {
				kappaVals[j] = m;
			}
			if (kappaVals[j] > M) {
				kappaVals[j] = M;
			}
		}

		
		//printf("\n Gradients: %f, %f, %f\n", gE[0], gE[1], gE[2]);

		//check if eq is sufficiently close to 1
		double diff1 = 1-eq;
		if (diff1 < 0.01) {
			printf("Equilibrium Probability is within tolerance of 1\n");
			break;
		}

		//check if eq is sufficiently close to previous iteration
		double input_diff = abs(eq-prev);
		if (input_diff < input_tol) {
			printf("Equilibrium probability is same as previous iteration, up to tol\n");
			break;
		}

		//check if there was a significant downward jump in eq
		if ((eq - prev) < -0.1) {
			if (jump > 1 && abs(maxEq-prev) < limit_tol) {
				printf("Limit cycle detected. Terminating iteration.\n");
				break;
			}

			jump++; double jump_val = prev;
			printf("Jump %d detected. Max value is %f \n", jump, jump_val);
			if (jump_val > maxEq) {
				maxEq = jump_val;
			}

			if (jump == max_jump) {
				printf("Too many jumps. Terminating iteration.\n");
				break;
			}

		}
		prev = eq;
	}


	double r = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);
	double eq = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
	printf("This configuration has eq = %f, and rate = %f\n", eq, r);

	


	//free memory
	delete []particleTypes; delete []Tconst; delete []kappaVals;
	delete []gE; delete []gR;

}

double analyticEqGradient(Database* db, int* particleTypes, int numTypes, int numInteractions, 
												double* kappaVals, double* eq, double* gE, std::vector<int> targets) {
	//evaluate the gradient of the equilibrium probability analytically

	//get the system parameters
	int N = db->getN();
	int ns = db->getNumStates();
	double tol = 1e-6;

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;
	makeKappaMap(numTypes, kappaVals, kappa);
	
	//declare the index mapping
	std::map<std::pair<int,int>,int> index;
	int* indexVals = new int[numInteractions];
	for (int i = 0; i < numInteractions; i++) {
		indexVals[i] = i;
	}
	makeIndexMap(numTypes, indexVals, index);

	//create an array for the derivative of each component
	double* eqDeriv = new double[ns];

	//first, reweight the eq distribution
	reweight(N, ns, db, particleTypes, eq, kappa);

	//loop over all bond types
	for (int bond = 0; bond < numInteractions; bond++) {
		double S = 0; gE[bond] = 0;

		//sum the eq prob over all states, weighted by n_b
		for (int i = 0; i < ns; i++) {
			int b = searchForBonds(N, db, i, bond, particleTypes, index);
			S += eq[i] * float(b);
		}

		//fill in the derivative array
		for (int i = 0; i < ns; i++) {
			int b = searchForBonds(N, db, i, bond, particleTypes, index);
			eqDeriv[i] = eq[i] / kappaVals[bond] * (float(b) - S);
		}

		//sum eqDeriv over the target states
		for (int i = 0; i < targets.size(); i++) {
			gE[bond] += eqDeriv[targets[i]];
		}

		//this checks if we accidentally overshoot a k->inf critical point and go negative
		if (abs(gE[bond]) < tol && gE[bond] < 0) {
			//printf("Got here, bond %d\n", bond);
			//gE[bond] *= -1;
			gE[bond] = 0;
		}
	}

	//compute the eq prob too
	double eqProb = 0;
	for (int i = 0; i < targets.size(); i++) {
		eqProb += eq[targets[i]];
	}


	//free memory
	delete []eqDeriv; delete []indexVals;

	return eqProb;

}

void makeIndexMap(int numTypes, int* indexVals, 
									std::map<std::pair<int,int>,int>& index) {
	//fill the map with values in kappaVals

	int counter = 0; 

	for (int i = 0; i < numTypes; i++) {
		for (int j = i; j < numTypes; j++) {
			index[std::make_pair(i,j)] = indexVals[counter];
			index[std::make_pair(j,i)] = indexVals[counter];
			counter++;
		}
	}
}

int searchForBonds(int N, Database* db, int state, int type, int* particleTypes,
									 std::map<std::pair<int,int>,int>& index) {
	//determine the number of bonds of a given type for given state

	int bonds = 0;

	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( (*db)[state].isInteracting(i,j,N)) {
				int p1 = particleTypes[i];
				int p2 = particleTypes[j];

				int bondIndex = index[std::make_pair(p1,p2)];
				if (bondIndex == type) {
					bonds++;
				}
			}
		}
	}

	//int b_actual = (*db)[state].getBonds();
	//printf("Found %d bonds, actual number is %d\n", bonds, b_actual-5);


	return bonds;
}

/***********************************************************************/
/********* Plot time dependent trajectory in phase space ***************/
/***********************************************************************/

void phaseTrajectory(int N, Database* db, int initial, int target, bool useFile, std::string& p_file) {
	//determine the points in phase space (prob,rate) along the time dependent 
	//optimized parameter protocol

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqTime.txt");




	//get the parameter sequence from the supplied text file
	Eigen::VectorXd beta; beta.setZero(100);
	Eigen::MatrixXd    E; E.setZero(100,3);
	getNEQparameters(p_file, beta, E);

	//compute the dimensionaless product beta*E
	for (int i = 0; i < 3; i++) {
		E.col(i) = E.col(i).cwiseProduct(beta);
	}

	//convert the energies to sticky parameters
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 100; j++) {
			double energy = E(j,i);
			double k = exp(energy) / sqrt(2.0*30*30*energy);
			E(j,i) = k;
		}
	}



	//do eq and rate calculations
	for (int timestep = 0; timestep < 100; timestep++) {
		
		//set kappa
		kappaVals[0] = E(timestep,0);
		kappaVals[1] = E(timestep,1);
		kappaVals[2] = E(timestep,2);

		//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";

		//make the map
		makeKappaMap(numTypes, kappaVals, kappa);
		//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
		//do rewieght
		reweight(N, num_states, db, particleTypes, eq, kappa);
		//get eq prob
		double eqProb = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);
		//copy Tconst into T
		std::copy(Tconst, Tconst+num_states*num_states, T);
		//fill in transposed entries such that T satisfies detailed balance
		satisfyDB(T, num_states, db, eq);
		//fill in diagonal with negative sum of row entries
		fillDiag(T, num_states);
		//get the transition rate
		if (num_states > 400) { //sparse solve
			computeMFPTsSP(num_states, T, targets, m);
		}
		else { //dense solve
			computeMFPTs(num_states, T, targets, m);
		}
		
		double rate = 1/m[initial];

		//output to file
		ofile << eqProb << ' ' << rate <<  "\n";

	}

	//close file
	ofile.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
	delete []T; delete []Tconst; delete []m; 




}

void getNEQparameters(std::string& p_file, Eigen::VectorXd& beta, Eigen::MatrixXd& E) {
	//extract the data from the file into the eigen arrays

	//open the file
	std::ifstream in_str(p_file);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", p_file.c_str());
		return;
	}

	//store the entries
	double b; double ee;
	int timestep = 0;
	while (in_str >> b) {
		beta(timestep) = b;
		for (int i = 0; i < 3; i++) {
			in_str >> ee;
			E(timestep, i) = ee;
		}
		timestep++;
	}

}



/***********************************************************************/
/********* Test possible sampling measures on the MC model ***************/
/***********************************************************************/






void createRateMatrix(int N, Database* db, int numTypes, double* kappaVals,
											int* particleTypes, 
											Eigen::MatrixXd& R, const Eigen::MatrixXd& r_form) {
	//use the database to construct the rate matrix at the given values
	//no derivative, so this is just a matrix

	//declare equilibrium measure and perform reweight
	int num_states = db->getNumStates();
	double* eq = new double[num_states];           //equilibrium measure
	std::map<std::pair<int,int>,double> kappa;    //map for interaction type to kappa
	makeKappaMap(numTypes, kappaVals, kappa);
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//fill in forward and satisfy detailed balance wrt eq for backward
	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (r_form(i,j) > 0) {
				R(i,j) = r_form(i,j);
				R(j,i) = r_form(i,j) * (eq[i] / eq[j]);
			}
		}
	}

	//finally, must set the diagonal to be negative of the row sum
	for (int i = 0; i < num_states; i++) {
		double R_sum = 0; 

		for (int j = 0; j < num_states; j++) {
			R_sum += R(i,j);
		}

		R(i,i) = -R_sum; 
	}


	delete []eq;
}



double smallestBarrier(int N, Database* db, int target, double* E) {
	//determine the weakest bond in the target state under the given configuration

	double minE = 40;
	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( (*db)[target].isInteracting(i,j,N)) {
		    double energy = E[toIndex(i,j,N)];
		    if (energy < minE) {
		    	minE = energy;
		    	//printf("minE is %f\n", minE);
		    }
			}
		}
	}

	return minE;
}

double harmonicMeanBarrier(int N, Database* db, int target, double* E) {
	//determine harmonic mean of the bond strengths in target state

	double m = 0;

	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( (*db)[target].isInteracting(i,j,N)) {
		    double energy = E[toIndex(i,j,N)];
		    m += 1.0 / energy;
			}
		}
	}

	return 1.0 / m;
}

void buildAM(int N, Database* db, int state, int* M) {
	//construct adjacency matrix from db

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			M[bd::toIndex(i,j,N)] = (*db)[state].isInteracting(i,j,N);
		}
	}
}

void misfoldVector(int N, Database* db, int target, Eigen::VectorXd& misfold, double* E) {
	//construct the vector of misfolded energies

	//get number of states
	int num_states = db->getNumStates();

	//make adjacency matrix of target state and comp state
	int* M_target = new int[N*N]; for (int i = 0; i < N*N; i++) M_target[i] = 0;
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;
	buildAM(N, db, target, M_target);

	//loop over states
	for (int i = 0; i < num_states; i++) {
		double badEnergy = 0;
		buildAM(N, db, i, M);
		for (int p1 = 0; p1 < N; p1++) {
			for (int p2 = p1+2; p2 < N; p2++) {
				int inState = M[bd::toIndex(p1,p2,N)];
				int inTarget = M_target[bd::toIndex(p1,p2,N)];
				if (inState == 1 && inTarget == 0) {
					badEnergy += E[toIndex(p1,p2,N)];
				}
			}
		}
		misfold(i) = badEnergy;
	}

	delete []M; delete []M_target;

}

double averageMisfoldEnergy(int N, double T, int initial, double* kappaVals, Database* db, 
														int* particleTypes, int numTypes, double* Tconst,
														const Eigen::VectorXd& misfold, std::vector<int>& targets) {
	//get average misfolded energy from time 0 to T

	int num_states = db->getNumStates();

	//first we need to construct the generator
	Eigen::MatrixXd R; R.setZero(num_states, num_states);
	Eigen::MatrixXd r_form; r_form.setZero(num_states, num_states);
	for (int i = 0; i < num_states*num_states; i++) {
		r_form(i) = Tconst[i];
	}
	createRateMatrix(N, db, numTypes, kappaVals, particleTypes,  R, r_form);

	//now we need a matrix exponential of R
	Eigen::MatrixXd L; L.setZero(num_states,num_states);
	R = R*T; 
	L = R.exp() - Eigen::MatrixXd::Identity(num_states,num_states);
	R = R / T;
	Eigen::RowVectorXd r; r.setZero(num_states);
	r = L.row(initial);

	//need to compute r*R^-1. fR = r. 
	Eigen::RowVectorXd f; f.setZero(num_states);
	f = r * R.completeOrthogonalDecomposition().pseudoInverse();
	f = (f+f.cwiseAbs())/2.0;
	//std::cout << f << "\n";
	//std::cout << (f*R - r).norm() << "\n";

	//abort();


	//take the row corresponding to initial of product LF, divide by T
	//double prod = f.sum();
	double prod = f.dot(misfold);

	return 1.0 / T * prod;

}

double averageMisfoldEnergyEnd(int N, double T, int initial, double* kappaVals, Database* db, 
														int* particleTypes, int numTypes, double* Tconst,
														const Eigen::VectorXd& misfold, std::vector<int>& targets) {
	//get average misfolded energy from time 0 to T

	int num_states = db->getNumStates();

	//first we need to construct the generator
	Eigen::MatrixXd R; R.setZero(num_states, num_states);
	Eigen::MatrixXd r_form; r_form.setZero(num_states, num_states);
	for (int i = 0; i < num_states*num_states; i++) {
		r_form(i) = Tconst[i];
	}
	createRateMatrix(N, db, numTypes, kappaVals, particleTypes,  R, r_form);

	//now we need a matrix exponential of R
	Eigen::MatrixXd L; L.setZero(num_states,num_states);
	R = R*T; 
	L = R.exp();
	Eigen::RowVectorXd r; r.setZero(num_states);
	r = L.row(initial);

	/*
	for (int i = 0; i < num_states; i++) {
		printf("State %d, prob %f, leave rate %f\n", i, r(i), R(i,i));
	}
	*/

	double prod = r.dot(misfold);

	return prod;

}

double averageMisfoldEnergyTrap(int N, double T, int initial, double* kappaVals, Database* db, 
														int* particleTypes, int numTypes, double* Tconst,
														const Eigen::VectorXd& misfold, std::vector<int>& targets) {
	//get average misfolded energy at first trap time

	int num_states = db->getNumStates();

	//first we need to construct the generator
	Eigen::MatrixXd R; R.setZero(num_states, num_states);
	Eigen::MatrixXd r_form; r_form.setZero(num_states, num_states);
	for (int i = 0; i < num_states*num_states; i++) {
		r_form(i) = Tconst[i];
	}
	createRateMatrix(N, db, numTypes, kappaVals, particleTypes,  R, r_form);

	//next we determine the "trap states" via a threshold on diagonal entries
	double threshold = 2; double min_rate = 1e6;
	//get minimum rate
	for (int i = 2; i < num_states; i++) {
		if (fabs(R(i,i)) < min_rate) {
			min_rate = fabs(R(i,i));
		}
	}
	//determine all  states within threshold
	std::vector<int> traps;
	for (int i = 2; i < num_states; i++) {
		if (fabs(R(i,i)) < threshold*min_rate) {
			traps.push_back(i);
			printf("State %d, rate %f\n", i, R(i,i));
		}
	}

	//determine the first hitting probability for the trap states
	Eigen::MatrixXd U; U.setZero(num_states,num_states);
	Eigen::MatrixXd P; P.setZero(num_states,num_states);
	//make rate matrix into probability matrix
	for (int i = 0; i < num_states; i++) {
		double Z = -R(i,i);
		for (int j = 0; j < num_states; j++) {
			P(i,j) = R(i,j) / Z;
		}
		P(i,i) = 0.0;
	}
	computeHittingProbability(P, num_states, traps, U);
	Eigen::VectorXd probs = U.row(initial); 
	
	double prod = probs.dot(misfold);

	return prod;

}

double finalTargetProb(int N, double T, int initial, double* kappaVals, Database* db, 
														int* particleTypes, int numTypes, double* Tconst,
														std::vector<int>& targets) {
	//return probability of target states at final time

	int num_states = db->getNumStates();

	//first we need to construct the generator
	Eigen::MatrixXd R; R.setZero(num_states, num_states);
	Eigen::MatrixXd r_form; r_form.setZero(num_states, num_states);
	for (int i = 0; i < num_states*num_states; i++) {
		r_form(i) = Tconst[i];
	}
	createRateMatrix(N, db, numTypes, kappaVals, particleTypes,  R, r_form);

	//now we need a matrix exponential of R
	Eigen::MatrixXd L; L.setZero(num_states,num_states);
	R = R*T; 
	L = R.exp();
	Eigen::RowVectorXd r; r.setZero(num_states);
	r = L.row(initial);

	double prob = r(targets[0]);

	return prob;

}

double stayingProb(int N, double T, int initial, double* kappaVals, Database* db, 
														int* particleTypes, int numTypes, double* Tconst,
														std::vector<int>& targets) {
	//return probability of target states at final time when initialized there

	int num_states = db->getNumStates();

	//first we need to construct the generator
	Eigen::MatrixXd R; R.setZero(num_states, num_states);
	Eigen::MatrixXd r_form; r_form.setZero(num_states, num_states);
	for (int i = 0; i < num_states*num_states; i++) {
		r_form(i) = Tconst[i];
	}
	createRateMatrix(N, db, numTypes, kappaVals, particleTypes,  R, r_form);

	//now we need a matrix exponential of R
	Eigen::MatrixXd L; L.setZero(num_states,num_states);
	R = R*T; 
	L = R.exp();
	Eigen::RowVectorXd r; r.setZero(num_states);
	r = L.row(targets[0]);

	double prob = r(targets[0]);

	return prob;

}


void testMeasures(int N, Database* db, int initial, int target, bool useFile) {
	//test measures - output eq prob and mfpt, along with test meaasures

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	readKappaFile(numInteractions, kappaVals);

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	std::cout << "All ground states:\n";
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}
	std::cout << "\n";

	//find all target states consistent with input target
	std::vector<int> targets; targets.push_back(target);

	//set interaction matrices
	int* P = new int[N*N]; 
	double* E = new double[N*N]; 

	std::map<std::pair<int,int>,double> kappaMap;
	makeKappaMap(numTypes, kappaVals, kappaMap);
	fillP(N, particleTypes, P, E, kappaMap);


	Eigen::VectorXd misfold;
	double T = 500.0;



	//setup kappa values
	for (int i = 0; i < N*N; i++) P[0] = 0; for (int i = 0; i < N*N; i++) E[0] = 0;
	kappaMap.clear();
	makeKappaMap(numTypes, kappaVals, kappaMap);
	fillP(N, particleTypes, P, E, kappaMap);

	//get eq prob
	double eq = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);

	//get mfpt 
	double rate = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);

	//get measure 1
	/*
	double c1 = .25;
	double m1 = 1.0-exp(-c1*smallestBarrier(N, db, target, E));
	*/
	double m1 = stayingProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	 Tconst, targets);

	//get measure 2
	misfold.setZero(num_states);
	misfoldVector(N, db, target, misfold, E);
	//std::cout << misfold << "\n";
	double m2 = averageMisfoldEnergyTrap(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	 Tconst, misfold, targets);
	
	std::cout << m2 << "\n";
	double c2 = 0.5;
	m2 = exp(-c2*m2);

	/*
	double m2 = finalTargetProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	 Tconst, targets);
	*/

	


	printf("X axis: eq = %f, m1 = %f\n", eq, m1);
	printf("Y axis: mfpt = %f, m2 = %f\n", rate, m2);




	//free memory
	delete []P; delete []E; delete []Tconst; 
	delete []kappaVals; delete []particleTypes;
}


void testMeasuresRecord(int N, Database* db, int initial, int target, bool useFile) {
	//test measures - output eq prob and mfpt, along with test meaasures

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	readKappaFile(numInteractions, kappaVals);

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	std::cout << "All ground states:\n";
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}
	std::cout << "\n";

	//find all target states consistent with input target
	std::vector<int> targets; targets.push_back(target);

	//set interaction matrices
	int* P = new int[N*N]; 
	double* E = new double[N*N]; 

	std::map<std::pair<int,int>,double> kappaMap;
	makeKappaMap(numTypes, kappaVals, kappaMap);
	fillP(N, particleTypes, P, E, kappaMap);


	Eigen::VectorXd misfold;
	double T = 5.0;
	double M = 80;
	double mult = 1.05;

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoCorrelations2.txt");

	//loop over kappa values for triangle
	for (int i = 0; i < M; i++) {

		//setup kappa values
		for (int i = 0; i < N*N; i++) P[0] = 0; for (int i = 0; i < N*N; i++) E[0] = 0;
		kappaMap.clear();
		kappaVals[1] *= mult;
		makeKappaMap(numTypes, kappaVals, kappaMap);
		fillP(N, particleTypes, P, E, kappaMap);

		//get eq prob
		double eq = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);

		//get mfpt 
		double rate = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);

		//get measure 1
		/*
		double c1 = .25;
		double m1 = 1.0-exp(-c1*smallestBarrier(N, db, target, E));
		*/

		
		double m1 = stayingProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	 Tconst, targets);
		

		/*
		double m1 = 1.0 - exp(-harmonicMeanBarrier(N, db, target, E));
		*/

		//get measure 2
		
		misfold.setZero(num_states);
		misfoldVector(N, db, target, misfold, E);
		//std::cout << misfold << "\n";
		double m2 = averageMisfoldEnergyTrap(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																		 Tconst, misfold, targets);
		double c2 = 0.5;
		m2 = exp(-c2*m2);
		

		/*
		double m2 = finalTargetProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	 Tconst, targets);
		*/

		ofile << eq << ' ' << rate << ' ' << m1 << ' ' << m2 << "\n";
	}


	//printf("X axis: eq = %f, m1 = %f\n", eq, m1);
	//printf("Y axis: mfpt = %f, m2 = %f\n", rate, m2);

	ofile.close();



	//free memory
	delete []P; delete []E; delete []Tconst; 
	delete []kappaVals; delete []particleTypes;
}

void testMeasuresScatter(int N, Database* db, int initial, int target, bool useFile) {
	//test measures - output eq prob and mfpt, along with test meaasures

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	readKappaFile(numInteractions, kappaVals);

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	std::cout << "All ground states:\n";
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}
	std::cout << "\n";

	//find all target states consistent with input target
	std::vector<int> targets; targets.push_back(target);

	//set interaction matrices
	int* P = new int[N*N]; 
	double* E = new double[N*N]; 

	std::map<std::pair<int,int>,double> kappaMap;
	makeKappaMap(numTypes, kappaVals, kappaMap);
	fillP(N, particleTypes, P, E, kappaMap);


	Eigen::VectorXd misfold;
	double T = 5.0;
	double M = 21;
	double mult = 1.5;
	double base = 0.3;

	kappaVals[0] = kappaVals[1] = kappaVals[2] = base;

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoCorrelationsScatter.txt");

	//loop over kappa values for triangle
	for (int i = 0; i < M; i++) {
		kappaVals[0] *= mult;
		printf("Outer loop number %d\n", i);
		for (int j = 0; j < M; j++) {
			kappaVals[1] *= mult;
			for (int k = 0; k < M; k++) {

				//setup kappa values
				for (int i = 0; i < N*N; i++) P[0] = 0; for (int i = 0; i < N*N; i++) E[0] = 0;
				kappaMap.clear();
				kappaVals[2] *= mult;
				makeKappaMap(numTypes, kappaVals, kappaMap);
				fillP(N, particleTypes, P, E, kappaMap);

				//get eq prob
				double eq = getEqProb(initial, kappaVals, db, particleTypes, numTypes, targets);

				//get mfpt 
				double rate = getRate(initial, kappaVals, db, particleTypes, numTypes, Tconst, targets);

				//get measure 1
				/*
				double c1 = .25;
				double m1 = 1.0-exp(-c1*smallestBarrier(N, db, target, E));
				*/
				double m1 = stayingProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																			 Tconst, targets);

				//get measure 2
				
				misfold.setZero(num_states);
				misfoldVector(N, db, target, misfold, E);
				//std::cout << misfold << "\n";
				double m2 = averageMisfoldEnergyEnd(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																				 Tconst, misfold, targets);
				double c2 = 1;
				m2 = exp(-c2*m2);
				

				/*
				double m2 = finalTargetProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																			 Tconst, targets);
				*/

				ofile << eq << ' ' << rate << ' ' << m1 << ' ' << m2 << "\n";
			}
			kappaVals[2] = base;
		}
		kappaVals[1] = base;
	}


	//printf("X axis: eq = %f, m1 = %f\n", eq, m1);
	//printf("Y axis: mfpt = %f, m2 = %f\n", rate, m2);

	ofile.close();



	//free memory
	delete []P; delete []E; delete []Tconst; 
	delete []kappaVals; delete []particleTypes;
}



}