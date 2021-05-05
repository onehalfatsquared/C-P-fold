#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include "point.h"
#include "database.h"
#include "bDynamics.h"
#include "tpt.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "../defines.h"
namespace bd{

void computeMFPTsSP(int num_states, double* T, std::vector<int> targets, double* m) {
	//compute the mean first passage time of every state to targets
	//solves A*tau = -1, where is A is rate matrix with rows/cols of the targets
	//zero'd out (removed)

	//uses sparse solver

	//Eigen::setNbThreads(0);
	//omp_set_num_threads(0);

	//printf("Solving with %d threads\n", Eigen::nbThreads());

	double tol = 1e-7;

	//find number of non-zero'd entries
	int num_targets = targets.size();
	int M = num_states - num_targets;
	Eigen::SparseMatrix<double> QMs(M,M);
	Eigen::SparseLU<Eigen::SparseMatrix<double>>   solver;
	Eigen::VectorXd b(M); 
	Eigen::VectorXd tau(M);  

	//construct a boolean array to check if a row should be skipped
	bool* skip = new bool[num_states]; for (int i = 0; i < num_states; i++) skip[i] = false;
	for (int i = 0; i < num_targets; i++) skip[targets[i]] = true;

	//fill in transition matrix for all non-skipped rows/cols
	typedef Eigen::Triplet<double> Tr;
	std::vector<Tr> tripletList;
	tripletList.reserve(num_states*10);

	int i_index = 0;
	for (int i = 0; i < num_states; i++) {
		if (!skip[i]) {
			int j_index = 0;
			for (int j = 0; j < num_states; j++) {
				if (!skip[j]) {
					double value = T[toIndex(i,j,num_states)];
					if (abs(value) > tol) {
						tripletList.push_back(Tr(i_index,j_index,value));
					}
					j_index++;
				}
			}
			i_index++;
		}
	}
	QMs.setFromTriplets(tripletList.begin(), tripletList.end());
	//QMs.makeCompressed();

	//fill b with -1
	b.fill(-1.0);

	//solve for tau
	solver.analyzePattern(QMs);
	solver.factorize(QMs);
	//std::cout << solver.info() << "\n";

	if (solver.info() != Eigen::Success) { //could not compute sparse lu, do dense solve

		//perform a dense solve
		computeMFPTs(num_states, T, targets, m);
		delete []skip;
		return;
	}

	else {
		//do the sparse solve
		tau = solver.solve(b);

		//std::cout << (QMs * tau - b).norm() << "\n";

		//store solution in m - re-add lost zeros
		int lostIndex = 0;
		for (int i = 0; i < num_states; i++) {
			if (!skip[i]) {
				m[i] = abs(tau(lostIndex));
				lostIndex++;
			}
			else {
				m[i] = 0;
			}
		}

		//free memory
		delete []skip;
		return;
	}



}

void computeMFPTs(int num_states, double* T, std::vector<int> targets, double* m) {
	//compute the mean first passage time of every state to targets
	//solves A*tau = -1, where is A is rate matrix with rows/cols of the targets
	//zero'd out (removed)

	//Eigen::setNbThreads(0);
	//omp_set_num_threads(0);

	//printf("Solving with %d threads\n", Eigen::nbThreads());

	//find number of non-zero'd entries
	int num_targets = targets.size();
	int M = num_states - num_targets;
	Eigen::MatrixXd QM(M,M); 
	Eigen::VectorXd b(M); 
	Eigen::VectorXd tau(M);  

	//construct a boolean array to check if a row should be skipped
	bool* skip = new bool[num_states]; for (int i = 0; i < num_states; i++) skip[i] = false;
	for (int i = 0; i < num_targets; i++) skip[targets[i]] = true;

	//fill in transition matrix for all non-skipped rows/cols
	int i_index = 0;
	for (int i = 0; i < num_states; i++) {
		if (!skip[i]) {
			int j_index = 0;
			for (int j = 0; j < num_states; j++) {
				if (!skip[j]) {
					QM(i_index, j_index) = T[toIndex(i,j,num_states)];
					j_index++;
				}
			}
			i_index++;
		}
	}


	//std::cout << QM.row(543)<< "\n";
	//abort();

	//fill b with -1
	b.fill(-1.0);

	//solve for tau
	tau = QM.lu().solve(b);

	//std::cout << tau << "\n";
	//std::cout << (QM * tau - b).norm() << "\n";

	//store solution in m - re-add lost zeros
	int lostIndex = 0;
	for (int i = 0; i < num_states; i++) {
		if (!skip[i]) {
			m[i] = abs(tau(lostIndex));
			lostIndex++;
		}
		else {
			m[i] = 0;
		}
	}

	//free memory
	delete []skip;
}

double computeTransitionRateTPT(int num_states, double* q, double* T, double* eq) {
	//compute the average transition rate from A to B with committor fn

	double S = 0; 
	int i, j;

	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, i, j);
		if (i != j) {
			S += eq[i] * T[entry] * (q[j]-q[i]) * (q[j]-q[i]);
		}
	}

	S /= 2;

	return S;
}

void computeFlux(int num_states, double* q, double* T, double* eq, double* flux) {
	//compute the probability flux from generator, invariant measure, and committor

	int i, j;

	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, i, j);
		if (i != j) {
			flux[entry] = eq[i]*T[entry]*q[j]*(1-q[i]);
		}
	}
}

void computeCommittor(double* q, double* T, int num_states, int initial, std::vector<int> targets) {
	//set up and solve equation for the committor function

	//initialize the matrix and vectors
	Eigen::MatrixXd TM(num_states,num_states); 
	Eigen::VectorXd b(num_states); 
	Eigen::VectorXd x(num_states);      


	//fill the matrix and RHS appropriately
	//copy T into TM, fill b with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		TM(i) = T[i]; 
	}
	b.fill(0.0);

	//loop over initial and target states. edit corresponding rows.
	//initial
	for (int j = 0; j < num_states; j++) {
		if (j != initial) {
			TM(initial, j) = 0;
		}
		else {
			TM(initial, j) = 1;
		}
	}

	//target
	for (int ind = 0; ind < targets.size(); ind++) {
		int target = targets[ind];
		for (int j = 0; j < num_states; j++) {
			if (j != target) {
				TM(target, j) = 0;
			}
			else {
				TM(target, j) = 1;
			}
		}
		b(target) = 1;
	}

	//solve the linear system
	x = TM.lu().solve(b);

	//store solution in q
	for (int i = 0; i < x.size(); i++) {
		q[i] = abs(x(i));
	}
}

void computeHittingProbability(double* P, int num_states, std::vector<int> endStates, 
															 double* U) {
	//compute the hitting probabilities for the states in endStates

	//init and fill all the required matrices
	Eigen::MatrixXd R(num_states,num_states); Eigen::MatrixXd Q(num_states,num_states); 
	Eigen::MatrixXd Pm(num_states,num_states); Eigen::MatrixXd Um(num_states,num_states); 
	for (int i = 0; i < num_states*num_states; i++) {
		R(i) = 0; Um(i) = 0; 
		Q(i) = P[i]; Pm(i) = P[i];
	}


	//let R be matrix with only columns corresponding to end states, diagonal = 1
	//let Q be matrix with rows/cols corresponding to end state zero'd out
	for (int i = 0; i < endStates.size(); i++) {
		int state = endStates[i];
		R.col(state) = Pm.col(state); R(state,state) = 1;
		Q.row(state) = Um.row(state); Q.col(state) = Um.col(state); 
	}

	//solve for Um via (I-Q)Um = R
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(num_states,num_states) - Q;
	Um = D.lu().solve(R);

	//store back in U array
	for (int i = 0; i < num_states*num_states; i++) {
		U[i] = Um(i);
	}
}

void computeHittingProbability(Eigen::MatrixXd& P, int num_states, std::vector<int> endStates, 
															 Eigen::MatrixXd& U) {
	//compute the hitting probabilities for the states in endStates

	//init and fill all the required matrices
	Eigen::MatrixXd R(num_states,num_states); Eigen::MatrixXd Q(num_states,num_states); 
	Eigen::MatrixXd Pm(num_states,num_states); Eigen::MatrixXd Um(num_states,num_states); 
	for (int i = 0; i < num_states*num_states; i++) {
		R(i) = 0; Um(i) = 0; 
		Q(i) = P(i); Pm(i) = P(i);
	}


	//let R be matrix with only columns corresponding to end states, diagonal = 1
	//let Q be matrix with rows/cols corresponding to end state zero'd out
	for (int i = 0; i < endStates.size(); i++) {
		int state = endStates[i];
		R.col(state) = Pm.col(state); R(state,state) = 1;
		Q.row(state) = Um.row(state); Q.col(state) = Um.col(state); 
	}

	//solve for Um via (I-Q)Um = R
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(num_states,num_states) - Q;
	Um = D.lu().solve(R);

	//store back in U array
	for (int i = 0; i < num_states*num_states; i++) {
		U(i) = Um(i);
	}
}


void fillDiag(double* T, int num_states) {
	//fill in diagonal elements with negative sum of entries

	double S;
	for (int state = 0; state < num_states; state++) {
		S = 0;
		for (int target = 0; target < num_states; target++) {
			S += T[toIndex(state, target, num_states)];
		}
		if (S == 0) { //if state is inaccesible, set to 1 so system is solvable
			S = 1;
		}
		T[toIndex(state, state, num_states)] = -S;
		//printf("row %d, value %f\n", state, S);
	}
}


void satisfyDB(double* T, int num_states, Database* db, double* eq) {
	//make T satisfy detailed balance

	int row = 0; int column = 0; //row and column indices in matrix

	//if an entry is non-zero, use DB to fill its transpose
	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, row, column);
		if (T[entry] > 0) {
			T[toIndex(column, row, num_states)] = T[entry] * eq[row] / eq[column];
		}
	}
}

void computeFreeEnergy(int num_states, double* Z, double* F) {
	//compute the free energy from configurational partition function

	double minF = 1e10;

	//compute the helmholtz free energy of each cluster from Z
	for (int i = 0; i < num_states; i++) {
		F[i] = -1.0/BETA * log(Z[i]);
		if (F[i] < minF) {
			minF = F[i];
		}
	}

	//shift to make 0 the free energy minimum
	for (int i = 0; i < num_states; i++) {
		F[i] -= minF;
	}
}
void computePartitionFn(int num_states, Database* db, double* Z) {
	//compute the configurational part of the partition fn from each cluster

	createMeasure(num_states, db, Z, 1);
}

void createMeasure(int num_states, Database* db, double* eq, double kappa) {
	//create the equilibrium distribution for this problem

	double kap0 = KAP;                 //sticky parameter for initial measurement
	double stickyRatio = kappa / kap0; //ratio of current to old sticky parameters
	double beta = BETA;                //inverse temp

	double Z = 0;                     //normalizing constant for eq

	//loop over states, get the equilibrium measure entry
	for (int i = 0; i < num_states; i++) {
		//get initial eq prob and number of bonds
		int b = (*db)[i].getBonds();
		double prob = (*db)[i].getFrequency();

		//new eq prob is prob*(ratio)^(bonds)
		double new_prob = prob*pow(stickyRatio,b);

		//increment Z and add to array
		Z += new_prob;
		eq[i] = new_prob;
	}

	//re-normalize
	for (int i = 0; i < num_states; i++) eq[i] /= Z;
}

void createProbabilityMatrix(double* T, int num_states, double* P) {
	//create probability matrix from rate matrix

	int i,j; 

	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, i, j);
		if (i != j) {
			P[entry] = - T[entry] / T[toIndex(i,i,num_states)];
		}
	}
}

void createTransitionMatrix(double* T, int num_states, Database* db, 
														std::vector<int>& endStates) {
	//create rate matrix from data in DB - forward rates

	//declare storage for each state
	double mfpt; double S;

	//loop over all states information
	for (int state = 0; state < num_states; state++) {
		//get all relevant data
		mfpt = (*db)[state].getMFPT();
		S = (*db)[state].sumP(); //get normalizing constant for this row
		//printf("MFPT %f, sum %f\n", mfpt, S);

		//if S = 0, this is end state, add it to vector
		if (S == 0) {
			endStates.push_back(state);
		}

		//fill in value in transition matrix
		std::vector<Pair> P = (*db)[state].getP();
		for (int i = 0; i < P.size(); i++) {
			T[toIndex(state, P[i].index, num_states)] = (P[i].value / S) / mfpt;
		}
	}
}


void writeHittingProbabilityGS(double* kappa, double* data, std::vector<int> endStates, int M) {
	//output the results to a file

	std::ofstream ofile;
	ofile.open("hittingProb.txt");
	ofile << "k ";
	for (int i = 0; i < M; i++) {
		ofile << kappa[i] << ' ';
	}
	ofile << "\n";
	for (int end = 0; end < endStates.size(); end++) {
		ofile << endStates[end] << ' ';
		for (int i = 0; i < M; i++) {
			ofile << data[i*endStates.size()+end] << ' ';
		}
		ofile << "\n";
	}
	ofile.close();
}

void getHittingProbabilityGS(int initial, Database* db) {
	/*compute the hitting probability of each end state as fn of kappa
		output the results to a file */

	//parameters
	int N = db->getN();
	int num_states = db->getNumStates();
	std::vector<int> endStates;
	int M = 200; //number of kappa values to test

	//set the kappa values the test is being done for 
	double* kappa = new double[M];
	for (int i = 0; i < M/2; i++) kappa[i] = float(i)/10.0+0.1;
	for (int i = M/2; i < M; i++) kappa[i] = i-M/2+10.0;

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* P = new double[num_states*num_states]; //probability transition matrix
	double* U = new double[num_states*num_states]; //hitting probability matrix
	double* eq = new double[num_states];           //equilibrium measure
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	createTransitionMatrix(Tconst, num_states, db, endStates);

	//create a data array which will store the hitting probabilities
	double* data = new double[M*endStates.size()];

	//loop over kappa and do computation
	for (int k = 0; k < M; k++) {

		//get sticky parameter
		double sticky = kappa[k]; 

		//reset transition matrices
		for (int i = 0; i < num_states*num_states; i++) {
			P[i] = U[i] = 0;
		}

		//copy Tconst into T
		std::copy(Tconst, Tconst+num_states*num_states, T);

		//make array for equilibrium distribution
		createMeasure(num_states, db, eq, sticky);

		//fill in transposed entries such that T satisfies detailed balance
		satisfyDB(T, num_states, db, eq);

		//fill in diagonal with negative sum of row entries
		fillDiag(T, num_states);

		//use the filled rate matrix to compute probability matrix
		createProbabilityMatrix(T, num_states, P);

		//solve for hitting probabilities to endStates states
		computeHittingProbability(P, num_states, endStates, U);

		//store hitting probabilities in data
		for (int end = 0; end < endStates.size(); end++) {
			double hit_prob = U[toIndex(initial, endStates[end], num_states)];
			data[k*endStates.size() + end] = hit_prob;
		}
	}

	//write to file
	writeHittingProbabilityGS(kappa, data, endStates, M);

	//free the memory
	delete []T; delete []Tconst; delete []P; delete []U;
	delete []eq; delete []kappa; delete []data;
}

void performTPT(int initial, int target, Database* db, bool getIso) {
	//perform tpt calculations from initial to target states

	//parameters
	double kappa = 6.0;
	int N = db->getN();
	int num_states = db->getNumStates();
	std::vector<int> endStates;

	//init and compute configurational partition function and free energy
	double* Z = new double[num_states]; double* F = new double[num_states];
	for (int i = 0; i < num_states; i++) Z[i] = F[i] = 0;
	computePartitionFn(num_states, db, Z); 
	computeFreeEnergy(num_states, Z, F);

	//initialize the transition rate matrix
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] =  0;
	}

	//get bonds -> bonds+1 entries from mfpt estimates
	createTransitionMatrix(T, num_states, db, endStates);

	//init array for equilibrium distribution and compute it
	double* eq = new double[num_states];
	createMeasure(num_states, db, eq, kappa);
	for (int i = 0; i < num_states; i++) std::cout << eq[i] << "\n";

	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);

	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);

	//solve for the committor
	//initialize committor in q
	double* q = new double[num_states];
	for (int i = 0; i < num_states; i++) q[i] = 0;

	//build target state vector. check if including isomorphic states
	//Note: this does nothing if states are already lumped
	std::vector<int> targets; 
	if (getIso == 1) {
		//find all states isomorphic to target
		findIsomorphic(N, num_states, target, db, targets);
	}
	else if (getIso == 0) {
		targets.push_back(target);
	}
	
	//solve dirichlet problem for committor, q 
	computeCommittor(q, T, num_states, initial, targets);

	//init and compute the probability fluxes
	double* flux = new double[num_states*num_states]; 
	for (int i = 0; i < num_states*num_states; i++) flux[i] = 0;
	computeFlux(num_states, q, T, eq, flux);

	double R = computeTransitionRateTPT(num_states, q, T, eq);
	std::cout << "Trans rate " << R << "\n";
	double* m = new double[num_states];
	computeMFPTs(num_states, T, targets, m);
	std::cout << 1/m[1] << "\n";

	//make a graph structure of the database
	Graph* g = makeGraph(db);

	//print out graphviz file with tpt data
	printGraphRev(g, initial, F, flux, 1, 1, 0);

	//free the memory
	delete []T; delete []q; delete []eq; delete []flux;
	delete []Z; delete []F; 
	delete g; delete []m;

}

void getEqRowHitProbability(int num_states, double* eq, double* eqHit, Database* db) {
	//normalize the equilibrium probability over states with same number of bonds

	std::vector<int> same_bonds; 
	bool* found = new bool[num_states];
	for (int i = 0; i < num_states; i++) found[i] = false;


	for (int state = 0; state < num_states; state++) {

		//clear the vector of states with same num of bonds
		same_bonds.clear();
		double Z = 0;  //normalizer

		//check if state has been found previously
		if (!found[state]) {
			int bonds = (*db)[state].getBonds(); same_bonds.push_back(state);
			Z += eq[state]; found[state] = true;
			for (int next = state+1; next < num_states; next++) {
				int next_bonds = (*db)[next].getBonds(); 
				if (next_bonds == bonds) {
					same_bonds.push_back(next);
					Z += eq[next]; found[next] = true;
				}
			}
		}

		//now re-weight by Z and store in eqHit
		for (int entry = 0; entry < same_bonds.size(); entry++) {
			int entry_state = same_bonds[entry];
			eqHit[entry_state] = eq[entry_state] / Z;
		}
	}

	//free memory
	delete []found;

}

void makeProbCompareGraph(int initial, Database* db) {
	//compute the invariant measure in high kappa limit, normalize over rows,
	//make graph, pass it to graphviz code

	//parameters
	double kappa = 2000.0;
	int N = db->getN();
	int num_states = db->getNumStates();

	//init and compute configurational partition function and free energy
	double* Z = new double[num_states]; double* F = new double[num_states];
	for (int i = 0; i < num_states; i++) Z[i] = F[i] = 0;
	computePartitionFn(num_states, db, Z); 
	computeFreeEnergy(num_states, Z, F);

	double* eq = new double[num_states];           //equilibrium measure
	double* eqHit = new double[num_states];        // normalized equilibrium measure
	createMeasure(num_states, db, eq, kappa);
	getEqRowHitProbability(num_states, eq, eqHit, db);

	//make a graph structure of the database
	Graph* g = makeGraph(db);

	//call graphviz creator
	printGraphEqHit(g, initial, eqHit, F, 1);

	//free memory
	delete []Z; delete []F; delete []eq; delete []eqHit;
	delete g;


}

//TODO ROUND THE PROBS THEN DONE//




}