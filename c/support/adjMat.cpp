#include <math.h>
#include <stdio.h>
#include "nauty.h"
#include "database.h"
#include "adjacency.h"
#include <vector>
#include <eigen3/Eigen/Dense>
namespace bd{

class Database;

int toIndex(int r, int c, long m) {
  //map row and column number into index in 1d array. column indexed
  return m*c+r;
}

void index2ij(int index, int N, int& i, int& j) {
  // map a row index 1-d array into (i,j) 
  i = index % N;
  j = index / N;
}

int checkConnected(int* M, int N) {
	//check if subdiagonal sums to N-1
	int S = 0;
	for (int i = 0; i < N-1; i++) {
		S += M[(N+1)*i+1];
	}
	if (S == N-1) {
		return 1;
	}
	else {
		return 0;
	}
}

int checkSame(int* M1, int* M2, int N) {
	//check if adjacency matrices M1 and M2 are the same
	for (int i = 0; i < N*N; i++) {
		if (M1[i]-M2[i] != 0) {
			return 0;
		}
	}
	return 1;
}

void getAdj(double* X, int N, int* M) {
	//get the adjacnecy matrix from a cluster X
	int i, j; double tol = 1e-5;
	double* particles = new double[2*N]; double* Z = new double[2];
	c2p(X, particles, N);
	for (int index = 0; index < N*N; index++) {
		index2ij(index, N, i, j);
		if (j > i) {
			double d = euDist(particles, i, j, N, Z);
			if (d < 1.1 + tol) {
				M[index] = 1; M[toIndex(j,i,N)] = 1;
			}
		}
	}
	delete []particles; delete []Z;
}

void buildNautyGraph(int N, int M, int state, Database* db, graph* g) {
	//build nauty graph of given state

	//zero out any pre-existing graph
	EMPTYGRAPH(g, M, N);

	//add edges for every non-zero entry in adjacnecy matrix
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			if ((*db)[state].isInteracting(i,j,N)) {
				ADDONEEDGE(g, i, j, M);
			}
		}
	} 
}

bool checkIsomorphic(int N, int M, graph* g1, graph* g2) {
	//check if two canonically labeled graphs are isomorphic

	size_t k;
	for (k = 0; k < M*(size_t)N; ++k) {
			if (g1[k] != g2[k]) { //label is different, return 0
				return false;
			}
		}
		if (k == M*(size_t)N) { //graphs are isomorphic, return 1
			return true;
		}

}


void findIsomorphic(int N, int num_states, int state, Database* db, std::vector<int>& iso) {
	//finds all states isomorphic to given state in database. stored to iso.

	//set keywords for nauty
	int M = SETWORDSNEEDED(N);

	//initialize the nauty graphs and parameters
	graph g1[N*M]; graph g2[N*M];   //the starting graphs
	graph cg1[N*M]; graph cg2[N*M]; //graphs with canonical labeling
	int lab1[N], lab2[N];           //label arrays
	int ptn[N], orbits[N];          //needed to call functions
	static DEFAULTOPTIONS_GRAPH(options); //defualt options
	statsblk stats;                 //nauty statistics of graph
	options.getcanon = TRUE;        //get canonical labeling
	options.defaultptn = TRUE;      //ignore any coloring of graph

	//build nauty graph of the state, get canonical labeling
	EMPTYGRAPH(cg1, M, N);
	buildNautyGraph(N, M, state, db, g1);
	densenauty(g1, lab1, ptn, orbits, &options, &stats, M, N, cg1);

	//loop over all states, check if isomorphic
	for (int i = 0; i < num_states; i++) {
		//create graph
		EMPTYGRAPH(cg2, M, N);
		buildNautyGraph(N, M, i, db, g2);
		densenauty(g2, lab2, ptn, orbits, &options, &stats, M, N, cg2);

		/* Compare canonically labelled graphs */
		if (checkIsomorphic(N, M, cg1, cg2) == 1) { //graphs are isomorphic, add to iso
			iso.push_back(i);
		}
	} 
}


}