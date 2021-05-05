#pragma once
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/QR>
namespace bd { 
class Database; 

//general functions
//perform tpt on initial and target states
void performTPT(int initial, int target, Database* db, bool getIso);
// get hitting probabilities to each end state as fn of kappa
void getHittingProbabilityGS(int initial, Database* db);
//write the hitting probabilities to a file
void writeHittingProbabilityGS(double* kappa, double* data, std::vector<int> endStates, int M);
//normalize the invariant measure over state with same number of bonds
void getEqRowHitProbability(int num_states, double* eq, double* eqHit, Database* db); 
//make graph comparing eq with hitting prob
void makeProbCompareGraph(int initial, Database* db);

//transition matrix functions
//fill in rate matrix with data from mfpt estimates
void createTransitionMatrix(double* T, int num_states, Database* db, 
														std::vector<int>& endStates);
//create the invariant measure for the current problem
void createMeasure(int num_states, Database* db, double* eq, double kappa);
//create configurational partition fn for each cluster
void computePartitionFn(int num_states, Database* db, double* Z); 
//edit rate matrix to satisfy detailed balance
void satisfyDB(double* T, int num_states, Database* db, double* eq);
//fill in diagonal of rate matrix
void fillDiag(double* T, int num_states);
//create probability transition matrix
void createProbabilityMatrix(double* T, int num_states, double* P); 

//tpt functions
//set up and solve dirchlet problem for committor function
void computeCommittor(double* q, double* T, int num_states, int initial, std::vector<int>);
//compute the probability flux from committor
void computeFlux(int num_states, double* q, double* T, double* eq, double* flux);
//compute the average transition rate from A to B on reactive trajectories
double computeTransitionRateTPT(int num_states, double* q, double* T, double* eq);
//compute free energy of each cluster - from configurational partition fn
void computeFreeEnergy(int num_states, double* Z, double* F); 
//compute hitting probabilities for end states in reversible case
void computeHittingProbability(double* P, int num_states, std::vector<int> endStates, 
															 double* U);
void computeHittingProbability(Eigen::MatrixXd& P, int num_states, std::vector<int> endStates, 
															 Eigen::MatrixXd& U);
//compute the mfpt from worm to target states
void computeMFPTs(int num_states, double* T, std::vector<int> targets, double* m);
void computeMFPTsSP(int num_states, double* T, std::vector<int> targets, double* m);

























}