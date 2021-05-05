#pragma once
#include "point.h"
#include "bDynamics.h"
#include "tpt.h"
#include <vector>
#include <map>
#include <deque>
#include <random>
namespace bd { 
class Database; 

//general functions
void getBondTypes(int N, int* particleTypes, Database* db, std::vector<int> targets);
int setTypes(int N, int* particleTypes, int IC);
void setTypes(int N, int* particleTypes, std::deque<std::string> perms, int perm);
double getStickyProduct(int N, int state,  Database* db, int* particleTypes, 
												std::map<std::pair<int,int>,double> kappa);
void reweight(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa);
void reweight7(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa);
void getReweightMaps(std::map<std::pair<int,int>,double> kappa,
										 std::map<std::pair<int,int>,double>& a_new,
										 std::map<std::pair<int,int>,double>& gamma_new); 
double getPdet(int N, int state, Database* db, int* particleTypes, 
							 double a, std::map<std::pair<int,int>,double> a_new);
void initKappaVals(int numInteractions, double* kappaVals);
void allPerms(int N, std::deque<std::string>&);
void distinctPerms(int N, std::deque<std::string>&);
void checkPerm(std::string , std::deque<std::string>&);
void checkPositive(int numInteractions, double* kappaVals); 
void applyMax(int numInteractions, double* kappaVals, double M);

//generic plots
void constructSurfaceTOY(int N, Database* db, int initial, int target, bool useFile);
void constructScatterTOY(int N, Database* db, int initial, int target, bool useFile);
void constructScatterTOY1(int N, Database* db, int initial, int target);
void constructScatterTOYsensitivity(int N, Database* db, int initial, int target, bool useFile);
void perturbDB(Database* db, double freq_perturb_frac, double rate_perturb_frac);
void perturbDBconstP(Database* db);

//hitting probability optimization
void printHitProbability(int num_states, int initial, std::vector<int> targets, double* U);
double getHitProbability(int num_states, int initial, std::vector<int> targets, double* U);
double solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
								 					int numTypes, double* Tconst, std::vector<int> ground, 
								 					std::vector<int> targets);
void solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
								 					double* Tconst, std::vector<int> targets);

double computeGradHP(int initial, int numInteractions,  double* kappaVals, 
										 Database* db, int* particleTypes, int numTypes, 
										 double* Tconst, std::vector<int> ground, 
								 		 std::vector<int> targets, double* g);

void hittingProbMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void hittingProbMaxTOYperms(int N, Database* db, int initial, int target);

//equilibrium probability optimization
void eqProbMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void eqProbMaxTOYperms(int N, Database* db, int initial, int target);
double computeGradEQ(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, std::vector<int> targets, double* g);
double getEqProb(int initial, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, std::vector<int> targets);
double printEqProb(int initial, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, std::vector<int> targets);

//average transition rate optimization
void rateMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void rateMaxTOYperms(int N, Database* db, int initial, int target);
double computeGradRate(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, int numTypes, double* Tconst, std::vector<int> targets, 
								 double* g);
double getRateTPT(int initial, double* kappaVals, Database* db, int* particleTypes, int numTypes,
							 double* Tconst, std::vector<int> targets);
double getRate(int initial, double* kappaVals, Database* db, int* particleTypes, int numTypes,
							 double* Tconst, std::vector<int> targets);

//hitting probability optimization with rate constraints
void hittingProbMaxTOYc(int N, Database* db, int initial, int target, bool useFile);
void eqProbMaxTOYc(int N, Database* db, int initial, int target, bool useFile);
void rateMaxTOYc(int N, Database* db, double c, int initial, int target, bool useFile);
double lineSearchHit(int initial, int numInteractions, double* kappaVals, Database* db, 
										 int* particleTypes, int numTypes, double* Tconst, 
										 std::vector<int> ground, std::vector<int> targets, 
										 double c, double r, double H, double R, double* g, double& step);
double lineSearchEq(int initial, int numInteractions, double* kappaVals, Database* db, 
										int* particleTypes, int numTypes, double* Tconst, 
										std::vector<int> targets, double c, double r, double H, double R, 
										double* g, double& step);
double lineSearchRate(int initial, int numInteractions, double* kappaVals, Database* db, 
											int* particleTypes, int numTypes, double* Tconst, 
											std::vector<int> targets, double c, double r, double eq, double R, 
											double* g, double& step);
void applyRateConstraint(double R, double c, double r, int numInteractions, 
												 double* gH, double* gR); 



//sampling functions
int findMatrix(int N, double* X, Database* db);
void estimateHittingProbability(int N, Database* db, int target);

//testing
void evalStats(int N, Database* db, int initial, int target, bool useFile);
void computeParetoFront(int N, Database* db, int initial, int target, bool useFile);
double rateMaxPareto(double& c, int initial, Database* db, int num_states, int* particleTypes, 
										 int numInteractions, int numTypes, double* kappaVals, double* Tconst, 
										 std::vector<int> ground, std::vector<int> targets);
void computeParetoFrontGD(int N, Database* db, int initial, int target, bool useFile);
void miskinOptimize(int N, Database* db, int initial, int target, bool useFile);
double analyticEqGradient(Database* db, int* particleTypes, int numTypes, int numInteractions, 
												double* kappaVals, double* eq, double* gE, std::vector<int> targets);
void makeIndexMap(int numTypes, int* indexVals, 
									std::map<std::pair<int,int>,int>& index);
int searchForBonds(int N, Database* db, int state, int type, int* particleTypes,
									 std::map<std::pair<int,int>,int>& index);

/* designTPT functions. Intended as an interfacing layer between design code and 
  TPT code. TPT code was built before support for several particle types. 
  These functions will call TPT functions with the proper re-weighting. */

void performTPT(int N, Database* db, int initial, int target, bool useFile);
void lumpArrays(int old_states, int new_states, int* lumpMap, Database* db, double* F, 
								double* Fnew, double* flux, double* fluxNew, double* r, double* rNew, 
								double* eq, double* eqNew, double* q, double* qNew, double* Tnew, 
								int initial, int target);

/* irreversibleDesign functions. */
bool hasInvalidInteraction(int state, Database* db, int* particleTypes,
													 std::map<std::pair<int,int>,double> kappa);
void getKappaQuench(double* kappaVals, std::vector<Point>& kappaQuench);
void getColor(double* kappaVals, std::string& color);
void getPstring(int N, int* particleTypes, std::string& pString);

void killTransitions(double* Tnew, int new_states, double* T, int* lumpMap, Database* db, 
										 int* particleTypes, std::map<std::pair<int,int>,double> kappa);
void killTransitionsQuench(double* Tnew, int new_states, double* T, int* lumpMap, Database* db, 
										 int* particleTypes, std::map<std::pair<int,int>,double> kappa1, 
										 std::map<std::pair<int,int>,double> kappa2, int c1);

void irreversibleGraph(int N, Database* db, int initial, int target, bool useFile);
void generateAllGraphs(int N, Database* db, int initial);
bool generateGraph(int N, Database* db, int initial, int numTypes, int num_interactions, 
									 std::vector<int> ground, int* particleTypes, Point interaction, 
									 double* Tconst, double* kappaVals, int& max_node);
void generateQuenchedGraph(int N, Database* db, int initial, int numTypes, int numInteractions, 
									 std::vector<int> ground, int* particleTypes, double* Tconst, 
									 double* kappaVals1, double* kappaVals2, int c1, int count);
void testQuench(int N, Database* db, int initial, int numTypes, int numInteractions, 
								std::vector<int> ground, int* particleTypes, std::ofstream& ofile, 
								double* Tconst, double* kappaVals);
void findQuenches(int N, Database* db);
void graphQuenches(int N, Database* db, int initial);

void evolveProbability(int N, Database* db, int initial, bool useFile);
void testTransitionTimes(int N, Database* db, int initial, int scheme);


// Time dependent plotting stuff
void getNEQparameters(std::string& p_file, Eigen::VectorXd& beta, Eigen::MatrixXd& E);
void phaseTrajectory(int N, Database* db, int initial, int target, bool, std::string& p_file);


//test measures for sampling 
void testMeasures(int N, Database* db, int initial, int target, bool useFile);
void testMeasuresRecord(int N, Database* db, int initial, int target, bool useFile);
void testMeasuresScatter(int N, Database* db, int initial, int target, bool useFile);








}