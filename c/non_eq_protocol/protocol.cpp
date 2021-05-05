#include "design.h"
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <deque>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "point.h"
#include "pair.h"
#include "adjacency.h"
#include "database.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "protocol.h"
#include "../defines.h"
namespace bd{

#define LOWER 0.005
#define GAM exp(5.58)
#define WHICH 1


/****************************************************/
/*********** Optimization class    *******************/
/****************************************************/



OptInfo::OptInfo(int N_, int num_states_, int initial_, int target_, 
						double E_, double rho_, Database* db_) {
	problem_type = 1;
	N = N_; num_states = num_states_; initial = initial_; target = target_;
	numInteractions = 0; numTypes = 0;
	E = E_; rho = rho_;
	h1 = 0; h2 = 0;
	db = db_;
	particleTypes = NULL;
}

OptInfo::OptInfo(int N_, int initial_, int target_, 
						double h1_, double h2_) {
	problem_type = 2;
	N = N_; num_states = 2; initial = initial_; target = target_;
	numInteractions = 0; numTypes = 0;
	E = 0; rho = 0;
	h1 = h1_; h2 = h2_;
	db = NULL;
	particleTypes = NULL;
}

OptInfo::OptInfo(int N_, int num_states_, int initial_, int target_, 
						const Eigen::MatrixXd& h_) {
	problem_type = 3;
	N = N_; num_states = num_states_; initial = initial_; target = target_;
	numInteractions = 0; numTypes = 0;
	h = h_;
	E = 0; rho = 0;
	h1 = 0; h2 = 0;
	db = NULL;
	particleTypes = NULL;
}

OptInfo::OptInfo(int N_, int num_states_, int initial_, int target_, 
						double h1_, double h2_) {
	problem_type = 4;
	N = N_; num_states = num_states_; initial = initial_; target = target_;
	numInteractions = 0; numTypes = 0;
	E = 0; rho = 0;
	h1 = h1_; h2 = h2_;
	db = NULL;
	particleTypes = NULL;
}

void OptInfo::CreateRateMatrix(Eigen::MatrixXd& R, double beta) {
	//call the correct function

	if (problem_type == 1) {
		createRateMatrix(R, num_states, db, beta, E, rho);
	}
	else if (problem_type == 2) {
		createRateMatrix(R, beta, h1, h2);
	}
	else if (problem_type == 3) {
		createRateMatrix(R, num_states, beta, h);
	}
	else if (problem_type == 4) {
		createRateMatrix(R, num_states, beta, h1, h2);
	}
}

void OptInfo::CreateRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, double beta) {
	//call the correct function

	//printf("Prob type is %d\n", problem_type);

	if (problem_type == 1) {
		createRateMatrix(R, D, num_states, db, beta, E, rho);
	}
	else if (problem_type == 2) {
		createRateMatrix(R, D, beta, h1, h2);
	}
	else if (problem_type == 3) {
		createRateMatrix(R, D, num_states, beta, h);
	}
	else if (problem_type == 4) {
		createRateMatrix(R, D, num_states, beta, h1, h2);
	}
}

/****************************************************/
/*********** Rate Functions - SA    *******************/
/****************************************************/


void createRateMatrix(Eigen::MatrixXd& R, int num_states, Database* db, 
												double beta, double E, double rho) {
	//use the database to construct the rate matrix at the given values

	//pre-compute some constants that go in later expressions
	double gamma = exp(beta*E);
	double kappa = gamma / sqrt(2.0 * rho * rho * beta * E);
	double kap0  = KAP;
	double gamma0 = GAM;

	//loop over all states. compute forward entries in both matrices
	for (int state = 0; state < num_states; state++) {
		//determine the forward rate from mfpt and probability distribution
		double mfpt = (*db)[state].getMFPT();
		double S = (*db)[state].sumP(); //get normalizing constant for this row

		//get each pair of transitioning states and fill in the entry
		std::vector<Pair> P = (*db)[state].getP();
		for (int t_index = 0; t_index < P.size(); t_index++) {
			//set the forward rates
			int target = P[t_index].index;
			double c = ((double(P[t_index].value) / S) / mfpt);
			R(state, target) = c * (1.0 / beta) * (1.0-exp(-beta)) / (1.0-exp(-1.0)) ;

			//set the backward rates
			int b_target = (*db)[target].getBonds();
			int b_state  = (*db)[state].getBonds();
			int b_diff = b_state - b_target;

			double eq_target = (*db)[target].getFrequency();
			double eq_state  = (*db)[state].getFrequency();

			double reweight_bond = (kap0 / kappa);
			if (b_diff == -2) {
				reweight_bond *= (gamma0 / gamma);
			}
			R(target, state) = R(state, target) * (eq_state / eq_target) * reweight_bond;

			if (beta < 0.7) {
				//R(state,target) = 0; 
			}
			//printf("i %d, j %d, rate %f\n", target, state, R(target,state));
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
}

void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, int num_states, Database* db, 
												double beta, double E, double rho) {
	//use the database to construct the rate matrix and its beta derivative at the
	//given values

	//pre-compute some constants that go in later expressions
	double gamma = exp(beta*E);
	double kappa = gamma / sqrt(2.0 * rho * rho * beta * E);
	double kap0  = KAP;
	double gamma0 = GAM;

	//loop over all states. compute forward entries in both matrices
	#pragma omp parallel for
	for (int state = 0; state < num_states; state++) {
		//determine the forward rate from mfpt and probability distribution
		double mfpt = (*db)[state].getMFPT();
		double S = (*db)[state].sumP(); //get normalizing constant for this row

		//get each pair of transitioning states and fill in the entry
		std::vector<Pair> P = (*db)[state].getP();
		for (int t_index = 0; t_index < P.size(); t_index++) {
			//set the forward rates
			int target = P[t_index].index;
			double c = ((double(P[t_index].value) / S) / mfpt);
			R(state, target) = c * (1.0 / beta) * (1.0-exp(-beta)) / (1.0-exp(-1.0)) ;
			D(state, target) = - R(state,target) * (1.0/beta - exp(-beta)/(1.0-exp(-beta))) ;
			if (beta < 0.7) {
				R(state,target) = 0; D(state,target) = 0;
			}

			//set the backward rates
			int b_target = (*db)[target].getBonds();
			int b_state  = (*db)[state].getBonds();
			int b_diff = b_state - b_target;

			double eq_target = (*db)[target].getFrequency();
			double eq_state  = (*db)[state].getFrequency();

			double reweight_bond = (kap0 / kappa);
			if (b_diff == -2) {
				reweight_bond *= (gamma0 / gamma);
			}
			//printf("rweight term is %f\n", reweight_bond);
			R(target, state) = R(state, target) * (eq_state / eq_target) * reweight_bond;
			D(target, state) = R(target, state) * (1.0 / beta) * (b_diff * beta*E - 0.5 + beta*exp(-beta)/(1.0-exp(-beta)));
		
			if (beta < 0.7) {
				//R(state,target) = 0; D(state,target) = 0;
			}
		}

	}

	//finally, must set the diagonal to be negative of the row sum
	for (int i = 0; i < num_states; i++) {
		double R_sum = 0; double D_sum = 0;

		for (int j = 0; j < num_states; j++) {
			R_sum += R(i,j); D_sum += D(i,j);
		}

		R(i,i) = -R_sum; D(i,i) = -D_sum;
	}
}

void createRateMatrix(Eigen::MatrixXd& R, double beta, double h1, double h2) {
	//create rate matrix for 2 state hopping problem

	double F = exp(-beta*h1); double B = exp(-beta*h2);
	R(0,0) = -F; R(0,1) = F;
	R(1,0) = B; R(1,1) = -B;
}

void createRateMatrix(Eigen::MatrixXd& R, int num_states, double beta, 
											const Eigen::MatrixXd& h) {
	//create rate matrix for model system markov chain

	double rho = 1.0 / sqrt(2.0);
	double E = 6;
	double kappa = 1.0 / sqrt(2*rho*rho*beta*E) * exp(beta*E);

	//fill in transitions from barrier heights
	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (h(i,j) > 0) {
				if (WHICH == 0) {
					R(i,j) = exp(-beta*h(i,j));
				}
				else if (WHICH == 1) {
					if (i > j) {

						R(i,j) = h(i,j) / beta * (1.0 / kappa) * sinh(sqrt(beta*E));
					}
					else {
						R(i,j) = h(i,j) / beta * exp(-1.0/(beta*beta*beta*beta));
					}
				}
			}
		}
	}

	//fill in diagonal entries
	for (int i = 0; i < num_states; i++) {
		R(i,i) = -R.row(i).sum();
	}
}


void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, double beta, 
											double h1, double h2) {
	//create rate matrix for 2 state hopping problem

	double F = exp(-beta*h1); double B = exp(-beta*h2);
	R(0,0) = -F; R(0,1) = F;
	R(1,0) = B; R(1,1) = -B;

	D(0,0) = h1*F; D(0,1) = -h1*F;
	D(1,0) = -h2*B; D(1,1) = h2*B;
}

void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, int num_states, double beta, 
											const Eigen::MatrixXd& h) {
	//create rate matrix for model system markov chain and derivative

	//fill in transitions from barrier heights

	double rho = 1.0 / sqrt(2.0);
	double E = 6;
	double kappa = 1.0 / sqrt(2*rho*rho*beta*E) * exp(beta*E);

	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (h(i,j) > 0) {
				if (WHICH == 0) {
					R(i,j) = exp(-beta*h(i,j));
					D(i,j) = -h(i,j) * exp(-beta*h(i,j));
				}
				else if (WHICH == 1) {
					if (i > j) {
						R(i,j) = h(i,j) / beta * (1.0 / kappa) * sinh(sqrt(beta*E));
						D(i,j) = - h(i,j) * (1.0/(beta*beta)) * (1.0/kappa) * (beta*E+0.5-0.5*sqrt(beta)*cosh(sqrt(beta))) ;
					}
					else {
						R(i,j) = h(i,j) / beta * exp(-1.0/(beta*beta*beta*beta));
						//std::cout << R(i,j) << "\n";
						D(i,j) =  - R(i,j) / beta * (1-4.0/(beta*beta*beta*beta));
					}
				}
			}
		}
	}

	//fill in diagonal entries
	for (int i = 0; i < num_states; i++) {
		R(i,i) = -R.row(i).sum();
		D(i,i) = -D.row(i).sum();
	}
}

void createRateMatrix(Eigen::MatrixXd& R, int num_states, double beta, 
											double h1, double h2) {
	//rate matrix for a problem with 1 row of states to transition to

	double pi = 3.1415926;
	double c = 100.0;

	double down_rate = 1.0 / (double(num_states) - 1.0);
	double sig = 1.0 / pi * atan(c*(beta-1.0)) + 0.5;

	//fill in rates for all non-target states
	for (int i = 1; i < num_states-1; i++) {
		//uncomment for energy barriers
		/*
		R(0,i) = down_rate; 
		//R(i,0) = exp(-beta*h1); 
		//D(i,0) = -h1*exp(-beta*h1); D(i,i) = - D(i,0);
		*/

		//uncomment for step functions
		/*
		if (beta > 1.0) {
			R(0,i) = down_rate; 
		}
		if (beta < 1.0) {
			R(i,0) = 1.0;
		}
		*/

		//uncomment for continuous step - sigmoid
		R(0,i) = down_rate * sig;
		R(i,0) = (1.0-sig);

		R(i,i) = - R(i,0);
	}

	//do the target state
	int target = num_states-1;

	//uncomment for step fns
	/*
	if (beta > 1.0) {
		R(0,target) = down_rate;
	}
	if (beta < 1.0) {
		//R(target,0) = exp(-beta*h2); 
		//D(target,0) = -h2*exp(-beta*h2); 
		R(target,0) = 0.01;
	}
	*/

	//sigmoid
	double gamma = 0.01;
	R(0,target) = down_rate * sig;
	R(target,0) = (1.0-sig) * gamma;

	R(target,target) = - R(target,0);
	
	//fill in (0,0)
	R(0,0) = -R.row(0).sum();

}

void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, int num_states, double beta, 
											double h1, double h2) {
	//rate matrix for a problem with 1 row of states to transition to

	double pi = 3.1415926;
	double c = 100.0;

	double down_rate = 1.0 / (double(num_states) - 1.0);
	double sig = 1.0 / pi * atan(c*(beta-1.0)) + 0.5;
	double sigd = c / pi * 1.0 / (1.0 + c*c*(beta-1.0)*(beta-1.0));

	//fill in rates for all non-target states
	for (int i = 1; i < num_states-1; i++) {
		//uncomment for energy barriers
		/*
		R(0,i) = down_rate; 
		//R(i,0) = exp(-beta*h1); 
		//D(i,0) = -h1*exp(-beta*h1); D(i,i) = - D(i,0);
		*/

		//uncomment for step functions
		/*
		if (beta > 1.0) {
			R(0,i) = down_rate; 
		}
		if (beta < 1.0) {
			R(i,0) = 1.0;
		}
		*/

		//uncomment for continuous step - sigmoid
		R(0,i) = down_rate * sig;
		D(0,i) = down_rate * sigd;
		R(i,0) = (1.0-sig);
		D(i,0) = -sigd;

		R(i,i) = - R(i,0);
		D(i,i) = - D(i,0);
	}

	//do the target state
	int target = num_states-1;

	//uncomment for step fns
	/*
	if (beta > 1.0) {
		R(0,target) = down_rate;
	}
	if (beta < 1.0) {
		//R(target,0) = exp(-beta*h2); 
		//D(target,0) = -h2*exp(-beta*h2); 
		R(target,0) = 0.01;
	}
	*/

	//sigmoid
	double gamma = 0.01;
	R(0,target) = down_rate * sig;
	D(0,target) = down_rate * sigd;
	R(target,0) = (1.0-sig) * gamma;
	D(target,0) = -sigd * gamma;

	R(target,target) = - R(target,0);
	D(target,target) = - D(target,0);

	//fill in (0,0)
	R(0,0) = -R.row(0).sum();
	D(0,0) = -D.row(0).sum();


}

/****************************************************/
/*********** Transition Operators  *******************/
/****************************************************/

void createTransitionOperator(int N, int num_states, OptInfo* problem, 
															const Eigen::VectorXd& betaP, const Eigen::VectorXd& t_disc, 
															double t0, double tf, int interval, 
															Eigen::MatrixXd& trans_op) {
	//construct the transition operator from t0 to tf, given the protocol

	//begin by setting the transition operator equal to the identity
	trans_op = Eigen::MatrixXd::Identity(num_states,num_states);

	//define a storage for the rate matrices
	Eigen::MatrixXd rate; 

	//check if t0 is 0 and tf is not T. if yes, propogate until tf
	if (t0 == 0 && tf != t_disc(N)) {
		//the final index is interval
		int final = interval; 
		//printf("Final interval is %d\n", interval);
		//printf("Flag1\n");

		//loop over each interval that gets the full time step
		for (int i = 0; i < final; i++) {
			rate.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, betaP(i));
			double delta_t = t_disc(i+1) - t_disc(i);
			rate *= delta_t;
			trans_op *= rate.exp();
		}

		//the final interval gets a partial step
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, betaP(final));
		double delta_t = tf - t_disc(final);
		rate *= delta_t;
		trans_op *= rate.exp();

		return;
	}

	//check if tf is the final time and t0 is not 0
	if (tf == t_disc(N) && t0 != 0) {
		//int initial index is interval
		int first = interval;
		//printf("Flag2\n");
		//printf("First interval is %d\n", interval);

		//the first interval gets a partial step
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, betaP(first));
		double delta_t = t0 - t_disc(first);
		rate *= delta_t;
		trans_op *= rate.exp();

		//loop over the rest of the intervals
		for (int i = first+1; i < N; i++) {
			rate.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, betaP(i));
			double delta_t = t_disc(i+1) - t_disc(i);
			rate *= delta_t;
			trans_op *= rate.exp();
		}

		return;
	}

	//check if t0 == 0 and tf == T
	if (tf == t_disc(N) && t0 == 0) {

		//printf("Flag3\n");
		for (int i = 0; i < N; i++) {
			rate.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, betaP(i));
			double delta_t = t_disc(i+1) - t_disc(i);
			rate *= delta_t;
			trans_op *= rate.exp();
		}

		return;
	}

}

double computeFinalProb(int N, int num_states, OptInfo* problem,
												const Eigen::VectorXd& betaP, const Eigen::VectorXd& t_disc,  
												int initial, int target) {

	Eigen::MatrixXd trans_op = Eigen::MatrixXd::Identity(num_states,num_states);
	Eigen::MatrixXd rate; rate.setZero(num_states,num_states);

	for (int i = 0; i < N; i++) {
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, betaP(i));
		double delta_t = t_disc(i+1) - t_disc(i);
		rate *= delta_t;
		trans_op *= rate.exp();
	}

	return trans_op(initial,target);
}


/****************************************************/
/*********** Helper  Functions    *******************/
/****************************************************/


void applyLowerBound(Eigen::VectorXd& vec, int length, double bound) {
	//apply the lower bound to the given vector

	for (int i = 0; i < length; i++) {
		if (vec(i) < bound) {
			vec(i) = bound;
		}
	}
}


/****************************************************/
/*********** Optimization Functions    **************/
/****************************************************/

double lineSearch(int N, Eigen::VectorXd& betaP, int num_states, OptInfo* problem,
								const Eigen::VectorXd& t_disc, const Eigen::VectorXd& R,
								int initial, int target) {
	//armijo line search
	//starting with step 1, reduce by 0.5 until objective function is increased

	int maxReductions = 16;
	double initStep = 1.0;
	double alpha = 0.5;       //reduction factor
	double tol = 1e-16;

	//current value
	double f0 = computeFinalProb(N, num_states, problem, betaP, t_disc, initial, target);
	double f1;

	//loop over step sizes
	for (int i = 0; i < maxReductions; i++) {
		Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest = betaP + initStep * R;
		applyLowerBound(betaTest, N, 0.005);
		//std::cout << betaTest << "\n";
		f1 = computeFinalProb(N, num_states, problem, betaTest, t_disc, initial, target);

		//printf("diff is %.12f\n", f1-f0);

		if (f1 > f0 - tol) {
			//stop 
			betaP = betaTest;
			return f1;
		}


		initStep *= alpha;
	}

	//printf("LInesearch failed\n");

	return -1;
}

void evalGradient(int N, const Eigen::VectorXd& betaP, int num_states, 
							 OptInfo* problem, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, int target,
							 Eigen::VectorXd& R) {

	//loop over each component of R, compute the gradient
	for (int j = 0; j < N; j++) {
		//compute the derivative of the rate matrix at current beta val
		Eigen::MatrixXd rate; rate.setZero(num_states,num_states);
		Eigen::MatrixXd LD; LD.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, LD, betaP(j));
		//std::cout << rate << "\n";
		//std::cout << rate(7,7) << "\n";
		//Eigen::MatrixXd test(num_states,num_states); test = rate * 100; test = test.exp();
		//std::cout << test << "\n";
		//abort();

		//apply the transition operator from 0 to t_center(j), then LD, then to T
		Eigen::MatrixXd T1; T1.setZero(num_states,num_states);
		Eigen::MatrixXd T2; T2.setZero(num_states,num_states);
		Eigen::MatrixXd trans_op; trans_op.setZero(num_states,num_states);
		
		createTransitionOperator(N, num_states, problem, betaP, t_disc, 0, 
														t_centers(j), j, T1);
		createTransitionOperator(N, num_states, problem, betaP, t_disc, t_centers(j), 
														t_disc(N), j, T2);
		trans_op = T1 * LD * T2;
		//std::cout << LD.row(initial) << "\n";
		//abort();

		//std::cout << trans_op << "\n";


		//extract the relevant entry: row = initial, column = target
		R(j) = trans_op(initial, target);

		//check for a move that violates the lower bound
		if (betaP(j) <= LOWER && R(j) < 0) {
			R(j) = 0;
		}
	}
	//abort();
	//std::cout << R << "\n\n";

}

void performCG(int N, int max_iters, Eigen::VectorXd& betaP, int num_states, 
							 OptInfo* problem, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, int target, double step,
							 double lowerBound, double grad_tol, double input_tol, bool verbose) {
	//perform conjugate gradient on the optimization problem

	//declare a previous beta storage to check input tolerance
	Eigen::VectorXd prev_beta; prev_beta.setZero(N); prev_beta = betaP;

	//declare gradient and direction vectors
	Eigen::VectorXd R; R.setZero(N); //gradient storage current
	Eigen::VectorXd prev_R; prev_R.setZero(N);       //gradient storage past
	Eigen::VectorXd p; p.setZero(N); //direction storage current
	Eigen::VectorXd prev_p; prev_p.setZero(N);       //direction storage past

	//do the iterations
	for (int iter = 0; iter < max_iters; iter++) {

		//init the gradient vector and evaluate it
		R.setZero(N);
		evalGradient(N, betaP, num_states, problem, t_disc, t_centers, initial, target, R);

		//perform the CG update
		//get gradient norm
		double r_norm = R.norm();

		//compute the coefficient on prev_p - g_k^2/g_{k-1}^2
		//double c = (R.dot(R-prev_R)) / (prev_R.dot(prev_R));
		double c;
		if (iter == 0) {
			c = 0.0;
		}
		else {
			//pr form
			//c = (R.dot(R-prev_R)) / (prev_R.dot(prev_R));
			//c = max(c,0.0);

			//fr form
			c = (R.dot(R)) / (prev_R.dot(prev_R));
			c = 0;
		}
		printf("c is %f\n", c);

		//set the search direction and normalize it
		p = R + c * prev_p;
		double p_norm = p.norm();
		p /= p_norm;

		//do line search
		double obj = lineSearch(N, betaP, num_states, problem, t_disc, p, initial, target);
		if (obj < 0) {
			betaP +=  R;
			applyLowerBound(betaP, N, lowerBound);
		}
		//printf("Prob is %.12f\n", obj);
		
		//std::cout << betaP << "\n\n";

		//check for gradient tolerance
		if (r_norm < grad_tol) {
			printf("Met the gradient tol %f after %d iterations\n", r_norm, iter);
			break;
		}

		//check for input tolerance
		prev_beta -= betaP;
		double delta = prev_beta.norm();
		//printf("delta is %.12f\n", delta);
		if (delta < input_tol) {
			printf("Met the input tol %f after %d iterations\n", delta, iter);
			break;
		}
		prev_beta = betaP;

		//print status if requested
		if (verbose) {
			printf("Iteration %d completed. gres %f, ires %f\n", iter, r_norm, delta);
		}

		//set previous values
		prev_R = R; prev_p = p;
	}
}


void performGD(int N, int max_iters, Eigen::VectorXd& betaP, int num_states, 
							 OptInfo* problem, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, int target, double step,
							 double lowerBound, double grad_tol, double input_tol, bool verbose) {
	//perform gradient descent on the given optimization problem

	Eigen::VectorXd prev_beta; prev_beta.setZero(N); prev_beta = betaP;

	for (int iter = 0; iter < max_iters; iter++) {

		//init the gradient vector and evaluate it
		Eigen::VectorXd R; R.setZero(N);
		evalGradient(N, betaP, num_states, problem, t_disc, t_centers, initial, target, R);

		//perform the GD update
		double r_norm = R.norm();
		if (step == 0) {//do line search
			R /= r_norm;
			//std::cout << betaP << "\n";
			//std::cout << R << "\n";
			//std::cout << r_norm << " hi\n";
			double obj = lineSearch(N, betaP, num_states, problem, t_disc, R, initial, target);
			if (obj < 0) {
				betaP +=  R * r_norm;
				applyLowerBound(betaP, N, lowerBound);
			}
			//printf("Prob is %.12f\n", obj);
		}
		else {
			betaP += step * R;
			applyLowerBound(betaP, N, lowerBound);
		}
		//std::cout << betaP << "\n\n";

		//check for gradient tolerance
		if (r_norm < grad_tol) {
			printf("Met the gradient tol %f after %d iterations\n", r_norm, iter);
			break;
		}

		//check for input tolerance
		prev_beta -= betaP;
		double delta = prev_beta.norm();
		//printf("delta is %.12f\n", delta);
		if (delta < input_tol) {
			printf("Met the input tol %f after %d iterations\n", delta, iter);
			break;
		}
		prev_beta = betaP;

		if (verbose) {
			printf("Iteration %d completed. gres %f, ires %f\n", iter, r_norm, delta);
		}
	}
}


void findProtocol(int np, Database* db, int initial, int target) {
	//determine a temperature protocol to maximize prob of target state

	//set problem parameters
	int num_states = db->getNumStates();  //number of states in markov chain
	double T = 4;                            //final time
	double E = 15.0;                       //bond strengths
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 16;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//steepest descent parameters
	int max_iters = 500;             //number of iterations of gradient descent
	double step = 0.005;               //constant step size
	double grad_tol = 1e-8;          //termination tolerance on gradient norm
	double input_tol = 1e-8;         //termintion tolerance on dx norm
	double lowerBound = 1.0 / (2*E); //lower bound for beta (beta*E > 0.5)

	//init vectors for beta storage and construct initial guess
	Eigen::VectorXd betaP; betaP.setZero(N); betaP.setOnes(N); betaP*=0.5;//betaP.setLinSpaced(N,0.4,1.5);
	Eigen::VectorXd prev_beta; prev_beta.setZero(N); prev_beta += betaP;

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, E, rho, db);

	double test_val = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
	printf("The initial guess gives final prob %f\n", test_val);

	//perform steepest descent
	performGD(N, max_iters, betaP, num_states, &problem, t_disc, t_centers, initial, target,
						step, lowerBound, grad_tol, input_tol, true);

	//output the result
	std::cout << betaP << "\n";

	double final_val = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
	printf("The protocol gives final prob %f\n", final_val);


}

void findProtocol2state() {
	//test the functions on the 2 state problem

	//set problem parameters
	int num_states = 2;  //number of states in markov chain
	int T = 32;                            //final time
	double h1 = 4.0;                      //set barrier height1
	double h2 = 6.0;                      //barrier height 2
	int initial = 0;                      //starting state
	int target = 1;                       //final state

	//discretization parameters
	int N = 1;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;
	t_centers(0) = 0;

	//steepest descent parameters
	int max_iters = 500;             //number of iterations of gradient descent
	double step = 0.5;               //constant step size
	double grad_tol = 1e-17;          //termination tolerance on gradient norm
	double input_tol = 1e-14;         //termintion tolerance on dx norm
	double lowerBound = 0.005; 

	step = 0;

	//init vectors for beta storage and construct initial guess
	Eigen::VectorXd betaP; betaP.setZero(N); betaP.setOnes(N);
	Eigen::VectorXd prev_beta; prev_beta.setZero(N); prev_beta += betaP;

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, initial, target, h1, h2);

	double test_val = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
	printf("The initial guess gives final prob %f\n", test_val);

	//perform steepest descent
	performGD(N, max_iters, betaP, num_states, &problem, t_disc, t_centers, initial, target,
						step, lowerBound, grad_tol, input_tol, true);

	//output the result
	std::cout << betaP << "\n";

	double final_val = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
	printf("The protocol gives final prob %f\n", final_val);
}

void getGaps2state() {
	//determine the probability gaps between the optimal constant and the optimal protocol
	//for the 2 state system as a function of final time


	//set problem parameters
	int num_states = 2;  //number of states in markov chain
	double h1 = 4.0;                      //set barrier height1
	double h2 = 6.0;                      //barrier height 2
	int initial = 0;                      //starting state
	int target = 1;                       //final state

	//steepest descent parameters
	int max_iters = 5000;             //number of iterations of gradient descent
	double step = 0.95;               //constant step size
	double grad_tol = 1e-14;          //termination tolerance on gradient norm
	double input_tol = 1e-13;         //termintion tolerance on dx norm
	double lowerBound = 0.005;       //lower bound for beta, cant go negative

	//discretization parameters
	int N = 1;                            //number of intervals in termporal discretization
	long int maxT = 10000000000;                        //the biggest Tf to sample
	
	std::vector<double> op;
	std::vector<double> cp;
	std::vector<double> optBeta;

	step = 0; //do linesearch

	int mult = 2; //multiply by this to sample logarithmically
	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, initial, target, h1, h2);

	//loop over final times, get the gap between the optimal protocol and optimal constant
	//#pragma omp parallel for
	for (long int T = 1; T < maxT; T*=mult) {

		double dt = double(T) / double(N);    //size of time intervals

		Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
		t_disc.setLinSpaced(N+1,0,T);

		Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
		t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
		t_shift += t_disc;

		Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
		t_centers = t_shift.segment(0,N) / 2.0;

		//init vectors for beta storage and construct initial guess
		Eigen::VectorXd betaP; betaP.setZero(N); 
		Eigen::VectorXd prev_beta; prev_beta.setZero(N); 

		//determine the best constant
		int tests = 300;
		Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,5.5);
		double maxP = 0; double maxBeta = 0;
		for (int i = 0; i < tests; i++) {
			betaP.setOnes(N); betaP *= betaTest(i);
			double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
			if (testP > maxP) {
				maxP = testP;
				maxBeta = betaTest(i);
			}
		}

		//determine the best protocol
		betaP.setOnes(N); //betaP*= maxBeta;
		betaP *= 2.5;
		performGD(N, max_iters, betaP, num_states, &problem, t_disc, t_centers, initial, target,
							step, lowerBound, grad_tol, input_tol, false);
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);

		//store optimal results
		op.push_back(testP); cp.push_back(maxP); optBeta.push_back(maxBeta);
		std::cout << maxP << " " << maxBeta <<  "\n";
		printf("GD gives %f\n", testP);
		std::cout << betaP << "\n";
	}

	//create a file for output
	std::ofstream ofile;
	ofile.open("optimal2stateGaps.txt");

	//print the result
	int count = 0;
	for (long int T = 1; T < maxT; T*=mult) {
		ofile << T << " " << cp[count] << " " << op[count] << " " << optBeta[count] << "\n";
		count++;
	}

	//close outfile
	ofile.close();
	


}

void flowerTestConstantBeta(int np, Database* db) {

	//set problem parameters
	int num_states = db->getNumStates();  //number of states in markov chain
	double T = 0.25;                            //final time
	double E = 15.0;                       //bond strengths
	double rho = RANGE;                   //range parameter to morse potential

	int initial = 0;
	int target = 90;

	//discretization parameters
	int N = 32;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;


	//init vectors for beta storage
	Eigen::VectorXd betaP; betaP.setZero(N); 

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, E, rho, db);

	//create output file
	std::ofstream ofile;
	ofile.open("constBetaFlower.txt");

	
	//determine the best constant
	int tests = 100;
	Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,3);
	for (int i = 0; i < tests; i++) {
		betaP.setOnes(N); betaP *= betaTest(i);
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
		ofile << betaTest(i) << " " << testP << "\n";
	}

	//close file
	ofile.close();
}

/****************************************************/
/*********** Model System Functions    **************/
/****************************************************/

void setBarrierHeights(Eigen::MatrixXd& h) {
	//set the hand-written barrier heights for the model problem

	//exponential format
	if (WHICH == 0) {
		h(0,1) = 2.0; h(0,2) = 1.0; h(0,3) = 2.0;
		h(1,4) = 3.0; h(1,6) = 3.0; h(1,0) = 5.0;
		h(2,5) = 1.0; h(2,8) = 3.0; h(2,0) = 5.0;
		h(3,8) = 3.0; h(3,9) = 3.0;  h(3,0) = 5.0;
		h(4,10) = 3.0; h(4,1) = 4.0;
		h(5,10) = 4.0; h(5,2) = 2.0;
		h(6,10) = 3.0; h(6,13) = 7.0; h(6,1) = 4.0;
		h(7,13) = 5.0; h(7,11) = 4.0; h(7,2) = 11.0; h(7,3) = 10.0;
		h(8,11) = 4.0; h(8,12) = 3.0; h(8,2) = 5.0; h(8,3) = 4.0;
		h(9,12) = 3.0; h(9,3) = 6.0;
		h(10,4) = 8.0; h(10,5) = 5.0; h(10,6) = 8.0;
		h(11,7) = 6.0; h(11,8) = 8.0;
		h(12,8) = 6.0; h(12,9) = 7.0;
		h(13,6) = 11.0; h(13,7) = 10.0;

		//remove these to force a reverse in the folding
		//h(2,7) = 10.0; h(3,7) = 9.0
	}

	//linear in beta format
	else if (WHICH == 1) {
		h(0,1) = 1.0; h(0,2) = 1.0; h(0,3) = 1.0;
		h(1,4) = 1.0; h(1,6) = 1.0; h(1,0) = 1.0;
		h(2,5) = 1.0; h(2,7) = 1.0; h(2,0) = 1.0;
		h(3,7) = 1.0; h(3,8) = 1.0; h(3,0) = 1.0;
		h(4,9) = 1.0; h(4,1) = 1.0;
		h(5,9) = 1.0; h(5,2) = 1.0;
		h(6,9) = 1.0; h(6,11) = 1.0; h(6,1) = 1.0;
		h(7,10) = 1.0; h(7,11) = 1.0; h(7,2) = 1.0; h(7,3) = 1.0;
		h(8,10) = 1.0; h(8,3) = 1.0;
		h(9,4) = 1.0; h(9,5) = 1.0; h(9,6) = 1.0;
		h(10,7) = 1.0; h(10,8) = 1.0; 
		h(11,6) = 0.1; h(11,7) = 0.1;
	}

}

void modelFlowerSystem() {
	//set up and solve for the optimal protocol in a simple model system intended to mimic 
	//the structue of the N=7 disk system

	//set problem parameters
	int num_states = 12;                  //number of states in markov chain
	int initial = 0;                      //starting state
	int target = 11;                      //final state

	//steepest descent parameters
	int max_iters = 2000;             //number of iterations of gradient descent
	double step = 0.95;               //constant step size
	double grad_tol = 1e-14;          //termination tolerance on gradient norm
	double input_tol = 1e-14;         //termintion tolerance on dx norm
	double lowerBound = 0.005;       //lower bound for beta, cant go negative

	//discretization parameters
	int N = 125;                            //number of intervals in termporal discretization
	long int maxT = 20;                        //the biggest Tf to sample
	
	std::vector<double> op;
	std::vector<double> cp;

	step = 0; //do linesearch

	int mult = 2; //multiply by this to sample logarithmically

	//create the barrier heights matrix
	Eigen::MatrixXd h; h.setZero(num_states,num_states);
	setBarrierHeights(h); 

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, h);

	//loop over final times, get the gap between the optimal protocol and optimal constant
	//#pragma omp parallel for
	for (long int T = maxT-1; T < maxT; T*=mult) {

		double dt = double(T) / double(N);    //size of time intervals

		Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
		t_disc.setLinSpaced(N+1,0,T);

		Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
		t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
		t_shift += t_disc;

		Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
		t_centers = t_shift.segment(0,N) / 2.0;

		//init vectors for beta storage and construct initial guess
		Eigen::VectorXd betaP; betaP.setZero(N); 
		Eigen::VectorXd prev_beta; prev_beta.setZero(N); 

		//determine the best constant
		int tests = 300;
		Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,5.5);
		double maxP = 0; double maxBeta = 0;
		for (int i = 0; i < tests; i++) {
			betaP.setOnes(N); betaP *= betaTest(i);
			double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
			std::cout << betaTest(i) << " " << testP << "\n";
			if (testP > maxP) {
				maxP = testP;
				maxBeta = betaTest(i);
			}
		}

		//determine the best protocol
		betaP.setOnes(N); betaP*= maxBeta;
		//betaP.setLinSpaced(N,0.5,3);
		performCG(N, max_iters, betaP, num_states, &problem, t_disc, t_centers, initial, target,
							step, lowerBound, grad_tol, input_tol, true);
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);

		//store optimal results
		op.push_back(testP); cp.push_back(maxP);
		std::cout << maxP << " " << maxBeta <<  "\n";
		printf("GD gives %f\n", testP);
		std::cout << betaP << "\n";
	}

	//create a file for output
	std::ofstream ofile;
	ofile.open("flowerModel.txt");

	//print the result
	int count = 0;
	for (long int T = maxT-1; T < maxT; T*=mult) {
		ofile << T << " " << cp[count] << " " << op[count] << "\n";
		count++;
	}

	//close outfile
	ofile.close();
	

}



/****************************************************/
/*********** Testing Functions    *******************/
/****************************************************/

void testProtocol(int np, Database* db, int initial, int target) {
	//test a given protocol

	//set problem parameters
	int num_states = db->getNumStates();  //number of states in markov chain
	int T = 3;                            //final time
	double E = 15.0;                       //bond strengths
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 32;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;


	//init vectors for beta storage and construct initial guess
	Eigen::VectorXd betaP; betaP.setZero(N); betaP.setOnes(N); betaP*=0.4;
	//betaP.setLinSpaced(N, 0.1,2);

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, E, rho, db);
	
	//get the probability of target state with this protocol
	double test_val = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);

	printf("The protocol being used is:\n");
	std::cout << betaP << "\n\n";
	printf("The protocol gives final prob %f\n", test_val);

}

void printProbabilities(int N, int num_states, OptInfo* problem,
												const Eigen::VectorXd& betaP, const Eigen::VectorXd& t_disc,  
												int initial, int target, std::ostream& ofile) {

	//compute the probability of being in each state at each time and output to file

	//create the rate matrix
	Eigen::MatrixXd rate; rate.setZero(num_states,num_states);

	//create the probability row vector
	Eigen::MatrixXd probs; probs.setZero(1,num_states);
	probs(0,initial) = 1.0;

	//print the number of states, the protocol, and initial condition to file
	ofile << num_states << "\n";
	for (int i = 0; i < N; i++) {
		ofile << betaP(i) << " ";
	}
	ofile << "\n";
	ofile << probs << "\n";

	for (int i = 0; i < N; i++) {
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, betaP(i));
		double delta_t = t_disc(i+1) - t_disc(i);
		rate *= delta_t;
		probs *= rate.exp();

		ofile << probs << "\n";
	}

}

/****************************************************/
/*********** dlib opt attempt    *******************/
/****************************************************/


double dlib_objective_call(int N, int num_states, OptInfo* problem, column_vector& betaD,
												 const Eigen::VectorXd& t_disc, const Eigen::VectorXd& t_centers,
												 int initial, int target, double lowerBound) {
	//call the dlib optimizer using an objective delta stop condition

	double m = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
             dlib::objective_delta_stop_strategy(1e-14,1000).be_verbose(), 

             [&](const column_vector& b) {
             	Eigen::VectorXd betaP; betaP.setZero(N); 
             	for (int i = 0; i < N; i++) betaP(i) = b(i);
             	return -computeFinalProb(N, num_states, problem, 
             	betaP, t_disc, initial, target);}, 

             [&](const column_vector& b) {
             	//compute gradient 
             	Eigen::VectorXd betaP; betaP.setZero(N); 
             	for (int i = 0; i < N; i++) betaP(i) = b(i);
             	Eigen::VectorXd	R; R.setZero(N);
             	evalGradient(N, betaP, num_states, problem, t_disc, t_centers, initial, target, R);
             	column_vector G(N); for (int i = 0; i < N; i++) G(i) = -R(i);
             	return G;
             }, 
             betaD, lowerBound, 5.0);

	return m;

}

void findProtocolSA_dlib(int np, Database* db, int initial, int target) {
	//determine a temperature protocol to maximize prob of target state

	//set problem parameters
	int num_states = db->getNumStates();  //number of states in markov chain
	double T = 10;                            //final time
	double E = 3.0;                       //bond strengths
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 32;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	double lowerBound = 1.0 / (2*E); //lower bound for beta (beta*E > 0.5)

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, E, rho, db);

	//init vectors for beta storage and construct initial guess
	Eigen::VectorXd betaP; betaP.setZero(N); 
	Eigen::VectorXd prev_beta; prev_beta.setZero(N); 

	//determine the best constant
	int tests = 100;
	Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,3);
	double maxP = 0; double maxBeta = 0;
	for (int i = 0; i < tests; i++) {
		betaP.setOnes(N); betaP *= betaTest(i);
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
		std::cout << betaTest(i) << " " << testP << "\n";
		if (testP > maxP) {
			maxP = testP;
			maxBeta = betaTest(i);
		}
	}

	//determine the best protocol
	betaP.setOnes(N); betaP*= maxBeta;
	//betaP *= 0.8;
	//betaP.setLinSpaced(N,0.9,1.4);

	//set up dlib storage types
	column_vector betaD(N);
	for (int i = 0; i < N; i++) betaD(i) = betaP(i);

	//perform optimization - try to escape local minima
	bool same = false; double obj = 0; double tol = 1e-6; 
	int attempt = 0; int max_attempts = 1;
	while (!same) { //if two consecutive ptimization dont give same obj value, perturb
		double m = dlib_objective_call(N, num_states, &problem, betaD, t_disc, t_centers, 
																	 initial, target,lowerBound);
		if (fabs(obj-m) < tol) { //escape
			same = true;
		}
		else { //perturb
			betaD /= 1.001;
			obj = m;
		}
		attempt++;
		if (attempt >= max_attempts) {
			break;
		}
	}

	//go from column vec to eigen vec
	for (int i = 0; i < N; i++) betaP(i) = betaD(i);
	
	//evaluate the final prob with this protocol
	double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);

	//store optimal results
	std::cout << maxP << " " << maxBeta <<  "\n";
	printf("GD gives %f\n", testP);
	std::cout << betaP << "\n";



}

void modelFlowerSystem_dlib() {
	//use dlib to do the optimization

	//set problem parameters
	int num_states = 12;                  //number of states in markov chain
	int initial = 0;                      //starting state
	int target = 11;                      //final state

	//steepest descent parameters
	int max_iters = 2000;             //number of iterations of gradient descent
	double step = 0.95;               //constant step size
	double grad_tol = 1e-14;          //termination tolerance on gradient norm
	double input_tol = 1e-14;         //termintion tolerance on dx norm
	double lowerBound = 0.005;       //lower bound for beta, cant go negative

	//discretization parameters
	int N = 250;                            //number of intervals in termporal discretization
	double T = 50;                          //final time
	
	step = 0; //do linesearch


	//create the barrier heights matrix
	Eigen::MatrixXd h; h.setZero(num_states,num_states);
	setBarrierHeights(h); 

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, h);

	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//init vectors for beta storage and construct initial guess
	Eigen::VectorXd betaP; betaP.setZero(N); 
	Eigen::VectorXd prev_beta; prev_beta.setZero(N); 

	//determine the best constant
	int tests = 300;
	Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,5.5);
	double maxP = 0; double maxBeta = 0;
	for (int i = 0; i < tests; i++) {
		betaP.setOnes(N); betaP *= betaTest(i);
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
		std::cout << betaTest(i) << " " << testP << "\n";
		if (testP > maxP) {
			maxP = testP;
			maxBeta = betaTest(i);
		}
	}

	//determine the best protocol
	betaP.setOnes(N); betaP*= maxBeta;
	//betaP *= 1.5;
	//betaP.setLinSpaced(N,0.9,1.4);

	//set up dlib storage types
	column_vector betaD(N);
	for (int i = 0; i < N; i++) betaD(i) = betaP(i);

	//perform optimization - try to escape local minima
	bool same = false; double obj = 0; double tol = 1e-6; 
	int attempt = 0; int max_attempts = 10;
	while (!same) { //if two consecutive ptimization dont give same obj value, perturb
		double m = dlib_objective_call(N, num_states, &problem, betaD, t_disc, t_centers, 
																	 initial, target,lowerBound);
		if (fabs(obj-m) < tol) { //escape
			same = true;
		}
		else { //perturb
			betaD /= 1.001;
			obj = m;
		}
		attempt++;
		if (attempt > max_attempts) {
			break;
		}
	}

	//go from column vec to eigen vec
	for (int i = 0; i < N; i++) betaP(i) = betaD(i);
	
	//evaluate the final prob with this protocol
	double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);

	//store optimal results
	std::cout << maxP << " " << maxBeta <<  "\n";
	printf("GD gives %f\n", testP);
	std::cout << betaP << "\n";

	//create a file for output
	std::ofstream ofile;
	ofile.open("saModel.txt");

	//output
	printProbabilities(N, num_states, &problem, betaP, t_disc, initial, target, ofile);

	//close outfile
	ofile.close();

}

void modelRowSystem_dlib() {
	//use dlib to do the optimization

	//set problem parameters
	int num_states = 50;                  //number of states in markov chain
	int initial = 0;                      //starting state
	int target = num_states-1;                      //final state

	//steepest descent parameters
	int max_iters = 2000;             //number of iterations of gradient descent
	double step = 0.95;               //constant step size
	double grad_tol = 1e-14;          //termination tolerance on gradient norm
	double input_tol = 1e-14;         //termintion tolerance on dx norm
	double lowerBound = 0.05;       //lower bound for beta, cant go negative

	//discretization parameters
	int N = 100;                            //number of intervals in termporal discretization
	long int maxT = 101;                     //the biggest Tf to sample
	
	std::vector<double> op;
	std::vector<double> cp;

	step = 0; //do linesearch

	int mult = 2; //multiply by this to sample logarithmically

	//pick the barrier heights
	double h1 = 2.0; double h2 = 4.0;

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, h1, h2);

	//loop over final times, get the gap between the optimal protocol and optimal constant
	//#pragma omp parallel for
	for (long int T = maxT-1; T < maxT; T*=mult) {

		double dt = double(T) / double(N);    //size of time intervals

		Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
		t_disc.setLinSpaced(N+1,0,T);

		Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
		t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
		t_shift += t_disc;

		Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
		t_centers = t_shift.segment(0,N) / 2.0;
		//t_centers = t_disc.segment(1,N+1);

		//init vectors for beta storage and construct initial guess
		Eigen::VectorXd betaP; betaP.setZero(N); 
		Eigen::VectorXd prev_beta; prev_beta.setZero(N); 

		//determine the best constant
		int tests = 2;
		Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,5.5);
		double maxP = 0; double maxBeta = 0;
		for (int i = 0; i < tests; i++) {
			betaP.setOnes(N); betaP *= betaTest(i);
			double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
			std::cout << betaTest(i) << " " << testP << "\n";
			if (testP > maxP) {
				maxP = testP;
				maxBeta = betaTest(i);
			}
		}

		//determine the best protocol
		betaP.setOnes(N); //betaP*= maxBeta;
		//betaP *= 1.5;
		//betaP.setLinSpaced(N,0,2*3.14*20);
		//betaP = betaP.sin()+1.005*Eigen::VectorXd::Ones(N);

		//set up dlib storage types
		column_vector betaD(N);
		for (int i = 0; i < N; i++) betaD(i) = betaP(i);

		//perform optimization - try to escape local minima
		bool same = false; double obj = 0; double tol = 1e-3; 
		int attempt = 0; int max_attempts = 10;
		while (!same) { //if two consecutive ptimization dont give same obj value, perturb
			double m = dlib_objective_call(N, num_states, &problem, betaD, t_disc, t_centers, 
																		 initial, target, lowerBound);
			if (fabs(obj-m) < tol) { //escape
				same = true;
			}
			else { //perturb
				//betaD *= 1.001;
				obj = m;
			}
			attempt++;
			if (attempt > max_attempts) {
				break;
			}
		}

		//go from column vec to eigen vec
		for (int i = 0; i < N; i++) betaP(i) = betaD(i);
		
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);

		//store optimal results
		op.push_back(testP); cp.push_back(maxP);
		std::cout << maxP << " " << maxBeta <<  "\n";
		printf("GD gives %f\n", testP);
		std::cout << betaP << "\n";
	}

	//create a file for output
	std::ofstream ofile;
	ofile.open("rowModel.txt");

	//print the result
	int count = 0;
	for (long int T = maxT-1; T < maxT; T*=mult) {
		ofile << T << " " << cp[count] << " " << op[count] << "\n";
		count++;
	}

	//close outfile
	ofile.close();

}

void testRowSystemProtocols() {

	//set problem parameters
	int num_states = 51;                  //number of states in markov chain
	int initial = 0;                      //starting state
	int target = num_states-1;                      //final state

	//discretization parameters
	int N = 200;                            //number of intervals in termporal discretization
	int T = 100;                     //the biggest Tf to sample
	

	//pick the barrier heights
	double h1 = 2.0; double h2 = 7.0;

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(N, num_states, initial, target, h1, h2);

	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//init vectors for beta storage and construct initial guess
	Eigen::VectorXd betaP; betaP.setZero(N); 
	double hot = 0.5; double cold = 1.5;
	for (int i = 0; i < N/4; i++) {
		betaP(4*i) = cold; betaP(4*i+1) = cold; 
		betaP(4*i+2) = hot; betaP(4*i+3) = hot; 
	}
	
	//test the current protocol
	double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
	printf("Probability is %f\n", testP);

	//get the best constant
	int tests = 100;
	Eigen::VectorXd betaTest; betaTest.setZero(N); betaTest.setLinSpaced(tests,0.005,2);
	double maxP = 0; double maxBeta = 0;
	for (int i = 0; i < tests; i++) {
		betaP.setOnes(N); betaP *= betaTest(i);
		double testP = computeFinalProb(N, num_states, &problem, betaP, t_disc, initial, target);
		std::cout << betaTest(i) << " " << testP << "\n";
		if (testP > maxP) {
			maxP = testP;
			maxBeta = betaTest(i);
		}
	}

}








}