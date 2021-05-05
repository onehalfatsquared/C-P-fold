#pragma once
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <random>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include "pair.h"
#include "design.h"
#include "nauty.h"
#include "../defines.h"


class RandomNo{
	std::mt19937 generator; 
	std::uniform_real_distribution<double> uDist;
	std::normal_distribution<double> gDist;

	public:
		RandomNo(unsigned int seed = std::random_device{}())
        : generator{seed},uDist{0,1},gDist{0,1} {}

    double getU() {
    	return uDist(generator);
    }
    double getG() {
    	return gDist(generator);
    }
};


namespace ga {



class Person {
public:

	//con(de)structors
	Person();
	virtual ~Person();

	Person(int N, int num_interactions, int numTypes,  int* types, double* kappaVals);

	void copy(const Person& old);
	Person(const Person& old) {
		copy(old);
	}

	void destroy();

	Person& operator=(const Person& old) {
		if (this != &old) {
			destroy();
			copy(old);
		}
		return *this;
	}

	//parameters
	int num_interactions; int numTypes;
	int N;
	double Rate;
	double Eq;
	double fitness;
	int* types;
	double* kappaVals;
	std::map<std::pair<int,int>,double> kappa;

	//functions 
	void setKappa(int num_interactions, double* kappaVals);
	Person mate(Person partner, bool, RandomNo*);
	void evalStats(int N, bd::Database* db, int initial, std::vector<int> targets, 
								 double* eq, double* Tconst, double* T, double* m);
	void evalFitness(double eq, double rate);

	void evalStats(double Tf, int samples, int* M_target, double* X_target);
	void evalStats(double Tf, double ts, int samples, int* M_target);
	void applyBound(double lower, double upper);

private:


};

void readTargetFile(int* M_target, double* X_target);


double sampleKappa(int N, RandomNo* rngee);
int sampleType(int N, int numInteractions, RandomNo* rngee);
void sampleParameters(int N, int numInteractions, double* kappaVals, int* particleTypes, 
											int numTypes, bool useFile, RandomNo* rngee);

void non_dominated_set(int pop_size, double* xAll, double* yAll, std::vector<int>& nonDom);



void perform_evolution(int N, bd::Database* db, int initial, int target, bool useFile);
void perform_evolution_sampling(int N, bool useFile);

void printPopulation(std::vector<Person> population, int pop_size, std::ofstream& ofile );
void printTypes(std::vector<Person> population, int pop_size, int N);

double getMisfoldEnergy(int N, int* M, int* M_target, double* E);
double harmonicBarrier(int N,  int* M_target, double* E);























}