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
#include "genetics.h"
#include "../defines.h"


struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

// now defined in genetics.h
/*
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
*/

namespace lattice {


/* Particle struct stores lattice coordinates and type */
struct Particle{
	//null constructor
	Particle() {x = y = type = 0;}
		
	//declare variables
	int type;
	int x; int y;

	//non-null constructor
	Particle(int x_, int y_, int type_) :  x(x_), y(y_), type(type_) {}
	Particle(int x_, int y_) :  x(x_), y(y_), type(0) {}
};

/* Database and state classes for proteins */
class Database;

class State{
	public:
		State();
		~State();
		//copy constructors
		State(const State& old) {
			copy(old);
		}
		State& operator=(const State& old) {
			if (this != &old) {
				destroy();
				copy(old);
			}
			return *this;
		}

		friend Database;
		friend Database* readData(std::string& filename);
		friend void buildPDB(int N);
		friend void addState(int N, Particle* chain, int* AM, std::vector<State>& new_states);

		//accessor functions
		double getFrequency() const {return freq;}
		int getBonds() const {return bond;}
		const std::vector<int> getCoordinates() const;
		bool isInteracting(int i, int j) const;
		double getMFPT() const {return mfpt;}
		double getSigma() const {return sigma;}
		int sumP() const;
		std::vector<bd::Pair> getP() const {return P;}

		//quantities to update
		int num_neighbors;
		double freq, mfpt, sigma;
		//Pair* P; Pair* Z; Pair* Zerr;
		std::vector<bd::Pair> P;

		//print function
		std::ostream& print(std::ostream&, int) const;



	private:
		int N;
		int bond;
		bool* am; 
		int* coordinates; 

		void destroy();
		void copy(const State& old);
};

class Database {
	public:
		Database(int N_, int num_states_); 
		~Database();
		State& operator[](int index) {return states[index];}
		const State& operator[](int index) const {return states[index];}
		friend std::ostream& operator<<(std::ostream&, const Database&);
		friend void buildPDB(int N);

		//accessor functions
		int getN() const {return N;}
		int getNumStates() const {return num_states;}

	private:	
		int N; int num_states; State* states; 

		//copy constructors - restricts compiling when user tries to copy a database
		Database(const Database&) {
			throw 1;
		}
		Database& operator=(const Database&) {
			throw 1;
		}

};

Database* readData(std::string& filename);


class Person2 : public ga::Person {
	public:
		Person2();
		Person2(int N_, int num_interactions_, int numTypes_, int* t, double* kV);
		Person2(const Person2& old);

		

		void evalStats(int N, Database* db, int initial, std::vector<int> targets, 
								 double* eq, double* Tconst, double* T, double* m);
		void evalStatsSampling(double Tf, double ts, int samples, int* M_target);
		Person2 mate(Person2 partner, bool, RandomNo*);
};




//allow unordered map to exist between pair and Particle*
typedef std::unordered_map<std::pair<int,int>,Particle*,pair_hash> particleMap; 



//functions to set up a chain
void getTypes(int N, int* types, bool useFile);
void initChain(int N, Particle* chain, particleMap& cMap, 
							 bool useFile);
int toIndex(int r, int c, int m);
void index2ij(int index, int N, int& i, int& j);





//functions to determine moves
void checkRotation(int x, int y, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
									 particleMap& cMap);

void getRotations(int particle, int neighbor, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap);

void checkCorner(int x, int y, std::vector<std::pair<int,int>>& moves,
									 particleMap& cMap);

void getCorners(int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap);

void getMoves(int N, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap);


//energy functions
void getBonds(int N, Particle* chain, std::vector<std::pair<int,int>>& bonds);

//mcmc functions
int randomInteger(int N, RandomNo* rngee);
void acceptMove(int particle, int x_old, int y_old, Particle* chain, 
								particleMap& cMap);

void rejectMove(int particle, int x_old, int y_old, Particle* chain);
bool takeStep(int N, Particle* chain, particleMap& cMap,
							RandomNo* rngee, double eps, double& energy);
void runMCMC(int N, bool useFile);
void testSampling(int N);
void testYield(int N);

//sampling functions
void buildPDB(int N);
void addState(int N, Particle* chain, int* AM, std::vector<State>& new_states);
void updatePDB(int N, Database* db);

void estimateMFPT(int N, int state, Database* db);
void estimateEqProbs(int N, Database* db);

//design functions
void constructScatterTOYL(int N, Database* db, int initial, int target, bool useFile);
void HPscatter(int N, Database* db, int initial, int target);
void evalStatsL(int N, Database* db, int initial, int target, bool useFile);

//genetic algorithm
void performGAevolution(int N, Database* db, int initial, int target, bool useFile);
void performGAlattice_sampling(int N, bool useFile);
void testMeasuresRecord(int N, Database* db, int initial, int target, bool useFile);
}