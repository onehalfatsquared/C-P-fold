#include "latticeP.h"
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/Eigen/Dense>
#include <omp.h>



namespace lattice {

/******************************************************************************/
/***************** Class functions for L proteins *****************************/
/******************************************************************************/


//state constructor
State::State() {
	am = NULL; coordinates = NULL; 

	freq = 0; bond = 0; 
	mfpt = 0; 
	sigma = 0; num_neighbors = 0;
	N = 0;
}

//state deconstructor
State::~State() {
	destroy();
}

void State::destroy() {
	delete []am; delete []coordinates; 
}

void State::copy(const State& old) {
	freq = old.freq;
	bond = old.bond;
	mfpt = old.mfpt;
	sigma = old.sigma;
	num_neighbors = old.num_neighbors;
	N = old.N;

	am = new bool[N*N];
	for(int i = 0; i < N*N; ++i) {
		am[i] = old.am[i];
	}

	coordinates = new int[DIMENSION*N];
	for (int i = 0; i < DIMENSION*N; ++i) {
		coordinates[i] = old.coordinates[i];
	}

}

//database constructor
Database::Database(int N_, int num_states_) {
	N = N_; num_states = num_states_;
	states = new State[num_states];
	for (int i = 0; i < num_states; i++) {
		states[i].N = N;
	}
}

//database deconstructor
Database::~Database() {
	delete []states;
}

//sum the entries of s.P
int State::sumP() const{
	int S = 0;
	for (int i = 0; i < num_neighbors; i++) {
		S += P[i].value;
	}
	return S;
}

bool State::isInteracting(int i, int j) const {
	if (abs(i-j) >= 2) {
		return (am[j*N+i] || am[i*N+j]);
	}
	else if (abs(i-j) == 1) {
		return true;
	}
	else {
		return false;
	}
}


//pull a random set of coordinates from the available
const std::vector<int> State::getCoordinates() const {
	std::vector<int> cv;
	for (int i = 0; i < DIMENSION * N; i++) {
		cv.push_back(coordinates[i]);
	}
	return cv;
}


//function to read in the database and store in database class
Database* readData(std::string& filename) {
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return NULL;
	}

	//read first line, N = number of particles
	int N;
	in_str >> N;

	//read second line - number of states
	int num_lines;
	in_str >> num_lines;


	bool val; //check if there is another line
	int coord; //number of sample coordinates for a state
	int index = 0; //loop over states
	int x, y; //coordinates on lattice
	char extra; //flag for whether mfpt estimates are in file
	int v1; double v2; //storing index and info in a pair

	//call the database class constructor
	Database* database = new Database(N, num_lines);

	//fill the database state classes
	while (in_str >> val) {
		//create reference to state, database[index]
		State& s = (*database)[index];

		//fill in adjacency matrix
		s.am = new bool[N*N];
		s.am[0] = val;
		for (int i = 1; i < N*N; i++) {
			in_str >> s.am[i];
		}

		//fill in frequency, bonds, and coords
		in_str >> s.freq;
		in_str >> s.bond;

		//fill in the sample coordinates
		s.coordinates = new int[DIMENSION*N];
		for (int j = 0; j < DIMENSION*N; j++) {
			in_str >> x;
			s.coordinates[j] = x;
		}

		//check the extra value for existence of mfpt estimates
		in_str >> extra;
		if (extra == 'N') {//no mfpt estimates, initialize to 0
			s.mfpt = 0;
			s.sigma = 0;
			s.num_neighbors = 0;
		}
		else if (extra == 'Y') {//mfpt estimates exist, read in
			in_str >> s.mfpt;
			in_str >> s.sigma; 
			in_str >> s.num_neighbors;

			for (int i = 0; i < s.num_neighbors; i ++) {
				in_str >> v1; in_str >> v2;
				s.P.push_back(bd::Pair(v1,v2));
			}
		}

		//next state
		index++; 
	}
	in_str.close();
	return database;
}

//write functions to output the updated database to a file
std::ostream& State::print(std::ostream& out_str, int N) const {
	for (int i = 0; i < N*N; i++) {
		out_str << am[i] << ' ';
	}
	out_str << freq << ' ';
	out_str << bond << ' ';
	for (int j = 0; j < DIMENSION*N; j++) {
		out_str << coordinates[j] << ' ';
	}
	out_str << "Y" << ' ';
	out_str << mfpt << ' ';
	out_str << sigma << ' ';
	out_str << num_neighbors << ' ';
	for (int i = 0; i < num_neighbors; i ++) {
		out_str << P[i].index << ' ' << P[i].value << ' ';
	}
	
	out_str << '\n';
	
}

std::ostream& operator<<(std::ostream& out_str, const Database& db) {
	//overload << to print out the database


	//print the states
	out_str << db.N << '\n';
	out_str << db.num_states  << '\n';
	for (int i = 0; i < db.num_states; i++) {
		db[i].print(out_str, db.N);
	}
}


/******************************************************************************/
/***************** Setup and General Functions *******************************/
/******************************************************************************/


void getTypes(int N, int* types, bool useFile) {
	//fill in the types of particles

	if (useFile) {
		//read it from the file

	}
	else{
		//set to all zeros
		for (int i = 0; i < N; i++) {
			types[i] = 0;
		}
	}
	
}

void initChain(int N, Particle* chain, particleMap& cMap,
							 bool useFile) {
	//initialize a linear chain, with some type distribution

	//get the types of each particle
	int* types = new int[N];
	getTypes(N, types, useFile);

	//initialize the chain - linear on x axis
	for (int i = 0; i < N; i++) {
		chain[i] = Particle(i, 0, types[i]);
		cMap[std::make_pair(i, 0)] = &(chain[i]); 
	}

	//free type memory
	delete []types;
}

void initChain(int N, Particle* chain, particleMap& cMap,
							 int* types) {
	//initialize a linear chain, with some type distribution

	//initialize the chain - linear on x axis
	for (int i = 0; i < N; i++) {
		chain[i] = Particle(i, 0, types[i]);
		cMap[std::make_pair(i, 0)] = &(chain[i]); 
	}
}

void initChain(int N, int* X, Particle* chain, particleMap& cMap,
							 bool useFile) {
	//initialize a chain with coordinates X, with some type distribution

	//get the types of each particle
	int* types = new int[N];
	getTypes(N, types, useFile);

	//initialize the chain - linear on x axis
	for (int i = 0; i < N; i++) {
		int x = X[2*i]; int y = X[2*i+1];
		chain[i] = Particle(x,y, types[i]);
		cMap[std::make_pair(x,y)] = &(chain[i]); 
	}

	//free type memory
	delete []types;
}

int toIndex(int r, int c, int m) {
  //map row and column number into index in 1d array. column indexed
  return m*c+r;
}

void index2ij(int index, int N, int& i, int& j) {
  // map a row index 1-d array into (i,j) 
  i = index % N;
  j = index / N;
}

void printChain(int N, Particle* chain) {
	//print the chain to terminal

	//declare array to put positions of particles in
	int* positions = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		positions[i] = 0;
	}

	//determine the minimum x and y coord in chain, gets mapped to 0
	int min_x = 100; int min_y = 100;
	for (int i = 0; i < N; i++) {
		int x = chain[i].x; int y = chain[i].y;
		if (x <= min_x) {
			min_x = x;
		}
		if (y <= min_y) {
			min_y = y;
		}
	}

	for (int i = 0; i < N; i++) {
		int x = chain[i].x - min_x;
		int y = chain[i].y - min_y;
		positions[toIndex(x,y,N)] = i+1;
		std::cout << x << ' ' << y << "\n";
	}

	std::cout << "\n";
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			int value = positions[toIndex(j,i,N)];
			if (N < 10) {
				std::cout << value << ' ';
			}
			else {
				if (value < 10) {
					std::cout << value << "  ";
				}
				else {
					std::cout << value << ' ';
				}
			}
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	//free memory
	delete []positions;
}

/******************************************************************************/
/***************** Monte Carlo Moves ******************************************/
/******************************************************************************/

void checkRotation(int x, int y, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
									 particleMap& cMap) {
	//check if rotating to (x,y) violates any rules, add to moves if doesnt

	//first check if (x,y) is in the cMap
	particleMap::iterator occupied = cMap.find(std::make_pair(x,y));
	if (occupied != cMap.end()) { //position found in map
		return;
	}

	//check for 180 degree rotation
	int xP = chain[particle].x; int yP = chain[particle].y;
	if (abs(x-xP) == 2 || abs(y-yP) == 2) {
		return;
	}

	//if we reach here, this rotation can be added to moves
	moves.push_back(std::make_pair(x,y));
}

void getRotations(int particle, int neighbor, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap) {
	//compute all valid rotations of an end particle

	//get the coordinates to rotate about
	int xN = chain[neighbor].x; int yN = chain[neighbor].y;

	//check all rotations
	checkRotation(xN-1, yN, particle, chain, moves, cMap);
	checkRotation(xN+1, yN, particle, chain, moves, cMap);
	checkRotation(xN, yN-1, particle, chain, moves, cMap);
	checkRotation(xN, yN+1, particle, chain, moves, cMap);

}

void checkCorner(int x, int y, std::vector<std::pair<int,int>>& moves,
									 particleMap& cMap) {
	//check if doing a corner move to (x,y) violates any rules

	//first check if (x,y) is in the cMap
	particleMap::iterator occupied = cMap.find(std::make_pair(x,y));
	if (occupied != cMap.end()) { //position found in map
		return;
	}

	//if we reach here, this move can be added to moves
	moves.push_back(std::make_pair(x,y));
}

void getCorners(int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap) {
	//compute all valid corner moves

	//get coordinates of the neigboring particles
	int x = chain[particle].x;    int y = chain[particle].y;
	int xL = chain[particle-1].x; int yL = chain[particle-1].y; 
	int xR = chain[particle+1].x; int yR = chain[particle+1].y; 

	//check if this particle is at a corner - add the appropriate move
	if (abs(xL-xR) == 1 && abs(yL-yR) == 1) {
		int slope = (xL-xR) / (yL-yR);
		int pos = yL + slope * (x - xL);
		//printf("Slope %d\n", slope);
		if (slope > 0) {
			if (y > pos) {
				checkCorner(x+1, y-1, moves, cMap);
			}
			if (y < pos) {
				checkCorner(x-1, y+1, moves, cMap);
			}
		}
		if (slope < 0) {
			if (y > pos) {
				checkCorner(x-1, y-1, moves, cMap);
			}
			if (y < pos) {
				checkCorner(x+1, y+1, moves, cMap);
			}
		}
	}

}

void getMoves(int N, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap) {
	//get all moves for a given particle

	//check for end particles
	if (particle == 0) {
		getRotations(particle, particle+1, chain, moves, cMap);
	}
	else if (particle == N-1) {
		getRotations(particle, particle-1, chain, moves, cMap);
	}
	else { //this is an interior particle
		getCorners(particle, chain, moves, cMap);
	}

}

/******************************************************************************/
/***************** Energy Functions ******************************************/
/******************************************************************************/

void getBonds(int N, Particle* chain, std::vector<std::pair<int,int>>& bonds) {
	//determine the non-trivial list of bonds

	for (int i = 0; i < N; i++) {
		int xi = chain[i].x; int yi = chain[i].y;

		for (int j = i+2; j < N; j++) {
			int xj = chain[j].x; int yj = chain[j].y;

			//check if particles are distance 1 apart. 1-norm distance
			int dist = abs(xj-xi) + abs(yj-yi);
			if (dist == 1) {
				bonds.push_back(std::make_pair(i,j));
			}
		}
	}

}

double getEnergy(int N, Particle* chain, double eps) {
	//get energy, assuming uniform eps per bond

	std::vector<std::pair<int,int>> bonds;
	getBonds(N, chain, bonds);
	double e1 = -eps * bonds.size();

	return e1;
}

double getEnergy(int N, Particle* chain, double* E) {
	//get energy for non-uniform bonds strengths

	std::vector<std::pair<int,int>> bonds;
	getBonds(N, chain, bonds);

	double e = 0;
	for (int i = 0; i < bonds.size(); i++) {
		e += - E[toIndex(bonds[i].first, bonds[i].second, N)];
	}

	return e;

}


/******************************************************************************/
/***************** Monte Carlo Functions **************************************/
/******************************************************************************/

int randomInteger(int N, RandomNo* rngee) {
	//pick an integer uniformly from 0 to N-1

	double U = N * rngee->getU();
	return floor(U);

}

void acceptMove(int particle, int x_old, int y_old, Particle* chain, 
								particleMap& cMap) {
	//replace the old cMap entry with the new one
	cMap.erase(std::make_pair(x_old,y_old));

	int x_new = chain[particle].x; int y_new = chain[particle].y;
	cMap[std::make_pair(x_new, y_new)] = &(chain[particle]);

}

void rejectMove(int particle, int x_old, int y_old, Particle* chain) {
	//put the old coordinates back into the correct chain entry

	chain[particle].x = x_old; chain[particle].y = y_old;

}


bool takeStep(int N, Particle* chain, particleMap& cMap,
							RandomNo* rngee, double eps, double& energy) {
	//perform an MCMC step - return the energy of the returned state

	//std::cout << "hello\n";
	//set the initial energy
	double e0 = energy;

	//first we pick a random particle and get its coordinates
	int particle = randomInteger(N, rngee);
	int x_old = chain[particle].x; int y_old = chain[particle].y;
	//printf("Particle %d\n", particle);

	//std::cout << "hello\n";

	//next we generate the set of moves
	std::vector<std::pair<int,int>> moves;
	getMoves(N, particle, chain, moves, cMap);

	
	//std::cout << "hello\n";

	//pick a move to perform, or return if there are no choices
	int M = moves.size();
	if (M == 0) {
		return false;
	}
	int move = randomInteger(M, rngee);
	//printf("Picked move %d of %d\n", move+1, M);


	//update the position of particle in chain to whats in move
	chain[particle].x = moves[move].first; chain[particle].y = moves[move].second; 

	//get energy from chain
	double e1 = getEnergy(N, chain, eps);

	//get acc probability - do accept/reject step
	double a = std::min(1.0, exp(-(e1-e0)));
	//printf("e1 = %f, e0 = %f, Acc = %f\n", e1, e0, a);
	double U = rngee->getU();
	if (U <= a) {
		acceptMove(particle, x_old, y_old, chain, cMap);
		energy = e1;
		//printf("Move accept %f\n", e1);
		//printChain(N, chain);
		return true;
	}
	else {
		rejectMove(particle, x_old, y_old, chain);
		//printf("Move reject\n");
		return false;
	}


}

bool takeStep(int N, Particle* chain, particleMap& cMap,
							RandomNo* rngee, double* E, double& energy) {
	//perform an MCMC step - return the energy of the returned state

	//std::cout << "hello\n";
	//set the initial energy
	double e0 = energy;

	//first we pick a random particle and get its coordinates
	int particle = randomInteger(N, rngee);
	int x_old = chain[particle].x; int y_old = chain[particle].y;
	//printf("Particle %d\n", particle);

	//std::cout << "hello\n";

	//next we generate the set of moves
	std::vector<std::pair<int,int>> moves;
	getMoves(N, particle, chain, moves, cMap);

	
	//std::cout << "hello\n";

	//pick a move to perform, or return if there are no choices
	int M = moves.size();
	if (M == 0) {
		return false;
	}
	int move = randomInteger(M, rngee);
	//printf("Picked move %d of %d\n", move+1, M);


	//update the position of particle in chain to whats in move
	chain[particle].x = moves[move].first; chain[particle].y = moves[move].second; 

	//get energy from chain
	double e1 = getEnergy(N, chain, E);

	//get acc probability - do accept/reject step
	double a = std::min(1.0, exp(-(e1-e0)));
	//printf("e1 = %f, e0 = %f, Acc = %f\n", e1, e0, a);
	double U = rngee->getU();
	if (U <= a) {
		acceptMove(particle, x_old, y_old, chain, cMap);
		energy = e1;
		//printf("Move accept %f\n", e1);
		//printChain(N, chain);
		return true;
	}
	else {
		rejectMove(particle, x_old, y_old, chain);
		//printf("Move reject\n");
		return false;
	}


}

void runMCMC(int N, bool useFile) {
	//run the monte carlo simulation

	//construct the chain of particles and the lattice mapping
	Particle* chain = new Particle[N];
	particleMap cMap;

	//initialize as linear chain
	initChain(N, chain, cMap, useFile);

	//set parameters
	double energy = 0;
	double eps = 2.0;
	double max_it = 1000;

	//initialize rng
	RandomNo* rngee = new RandomNo();

	//do the monte carlo steps
	for (int i = 0; i < max_it; i++) {
		takeStep(N, chain, cMap, rngee, eps, energy);
		//printf("Iteration: %d, energy %f\n", i, energy);
		std::cout << "Iteration " << i << "\n"; 
	}


	//delete memory
	delete []rngee; delete []chain;
}

/******************************************************************************/
/******************** Sampling Functions **************************************/
/******************************************************************************/

void buildPDB(int N) {
	//construct a new database file that only has the linear state in it

	//make the db with only 1 state - linear chain
	Database* db = new Database(N, 1);

	//create reference to that state
	State& s = (*db)[0];

	//fill in the state info

	//adjacency matrix
	s.am = new bool[N*N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j == i+1 || j == i-1) {
				s.am[toIndex(i,j,N)] = 1;
			}
			else {
				s.am[toIndex(i,j,N)] = 0;
			}
		}
	}

	//bonds
	s.bond = N-1;

	//coordinates
	s.coordinates = new int[DIMENSION*N];
	for (int i = 0; i < N; i++) {
		s.coordinates[2*i] = i; s.coordinates[2*i+1] = 0; 
	}

	//print out the new database
	std::string out = "N" + std::to_string(N) + "DB.txt";
	std::ofstream out_str(out);
	out_str << *db;
	delete db;
}

void getAM(int N, Particle* chain, int* AM) {
	//determine the non-trivial list of bonds

	for (int i = 0; i < N-1; i++) {
		int xi = chain[i].x; int yi = chain[i].y;
		AM[toIndex(i,i+1,N)] = 1; AM[toIndex(i+1,i,N)] = 1;

		for (int j = i+2; j < N; j++) {
			int xj = chain[j].x; int yj = chain[j].y;

			//check if particles are distance 1 apart. 1-norm distance
			int dist = abs(xj-xi) + abs(yj-yi);
			//printf("Particles %d and %d, dist %d\n", i,j,dist);
			if (dist == 1) {
				AM[toIndex(i,j,N)] = 1; AM[toIndex(j,i,N)] = 1;
			}
			else {
				AM[toIndex(i,j,N)] = 0; AM[toIndex(j,i,N)] = 0;
			}
		}
	}
}

void printAM(int N, int* AM) {
	//print out the adj matrix 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%d ", AM[toIndex(i,j,N)]);
		}
		printf("\n");
	}
}

bool checkSame(int N, int* AM, State& s) {
	//check if states are same by adjacency matrix

	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if (s.isInteracting(i,j) != AM[toIndex(i,j,N)]) {
				return false;
			}
		}
	}

	return true;
}

int searchDB(int N, Database* db, std::vector<State> new_states, int* AM) {
	//check if the current adj matrix is in db. if yes, return state #. if not, return -1.

	int num_states = db->getNumStates();

	for (int state = 0; state < num_states; state++) {
		State& s = (*db)[state];
		bool same = checkSame(N, AM, s);
		if (same)
			return state;
	}

	//std::cout << new_states.size() << "\n\n\n\n\n";
	for (int state = 0; state < new_states.size(); state++) {
		State& s = new_states[state];
		bool same = checkSame(N, AM, s);
		if (same)
			return state+num_states;
	}

	return -1;
}

void addState(int N, Particle* chain, int* AM, std::vector<State>& new_states) {
	//add a new state to a vector

	State s = State();

	int b = 0;
	s.N = N;

	s.am = new bool[N*N];
	//construct AM
	for (int i = 0; i < N*N; i++) {
		s.am[i] = AM[i];
		if (AM[i] == 1) {
			b++;
			//std::cout <<"hello" << s.isInteracting(i / N, i % N) << "\n";
		}
	}

	//add num bonds
	b += N-1;
	s.bond = b;

	//add configuration
	s.coordinates = new int[DIMENSION*N];
	for (int i = 0; i < N; i++) {
		s.coordinates[2*i] = chain[i].x; s.coordinates[2*i+1] = chain[i].y;
	}

	//push to vector
	new_states.push_back(s);
}

void updatePDB(int N, Database* db) {
	//updates a database with new states

	//construct the chain of particles and the lattice mapping
	Particle* chain = new Particle[N];
	particleMap cMap;

	//initialize as linear chain
	initChain(N, chain, cMap, false);

	//set parameters
	double energy = 0;
	int max_it = 1e6;
	double eps = 1.0;

	//initialize rng
	RandomNo* rngee = new RandomNo();

	//set up an adjacency matrix
	int* AM = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		AM[i] = 0;
	}

	//create vector of new states
	std::vector<State> new_states;
	int count = 1;

	//do the monte carlo steps
	for (int i = 0; i < max_it; i++) {
		//take a step
		bool accept = takeStep(N, chain, cMap, rngee, eps, energy);
		//std::cout << accept << "\n";

		//check if the state changed from previous step
		if (accept) {
			//get adjacency matrix for current state
			getAM(N, chain, AM);
			//printAM(N, AM);

			//check if this state has been seen before
			int state = searchDB(N, db, new_states, AM);
			//std::cout << state << "\n";
			if (state == -1) {
				//add the new state to vector
				addState(N, chain, AM, new_states);
				printf("Found new state. Total found this run: %d\n", count);
				printChain(N, chain);
				count++;
			}
		}

		if (i % (max_it / 100) == 0) {
			printf("Finished step %d of %d\n", i, max_it);
		}
	}

	//construct a new database
	int num_states_old = db->getNumStates();
	int num_states_new = new_states.size();
	int total = num_states_old + num_states_new;

	Database* newDB = new Database(N, total);

	//copy the first states from the old db
	for (int i = 0; i < num_states_old; i++) {
		(*newDB)[i] = (*db)[i];
	}

	//copy the rest from the vector
	for (int i = 0; i < num_states_new; i++) {
		(*newDB)[i+num_states_old] = new_states[i];
	}

	//print the new db
	std::string out = "N" + std::to_string(N) + "DBupdate.txt";
	std::ofstream out_str(out);
	out_str << *newDB;

	//delete memory
	delete rngee; 
	delete []chain; 
	delete []AM;
	delete newDB;

}

void sampleStats(double* X, int N, double& M, double& V) {
	//return the sample mean of the array X with N elements

	//init the mean and variance at 0, get sample size
	M = 0; V = 0;

	//compute the mean
	for (int i = 0; i < N; i++) {
		M += X[i];
	}
	M /= float(N);

	//compute the variance 
	for (int i = 0; i < N; i++) {
		V += (X[i]-M) * (X[i]-M);
	}
	V /= float(N-1.0);

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

void minVarEstimate(int sampleSize, double* means, double* variances, double& M, double& V) {
	//combine the estimated means in a linear combination to minimize variance
	//formula: M = 1/S sum(V_jm_j), S = 1/(sum(1/V_j))

	//init the output mean and variance, and normalizer
	M = 0; V = 0;
	double S = 0;

	//compute the normalizing term, S
	for (int i = 0; i < sampleSize; i++) {
		S += 1.0/variances[i];
	}

	//compute the minimizing variance and the corresponding mean
	for (int i = 0; i < sampleSize; i++) {
		M += means[i]/variances[i];
		V += 1.0/variances[i];
	}
	M /= S;
	V /= S*S;
}


void updatePM(int new_state, std::vector<bd::Pair>& PM) {
	//find pair with index = state and increment value by 1

	int i;

	for (i = 0; i < PM.size(); i++) {
		if (PM[i].index == new_state) {
			PM[i].value += 1;
			break;
		}
	}

	if (i == PM.size()) {//this state is being hit for the first time
		PM.push_back(bd::Pair(new_state,1));
	}

}

void equilibrate(Particle* chain, particleMap& cMap, Database* db, int state, int N,
								 RandomNo* rngee) {
	//get a sample of the mean first passage time, record the state that gets visited

	//parameters for the estimator and bond checking
	int max_it = 500;             //cut off if not done after max_it samples
	int new_state = state;            //new state id
	bool accepted;                    //flag to check if MCMC accepted proposal
	int b = (*db)[state].getBonds();  //num bonds in starting cluster
	double timer = 0;                //init the timer for the mfpt
	bool reset = false;              //if this becomes true, reset with no sample
	double eps = 1000;               //energy per bond, so chain does not break
	double energy0 = getEnergy(N, chain, eps);  //energy for initial state
	double energy;                   //dummy energy

	//make another chain and cMap for reversions
	Particle* prev_chain = new Particle[N];
	particleMap prevMap;
	for (int i = 0;  i < N; i++) {
		prev_chain[i] = chain[i];
	}
	prevMap = cMap;

	//time per step is constant, set it
	double dt = 1;

	//set up an adjacency matrix
	int* AM = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		AM[i] = 0;
	}

	//generate samples using MCMC until max_its or num samples is reached
	for (int i = 0; i < max_it; i++) {
		//get the sample
		energy = energy0;
		accepted = takeStep(N, chain, cMap, rngee, eps, energy);

		//if the position has changed, update timer and check for bond formation
		if (accepted) {

			//increment the timer
			timer += 1;

			//check if state changed
			//get adjacency matrix for current state
			getAM(N, chain, AM);

			//search database for current state
			std::vector<State> useless;
			new_state = searchDB(N, db, useless, AM);

			// if reset is true, this sample is invalid. reset clock and config
			if (reset) { 
				timer = 0; 
				for (int i = 0;  i < N; i++) {
					chain[i] = prev_chain[i];
				}
				cMap = prevMap;
			}

			//if the new_state is different from state update estimates
			if (state != new_state) { 
				//reset timer and reflect the state
				timer = 0;  new_state = state;
				for (int i = 0;  i < N; i++) {
					chain[i] = prev_chain[i];
				}
				cMap = prevMap;
			}

			//update the previous config
			for (int i = 0;  i < N; i++) {
				prev_chain[i] = chain[i];
			}
			prevMap = cMap;
		}
	}

	//free memory 
	delete []prev_chain; delete []AM;
}

void getSamplesMFPT(Particle* chain, particleMap& cMap, Database* db, int state, int N,
	std::vector<bd::Pair>& PM, std::vector<double>& mfptVec, RandomNo* rngee) {
	//get a sample of the mean first passage time, record the state that gets visited

	//parameters for the estimator and bond checking
	int max_it = 1e6;             //cut off if not done after max_it samples
	int new_state = state;            //new state id
	bool accepted;                    //flag to check if MCMC accepted proposal
	int b = (*db)[state].getBonds();  //num bonds in starting cluster
	double timer = 0;                //init the timer for the mfpt
	bool reset = false;              //if this becomes true, reset with no sample
	int hits = 0;                    //counter for number of samples
	double eps = 1000;               //energy per bond, so chain does not break
	double energy0 = getEnergy(N, chain, eps);  //energy for initial state
	double energy;                   //dummy energy
	int samples = 50000;              //number of samples per walker

	//make another chain and cMap for reversions
	Particle* prev_chain = new Particle[N];
	particleMap prevMap;
	for (int i = 0;  i < N; i++) {
		prev_chain[i] = chain[i];
	}
	prevMap = cMap;

	//time per step is constant, set it
	double dt = 1;

	//set up an adjacency matrix
	int* AM = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		AM[i] = 0;
	}

	//generate samples using MCMC until max_its or num samples is reached
	for (int i = 0; i < max_it; i++) {
		//get the sample
		energy = energy0;
		accepted = takeStep(N, chain, cMap, rngee, eps, energy);
		//increment the timer
		timer += 1;

		//if the position has changed, update timer and check for bond formation
		if (accepted) {

			//timer += 1;
			//check if state changed
			//get adjacency matrix for current state
			getAM(N, chain, AM);

			//search database for current state
			std::vector<State> useless;
			new_state = searchDB(N, db, useless, AM);

			// if reset is true, this sample is invalid. reset clock and config
			if (reset) { 
				timer = 0; 
				for (int i = 0;  i < N; i++) {
					chain[i] = prev_chain[i];
				}
				cMap = prevMap;
			}

			//if the new_state is different from state update estimates
			if (state != new_state) { 
				//record which state is hit
				updatePM(new_state, PM); hits++;

				//get an mfpt estimate, add to vector
				double tau = (timer+1.0)/2.0;
				mfptVec.push_back(tau*dt);
				//printf("%f\n", tau*dt);

				if (state == 4 && new_state == 0) {
					printChain(N, prev_chain);
					printChain(N, chain);
					std::vector<std::pair<int,int>> bonds; std::vector<std::pair<int,int>> bonds_prev;
					getBonds(N, chain, bonds); getBonds(N, prev_chain, bonds_prev);
					printf("Previous Bonds = %lu, Current bonds = %lu\n", bonds_prev.size(), bonds.size());
					abort();
				}

				//reset timer and reflect the state
				timer = 0;  new_state = state;
				for (int i = 0;  i < N; i++) {
					chain[i] = prev_chain[i];
				}
				cMap = prevMap;
			}

			//check if num_samples has been reached
			if (hits == samples) {
				break;
			}

			//update the previous config
			for (int i = 0;  i < N; i++) {
				prev_chain[i] = chain[i];
			}
			prevMap = cMap;
		}
	}

	//free memory 
	delete []prev_chain; delete []AM;
}

void estimateMFPT(int N, int state, Database* db) {
	/*estimate mean first passage time starting in state and going to state with
	one additional bond. Uses parallel implementations of a single walker with
	long trajectory.*/

	//set parameters
	int num_states = db->getNumStates(); //total number of states

	//check if this state has max number of bonds
	int maxB = 0;
	for (int i = 0; i < num_states; i++) {
		int b = (*db)[i].getBonds();
		if (b > maxB) {
			maxB = b;
		}
	}
	if ((*db)[state].getBonds() == maxB) {
		printf("This state is a ground state. No MPFT estimation necessary.\n");
		return;
	}

	//quantities to estimate
	std::vector<bd::Pair> PM; std::vector<bd::Pair> PMshare;
	double mfpt = 0;
	double sigma2 = 0;

	//store mfpt estimates on each thread to get standard deviation
	double* mfptSamples; double* mfptVar; int num_threads;

	//debug line
	for (int i = 0; i < PM.size(); i++) printf("0 thread has:\n %d, %f\n", PM[i].index, PM[i].value);

	//open parallel region
	#pragma omp parallel private(PM) shared(PMshare) 
	{
		//initialize final samples storage - only on one processor - then barrier
		num_threads = omp_get_num_threads();
		#pragma omp single
		{
			mfptSamples = new double[num_threads];
			mfptVar     = new double[num_threads];
			for (int i = 0; i < num_threads; i++) {
				mfptSamples[i] = 0; mfptVar[i] = 0;
			}
		}
		#pragma omp barrier

		//init the private mfpt sample storage
		std::vector<double> mfptVec;

		//init the random number generator 
		RandomNo* rngee = new RandomNo(); 
		//printf("Thread: %d, number %f\n", omp_get_thread_num(), rngee->getU());

		//get starting coordinates randomly from the database
		const std::vector<int> c = (*db)[state].getCoordinates();

		//make array with IC
		int* X = new int[DIMENSION*N];
		for (int i = 0; i < DIMENSION*N; i++) {
			X[i] = c[i];
		}

		//construct the chain of particles and the lattice mapping
		Particle* chain = new Particle[N];
		particleMap cMap;

		//initialize as linear chain
		initChain(N, X, chain, cMap, false);

		//equilibrate the trajectories
		equilibrate(chain, cMap, db, state, N, rngee);

		//get samples - has to update PM
		getSamplesMFPT(chain, cMap, db, state, N, PM, mfptVec, rngee);

		//get sample means and variances
		double M; double V;
		sampleStats(mfptVec, M, V);
		if (mfptVec.size() == 0) {
			M = 0; V = 0;
		}
		mfptSamples[omp_get_thread_num()] = M;
		mfptVar[omp_get_thread_num()] = V;

		//do update on PM vectors - need barrier - one at a time
		if (omp_get_thread_num() == 0) {
			PMshare = PM;
		}
		#pragma omp barrier
		#pragma omp critical
		{
			if (omp_get_thread_num() != 0) {
				combinePairs(PMshare, PM);
			}
		}

		//free memory
		delete []X; delete rngee; delete []chain;

		//end parallel region
	}

	//update estimates
	//combine the mfptSamples to get mean and variance
	sampleStats(mfptSamples, num_threads, mfpt, sigma2);
	double sigma = sqrt(sigma2);

	//update database
	(*db)[state].mfpt = mfpt;
	(*db)[state].num_neighbors = PMshare.size();
	(*db)[state].P = PMshare;
	(*db)[state].sigma = sigma;

	//print out final estimates - debug
	///*
	double sum = 0;
	for (int i = 0; i < PMshare.size(); i++) {
		printf("State = %d, visits = %f\n", PMshare[i].index, PMshare[i].value);
		sum +=PMshare[i].value;
	}
	printf("sum of hits = %f\n", sum);
	printf("Total Estimate = %f +- %f\n", mfpt, sigma);
	for (int i = 0; i < num_threads; i++) printf("MFPT estimate %d = %f +- %f\n", i, mfptSamples[i], sqrt(mfptVar[i]));
	//*/


	//free memory
	delete []mfptSamples; delete []mfptVar;
}


void estimateEqProbs(int N, Database* db) { 
	//estimate the equilibrium probabilities for each state
	//use MCMC estimator, use every c moves
	//seperate across threads to get better estimate

	//set parameters
	int num_states = db->getNumStates(); //total number of states
	int max_its = 5e7;                   //max number of MCMC steps to take
	int eq_its = 500;                    //number of steps to equilibrate for
	int cut = 5;                         //keep a sample every cut iterations
	double eps = EPS;

	//store freq estimates on each thread to get standard deviation
	double* eqShared = new double[num_states];
	double* eqVar = new double[num_states];
	for (int i = 0; i < num_states; i++) {
		eqShared[i] = eqVar[i] = 0;
	}

	int num_threads;
	double* X;                          //store estimates from each thread

	//open parallel region
	#pragma omp parallel shared(eqShared, eqVar, X) 
	{
		//init the private freq storage
		int* freqPrivate = new int[num_states]; 
		double* eqPrivate = new double[num_states];
		for (int i = 0; i < num_states; i++) {
			freqPrivate[i] = 0; eqPrivate[i] = 0.0;
		}

		//get number of threads - init estimate storage
		num_threads = omp_get_num_threads();
		X = new double[num_threads];

		//init the random number generator 
		RandomNo* rngee = new RandomNo(); 

		//init the adjacency matrix
		int* AM = new int[N*N]; for (int i = 0; i < N*N; i++) AM[i] = 0;

		//construct the chain of particles and the lattice mapping
		Particle* chain = new Particle[N];
		particleMap cMap;

		//initialize as linear chain
		initChain(N, chain, cMap, false);
		double energy = 0;

		//equilibrate the trajectories
		for (int i = 0; i < eq_its; i++) {
			takeStep(N, chain, cMap, rngee, eps, energy);
		}
		#pragma omp barrier

		//do MCMC, store every cut iterations
		for (int i = 0; i < max_its; i++) {
			//do mcmc step
			takeStep(N, chain, cMap, rngee, eps, energy);

			//check if we store samples this iteration
			if (i % cut == 0) {
				//check which state we are in
				getAM(N, chain, AM);

				//search database for current state
				std::vector<State> useless;
				int new_state = searchDB(N, db, useless, AM);
				//printf("New state is %d of %d\n", new_state, num_states);
				freqPrivate[new_state] += 1;
			}
		}
		#pragma omp barrier
	
		//normalize each of the freq arrays into the private eq array
		int Z = 1 + ((max_its-1) / cut);
		for (int i = 0; i < num_states; i++) {
			int count = freqPrivate[i];
			eqPrivate[i] = double(count) / double(Z);
		}

		//barrier so all threads finish constructing eqPrivate
		#pragma omp barrier
		
		double M, V;
		//for each state, do a critical step to get mean and variance over threads
		for (int i = 0; i < num_states; i++) {
			X[omp_get_thread_num()] = eqPrivate[i];
			//printf("Thread: %d, i: %d, value %f\n", omp_get_thread_num(), i, eqPrivate[i]);
			#pragma omp barrier
			#pragma omp single
			{
				sampleStats(X, num_threads, M, V);
				eqShared[i] = M;
				eqVar[i] = V;
			}
		}

		//free memory
		delete rngee; delete []chain; delete []AM; delete []X; 
		delete []freqPrivate; delete []eqPrivate;

		//end parallel region
	}

	//update estimates in db
	for (int i = 0; i < num_states; i++) {
		double freqU = eqShared[i];
		(*db)[i].freq = freqU;
	}
	

	//print out final estimates - test & debug
	for (int i = 0; i < num_states; i++) {
		printf("State %d: Eq = %f, Std = %f\n", i, eqShared[i], sqrt(eqVar[i]));
	}
	
	
	//free memory
	delete []eqShared; delete []eqVar;
}


/******************************************************************************/
/***************** Design Functions  ******************************************/
/******************************************************************************/

double getStickyProductL(int N, int state,  Database* db, int* particleTypes, 
												const std::map<std::pair<int,int>,double>& kappa) {
	//return the product of sticky param ^ num bonds for each type

	double stickyProd = 1.0;

	//loop over the adjacency matrix, determine bond types, get factors
	//todo - the isInteracting might fail here b/c AM is only upper triangular?
	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ((*db)[state].isInteracting(i,j)) {
				int p1 = particleTypes[i]; int p2 = particleTypes[j];
				stickyProd *= kappa.find({p1,p2})->second;
				//printf("state: %d, i %d, j %d\n",state, i, j);
			}
		}
	}

	return stickyProd;
}

void reweightL(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  const std::map<std::pair<int,int>,double>& kappa) {
	//perform the re-weighting of the eq measure for new sticky params

	double kap0 = exp(EPS);            //sticky parameter for initial measurement

	double Z = 0;                     //normalizing constant for eq

	//loop over states, get the equilibrium measure entry
	for (int i = 0; i < num_states; i++) {
		//get initial eq prob and number of bonds
		int b = (*db)[i].getBonds();
		int ntBonds = b - (N-1);               //number of non-trivial bonds
		double prob = (*db)[i].getFrequency();

		//get the denominator of reweight factor - KAP^b_i
		double denom = pow(kap0, ntBonds);

		//get the numerator of the re-weight factor
		double num = getStickyProductL(N, i, db, particleTypes, kappa); 

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
	
void checkPositiveL(int numInteractions, double* kappaVals) {
	//check if each kappa is positive, if not set to 0.01

	for (int i = 0; i < numInteractions; i++) {
		if (kappaVals[i] < 1.0) {
			kappaVals[i] = 1.001;
		}
	}
}

void createTransitionMatrix(double* T, int num_states, Database* db, 
														std::vector<int>& endStates) {
	//create rate matrix from data in DB - forward rates

	//declare storage for each state
	double mfpt; double S;
	double mult = 300;



	//loop over all states to get max bonds number
	int maxB = 0;
	for (int i = 0; i < num_states; i++) {
		int b = (*db)[i].getBonds();
		if (b > maxB) {
			maxB = b;
		}
	}

	//loop over all states information
	for (int state = 0; state < num_states; state++) {
		//get all relevant data
		mfpt = (*db)[state].getMFPT(); 
		S = (*db)[state].sumP(); //get normalizing constant for this row

		//if S = 0, this is end state, add it to vector
		if ((*db)[state].getBonds() == maxB) {
			endStates.push_back(state);
		}

		//fill in value in transition matrix
		std::vector<bd::Pair> P = (*db)[state].getP();
		for (int i = 0; i < P.size(); i++) {
			T[toIndex(state, P[i].index, num_states)] = mult * (P[i].value / S) / mfpt;
		}
	}
}

void buildNautyGraph(int N, int M, int state, Database* db, graph* g) {
	//build nauty graph of given state

	//zero out any pre-existing graph
	EMPTYGRAPH(g, M, N);

	//add edges for every non-zero entry in adjacnecy matrix
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			//std::cout << i << ' ' << j << "\n";
			if ((*db)[state].isInteracting(i,j)) {
				ADDONEEDGE(g, i, j, M);
			}
		}
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
		if (bd::checkIsomorphic(N, M, cg1, cg2) == 1) { //graphs are isomorphic, add to iso
			iso.push_back(i);
		}
	} 
}

double getEqProb(double* eq, std::vector<int> targets) {
	//sum the equilibrium probabilities of the target states
	//if cond == true, return eq prob conditional on being in a ground state

	//get the new equilibrium probability of target state
	double prob = 0;
	for (int i = 0;  i< targets.size(); i++) {
		prob += eq[targets[i]];
	}

	return prob;
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

void constructScatterTOYL(int N, Database* db, int initial, int target, bool useFile) {
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
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = bd::setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	bd::initKappaVals(numInteractions, kappaVals);

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
	std::vector<int> targets; targets.push_back(target);

	//the following doesnt work for distinguishing lattice ground states
	/*
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	*/

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
	M = 33;
	double base = 1.0; //smallest kappa to use - 1
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
		for (int y = 0; y < 1; y++) {
			for (int z = 0; z < 1; z++) {
				//set kappa
				kappaVals[0] = Ks[0];
				kappaVals[1] = Ks[x];
				kappaVals[2] = Ks[0];

				//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << "\n";

				//make the map
				bd::makeKappaMap(2, kappaVals, kappa);
				//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
				//do rewieght
				reweightL(N, num_states, db, particleTypes, eq, kappa);
				//get eq prob
				double eqProb = getEqProb(eq, targets);
				//copy Tconst into T
				std::copy(Tconst, Tconst+num_states*num_states, T);
				//fill in transposed entries such that T satisfies detailed balance
				satisfyDB(T, num_states, db, eq);
				//fill in diagonal with negative sum of row entries
				bd::fillDiag(T, num_states);
				//get the transition rate
				bd::computeMFPTsSP(num_states, T, targets, m);
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

void getBondTypesL(int N, int* particleTypes, Database* db, std::vector<int> targets) {
	//for each state in targets, check how many of each bond type each state has
	//hard-coded for 2 types

	for (int index = 0; index < targets.size(); index++) {
		int AA = 0; int AB = 0; int BB = 0;
		int state = targets[index];
		for (int i = 0; i < N; i++) {
			for (int j = i+2; j < N; j++) {
				if ( (*db)[state].isInteracting(i,j)) {
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
		printf("State: %d, AA: %d, AB %d, BB %d\n", state, AA, AB, BB);
	}
}

void solveFKE(int N, double* T, int num_states, double tf, Eigen::VectorXd& sol) {
	//evaluate exponential of the generator

	Eigen::MatrixXd Q; Q.setZero(num_states, num_states);
	for (int i = 0; i < num_states * num_states; i++) {
		Q(i) = T[i];
	}
	Q = Q * tf;

	Eigen::MatrixXd E; E.setZero(num_states, num_states);
	E = Q.exp();

	sol = E.row(0);
}

void evalStatsL(int N, Database* db, int initial, int target, bool useFile) {
	//evaluate the eq prob, hitting prob, rate, and configurations for the
	//properties in the input files

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = bd::setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	bd::initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0; T[i] = 0;
	}

	//init eq and m with zeros
	for (int i = 0; i < num_states; i++) {
		eq[i] = m[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; targets.push_back(target);

	//the following doesnt work for distinguishing lattice ground states
	/*
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	*/

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	bd::readKappaFile(numInteractions, kappaVals);

	bd::makeKappaMap(numTypes, kappaVals, kappa);
	//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
	//do rewieght
	reweightL(N, num_states, db, particleTypes, eq, kappa);
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	bd::fillDiag(T, num_states);
	//get the transition rate
	bd::computeMFPTsSP(num_states, T, targets, m);
	double rate = 1.0 / m[initial];

	for (int i = 0; i < 300; i++) {
		Eigen::VectorXd probs; probs.setZero(num_states);
		solveFKE(N, T, num_states, (10.0*i)/300.0, probs);
		std::cout << probs(target) << "\n";
	}
	Eigen::VectorXd probs; probs.setZero(num_states);
	solveFKE(N, T, num_states, 100.0/300.0, probs);
	printf("Final Time Occupation Probabilities:\n");
	for (int i = 0; i < num_states; i++) {
		printf("State: %d, Prob %f, mfpt %f\n", i, probs(i), m[i]);
	}

	//print results
	printf("Rate to target: %f\n", rate);
	printf("Ground State Eq Probabilities:\n");
	for (int i = 0; i < ground.size(); i++) {
		printf("State: %d, Eq: %f\n", ground[i], eq[ground[i]]);
	}

	//print the bond types
	std::cout << "Bond distributions:\n";
	getBondTypesL(N, particleTypes, db, ground);
	std::cout << "\n";

	//free memory
	delete []particleTypes; delete []Tconst; delete []kappaVals;
	delete []T; delete []m; delete []eq;


}




void HPscatter(int N, Database* db, int initial, int target) {
	/*Construct a scatter plot of avg rate vs equilibrium probability using 
	  every combination of A and B particles under the HP model
	  Here H = A, P = B, and AA = -1, AB = 0, BB = 0 */

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	bd::readKappaFile(numInteractions, kappaVals);

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
	std::vector<int> targets; targets.push_back(target);

	//the following doesnt work for distinguishing lattice ground states

	/*
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	*/

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;
	bd::makeKappaMap(numTypes, kappaVals, kappa);

	//declare outfile
	std::ofstream ofile;
	ofile.open("rateEqScatter.txt");

	//loop over all permutations
	std::deque<std::string> perms; 
	bd::allPerms(N,perms);
	//bd::distinctPerms(N, perms);
	int num_perms = perms.size();
	double* permProb = new double[num_perms];
	double* permRate = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		bd::setTypes(N, particleTypes, perms, p);

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//do rewieght
		reweightL(N, num_states, db, particleTypes, eq, kappa);
		//get eq prob - conditional on being in a ground state
		double eqProb = getEqProb(eq, targets);
		double groundProb = getEqProb(eq, ground);
		//printf("Eq %f, ground %f\n", eqProb, groundProb);
		eqProb /= groundProb;
		//copy Tconst into T
		std::copy(Tconst, Tconst+num_states*num_states, T);
		//fill in transposed entries such that T satisfies detailed balance
		//for (int q = 0; q < num_states; q++) eq[q] = 1.0;
		//for (int q = 0; q < num_states*num_states; q++) T[q] = 100.0;
		satisfyDB(T, num_states, db, eq);

		//fill in diagonal with negative sum of row entries
		bd::fillDiag(T, num_states);
		//get the transition rate
		bd::computeMFPTsSP(num_states, T, targets, m);
		double rate = 1/m[initial];

		permProb[p] = eqProb; permRate[p] = rate;

		//output to file
		ofile << eqProb << ' ' << rate << "\n";

		
	}

	//print out the probability and rate with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Eq Prob: %f, Rate: %f, ID: %s\n", permProb[p], permRate[p], perms[p].c_str());
	}


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []T; delete []eq; delete []m;
	delete []permProb; delete []permRate;

}

/**********************************************************************/
/***************** Measure testing on a lattice  **********************/
/**********************************************************************/

void buildAM(int N, Database* db, int target, int* M) {
	//build adj mat from database

	//zero out the matrix from previous searches
	for (int i = 0; i < N*N; i++) M[i] = 0;

	//fill the matrix with db data
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			M[i*N+j] = (*db)[target].isInteracting(i,j);
		}
	}

}

double harmonicMeanBarrier(int N, Database* db, int target, double* E) {
	//determine harmonic mean of the bond strengths in target state

	double m = 0;

	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( (*db)[target].isInteracting(i,j)) {
		    double energy = E[toIndex(i,j,N)];
		    printf("energy is %f\n", energy);
		    m += 1.0 / energy;
			}
		}
	}

	return 1.0 / m;
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

void createRateMatrix(int N, Database* db, int numTypes, double* kappaVals,
											int* particleTypes, 
											Eigen::MatrixXd& R, const Eigen::MatrixXd& r_form) {
	//use the database to construct the rate matrix at the given values
	//no derivative, so this is just a matrix

	//declare equilibrium measure and perform reweight
	int num_states = db->getNumStates();
	double* eq = new double[num_states];           //equilibrium measure
	std::map<std::pair<int,int>,double> kappa;    //map for interaction type to kappa
	bd::makeKappaMap(numTypes, kappaVals, kappa);
	reweightL(N, num_states, db, particleTypes, eq, kappa);

	double mult = 300.0;

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

	R *= mult;


	delete []eq;
}

double averageMisfoldEnergyTrap(int N, double Tf, int initial, double* kappaVals, Database* db, 
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

	//next we determine the "trap states" by sorting a list of rates
	std::vector<std::pair<double,int>> exit_rates;
	for (int i = 0; i < num_states; i++) {
		exit_rates.push_back(std::make_pair(fabs(R(i,i)),i));
	}

	//sort vector
	std::sort(exit_rates.begin(), exit_rates.end());

	//extract the top m elements and add to trap states
	int m = 6;
	std::vector<int> traps;
	for (int i = 0; i < m; i++) {
		traps.push_back(exit_rates[i].second);
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
	bd::computeHittingProbability(P, num_states, traps, U);
	Eigen::VectorXd probs = U.row(initial); 

	for (int i = 0; i < num_states; i++) {
		printf("State %d, prob %f\n", i, probs(i));
	}
	
	double prod = probs.dot(misfold);

	return prod;

}

void testMeasuresRecord(int N, Database* db, int initial, int target, bool useFile) {
	//test measures - output eq prob and mfpt, along with test meaasures

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = bd::setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	bd::readKappaFile(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0; T[i] = 0;
	}

	//init eq and m with zeros
	for (int i = 0; i < num_states; i++) {
		eq[i] = m[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; targets.push_back(target);


	//set interaction matrices
	int* P = new int[N*N]; 
	double* E = new double[N*N]; 

	std::map<std::pair<int,int>,double> kappaMap;
	bd::makeKappaMap(numTypes, kappaVals, kappaMap);
	bd::fillP(N, particleTypes, P, E, kappaMap);


	Eigen::VectorXd misfold;
	double Tf = 5.0;
	double M = 110;
	double mult = 1.08;

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoMap_lattice.txt");

	//loop over kappa values for triangle
	for (int i = 0; i < M; i++) {

		//setup kappa values
		for (int i = 0; i < N*N; i++) P[0] = 0; for (int i = 0; i < N*N; i++) E[0] = 0;
		kappaMap.clear();
		kappaVals[0] *= mult;
		bd::makeKappaMap(numTypes, kappaVals, kappaMap);
		bd::fillP(N, particleTypes, P, E, kappaMap);

		//do rewieght
		reweightL(N, num_states, db, particleTypes, eq, kappaMap);
		//get eq prob - conditional on being in a ground state
		double eqProb = getEqProb(eq, targets);
		for (int j = 0; j < num_states; j++) {
			printf("State %d, eq %f\n", j, eq[j]);
		}
		//copy Tconst into T
		std::copy(Tconst, Tconst+num_states*num_states, T);
		//fill in transposed entries such that T satisfies detailed balance
		//for (int q = 0; q < num_states; q++) eq[q] = 1.0;
		//for (int q = 0; q < num_states*num_states; q++) T[q] = 100.0;
		satisfyDB(T, num_states, db, eq);

		//fill in diagonal with negative sum of row entries
		bd::fillDiag(T, num_states);
		//get the transition rate
		bd::computeMFPTsSP(num_states, T, targets, m);
		double rate = 1/m[initial];

		//get measure 1
		/*
		double c1 = .25;
		double m1 = 1.0-exp(-c1*smallestBarrier(N, db, target, E));
		*/

		
		//double m1 = stayingProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	// Tconst, targets);
		

		double m1 = 1.0 - exp(-harmonicMeanBarrier(N, db, target, E));

		//get measure 2
		
		misfold.setZero(num_states);
		misfoldVector(N, db, target, misfold, E);
		for (int j = 0; j < num_states; j++) {
			printf("State %d, misfold %f\n", j, misfold(j));
		}
		double m2 = averageMisfoldEnergyTrap(N, Tf, initial, kappaVals, db, particleTypes, numTypes, 
																		 Tconst, misfold, targets);
		double c2 = 1.0;
		m2 = exp(-c2*m2);
		

		/*
		double m2 = finalTargetProb(N, T, initial, kappaVals, db, particleTypes, numTypes, 
																	 Tconst, targets);
		*/

		ofile << eqProb << ' ' << rate << ' ' << m1 << ' ' << m2 << "\n";
		std::cout  << eqProb << ' ' << rate << ' ' << m1 << ' ' << m2 << "\n";

	}


	ofile.close();



	//free memory
	delete []P; delete []E; delete []Tconst; delete []T; delete []m; delete []eq;
	delete []kappaVals; delete []particleTypes;
}

/**********************************************************************/
/***************** Genetic Algorithm Stuff  ***************************/
/**********************************************************************/

Person2::Person2() : ga::Person()  {}

Person2::Person2(int N_, int num_interactions_, int numTypes_, int* t, double* kV) : 
ga::Person(N_, num_interactions_, numTypes_, t, kV) {}

Person2::Person2(const Person2& old) : ga::Person(old) {}

void Person2::evalStats(int N, Database* db, int initial, std::vector<int> targets, 
								 double* eq, double* Tconst, double* T, double* m) {
	//evaluate stats for a given person

	//get the number of states
	int num_states = db->getNumStates();

	//zero out the arrays before setting
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] = 0;
	}
	for (int i = 0; i < num_states; i++) {
		eq[i] = m[i] = 0;
	}

	//do rewieght
	reweightL(N, num_states, db, types, eq, kappa);
	//get eqProb
	double eqProb = eq[targets[0]];
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	bd::fillDiag(T, num_states);
	//get the transition rate
	if (num_states > 300) { //do a sparse solve
		bd::computeMFPTsSP(num_states, T, targets, m);
	}
	else { //do a dense solve
		bd::computeMFPTs(num_states, T, targets, m);
	}
	double rate = 1.0 / m[initial];

	Rate = rate; Eq = eqProb;
}

Person2 Person2::mate(Person2 partner, bool useFile, RandomNo* rngee) {
	//create offspring from two parents

	//create new trait array
	double* kV = new double[num_interactions];
	int*    t  = new int[N];

	//take each trait w/ 50/50 chance from each parent
	for (int i = 0; i < num_interactions; i++) {
		double p = rngee->getU();
		
		if (p < 0.45) { //parent 1
			kV[i] = kappaVals[i];
		}
		else if (p > 0.45 && p < 0.9) { //parent 2
			kV[i] = partner.kappaVals[i];
		}
		else { //mutate
			kV[i] = ga::sampleKappa(N, rngee);
		}

		//include a chance to make a large parameter even larger
		if (kV[i] > 80 && kV[i] < 1e4) {
			double p2 = rngee->getU();
			if (p2 < 0.45) {
				kV[i] *= 2;
			}
		}

		//try setting a large parameter to a maximum value
		if (N == 7) {
			if (kV[i] > 1e2) {
				kV[i] = 1500;
			}
		}
		else {
			if (kV[i] > 1e3) {
				//kV[i] = 1e5;
			}
		}
	}

	if (useFile) {
		for (int i = 0; i < N; i++) {
			t[i] = types[i];
		}
	}
	else {
		//double p = rngee->getU();     //use this to keep the parents whole string
		for (int i = 0; i < N; i++) {
			double p = rngee->getU();  //use this to get particles individually
			
			if (p < 0.45) { //parent 1
				t[i] = types[i];
			}
			else if (p > 0.45 && p < 0.9) { //parent 2
				t[i] = partner.types[i];
			}
			else { //mutate
				t[i] = ga::sampleType(N, numTypes, rngee);
			}
		}
	}

	//create the offspring
	Person2 kid = Person2(N, num_interactions, numTypes, t, kV);

	//free memory
	delete []kV; delete []t;

	//birth the child
	return kid;
}

void printPopulation(std::vector<Person2> population, int pop_size, std::ofstream& ofile ) {
	//print the current population and stats to file
	int nt = population[0].num_interactions;

	for (int i = 0; i < pop_size; i++) {
		ofile << population[i].Eq << ' ' << population[i].Rate << ' ';
		for (int j = 0; j < nt; j++) {
			ofile << population[i].kappaVals[j] << ' ';
		}
		ofile << "\n";
	}
}

void printTypes(std::vector<Person2> population, int pop_size, int N) {
	//print the final types of each particle
	for (int i = 0; i < pop_size; i++) {
		std::cout << "Person " << i << ":";
		for (int j = 0; j < N; j++) {
			std::cout << population[i].types[j];
		}
		std::cout << ", EQ: " << population[i].Eq;
		std::cout << ", Rate: " << population[i].Rate;

		//for rate debugging
		/*
		std::cout <<" Rate: " << population[i].Rate << " ";
		std::cout << population[i].kappaVals[0] << ' ';
		std::cout << population[i].kappaVals[1] << ' ';
		std::cout << population[i].kappaVals[2] << ' ';
		*/
		
		std::cout << "\n";
	}
}


void performGAevolution(int N, Database* db, int initial, int target, bool useFile) {
	//run the genetic algorithm on the colloid database to find the pareto front
	//for the target state. Uses fixed particle types from file. 

	//get rid of the eigen parallelism
	Eigen::setNbThreads(0);

	//parameters to the genetic algorithm
	int generations = 1000;
	int pop_size    = 2000;
	double elite_p  = 0.1;
	double mates_p  = 0.4;
	bool printAll   = false;         //set true to make movie of output

	//if we have prior estimates of the max rate and eqProb, set here.
	//otherwise, this will update, adaptively. 
	double rateMax = 0.1; double eqMax = 0.1;

	//get database info - perturb if desired
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		//int IC = 1; 
		//numTypes = bd::setTypes(N, particleTypes, IC);
		numTypes = 3;
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix - only forward entries
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
	std::vector<int> targets; targets.push_back(target);
	/*
	bd::findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	*/


	//construct the initial population
	printf("Generating the initial population\n");
	Person2* pop_array = new Person2[pop_size];
	std::vector<Person2> population;
	#pragma omp parallel 
	{
	//declare all arrays we need to do calculations
	double* T = new double[num_states*num_states]; //rate matrix
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 
	RandomNo* rngee = new RandomNo();              //random number generator

	//loop over population
	#pragma omp for
	for (int i = 0; i < pop_size; i++) {
		std::cout << i << "\n";
		//draw parameters from a distribution
		ga::sampleParameters(N, numInteractions, kappaVals, particleTypes, numTypes, useFile, rngee);
		
		//create a person, evaluate their stats
		Person2 p = Person2(N, numInteractions, numTypes, particleTypes, kappaVals);
		//need new evalStats
		p.evalStats(N, db, initial, targets, eq, Tconst, T, m);
		p.evalFitness(eqMax, rateMax);
		pop_array[i] = p;
		//printf("e %f, r %f, f %f\n", pop_array[i].Eq, pop_array[i].Rate, pop_array[i].fitness);
		printf("Finsihing sample %d on thread %d\n", i, omp_get_thread_num());
	}
	//free memory
	delete []T; delete []eq; delete []m;
	delete rngee;
	//end parallel region
	}

	//move from array to vector
	for (int i = 0; i < pop_size; i++) {
		population.push_back(pop_array[i]);
	}
	delete []pop_array;

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoGAlattice.txt");

	//print the initial population
	if (printAll)
		printPopulation(population, pop_size, ofile);

	//init storage for storing eq and rate
	double* popEq = new double[pop_size];
	double* popRate = new double[pop_size];
	double* popFitness = new double[pop_size];

	//loop over generations
	for (int gen = 0; gen < generations; gen++) {
		//get the eq, rate, and fitness of each person

		printf("Spawning generation %d of %d\n", gen+1,generations);
		for (int i = 0; i < pop_size; i++) {
			popEq[i] = population[i].Eq; 
			popRate[i] = population[i].Rate; 
			popFitness[i] = population[i].fitness;

			//printf("Person %d, eq %f, rate %f\n", i, population[i].Eq, population[i].Rate);
			//printf("%f, %f, %f\n", population[i].kappaVals[0],population[i].kappaVals[1],population[i].kappaVals[2]);
		}

		//get the current max eq and rate to update the scalings
		double eM = *std::max_element(popEq, popEq + pop_size); 
		double rM = *std::max_element(popRate, popRate + pop_size); 
		if (eM > eqMax) {
			eqMax = eM;
		}
		if (rM > rateMax) {
			rateMax = rM;
		}

		//printf("Gen %d, em %f, rm %f\n", gen, eqMax, rateMax);

		//next, we sort the population by fitness, high to low
		std::vector<int> p(pop_size);
		std::iota(p.begin(), p.end(), 0);
		//std::sort(p.begin(), p.end(), [&](int i1, int i2) { return popFitness[i1] > popFitness[i2]; });
		
		double f = population[p[0]].fitness;
		std::cout << "max f " << f << "\n";

		//determine the non-dominated points
		std::vector<int> nonDom;
		ga::non_dominated_set(pop_size, popEq, popRate, nonDom);
		
		//create the new generation, perform elitism step
		std::vector<Person2> new_generation;
		int elites = nonDom.size(); 
		printf("Gen %d, num elites %d\n", gen, elites);
		for (int i = 0; i < elites; i++) {
			new_generation.push_back(population[nonDom[i]]);
		}

		//fill rest by mating the top percent of the present generation
		int rest = pop_size - elites;
		int top = mates_p * pop_size;
		Person2* pop_array = new Person2[rest];
		#pragma omp parallel 
		{
		//init the arrays
		double* T = new double[num_states*num_states]; //rate matrix
		double* eq = new double[num_states];           //equilibrium measure
		double* m = new double[num_states];            //mfpts 
		RandomNo* rngee = new RandomNo(); 
		//loop over population
		#pragma omp for
		for (int i = 0; i < rest; i++) {
			int r1 = p[floor(rngee->getU()*top)];
			int r2 = p[floor(rngee->getU()*top)];
			Person2 p1 = population[r1];
			Person2 p2 = population[r2];
			Person2 kid = p1.mate(p2, useFile, rngee);
			//need new eval stats
			kid.evalStats(N, db, initial, targets, eq, Tconst, T, m);
			kid.evalFitness(eqMax, rateMax);
			pop_array[i] = kid;
		}
		//end parallel region / free memory
		delete []T; delete []eq; delete []m;
		delete rngee;
		}

		//move from array to vector
		for (int i = 0; i < rest; i++) {
			new_generation.push_back(pop_array[i]);
		}
		delete []pop_array;

		//set the population equal to the newly generated one
		population = new_generation;

		//print if required
		if (printAll) {
			printPopulation(population, pop_size, ofile);
		}

		if (elites == pop_size) {
			break;
		}


	}

	//output the final results
	if (!printAll)
		printPopulation(population, pop_size, ofile);

	ofile.close();

	//print the particle types
	printTypes(population, pop_size, N);

	//free memory
	delete []particleTypes; delete []kappaVals; delete []Tconst;
	delete []popRate; delete []popEq; delete []popFitness;


}



/**********************************************************************/
/***************** Genetic Algorithm Sampling Stuff  ******************/
/**********************************************************************/

int yieldSample(int N, int Tf, int ts, double* E, int* M_target, int* types) {
	//return 1 if ground state forms

	//construct the chain of particles and the lattice mapping
	Particle* chain = new Particle[N];
	particleMap cMap;

	//initialize as linear chain
	initChain(N, chain, cMap, types);
	double energy0 = 0;

	//make another chain and cMap for reversions
	Particle* prev_chain = new Particle[N];
	particleMap prevMap;
	for (int i = 0;  i < N; i++) {
		prev_chain[i] = chain[i];
	}
	prevMap = cMap;

	//create a random number generator
	RandomNo* rngee = new RandomNo(); 

	//time per step is constant, set it
	double dt = 1;

	//set up an adjacency matrix for current and previous state
	int* M = new int[N*N];
	int* M_prev = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		M[i] = M_prev[i] = 0;
	}
	getAM(N, chain, M);
	getAM(N, chain, M_prev);

	double energy = energy0;         
	double timer = 0;                //init the timer for the mfpt
	bool accepted;
	double misfoldE = 0;
	bool misfoldFlag = false;
	double lastEvent = 0;

	//generate samples using MCMC until trapped or Tf steps
	for (int i = 0; i < Tf; i++) {
		//get the sample
		accepted = takeStep(N, chain, cMap, rngee, E, energy);
	}

	int formed = 0;
	getAM(N, chain, M);
	//printAM(N,M);
	//abort();
	//printAM(N,M_target);
	bool same = bd::checkSame(M, M_target, N);
	if (same) {
		formed = 1;
		//printf("hello\n");
	}


	//free memory 
	delete []chain; delete []prev_chain; delete []M; delete []M_prev;
	delete rngee;

	return formed;
}

double rateSampleTrap(int N, int Tf, int ts, double* E, int* M_target, int* types) {
	//sample a lattice protein trajectory for Tf steps. If stuck in trap for ts steps,
	//return the misfolded energy

	//construct the chain of particles and the lattice mapping
	Particle* chain = new Particle[N];
	particleMap cMap;

	//initialize as linear chain
	initChain(N, chain, cMap, types);
	double energy0 = 0;

	//make another chain and cMap for reversions
	Particle* prev_chain = new Particle[N];
	particleMap prevMap;
	for (int i = 0;  i < N; i++) {
		prev_chain[i] = chain[i];
	}
	prevMap = cMap;

	//create a random number generator
	RandomNo* rngee = new RandomNo(); 

	//time per step is constant, set it
	double dt = 1;

	//set up an adjacency matrix for current and previous state
	int* M = new int[N*N];
	int* M_prev = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		M[i] = M_prev[i] = 0;
	}
	getAM(N, chain, M);
	getAM(N, chain, M_prev);

	double energy = energy0;         
	double timer = 0;                //init the timer for the mfpt
	bool accepted;
	double misfoldE = 0;
	bool misfoldFlag = false;
	double lastEvent = 0;

	//generate samples using MCMC until trapped or Tf steps
	for (int i = 0; i < Tf; i++) {
		//get the sample
		accepted = takeStep(N, chain, cMap, rngee, E, energy);
		//printf("Energy is %f\n", energy);
		//printChain(N,chain);
		timer += 1;

		//if the position has changed, update timer and check for bond formation
		if (accepted) {
			//printf("Hello at time %f\n", timer);

			//get adjacency matrix for current state
			getAM(N, chain, M);
			//printAM(N,M);

			//compare to target matrix
			bool same = bd::checkSame(M, M_prev, N);
			//if different, log new last event time and new adj matrix
			if (!same) {
				lastEvent = timer;
				//printf("Action at %f\n", timer);
				//printChain(N,chain);
				for (int i = 0; i < N*N; i++) M_prev[i] = M[i];
			}


			//update the previous config
			for (int i = 0;  i < N; i++) {
				prev_chain[i] = chain[i];
			}
			prevMap = cMap;
		}

		//check for traps
		if (timer - lastEvent > ts) {
				//hit a trap, compute misfold energy
				misfoldE = ga::getMisfoldEnergy(N, M, M_target, E);
				misfoldFlag = true;
				//printf("Stuck in state for %f, misfold E = %f\n",timer-lastEvent, misfoldE);
				//printChain(N, chain);
				break;
			}

	}

	if (!misfoldFlag) {
		misfoldE = ga::getMisfoldEnergy(N, M, M_target, E);
		misfoldE = -1;
	}

	//free memory 
	delete []chain; delete []prev_chain; delete []M; delete []M_prev;
	delete rngee;

	return misfoldE;
}

void Person2::evalStatsSampling(double Tf, double ts, int samples, int* M_target) {
	//evaluate the stats of given member by sampling

	//evaluate the eq first
	int* P = new int[N*N]; for (int i = 0; i < N*N; i++) P[i] = 0;
	double* E = new double[N*N]; for (int i = 0; i < N*N; i++) E[i] = 0;
	bd::fillP(N, types, P, E, kappa);
	double c1 = 1.0;

	//evaluate the eq probility measure
	double hBarrier = ga::harmonicBarrier(N, M_target, E);
	//std::cout << hBarrier << "\n";
	//std::cout << E[2] << ' ' << E[3] << "\n";
	Eq = 1.0 - exp(-c1*hBarrier);

	//sample for the rate
	double rateEst = 0; int trapped = 0;
	for (int i = 0; i < samples; i++) {
		//rateEst += rateSampleBarrier(N, Tf, DT, rho, beta, E, P, method, pot, M_target);
		//rateEst += rateSampleTrap(N, Tf, DT, ts, rho, beta, E, P, method, pot, M_target);
		/*
		double result = rateSampleTrap(N, Tf, ts, E, M_target, types);
		if (result >= 0) {
			rateEst += result;
			trapped++;
		}
		*/

		rateEst += yieldSample(N, Tf, ts, E, M_target, types);
		trapped = samples;
		//std::cout << rateEst << "\n";
	}

	if (trapped == 0) {
		rateEst = 0;
	}
	else {
		rateEst /= trapped;
	}
	double c2 = 1.0;
	//Rate = exp(-c2*rateEst);
	Rate = rateEst;

	delete []E; delete []P;

}













void performGAlattice_sampling(int N, bool useFile) {
	//run the genetic algorithm on the colloid database to find the pareto front
	//for the target state. Uses fixed particle types from file. 

	//get rid of the eigen parallelism
	Eigen::setNbThreads(0);

	//parameters to the genetic algorithm
	int generations = 30;
	int pop_size    = 50;
	double elite_p  = 0.1;
	double mates_p  = 0.1;
	bool printAll   = false;         //set true to make movie of output

	//sampling parameters
	double Tf = 5000.0;           //N=8 Tf = 500
	double ts = 1500.0;           //N=8 ts = 300
	int samples = 500;
	double u_bound = 10000;

	//if we have prior estimates of the max rate and eqProb, set here.
	//otherwise, this will update, adaptively. 
	double rateMax = 1.0; double eqMax = 1.0;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		//int IC = 1; 
		//numTypes = bd::setTypes(N, particleTypes, IC);
		numTypes = 8;
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//set up target
	int* M_target = new int[N*N]; for (int i = 0; i < N*N; i++) M_target[i] = 0;
	double* X_target = new double[N*DIMENSION];
	for (int i = 0; i < N*DIMENSION; i++) X_target[i] = 0;
	ga::readTargetFile(M_target, X_target);
	//printAM(N,M_target);
	//abort();
	
	


	//construct the initial population
	printf("Generating the initial population\n");
	Person2* pop_array = new Person2[pop_size];
	std::vector<Person2> population;
	#pragma omp parallel 
	{
	//declare all arrays we need to do calculations
	RandomNo* rngee = new RandomNo();              //random number generator

	//loop over population
	#pragma omp for schedule(dynamic)
	for (int i = 0; i < pop_size; i++) {
		//draw parameters from a distribution
		ga::sampleParameters(N, numInteractions, kappaVals, particleTypes, numTypes, useFile, rngee);
		
		//create a person, evaluate their stats
		Person2 p = Person2(N, numInteractions, numTypes, particleTypes, kappaVals);
		p.applyBound(0.1, u_bound);
		p.evalStatsSampling(Tf, ts, samples, M_target);
		p.evalFitness(eqMax, rateMax);
		pop_array[i] = p;
		//printf("e %f, r %f, f %f\n", pop_array[i].Eq, pop_array[i].Rate, pop_array[i].fitness);
		printf("Finsihing sample %d on thread %d\n", i, omp_get_thread_num());
	}
	//free memory
	delete rngee;
	//end parallel region
	}

	//move from array to vector
	for (int i = 0; i < pop_size; i++) {
		population.push_back(pop_array[i]);
	}
	delete []pop_array;

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoGAlattice_sampling.txt");

	//print the initial population
	if (printAll)
		printPopulation(population, pop_size, ofile);

	//init storage for storing eq and rate
	double* popEq = new double[pop_size];
	double* popRate = new double[pop_size];
	double* popFitness = new double[pop_size];

	//loop over generations
	for (int gen = 0; gen < generations; gen++) {
		//get the eq, rate, and fitness of each person

		printf("Spawning generation %d of %d\n", gen+1,generations);
		for (int i = 0; i < pop_size; i++) {
			popEq[i] = population[i].Eq; 
			popRate[i] = population[i].Rate; 
			popFitness[i] = population[i].fitness;

			//printf("Person %d, eq %f, rate %f\n", i, population[i].Eq, population[i].Rate);
			//printf("%f, %f, %f\n", population[i].kappaVals[0],population[i].kappaVals[1],population[i].kappaVals[2]);
		}

		//get the current max eq and rate to update the scalings
		double eM = *std::max_element(popEq, popEq + pop_size); 
		double rM = *std::max_element(popRate, popRate + pop_size); 
		if (eM > eqMax) {
			eqMax = eM;
		}
		if (rM > rateMax) {
			rateMax = rM;
		}

		//printf("Gen %d, em %f, rm %f\n", gen, eqMax, rateMax);

		//next, we sort the population by fitness, high to low
		std::vector<int> p(pop_size);
		std::iota(p.begin(), p.end(), 0);
		//std::sort(p.begin(), p.end(), [&](int i1, int i2) { return popFitness[i1] > popFitness[i2]; });
		
		double f = population[p[0]].fitness;
		std::cout << "max f " << f << "\n";

		//determine the non-dominated points
		std::vector<int> nonDom;
		ga::non_dominated_set(pop_size, popEq, popRate, nonDom);
		
		//create the new generation, perform elitism step
		std::vector<Person2> new_generation;
		int elites = nonDom.size(); 
		printf("Gen %d, num elites %d\n", gen, elites);
		for (int i = 0; i < elites; i++) {
			new_generation.push_back(population[nonDom[i]]);
		}

		//fill rest by mating the top percent of the present generation
		int rest = pop_size - elites;
		int top = mates_p * pop_size;
		Person2* pop_array = new Person2[rest];
		#pragma omp parallel 
		{
		//init the arrays
		RandomNo* rngee = new RandomNo(); 
		//loop over population
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < rest; i++) {
			int r1 = p[floor(rngee->getU()*top)];
			int r2 = p[floor(rngee->getU()*top)];
			Person2 p1 = population[r1];
			Person2 p2 = population[r2];
			Person2 kid = p1.mate(p2, useFile, rngee);
			kid.applyBound(0.1, u_bound);
			kid.evalStatsSampling(Tf, ts, samples, M_target);
			kid.evalFitness(eqMax, rateMax);
			pop_array[i] = kid;
		}
		//end parallel region / free memory
		delete rngee;
		}

		//move from array to vector
		for (int i = 0; i < rest; i++) {
			new_generation.push_back(pop_array[i]);
		}
		delete []pop_array;

		//set the population equal to the newly generated one
		population = new_generation;

		//print if required
		if (printAll) {
			printPopulation(population, pop_size, ofile);
		}


	}

	//output the final results
	if (!printAll)
		printPopulation(population, pop_size, ofile);

	ofile.close();

	//print the particle types
	printTypes(population, pop_size, N);

	//free memory
	delete []particleTypes; delete []kappaVals; delete []M_target; delete []X_target;
	delete []popRate; delete []popEq; delete []popFitness;
}


void testSampling(int N) {
	//test using the MCMC sampler for a given setup

	int* particleTypes = new int[N];
	int numTypes = bd::readDesignFile(N, particleTypes);
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	for (int i = 0; i < numInteractions; i++) kappaVals[i] = 0;
	bd::readKappaFile(numInteractions, kappaVals);
	std::map<std::pair<int,int>,double> kappa;
	bd::makeKappaMap(numTypes, kappaVals, kappa);

	int* P = new int[N*N]; for (int i = 0; i < N*N; i++) P[i] = 0;
	double* E = new double[N*N]; for (int i = 0; i < N*N; i++) E[i] = 0;
	bd::fillP(N, particleTypes, P, E, kappa);


	//construct the chain of particles and the lattice mapping
	Particle* chain = new Particle[N];
	particleMap cMap;

	//initialize as linear chain
	initChain(N, chain, cMap, particleTypes);
	double energy0 = 0;

	//make another chain and cMap for reversions
	Particle* prev_chain = new Particle[N];
	particleMap prevMap;
	for (int i = 0;  i < N; i++) {
		prev_chain[i] = chain[i];
	}
	prevMap = cMap;

	//create a random number generator
	RandomNo* rngee = new RandomNo(); 

	//time per step is constant, set it
	double dt = 1;

	//set up an adjacency matrix for current and previous state
	int* M = new int[N*N];
	int* M_prev = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		M[i] = M_prev[i] = 0;
	}
	getAM(N, chain, M);
	getAM(N, chain, M_prev);

	double energy = energy0;         
	double timer = 0;                //init the timer for the mfpt
	bool accepted;
	double lastEvent = 0;
	int Tf = 50000;

	//generate samples using MCMC until trapped or Tf steps
	for (int i = 0; i < Tf; i++) {
		//get the sample
		accepted = takeStep(N, chain, cMap, rngee, E, energy);
		//printf("Energy is %f\n", energy);
		//printChain(N,chain);
		timer += dt;

		//if the position has changed, update timer and check for bond formation
		if (accepted) {
			//printf("Hello at time %f\n", timer);

			//get adjacency matrix for current state
			getAM(N, chain, M);
			//printAM(N,M);

			//compare to target matrix
			bool same = bd::checkSame(M, M_prev, N);
			//if different, log new last event time and new adj matrix
			if (!same) {
				lastEvent = timer;
				printf("Action at %f\n", timer);
				printChain(N,chain);
				for (int i = 0; i < N*N; i++) M_prev[i] = M[i];
			}


			//update the previous config
			for (int i = 0;  i < N; i++) {
				prev_chain[i] = chain[i];
			}
			prevMap = cMap;
		}


	}


	//free memory 
	delete []chain; delete []prev_chain; delete []M; delete []M_prev;
	delete rngee;
	delete []particleTypes; delete []kappaVals; delete []P; delete []E;
}

void testYield(int N) {
	//test yields generated by Pareto parameters by sampling

	int samples = 5000;

	int* particleTypes = new int[N];
	int numTypes = bd::readDesignFile(N, particleTypes);
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	for (int i = 0; i < numInteractions; i++) kappaVals[i] = 0;
	bd::readKappaFile(numInteractions, kappaVals);
	std::map<std::pair<int,int>,double> kappa;
	bd::makeKappaMap(numTypes, kappaVals, kappa);

	int* P = new int[N*N]; for (int i = 0; i < N*N; i++) P[i] = 0;
	double* E = new double[N*N]; for (int i = 0; i < N*N; i++) E[i] = 0;
	bd::fillP(N, particleTypes, P, E, kappa);

	//set up target
	int* M_target = new int[N*N]; for (int i = 0; i < N*N; i++) M_target[i] = 0;
	double* X_target = new double[N*DIMENSION];
	for (int i = 0; i < N*DIMENSION; i++) X_target[i] = 0;
	ga::readTargetFile(M_target, X_target);
	//printAM(N,M_target);
	//abort();

	int folded = 0;
	#pragma omp parallel for reduction(+: folded)
	for (int sample = 0; sample < samples; sample++) {

		//construct the chain of particles and the lattice mapping
		Particle* chain = new Particle[N];
		particleMap cMap;

		//initialize as linear chain
		initChain(N, chain, cMap, particleTypes);
		double energy0 = 0;

		//create a random number generator
		RandomNo* rngee = new RandomNo(); 

		//time per step is constant, set it
		double dt = 1;

		//set up an adjacency matrix for current and previous state
		int* M = new int[N*N];
		for (int i = 0; i < N*N; i++) {
			M[i] = 0;
		}

		double energy = energy0;         
		double timer = 0;                //init the timer for the mfpt
		bool accepted;
		double lastEvent = 0;
		int Tf = 50000;

		//generate samples using MCMC until Tf steps
		for (int i = 0; i < Tf; i++) {
			//get the sample
			accepted = takeStep(N, chain, cMap, rngee, E, energy);
		}

		getAM(N, chain, M);

		//check if in target
		bool same = bd::checkSame(M, M_target, N);
		if (same) {
			folded++;
			//printf("hello\n");
		}

		//free memory
		delete []chain;  delete []M; 
		delete rngee;
		//end parallel region
	}

	double yield = float(folded) / float(samples);
	printf("Yield is %f\n", yield);


	//free memory 
	delete []M_target; delete []X_target;
	delete []particleTypes; delete []kappaVals; delete []P; delete []E;
}



}