#include "latticeP.h"



namespace lattice {

/******************************************************************************/
/***************** Class functions for proteins *******************************/
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

int toIndex(int r, int c, int m) {
  //map row and column number into index in 1d array. column indexed
  return m*c+r;
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
			std::cout << positions[toIndex(j,i,N)] << ' ';
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
							RandomNo* rngee, double& energy) {
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
	std::vector<std::pair<int,int>> bonds;
	getBonds(N, chain, bonds);
	double e1 = -0.5 * bonds.size();
	//printf("Energy %f, %lu \n", e1, bonds.size());


	//get acc probability - do accept/reject step
	double a = std::min(1.0, exp(-(e1-e0)));
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
	double max_it = 10000;

	//initialize rng
	RandomNo* rngee = new RandomNo();

	//do the monte carlo steps
	for (int i = 0; i < max_it; i++) {
		takeStep(N, chain, cMap, rngee, energy);
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

	for (int i = 0; i < N; i++) {
		int xi = chain[i].x; int yi = chain[i].y;

		for (int j = i+2; j < N; j++) {
			int xj = chain[j].x; int yj = chain[j].y;

			//check if particles are distance 1 apart. 1-norm distance
			int dist = abs(xj-xi) + abs(yj-yi);
			if (dist == 1) {
				AM[toIndex(i,j,N)] = 1;
			}
			else {
				AM[toIndex(i,j,N)] = 0;
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
	int max_it = 1e4;

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
		bool accept = takeStep(N, chain, cMap, rngee, energy);
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



}