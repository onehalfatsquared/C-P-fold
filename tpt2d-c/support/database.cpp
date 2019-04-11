#include "database.h"
#include "nauty.h"
#include <../Eigen/Dense>
#include <fstream>
#include <cstdlib>
namespace bd {


//state constructor
State::State() {
	am = NULL; coordinates = NULL; P = NULL;
	freq = 0; bond = 0; num = 0; denom = 0;
	num_coords = 0; mfpt = 0; Z = 0; Zerr = 0;
	sigma = 0;
}

//state deconstructor
State::~State() {
	delete []am; delete []coordinates; delete []P;
}

//database constructor
Database::Database(int N_, int num_states_) {
	N = N_; num_states = num_states_;
	states = new State[num_states];
}

//database deconstructor
Database::~Database() {
	delete []states;
}

//sum the entries of s.P
int State::sumP(int num_states) const{
	int S = 0;
	for (int i = 0; i < num_states; i++) {
		S += P[i];
	}
	return S;
}

//pull a random set of coordinates from the available
const Cluster& State::getRandomIC() const {//todo fix rng
	int rand_state = rand() % num_coords;
	return coordinates[rand_state];
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
	int index; //loop over states
	double x, y; //coordinates in a point
	char extra; //flag for whether mfpt estimates are in file

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
		in_str >> coord;
		s.num_coords = coord;

		//fill in the sample coordinates
		s.coordinates = new Cluster[coord];
		for (int i = 0; i < coord; i++) {
			s.coordinates[i].setNumPoints(N);
			for (int j = 0; j < N; j++) {
				in_str >> x >> y;
				s.coordinates[i][j] = Point(x,y);
			}
		}

		//check the extra value for existence of mfpt estimates
		in_str >> extra;
		if (extra == 'N') {//no mfpt estimates, initialize to 0
			s.num = 0;
			s.denom = 0;
			s.P = new int[num_lines];
			for (int i = 0; i < num_lines; i ++) {
				s.P[i] = 0;
			}
			s.mfpt = 0;
		}
		else if (extra == 'Y') {//mfpt estimates exist, read in
			in_str >> s.num;
			in_str >> s.denom;
			s.P = new int[num_lines];
			for (int i = 0; i < num_lines; i ++) {
				in_str >> s.P[i];
			}
			in_str >> s.mfpt;
			in_str >> s.Z; 
			in_str >> s.Zerr;
		}

		//next state
		index++;
	}
	in_str.close();
	return database;
}

//write functions to output the updated database to a file
std::ostream& State::print(std::ostream& out_str, int N, int num_states) const {
	for (int i = 0; i < N*N; i++) {
		out_str << am[i] << ' ';
	}
	out_str << freq << ' ';
	out_str << bond << ' ';
	out_str << num_coords << ' ';
	for (int i = 0; i < num_coords; i++) {
		for (int j = 0; j < N; j++) {
			out_str << coordinates[i][j].x << ' ' << coordinates[i][j].y << ' ';
		}
	}
	out_str << "Y" << ' ';
	out_str << num << ' ';
	out_str << denom << ' ';
	for (int i = 0; i < num_states; i++) {
		out_str << P[i] << ' ';
	}
	out_str << mfpt << ' ';
	out_str << Z << ' ';
	out_str << Zerr << '\n';
}

std::ostream& operator<<(std::ostream& out_str, const Database& db) {
	out_str << db.N << '\n';
	out_str << db.num_states << '\n';
	for (int i = 0; i < db.num_states; i++) {
		db[i].print(out_str, db.N, db.num_states);
	}
}







void makeNM(int N, int state, Database* db, Eigen::VectorXd x, Eigen::MatrixXd& J, Eigen::VectorXd& F) {
	//make the matrix and vector to solve for Newtons method

	double XD, YD;
	int count = 0;

	//loop over and construct system
	for (int i = 0; i <N; i ++) {
		for (int j = i+1; j < N; j++) {
			if ((*db)[state].isInteracting(i,j,N)) {
				XD = x(2*i) - x(2*j);
				YD = x(2*i+1) - x(2*j+1);
				F(count) = XD*XD + YD*YD -1.0;
				J(count, 2*i) = 2*XD; J(count, 2*j) = -2*XD;
				J(count, 2*i+1) = 2*YD; J(count, 2*j+1) = -2*YD;
				count +=1;
			}
		}
	}
}


bool checkPhysicalState(int N, int state, Database* db) {
	//check if a state is physical - use as start point in Newtons method

	//check if there are any example coordinates
	int coord = (*db)[state].getNumCoords();
	if (coord == 0) {
		printf("No coordinates\n");
		return false;
	}

	//set parameters to newtons method
	int max_iter = 50;
	double tol = 1e-8;

	//get the number of bonds
	int b = (*db)[state].getBonds(); 

	//initialize the matrix and vectors
	Eigen::MatrixXd J(b,2*N); 
	Eigen::VectorXd F(b); 
	Eigen::VectorXd dx(2*N);  
	Eigen::VectorXd x(2*N); 

	//initialize x. fill others with zeros.
	Cluster xc = (*db)[state].getRandomIC();
	for (int i = 0; i < N; i++) {
		x(2*i) = xc.points[i].x; x(2*i+1) = xc.points[i].y;
	}
	dx.fill(0.0); F.fill(0.0); J.fill(0.0);

	//output starting cluster - debug
	for (int i =0; i < 2*N; i++) printf("%f\n", x(i));
	printf("\n");


	//fill F and J. do solve with svd decomp.
	int iter;
	for (iter = 0; iter < max_iter; iter++) {
		makeNM(N, state, db, x, J, F);
		dx = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-F);
		x = x+dx;
		if (dx.norm() < tol) {
			break;
		}
	}

	//print message if no convergence
	if (iter == max_iter) {
		printf("State %d, Newton's Method terminated with residual %f\n", state, dx.norm());
	}

	printf("%d, %d\n", state, iter);

	//output final cluster - debug
	for (int i =0; i < 2*N; i++) printf("%f\n", x(i));

	//compare initial and post newton adjacency matrix

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

	//build nauty graph of the old state, get canonical labeling
	buildNautyGraph(N, M, state, db, g1);
	EMPTYGRAPH(cg1, M, N);
	densenauty(g1, lab1, ptn, orbits, &options, &stats, M, N, cg1);

	//get adjacnecy matrix of post newton state
	int* AM = new int[N*N]; 
	double* X = new double[2*N]; for (int i = 0; i < 2*N; i++) X[i]=x(i);
	getAdj(X, N, AM);

	//make graph
	//zero out any pre-existing graph
	EMPTYGRAPH(g2, M, N); EMPTYGRAPH(cg2, M, N);

	//add edges for every non-zero entry in adjacnecy matrix
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			if (AM[toIndex(i,j,N)] == 1) {
				ADDONEEDGE(g2, i, j, M);
			}
		}
	} 

	//free memory
	delete []X; delete []AM;

	//build nauty graph of the old state, get canonical labeling
	EMPTYGRAPH(cg2, M, N);
	densenauty(g2, lab2, ptn, orbits, &options, &stats, M, N, cg2);

	return checkIsomorphic(N, M, cg1, cg2);
}


void getPurgeStates(Database* db, std::vector<int>& toPurge) {
	//loop over all states. if unphysical, add to a list to purge.

	int N = db->getN(); int ns = db->getNumStates();

	for (int i = 5; i < 10; i++) {
		//printf("test says %d\n", checkPhysicalState(N, i, db));
		if (!checkPhysicalState(N, i, db)) {
			//state is unphysical, add to purge array.
			toPurge.push_back(i);
		}
	}
	printf("loop done\n");
}

void purgeUnphysical() {
	//write to new file, skip purge states
}



}