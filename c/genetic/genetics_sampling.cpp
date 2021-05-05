#include "genetics.h"
#include "bDynamics.h"

#define PI 3.1415926
#define CUT 1.2

namespace ga {

void printChain(double* X, int N) {
	//print out the chain

	for (int i = 0; i < N; i++) {
		printf("%f\n%f\n", X[DIMENSION*i], X[DIMENSION*i+1]);
	}
	printf("\n");

}

void readTargetFile(int* M_target, double* X_target) {
	//read in bonds present in the target state

	//file location
	std::string filename = "input/design/target.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return;
	}

	//get the top of file data
	int N; int c; //number of particles and number of constraints
	in_str >> N; in_str >> c;

	//store the adjacency matrix
	int p1; int p2;
	int count = 0;
	while (count < c) {
		//read line by line, add entry into AM
		in_str >> p1; in_str >> p2;

		M_target[bd::toIndex(p1,p2,N)] = 1;
		M_target[bd::toIndex(p2,p1,N)] = 1;

		count++;
	}

	//store the coordinates
	double x; double y; double z; count = 0;
	while (in_str >> x) {
		in_str >> y;

		X_target[count*DIMENSION] = x;
		X_target[count*DIMENSION+1] = y;

#if (DIMENSION == 3)
		in_str >> z;
		X_target[count*DIMENSION+2] = z;
#endif

		count++;
	}

	//fill in super and sub diagonal
	for (int i = 0; i < N-1; i++) {
		M_target[bd::toIndex(i,i+1,N)] = 1;
		M_target[bd::toIndex(i+1,i,N)] = 1;
	}
}

double eqSample(int N, double Tf, double DT, double rho, double beta, double* E, int* P,
								int method, int pot, int* M_target, double* X_target) {
	//start in the target configuration and see what fraction of sim time is spent there

	//set the target config
	double* X = new double[N*DIMENSION];
	for (int i = 0; i < N*DIMENSION; i++) {
		X[i] = X_target[i];
	}

	//set the adjacncy matrix
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;

	//solve sde, increment time in target state
	double t = 0; double eqTime = 0;
	while (t < Tf) {
		bd::solveSDE(X, N, DT, rho, beta, E, P, method, pot);
		t += DT;

		bd::getAdjCut(X, N, M, CUT);
		int same = bd::checkSame(M, M_target, N);

		if (same) {
			eqTime += DT;
		}
	}

	//free memory
	delete []X; delete []M;

	//return sample
	return eqTime / Tf;
}

double rateSample(int N, double Tf, double DT, double rho, double beta, double* E, int* P,
								int method, int pot, int* M_target) {
	//start in chain. determine minimal time to form target. 
	//if target does not form, return -1

	//set the target config
	double* X = new double[N*DIMENSION];
	bd::setupChain(X,N);

	//set the adjacncy matrix
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;

	//solve sde, check if target state is hit
	double t = 0; double hitTime = -1; 
	while (t < Tf) {
		bd::solveSDE(X, N, DT, rho, beta, E, P, method, pot);
		t += DT;

		bd::getAdjCut(X, N, M, CUT);
		int same = bd::checkSame(M, M_target, N);

		if (same) {
			hitTime = t;
			break;
		}
	}

	//free memory
	delete []X; delete []M;

	//return sample
	return hitTime;
}

double getMisfoldEnergy(int N, int* M, int* M_target, double* E) {
	//compute bond energy for bonds not present in target

	double badVibes = 0;
	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			int inSim = M[bd::toIndex(i,j,N)];
			int inTarget = M_target[bd::toIndex(i,j,N)];
			if (inSim == 1 && inTarget == 0) {
				badVibes += E[bd::toIndex(i,j,N)];
			}
		}
	}

	return badVibes;
}

double rateSampleBarrier(int N, double Tf, double DT, double rho, double beta, double* E, int* P,
								         int method, int pot, int* M_target) {
	//run trajectory until Tf. return misfolded energy barrier


	//set the target config
	double* X = new double[N*DIMENSION];
	bd::setupChain(X,N);

	//set the adjacncy matrix
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;

	//solve sde
	bd::solveSDE(X, N, Tf, rho, beta, E, P, method, pot);

	//get current adj mat
	bd::getAdjCut(X, N, M, CUT);

	//determine the energy of inccorect bonds
	double misfoldE = getMisfoldEnergy(N, M, M_target, E);
	//std::cout << "Bad energy " << misfoldE << "\n";

	//free memory
	delete []X; delete []M;

	//return sample
	return misfoldE;
}

double rateSampleTrap(int N, double Tf, double DT, double ts, double rho, double beta, double* E, int* P,
								         int method, int pot, int* M_target) {
	//run trajectory until stuck in a state for > t_s. return misfolded energy barrier

	//set the target config
	double* X = new double[N*DIMENSION];
	bd::setupChain(X,N);

	//set the adjacency matrix
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;
	int* M_prev = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;
	bd::getAdjCut(X, N, M, CUT);
	bd::getAdjCut(X, N, M_prev, CUT);

	//solve sde, check for traps
	double t = 0; double lastEvent = 0;
	double misfoldE = 0;
	while (t < Tf) {
		//solve sde
		bd::solveSDE(X, N, DT, rho, beta, E, P, method, pot);
		t += DT;

		//get adj matrix, compare to previous
		bd::getAdjCut(X, N, M, CUT);
		bool same = bd::checkSame(M, M_prev, N);
		//if different, log new last event time and new adj matrix
		if (!same) {
			lastEvent = t;
			//printf("Action at %f\n", t);
			//printChain(X,N);
			for (int i = 0; i < N*N; i++) M_prev[i] = M[i];
		}

		//check for traps
		if (t - lastEvent > ts) {
			//hit a trap, compute misfold energy
			misfoldE = getMisfoldEnergy(N, M, M_target, E);
			//printf("Stuck in state for %f, misfold E = %f\n",t-lastEvent, misfoldE);
			//printChain(X,N);
			break;
		}
	}

	//printChain(X,N);
	misfoldE = getMisfoldEnergy(N, M, M_target, E);

	//free memory
	delete []X; delete []M; delete []M_prev;

	//return sample
	return misfoldE;
}


void doSampling(int N, double Tf, int* types, std::map<std::pair<int,int>,double> kappa,
							  double& eq, double& time, int* M_target, double* X_target) {
	//sample a trajectory and compute the fraction of time spent in target as well
	//as the minimum time the target was hit

	//set parameters for sde solve
	double DT = 0.01;
	double beta = 1.0;  double rho = 10;
	double method = 1; double pot = 0;

	int* P = new int[N*N];
	double* E = new double[N*N];
	bd::fillP(N, types, P, E, kappa);

	eq = eqSample(N, Tf, DT, rho, beta, E, P,method, pot, M_target, X_target);
	//time = rateSample(N, Tf, DT, rho, beta, E, P,method, pot, M_target);
	time = exp(-rateSampleBarrier(N, Tf, DT, rho, beta, E, P, method, pot, M_target));

}

void Person::evalStats(double Tf, int samples, int* M_target, double* X_target) {
	//evaluate a rate and eq prob surrogate by sampling
	//eq - time spent in target after being initialize there
	//rate - average of exponential of misfolded energy

	//estimate objectives from given number of sample trajectories
	double eqEstimate = 0; double rateEstimate = 0;
	double hits = 0;

	//get the samples
	for (int i = 0; i < samples; i++) {
		//get sample
		double eqSample = 0; double timeSample = 0;
		doSampling(N, Tf, types, kappa, eqSample, timeSample, M_target, X_target);

		//update eq time estimate
		eqEstimate += eqSample; 

		//update rate estimate. gets 0 if structure did not form
		/*
		if (timeSample > 0.001) {
			rateEstimate += 1.0 / timeSample;
			hits++;
		}
		*/
		rateEstimate += timeSample;
	}

	//weight by num samples
	eqEstimate /= double(samples);
	/*
	if (hits > 0) {
		rateEstimate /= double(hits);
	}
	*/
	rateEstimate /= double(samples);

	//set the class variables
	Eq = eqEstimate;
	Rate = rateEstimate;

} 

double smallestBarrier(int N,  int* M_target, double* E) {
	//determine the weakest bond in the target state under the given configuration

	double minE = 40;
	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( M_target[bd::toIndex(i,j,N)] == 1 ) {
		    double energy = E[bd::toIndex(i,j,N)];
		    if (energy < minE) {
		    	minE = energy;
		    	//printf("minE is %f\n", minE);
		    }
			}
		}
	}

	return minE;
}

double harmonicBarrier(int N,  int* M_target, double* E) {
	//estimate a barrier to break first bond by harmonic average of bond strengths in target

	double recip = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( M_target[bd::toIndex(i,j,N)] == 1 ) {
		    double energy = E[bd::toIndex(i,j,N)];
		    //printf("i %d, j %d, E %f\n", i,j,energy);
		    recip += 1.0 /energy;
			}
		}
	}

	return 1.0 / recip ;
}

void Person::evalStats(double Tf, double ts, int samples, int* M_target) {
	//eval a eq and rate surrogate
	//eq - 1-exp(smallest barrier) -> no sampling needed
	//rate - expoential of average misfolded trap energy

	//sampling parameters
	double DT = 0.01;
	double beta = 1.0;  double rho = 20;
	double method = 1; double pot = 0;

	//evaluate the eq first
	int* P = new int[N*N]; for (int i = 0; i < N*N; i++) P[i] = 0;
	double* E = new double[N*N]; for (int i = 0; i < N*N; i++) E[i] = 0;
	bd::fillP(N, types, P, E, kappa);
	double c1 = 1.0;
	//Eq = 1.0 - exp(-c1*smallestBarrier(N, M_target, E));
	Eq = 1.0 - exp(-c1*harmonicBarrier(N, M_target, E));

	//sample for the rate
	double rateEst = 0;
	for (int i = 0; i < samples; i++) {
		//rateEst += rateSampleBarrier(N, Tf, DT, rho, beta, E, P, method, pot, M_target);
		rateEst += rateSampleTrap(N, Tf, DT, ts, rho, beta, E, P, method, pot, M_target);
	}
	rateEst /= samples;
	double c2 = 1.0;
	Rate = exp(-c2*rateEst);

	delete []E; delete []P;

}

void Person::applyBound(double lower, double upper) {
	//apply bounds to kappa values

	for (int i = 0; i < num_interactions; i++) {
		if (kappaVals[i] > upper) {
			kappaVals[i] = upper;
		}
		else if (kappaVals[i] < lower) {
			kappaVals[i] = lower;
		}
	}

	bd::makeKappaMap(numTypes, kappaVals, kappa);
}






void perform_evolution_sampling(int N, bool useFile) {
	//run the genetic algorithm on the colloid database to find the pareto front
	//for the target state. Uses fixed particle types from file. 

	//get rid of the eigen parallelism
	Eigen::setNbThreads(0);

	//parameters to the genetic algorithm
	int generations = 30;
	int pop_size    = 100;
	double elite_p  = 0.1;
	double mates_p  = 0.4;
	bool printAll   = false;         //set true to make movie of output

	//sampling parameters
	double Tf = 10.0;
	double ts = 4.0;
	int samples = 30;

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
		numTypes = 2;
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//set up target
	int* M_target = new int[N*N]; for (int i = 0; i < N*N; i++) M_target[i] = 0;
	double* X_target = new double[N*DIMENSION];
	for (int i = 0; i < N*DIMENSION; i++) X_target[i] = 0;
	readTargetFile(M_target, X_target);
	


	//construct the initial population
	printf("Generating the initial population\n");
	Person* pop_array = new Person[pop_size];
	std::vector<Person> population;
	#pragma omp parallel 
	{
	//declare all arrays we need to do calculations
	RandomNo* rngee = new RandomNo();              //random number generator

	//loop over population
	#pragma omp for schedule(dynamic)
	for (int i = 0; i < pop_size; i++) {
		//draw parameters from a distribution
		sampleParameters(N, numInteractions, kappaVals, particleTypes, numTypes, useFile, rngee);
		
		//create a person, evaluate their stats
		Person p = Person(N, numInteractions, numTypes, particleTypes, kappaVals);
		p.applyBound(0.1, 1000);
		//p.evalStats(Tf, samples, M_target, X_target);
		p.evalStats(Tf, ts, samples, M_target);
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
	ofile.open("paretoGAsampling.txt");

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
		non_dominated_set(pop_size, popEq, popRate, nonDom);
		
		//create the new generation, perform elitism step
		std::vector<Person> new_generation;
		int elites = nonDom.size(); 
		printf("Gen %d, num elites %d\n", gen, elites);
		for (int i = 0; i < elites; i++) {
			new_generation.push_back(population[nonDom[i]]);
		}

		//fill rest by mating the top percent of the present generation
		int rest = pop_size - elites;
		int top = mates_p * pop_size;
		Person* pop_array = new Person[rest];
		#pragma omp parallel 
		{
		//init the arrays
		RandomNo* rngee = new RandomNo(); 
		//loop over population
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < rest; i++) {
			int r1 = p[floor(rngee->getU()*top)];
			int r2 = p[floor(rngee->getU()*top)];
			Person p1 = population[r1];
			Person p2 = population[r2];
			Person kid = p1.mate(p2, useFile, rngee);
			kid.applyBound(0.1, 1000);
			//kid.evalStats(Tf, samples, M_target, X_target);
			kid.evalStats(Tf, ts, samples, M_target);
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








































}