#include "genetics.h"

#define PI 3.1415926

namespace ga {


Person::Person() {
	Rate = Eq = 0.0; fitness = 0.0;
	N = 0;
	num_interactions = 0;
	numTypes = 0;
	kappaVals = NULL;
	types = NULL;
}

Person::~Person() {
	destroy();
}

void Person::copy(const Person& old) {
	N = old.N;
	Rate = old.Rate;
	Eq = old.Eq;
	num_interactions = old.num_interactions;
	numTypes = old.numTypes;
	fitness = old.fitness;

	kappaVals = new double[num_interactions];
	for (int i = 0; i < num_interactions; i++) {
		kappaVals[i] = old.kappaVals[i];
	}

	types = new int[N];
	for (int i = 0; i < N; i++) {
		types[i] = old.types[i];
	}

	kappa = old.kappa;
}

void Person::destroy() {
	delete []kappaVals;
	delete []types;
}

Person::Person(int N_, int num_interactions_, int numTypes_, int* t, double* kV) {
	Rate = Eq = 0.0; fitness = 0.0;
	N = N_;
	num_interactions = num_interactions_;
	numTypes = numTypes_;
	kappaVals = new double[num_interactions];
	types     = new int[N];
	for (int i = 0; i < num_interactions; i++) {
		kappaVals[i] = kV[i];
	}
	for (int i = 0; i < N; i++) {
		types[i] = t[i];
	}
	bd::makeKappaMap(numTypes, kappaVals, kappa);
}

void Person::setKappa(int num_interactions_, double* kV) {
	//construct the array and set the kappa values 
	num_interactions = num_interactions_;
	kappaVals = new double[num_interactions];
	for (int i = 0; i < num_interactions; i++) {
		kappaVals[i] = kV[i];
	}
	bd::makeKappaMap(numTypes, kappaVals, kappa);
}

Person Person::mate(Person partner, bool useFile, RandomNo* rngee) {
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
			kV[i] = sampleKappa(N, rngee);
		}

		//include a chance to make a large parameter even larger
		if (kV[i] > 10 && kV[i] < 1e4) {
			double p2 = rngee->getU();
			if (p2 < 0.45) {
				kV[i] *= 10;
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
				kV[i] = 1e5;
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
				t[i] = sampleType(N, numTypes, rngee);
			}
		}
	}

	//create the offspring
	Person kid = Person(N, num_interactions, numTypes, t, kV);

	//free memory
	delete []kV; delete []t;

	//birth the child
	return kid;
}

void Person::evalStats(int N, bd::Database* db, int initial, std::vector<int> targets, 
											 double* eq, double* Tconst, double* T, double* m) {
	//compute the eq and rate for this person

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
	bd::reweight(N, num_states, db, types, eq, kappa);
	//get eq prob
	double eqProb = bd::getEqProb(initial, kappaVals, db, types, numTypes, targets);
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	bd::satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	bd::fillDiag(T, num_states);
	//get the transition rate
	//printf("%f, %f, %f\n", kappaVals[0], kappaVals[1], kappaVals[2]);
	if (num_states > 300) { //do a sparse solve
		bd::computeMFPTsSP(num_states, T, targets, m);
	}
	else { //do a dense solve
		bd::computeMFPTs(num_states, T, targets, m);
	}
	double rate = 1/m[initial];

	Rate = rate; Eq = eqProb;
}

void Person::evalFitness(double eqMax, double rateMax) {
	//evulate the fitness as weighted sum of eq and rate

	fitness = Eq / eqMax + Rate / rateMax;
}

double sampleKappa(int N, RandomNo* rngee) {
	//sample a kappa value


	double u = rngee->getU();
	/*
	if (N != 7) { //get any number from 0 to pi/2
		u *= PI * 0.5;
	}
	else if (N == 7) { //get a number from 0.1 to 1.57
		u = u * 1.47 + 0.1;
	}

	return tan(u);
	*/

	double E = 11.0*u + 0.8;
	return exp(E) / sqrt(2.0*40*40*E);
}

int sampleType(int N, int numTypes, RandomNo* rngee) {
	//sample a particle type

	double u = rngee->getU();
	return floor(numTypes*u);
}

void sampleParameters(int N, int numInteractions, double* kappaVals, int* particleTypes, 
											int numTypes, bool useFile, RandomNo* rngee) {
	//sample the parameters from some distribution

	//sample the kappas as tan(uniformly random)
	for (int j = 0; j < numInteractions; j++) {
		kappaVals[j] = sampleKappa(N, rngee);
	}

	//sample the types as bernoulli 50/50 - if not constant
	if (!useFile) {
		for (int j = 0; j < N; j++) {
			particleTypes[j] = sampleType(N, numTypes, rngee);
		}
	}

}


bool is_dominated(int pop_size, int me, double x, double y, double* xAll, double* yAll) {
	//determine of the points x and y are dominated by any others in the set
	double tol = 1e-12;

	for (int i = 0; i < pop_size; i++) {
		if (i != me) {
			double d1 = x - xAll[i]; double d2 = y - yAll[i];
			bool s1 = (d1 > tol); bool s2 = (d2 > tol);
			double S = fabs(d1) + fabs(d2);
			//printf("Diffs %f and %f, bools %d and %d\n", d1, d2, s1, s2);
			if (!s1 && !s2 && S > tol) { 
				return true;
			}
		}
	}

	return false;
}

void non_dominated_set(int pop_size, double* xAll, double* yAll, std::vector<int>& nonDom) {
	//find the elements of xAll and yAll that are non-dominated by anything else in the set

	for (int i = 0; i < pop_size; i++) {
		//printf("DOing person %d\n\n\n", i);
		double x = xAll[i]; double y = yAll[i];
		bool D = is_dominated(pop_size, i, x, y, xAll, yAll);
		if (!D) {
			nonDom.push_back(i);
		}
	}
}



void printPopulation(std::vector<Person> population, int pop_size, std::ofstream& ofile ) {
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

void printTypes(std::vector<Person> population, int pop_size, int N) {
	//print the final types of each particle
	for (int i = 0; i < pop_size; i++) {
		std::cout << "Person " << i << ":";
		for (int j = 0; j < N; j++) {
			std::cout << population[i].types[j];
		}

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



void perform_evolution(int N, bd::Database* db, int initial, int target, bool useFile) {
	//run the genetic algorithm on the colloid database to find the pareto front
	//for the target state. Uses fixed particle types from file. 

	//get rid of the eigen parallelism
	Eigen::setNbThreads(0);

	//parameters to the genetic algorithm
	int generations = 300;
	int pop_size    = 750;
	double elite_p  = 0.1;
	double mates_p  = 0.4;
	bool printAll   = false;         //set true to make movie of output
	bool perturb    = true;          //set true for sensitivity testing

	//if we have prior estimates of the max rate and eqProb, set here.
	//otherwise, this will update, adaptively. 
	double rateMax = 0.1; double eqMax = 0.1;

	//get database info - perturb if desired
	int num_states = db->getNumStates(); 
	if (perturb) {
		double freq_perturb_frac = 0.0; //perturb the eq prob data by up to this fraction
		double rate_perturb_frac = 0.75; //perturb the rate data by up to this fraction
		perturbDB(db, freq_perturb_frac, rate_perturb_frac);
		//perturbDBconstP(db);
	}

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

	//declare rate matrix - only forward entries
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	bd::createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; //targets.push_back(target);
	
	bd::findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}
	


	//construct the initial population
	printf("Generating the initial population\n");
	Person* pop_array = new Person[pop_size];
	std::vector<Person> population;
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
		//draw parameters from a distribution
		sampleParameters(N, numInteractions, kappaVals, particleTypes, numTypes, useFile, rngee);
		
		//create a person, evaluate their stats
		Person p = Person(N, numInteractions, numTypes, particleTypes, kappaVals);
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
	ofile.open("paretoGA.txt");

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
		double* T = new double[num_states*num_states]; //rate matrix
		double* eq = new double[num_states];           //equilibrium measure
		double* m = new double[num_states];            //mfpts 
		RandomNo* rngee = new RandomNo(); 
		//loop over population
		#pragma omp for
		for (int i = 0; i < rest; i++) {
			int r1 = p[floor(rngee->getU()*top)];
			int r2 = p[floor(rngee->getU()*top)];
			Person p1 = population[r1];
			Person p2 = population[r2];
			Person kid = p1.mate(p2, useFile, rngee);
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















}