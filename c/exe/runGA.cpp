#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "genetics.h"

/* Takes in a database with hyperstatic clusters having a 0 equilibrium probability 
   and fills in the values using theory results.  */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 5) {
		fprintf(stderr, "Usage: <Database File> <initial state> <target state> <usefile>  %s\n", argv[0]);
		return 1;
	}
	std::string infile1 (argv[1]);
	int initial = atoi(argv[2]);
	int target = atoi(argv[3]);
	bool useFile = atoi(argv[4]);

	//get the database here
	bd::Database* db = bd::readData(infile1);

	//get num of particles
	int N = db->getN();

	//call the genetic algorithm
	ga::perform_evolution(N, db, initial, target, useFile);
	//ga::perform_evolution_sampling(N, useFile);

	//free memory 
	delete db; 

	return 0;
}