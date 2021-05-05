#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "database.h"
#include "graph.h"

/* Input Description:
	 Outside File - location of file containing outside clusters
	 Databae File - location of file containing my database - NxLump.txt
*/

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: <Outside File> <Database File> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	std::string infileD (argv[2]);

	//get the database here
	bd::Database* db = bd::readData(infileD);
	printf("Database of states has been read.\n");

	//make array to store mapping, mine to andreas
	int ns = db->getNumStates();
	int* M2A = new int[ns];
	for (int i = 0; i < ns; i++) {
		M2A[i] = -1; //any state not found will be -1
	}

	//get the mapping between structures
	int NC;
	//bd::structureMap(infile, db, M2A, NC);
	printf("Structures have been mapped\n");

	//make the transition matrix
	//bd::printTM(db, M2A, NC);
	bd::printTM(db); //no mapping
	printf("Transition matrix has been printed\n");



	
	
	//free memory - delete database
	delete db; delete []M2A;

	return 0;
}