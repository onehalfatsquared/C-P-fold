#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "database.h"
#include "tpt.h"
#include "nauty.h"
#include "design.h"
#include "protocol.h"

/* Initial conditions for target
	N = 6. Trapezoid: 7
				 Chevron  : 
				 Triangle : 

	N = 7. Boat:   41
				 Turtle: 48
				 Tree:   80
				 Flower: 90
*/



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 4) {
		fprintf(stderr, "Usage: <Database File> <Initial State>"
		" <Target State>  %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	int initial = atoi(argv[2]);
	int target = atoi(argv[3]);


	//build database from both files
	bd::Database* db = bd::readData(infile);
	int N = db->getN();

	//run code
	//bd::findProtocol(N, db, initial, target);
	//bd::testProtocol(N, db, initial, target);
	//bd::findProtocol2state();
	//bd::getGaps2state();
	//bd::flowerTestConstantBeta(N, db);

	//bd::modelFlowerSystem();
	bd::modelFlowerSystem_dlib();
	//bd::modelRowSystem_dlib();
	//bd::findProtocolSA_dlib(N, db, initial, target);
	//bd::findProtocolSA_dlib_E(N, db, initial, target, true);
	//bd::testProtocol_E(N, db, initial, target, true);
	//bd::testEqDeriv(N, db);
	//bd::findProtocolSA_adjoint(N, db, initial, target, true);
	//bd::testProbGradient(N, db);

	//bd::testRowSystemProtocols();

	//free memory - delete database
	//delete db; 

	return 0;
}