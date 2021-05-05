#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "database.h"
#include "tpt.h"
#include "nauty.h"
#include "design.h"

/* Feed in the file with the un-lumped states, the initial, and target state.
	 Specify whether the particle identities come from file or will be manually set. 
	 The code will construct the set B of permutations of the target, and 
	 compute the hitting probabilities given a set of identities for the particles
*/

/* Initial conditions for target
	N = 6. Trapezoid: 10
				 Chevron  : 65
				 Triangle : 144
*/



int main(int argc, char* argv[]) {

	//handle input
	if (argc < 5 || argc > 6) {
		fprintf(stderr, "Usage: <Unlumped File> <Initial State>"
		" <Target State> <use file> <optional file> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	int initial = atoi(argv[2]);
	int target = atoi(argv[3]);
	bool useFile = atoi(argv[4]);
	std::string opt_file;
	if (argc == 6) {
		opt_file = (argv[5]);
	}

	//build database from both files
	bd::Database* db = bd::readData(infile);
	int N = db->getN();

	//do stuff

	//make surface of hitting probs w/ fixed AB interactions
	//bd::constructSurfaceTOY(N, db, initial, target, useFile);
	//bd::constructScatterTOY(N, db, initial, target, useFile);
	//bd::constructScatterTOY1(N, db, initial, target);
	//bd::constructScatterTOYsensitivity(N, db, initial, target, useFile);

	//do a maximization using the particle labels in file
	//bd::hittingProbMaxTOY(N, db, initial, target, useFile);
	//bd::eqProbMaxTOY(N, db, initial, target, useFile);
	//bd::rateMaxTOY(N, db, initial, target, useFile);

	//do maximization over all distinct permutations of particles
	//bd::hittingProbMaxTOYperms(N, db, initial, target);
	//bd::eqProbMaxTOYperms(N, db, initial, target);
	//bd::rateMaxTOYperms(N, db, initial, target);

	//constrained maximization
	//bd::hittingProbMaxTOYc(N, db, initial, target, useFile);
	//bd::eqProbMaxTOYc(N, db, initial, target, useFile);
	//bd::rateMaxTOYc(N, db, 0.9, initial, target, useFile);

	//do some sampling
	//bd::estimateHittingProbability(N, db, target);

	//evaluate the eq prob and rate at kappa in file
	bd::evalStats(N, db, initial, target, useFile);
	//bd::miskinOptimize(N, db, initial, target, useFile);
	//bd::computeParetoFront(N, db, initial, target, useFile);
	//bd::computeParetoFrontGD(N, db, initial, target, useFile);

	//tpt stuff
	//performTPT(N, db, initial, target, useFile);

	//irreversible stuff
	//irreversibleGraph(N, db, initial, target, useFile); //graphs for aleks 
	//generateAllGraphs(N, db, initial);
	//findQuenches(N, db);
	//graphQuenches(N, db, initial);

	//bd::evolveProbability(N, db, initial, useFile);
	//bd::testTransitionTimes(N, db, initial, target);

	//time dependent stuff
	//phaseTrajectory(N, db, initial, target, useFile, opt_file);

	//measure testing
	//bd::testMeasures(N, db, initial, target, useFile);
	//bd::testMeasuresRecord(N, db, initial, target, useFile);
	//bd::testMeasuresScatter(N, db, initial, target, useFile);






//

	//free memory - delete database
	delete db; 

	return 0;
}