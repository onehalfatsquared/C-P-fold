#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"
#include "sampling.h"

//MAKE SURE TO RUN THIS ON A PURGED DATA SET EX input/N7Purge.txt

// Run Type: -1 to build a database
//					 0 for full run
//					 1 for chain estimate, with resets
//						2 for getttng trajectories

// Chain State: index of the linear chain in the DB. 1 for 6, 0 for 7.

int main(int argc, char* argv[]) {

	//handle input
	if (argc < 3) {
		fprintf(stderr, "Usage: <Input File> <Run Type> <Chain State> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	int rType = atoi(argv[2]); 
	int source;
	if (rType > 0) {
		if (argc != 4) {
			fprintf(stderr, "Usage: <Input File> <Run Type> <Chain State> %s\n", argv[0]);
			return 1;
		}
		source = atoi(argv[3]); 
	} 

	bool created = true;
	//build a database if one does not exist
	if (rType == -1) {
		if (created) { //the database file has been created. fill it
			//get the database here
			bd::Database* db = bd::readData(infile);
			//set database variables
			int N = db->getN(); int num_states = db->getNumStates();
			bd::addToDB(N, db);
		}
		else { //no file exists 
			bd::buildEmptyDB(6); //this line will construct a new database file to fill
		}
	}



	//otherwise, do estimation with the database
	//get the database here
	bd::Database* db = bd::readData(infile);

	//set database variables
	int N = db->getN(); int num_states = db->getNumStates();

	//tell user information
	printf("Database of states has been read.\n");


	if (rType == 0) {
		//do a full run over all states in DB

		printf("Mean first passage time estimator beginning.\n");
		//call estimator over every state, i. if i has 11 bonds (N=7), do nothing.
		for (int i = 0; i < num_states; i++) {
			if (((*db)[i].getBonds() >= 11 && N == 7) || ((*db)[i].getBonds() >= 9 && N == 6)) {//these states are rigid
				//do nothing
			}
			else{
				//call the estimator
				bd::estimateMFPT(N, i, db);
			}
		}

		//output stuff - for debug
		//std::cout << *db;

		//output the mfpt results to file
		std::string out = infile.substr(6,2);
		out = out + "mfpt.txt";
		std::ofstream out_str(out);
		out_str << *db; 

	}

	else if (rType == 1) {
		//just perform test on linear chain with resets

		/*
		bd::estimateChain(N, source, db);

		//output the mfpt results to file
		std::string out = infile.substr(6,2);
		out = out + "mfptChain.txt";
		std::ofstream out_str(out);
		out_str << *db; 
		*/

		//estimate a quantity at the first hitting time
		bd::sampleFirstExit(N, source, db);
		//bd::sampleFirstExitR(N, source, db);
		//bd::sampleSecondExit(N, source, db);
		//bd::sampleSecondExit(N, db);

		//bd::sampleQSD(N, db);
	}

		else if (rType == 2) {
		//get some trajectories
		sampleTrajectories(N, db, source);
	}

	//free memory - delete database
	delete db;

	return 0;
}