#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"
#include "hydro.h"
#include "../defines.h"

/* This executable asks for how many simulation data sets to consider and which data set
	 It then constructs a sequence of databases to combine to get an overall db to print */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 4) {
		fprintf(stderr, "Usage: <Hydro Yes/No> <Num files> <empty db file>  %s\n", argv[0]);
		return 1;
	}

	//read input and set parameters
	int hydro = atoi(argv[1]);
	int num_files = atoi(argv[2]);
	std::string db_file (argv[3]);
	int maxT = 5000;

	int runType = 2;
	std::ofstream ofile;
	ofile.open("N7transitions.txt");

	std::string base = "input/hydro/N7/chain/"; 
	double dt; int n_save;
	//check if hydrodynamics were on or off, set parameters accordingly
	if (hydro == 0) { //hydro is off
		dt = 0.01;
		n_save = 40;
	  base += "trueNOHD/trueNOHD";
	}
	else if (hydro == 1) { //hydro is on
		dt = 0.1;
		n_save = 40;
		base += "trueHD/trueHD";
	}
	else if(hydro == 2) { //other test
		dt = 0.1;
		n_save = 40;
		base += "hd50/hd50_test";
	}

	//get timestep info
	double tps = dt*n_save;  //elapsed time per timestep
	
	//first create 1 database to be the master file
	//construct empty database from file with no mfpt data
	bd::Database* dbMaster = bd::readData(db_file);

	//get num of particles from db
	int N = dbMaster->getN();


	//fill HCC data structure with hydrodynamics data
	std::string file1 = base + "1.config";
	bd::HCC* hc = bd::extractData(file1, N, maxT);
	//use HD data to fill a database 
	if (runType == 0) {
		bd::determineTransitions(hc, dbMaster, tps);
	}
	else if (runType == 1) {
		lumpPerms(dbMaster);
		bd::determineTransitionStates(hc, dbMaster, tps, ofile);
	}
	else if (runType == 2) {
		bd::determineTransitionTimes(hc, dbMaster, tps, ofile);
	}

	//delete hc, loop over the rest of the files, creating hc and db, and combine db
	delete hc;
	for (int i = 2; i <= num_files; i++) {
		//create the i-th filename
		std::stringstream ss;
		ss << i;
		std::string file = base + ss.str() +".config";

		//print that we are now analyzing i-th file
		printf("Now analyzing file number %d\n", i);

		//construct empty db, get hd data, fill the db
		bd::Database* db = bd::readData(db_file);
		lumpPerms(db);
		bd::HCC* hc = bd::extractData(file, N, maxT);
		if (runType == 0) {
			bd::determineTransitions(hc, db, tps);
			bd::combineHittingData(dbMaster, db);
		}
		else if (runType == 1) {
			bd::determineTransitionStates(hc, db, tps, ofile);
		}
		else if (runType == 2) {
			bd::determineTransitionTimes(hc, db, tps, ofile);
		}

		delete hc; delete db;
	}

	//output the database to file
	std::string out = base + "DB.txt";
	std::ofstream out_str(out);
	out_str << *dbMaster; 

	out_str.close();
	ofile.close();

	//free the memory - just delete the class objects
	delete dbMaster;
	return 0;
}