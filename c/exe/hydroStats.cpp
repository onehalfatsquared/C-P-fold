#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"
#include "hydro.h"
#include "../defines.h"

/* This executable asks for how many simulation data sets to consider and which data set
	 It then computes statistics on the trajectories, ex. time avg r_g, 1st hitting distribution */

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

	std::string base = "input/hydro/N7/chain/"; 
	double dt; int n_save;
	//check if hydrodynamics were on or off, set parameters accordingly
	if (hydro == 0) { //hydro is off
		dt = 0.001;
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
		base += "hd45/hd45_";
	}

	//get timestep info
	double tps = dt*n_save;  //elapsed time per timestep

	//set the quantity to get stats on
	/* 0 -> radius of gyration 
	   1 -> 2-d bond orientational order parameter
	   2 -> end to end distance of chain
	   */

	int which = 0;

	//set the type of stats - 0 -> time average until FHT, 1 -> distribution at FHT
	//                        2 -> clusters at FHT
	//												3 -> dist. at second hitting time
	int stats_type = 1;

	//first create 1 database to be the master file
	//construct empty database from file with no mfpt data
	bd::Database* db = bd::readData(db_file);

	//get num of particles from db
	int N = db->getN();

	//create a vector to store samples
	std::vector<double> q;

	//declare an outfile
	std::ofstream ofile;

	//fill HCC data structure with hydrodynamics data
	std::string file1 = base + "1.config";
	bd::HCC* hc = bd::extractData(file1, N, maxT);
	//use HD data to compute stuff
	if (stats_type == 0) {
		timeAverageFHT(hc, db, q, which);
	}
	else if (stats_type == 1) {
		distributionFHT(hc, db, q, which);
	}
	else if (stats_type == 2) {
		ofile.open("fhtClusters.txt");
		clustersFHT(hc, db, ofile);
	}
	else if (stats_type == 3) {
		distributionFHT2(hc, db, q, which);
	}

	

	//delete hc, loop over the rest of the files, get samples
	delete hc; 
	for (int i = 2; i <= num_files; i++) {
		//create the i-th filename
		std::stringstream ss;
		ss << i;
		std::string file = base + ss.str() +".config";

		//print that we are now analyzing i-th file
		printf("Now analyzing file number %d\n", i);

		//construct empty db, get hd data, fill the db
		bd::HCC* hc = bd::extractData(file, N, maxT);
		//compute stats
		if (stats_type == 0) {
		timeAverageFHT(hc, db, q, which);
		}
		else if (stats_type == 1) {
			distributionFHT(hc, db, q, which);
		}
		else if (stats_type == 2) {
			clustersFHT(hc, db, ofile);
			if (i == num_files) {
				ofile.close();
			}
		}
		else if (stats_type == 3) {
			distributionFHT2(hc, db, q, which);
		}


		delete hc; 
	}


	//do stuff with the data in q
	if (stats_type == 0) {
		//just get mean and variance from time averages

		double M; double V;
		bd::sampleStats(q, M, V);
		printf("Mean %f, Std deviation %f\n", M, sqrt(V));
	}
	else if (stats_type == 1 || stats_type == 3) {
		//output a file with samples

		//declare outfile
		std::ofstream ofile;
		ofile.open("fhtDistribution.txt");
		for (int i = 0; i < q.size(); i++) {
			ofile << q[i] << "\n";
		}

		ofile.close();
	}


	//free memory
	delete db;
	return 0;
}