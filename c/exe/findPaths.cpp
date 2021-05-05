#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "database.h"
#include "graph.h"
#include "graphviz.h"
#include "tpt.h"

/* Input Description:
	 Input File - location of file containing mfpt data. Likely NxLump.txt
	 source - index of the worm state. 1 for 6, 0 for 7.
	 target - for a subgraph that ends at particular target
	 draw - if 1, draw cluster images on the graphviz output. if 0, draw circles.
	 clean - if 1, suppress output of data on edges. if 0, output rates
	 reduce - if 1, remove nodes with less than 7% probability.
*/

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 7) {
		fprintf(stderr, "Usage: <Input File> <source> <target> <draw> <clean> <reduce> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	int source = atoi(argv[2]);
	int target = atoi(argv[3]);
	int draw = atoi(argv[4]);
	int clean = atoi(argv[5]);
	int reduce = atoi(argv[6]);

	//get the database here
	bd::Database* db = bd::readData(infile);
	printf("Database of states has been read.\n");
	
	//make the graph
	bd::Graph* g = bd::makeGraph(db);
	printf("Database has been stored as a graph.\n");

	//get the geometric partition functions
	int ns = db->getNumStates();
	double* Z = new double[ns];
	computePartitionFn(ns, db, Z); 

	//print out the ending distribution
	bd::endDistribution(g);

	//print out most probable path
	//bd::MPP(g, source);

	//print the graph
	bd::printGraphPF(g, source, Z, draw, clean, reduce);

	//print out path with fastest rate
	if (target == -1) {
		bd::QP(g, source);
	}
	else {
		bd::QP(g, source, target);
	}

	//construct the subgraph that ends at target
	//bd::Graph* sub = bd::targetSubgraph(g, source, target);
	//printGraph(sub, source, draw, clean, reduce);

	//print most probable path ending at target
	//bd::MPP(sub,source);

	//construct subgraph that starts at some node
	//bd::Graph* sub2 = bd::sourceSubgraph(g, 17);
	//bd::printGraph(sub2, 17, draw, clean, reduce);
	//bd::endDistribution(sub2);

	//check end distr calc
	//bd::findConditionalEnd(g, source);

	//get graph with end state probabilities
	//bd::printGraphEndDistribution(g, source, reduce);


	//free memory - delete database
	delete db; delete g; //delete sub; delete sub2;
	delete []Z;

	return 0;
}