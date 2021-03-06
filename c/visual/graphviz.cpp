#include "database.h"
#include "graphviz.h"
#include "graph.h"
#include "pair.h"
#include <vector>
#include <deque>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include "../defines.h"

#define SCALE_REV 25
namespace bd{

void printCluster(std::ofstream& out_str, int index, int draw) {
	//print a node creation in graphviz

	if (draw == 0) {
		out_str << "\"" + std::to_string(index) + "\" [label= \"" + std::to_string(index)
		+ "\" , shape=circle, width = 2, regular = 1, style = filled, fillcolor=white]; \n";
	}
	else if (draw == 1) {
		out_str << "\"" + std::to_string(index) + "\" [label=\"\", shape=circle, width = 1, regular = 1," +
		" style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";
	}
} 

void sameRank(std::ofstream& out_str, std::vector<int> states) {
	//print a rank = same statement for all elements of states

	out_str << "{rank = same; ";
	for (int i = 0; i < states.size(); i++) {
		out_str << "\"" + std::to_string(states[i]) + "\";";
	}
	out_str << "}\n";
}

void makeEdgeClean(std::ofstream& out_str, int source, int target, double edgeWidth) {
	//draw an edge from source to target - no labels

	out_str << "\"" + std::to_string(source) + "\" -> \"" + 
	std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + "]\n";
}

void makeEdge(std::ofstream& out_str, int source, int target, double edgeWidth, double rate) {
	//draw an edge from source to target

	std::string s = std::to_string(rate);
	s.erase(s.end()-3,s.end());

	out_str << "\"" + std::to_string(source) + "\" -- \"" + 
	std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + 
	", label = " + s + "]\n";
}

void printGraph(Graph* g, int source, int draw, int clean, int reduce) {
	//construct a graphviz file to print the graph

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphviz.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "graph bd {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//print the source node
	printCluster(out_str, source, draw);
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double rate = (*g)[source].getEdgeRate(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (reduce == 0 || prob > PTOL) {
			if (!drawn[T]) {
				printCluster(out_str, T, draw);
				drawn[T] = 1; rankVec.push_back(T);
			}
			double edgeWidth = 40*node_prob*prob;
			if (clean == 0) {
				makeEdge(out_str, source, T, edgeWidth, rate);
			}
			else if (clean == 1) {
				makeEdgeClean(out_str, source, T, edgeWidth);
			}
			toVisit.push_back(T);
		}
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (rate > 0.01 && (reduce == 0 || prob > PTOL)) {
				if (!drawn[T]) {
					printCluster(out_str, T, draw);
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidth = 40*node_prob*prob;
				if (clean == 0) {
					makeEdge(out_str, node, T, edgeWidth, rate);
				}
				else if (clean == 1) {
					makeEdgeClean(out_str, node, T, edgeWidth);
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}

	//end reached - put a curly to end
	out_str << "}";	

	//free memory
	out_str.close();
	delete []drawn;
}

void printCluster(std::ofstream& out_str, int index, std::vector<double> end) {
	//print a node creation in graphviz. Include probability for each end state.

	//convert the vector of probs into a string
	std::string P = "(";
	for (int i = 0; i < end.size()-1; i++) {
		int val = end[i]*100;
		std::string s = std::to_string(val);
		//s.erase(s.end()-4,s.end());
		P = P + s + ",";
	}
	int val = end[end.size()-1]*100;
	std::string s = std::to_string(val);
	//s.erase(s.end()-4,s.end());
	P = P + s + ")";

	//output the node creation command
	//out_str << "\"" + std::to_string(index) + "\" [xlabel=\"" + P +" \", shape=circle, width = 1, regular = 1," +
	//" style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";

	//out_str << "\"" + std::to_string(index) + "\" [label=\"|" + P +" \", shape=record, width = 1, regular = 1," +
	//" style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";

	out_str << "\"" + std::to_string(index) + "\" [label=\"\n\n\n\n\n\n" + P +" \", shape=circle, width = 1, regular = 1," +
	" fontsize=50, style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";
} 

void printGraphEndDistribution(Graph* g, int source, int reduce) {
	//construct a graphviz file to print the graph - include probability for each end state
	//always draw, always clean, sometimes reduce

	//set usual input variables
	int clean = 1;

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphvizEND.txt";
	std::ofstream out_str(out);

	//probability tolerance for reduced graphs

	//write the graphviz header
	out_str << "graph bd {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//print the source node
	printCluster(out_str, source, (*g)[source].endDistr);
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double rate = (*g)[source].getEdgeRate(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (reduce == 0 || prob > PTOL) {
			if (!drawn[T]) {
				printCluster(out_str, T, (*g)[T].endDistr);
				drawn[T] = 1; rankVec.push_back(T);
			}
			double edgeWidth = 40*node_prob*prob;
			if (clean == 0) {
				makeEdge(out_str, source, T, edgeWidth, rate);
			}
			else if (clean == 1) {
				makeEdgeClean(out_str, source, T, edgeWidth);
			}
			toVisit.push_back(T);
		}
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (rate > 0.01 && (reduce == 0 || prob > PTOL)) {
				if (!drawn[T]) {
					printCluster(out_str, T, (*g)[T].endDistr);
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidth = 40*node_prob*prob;
				if (clean == 0) {
					makeEdge(out_str, node, T, edgeWidth, rate);
				}
				else if (clean == 1) {
					makeEdgeClean(out_str, node, T, edgeWidth);
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}

	//end reached - put a curly to end
	out_str << "}";	

	//free memory
	out_str.close();
	delete []drawn;
}

void printPath(std::vector<int> path, std::vector<double> val, std::string s) {
	//makes a graphviz file that prints the given path - with images, thats the point

	//declare a file to write output to
	std::string out;
	out = s + "graphviz.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "graph " + s + " {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//draw all the nodes
	for (int i = 0; i < path.size(); i++) {
		int index = path[i];
		out_str << "\"" + std::to_string(index) + "\" [label=\"\", shape=circle, width = 1, regular = 1," +
			" style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";
	}

	//draw all the edges
	for (int i = 0; i < path.size()-1; i++) {
		int source = path[i]; int target = path[i+1];
		std::string s = std::to_string(val[i]);
		s.erase(s.end()-3,s.end());
		int edgeWidth = 5;
		out_str << "\"" + std::to_string(source) + "\" -- \"" + 
			std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + 
			", label = " + s + "]\n";
	}

	//make all nodes have the same rank
	out_str << "{rank = same; ";
	for (int i = 0; i < path.size(); i++) {
		out_str << "\"" + std::to_string(path[i]) + "\";";
	}
	out_str << "}\n";

	//print ending curly
	out_str << "}";
}

void QP(Graph* g, int source) {
	//determine the quickest folding path (largest rate)

	//get number of nodes
	int ns = g->getN();

	//make arrays that will store max sums and paths
	double* M = new double[ns];
	std::vector<std::vector<int>> paths(ns);
	for (int i = 0; i < ns; i++) M[i] = 0;
	M[source] = 0; paths[source].push_back(source);

	//compute a topological ordering of the graph
	int* T = new int[ns];
	tOrder(g, T, source);

	//loop over nodes, go down edges, check for largest path
	for (int i = 0; i < ns; i++) {
		int node = T[i]; 
		int E = (*g)[node].getNumEdges();
		for (int edge = 0; edge < E; edge++) {
			int target = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double pathTotal = M[node] + rate;
			if (pathTotal > M[target]) {
				M[target] = pathTotal;
				std::vector<int> prevPath = paths[node];
				prevPath.push_back(target);
				paths[target] = prevPath;
			}
		}
	}
	
	//find the max of the array M, get corresponding path
	int ind; double maxVal = 0;
	for (int i = 0; i < ns; i++) {
		if (M[i] > maxVal) {
			maxVal = M[i];
			ind = i;
		}
	}
	std::vector<int> maxPath = paths[ind];

	//print to user
	printf("Quickest Path: ");
	for (int i = 0; i < maxPath.size()-1; i++) {
		printf("%d, ", maxPath[i]);
	}
	printf("%d. ", maxPath[maxPath.size()-1]);
	printf("Total Rate = %f\n", M[ind]);

	//fill the rates on the path
	std::vector<double> rates;
	fillRates(g, maxPath, rates);

	//get graphviz code to show path
	printPath(maxPath, rates, "QP");

	//free memory
	delete []T; delete []M;
}

void QP(Graph* g, int source, int target) {
	//determine the quickest folding path (largest rate) from source to target

	//get number of nodes
	int ns = g->getN();

	//make arrays that will store max sums and paths
	double* M = new double[ns];
	std::vector<std::vector<int>> paths(ns);
	for (int i = 0; i < ns; i++) M[i] = 0;
	M[source] = 0; paths[source].push_back(source);

	//compute a topological ordering of the graph
	int* T = new int[ns];
	tOrder(g, T, source);

	//loop over nodes, go down edges, check for largest path
	for (int i = 0; i < ns; i++) {
		int node = T[i]; 
		int E = (*g)[node].getNumEdges();
		for (int edge = 0; edge < E; edge++) {
			int target = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double pathTotal = M[node] + rate;
			if (pathTotal > M[target]) {
				M[target] = pathTotal;
				std::vector<int> prevPath = paths[node];
				prevPath.push_back(target);
				paths[target] = prevPath;
			}
		}
	}
	
	//Get the desired target max
	int ind = target; double maxVal = M[target];
	std::vector<int> maxPath = paths[ind];

	//print to user
	printf("Quickest Path ending at %d: ", target);
	for (int i = 0; i < maxPath.size()-1; i++) {
		printf("%d, ", maxPath[i]);
	}
	printf("%d. ", maxPath[maxPath.size()-1]);
	printf("Total Rate = %f\n", M[ind]);

	//fill the rates on the path
	std::vector<double> rates;
	fillRates(g, maxPath, rates);

	//get graphviz code to show path
	printPath(maxPath, rates, "QP");

	//free memory
	delete []T; delete []M;
}

void MPP(Graph* g, int source) {
	//determine the most probable path down the graph from given source

	//find the path by traveling from top to bottom
	std::vector<int> path; std::vector<double> prob;
	path.push_back(source);
	int E = (*g)[source].getNumEdges();
	while (E > 0) {
		findMaxProb(g, source, E, path, prob);
		E = (*g)[source].getNumEdges();
	}

	//print to user
	printf("Most Probable Path: ");
	for (int i = 0; i < path.size()-1; i++) {
		printf("%d, ", path[i]);
	}
	printf("%d\n", path[path.size()-1]);

	//make a graphviz file that shows the path
	printPath(path, prob, "MPP");
}

void printClusterRev(std::ofstream& out_str, int index, double pw) {
	//print a node creation in graphviz. size depends on free energy

	//output the node creation command
	out_str << "\"" + std::to_string(index) + "\" [label=\"\", shape=circle, width = 1, " +
	"penwidth = " + std::to_string(pw) + ", regular = 1," + " style = filled, fillcolor=white," 
	+ "image=\"c" + std::to_string(index) + ".png\"]; \n";
} 

void makeEdgeCleanRev(std::ofstream& out_str, int source, int target, double edgeWidth, 
											std::string color) {
	//draw an edge from source to target - no labels

	if (edgeWidth > 0.1) {
		out_str << "\"" + std::to_string(source) + "\" -> \"" + 
		std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) +
		", color= " + color + "]\n";
	}
}

void makeEdgeRev(std::ofstream& out_str, int source, int target, double edgeWidth, 
								 double rate, std::string color) {
	//draw an edge from source to target

	if (edgeWidth > 0.1) {
		std::string s = std::to_string(rate);
		s.erase(s.end()-3,s.end());

		out_str << "\"" + std::to_string(source) + "\" -> \"" + 
		std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + 
		", color= " + color + ", label = " + s + "]\n";
	}
}

double getPW(double f) {
	//return pen width of node for given free energy

	double pw = 15.0/(f + 1);
	return pw; 
}

double getMax(int num_states, double* array) {
	//get the max value of the array to scale by
	
	double M = 0;
	for (int i = 0; i < num_states*num_states; i++) {
		if (array[i] > M) {
			M = array[i];
		}
	}

	return M;
}

void printGraphRev(Graph* g, int source, double* F, double* flux, int draw, int clean, 
									 int reduce) {
	//construct a graphviz file to print the graph - directed for reversible case

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphvizRev.txt";
	std::ofstream out_str(out);

	//get the scaling for edgeWidth
	double fluxMax = getMax(ns, flux);
	double scaling = SCALE_REV / fluxMax;

	//write the graphviz header
	out_str << "digraph bd {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//print the source node
	double pw = getPW(F[source]);
	printClusterRev(out_str, source, pw);
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (reduce == 0 || prob > PTOL) {
			if (!drawn[T]) {
				double pw = getPW(F[T]);
				printClusterRev(out_str, T, pw);
				drawn[T] = 1; rankVec.push_back(T);
			}
			double edgeWidthF = flux[toIndex(source, T, ns)]*scaling;
			double edgeWidthR = flux[toIndex(T, source, ns)]*scaling;
			if (clean == 0) {
				makeEdgeRev(out_str, source, T, edgeWidthF, edgeWidthF/scaling, "green");
				makeEdgeRev(out_str, T, source, edgeWidthR, edgeWidthR/scaling, "red");
			}
			else if (clean == 1) {
				makeEdgeCleanRev(out_str, source, T, edgeWidthF, "green");
				makeEdgeCleanRev(out_str, T, source, edgeWidthR, "red");
			}
			toVisit.push_back(T);
		}
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (reduce == 0 || prob > PTOL) {
				if (!drawn[T]) {
					double pw = getPW(F[T]);
					printClusterRev(out_str, T, pw);
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidthF = flux[toIndex(node, T, ns)]*scaling;
				double edgeWidthR = flux[toIndex(T, node, ns)]*scaling;
				if (clean == 0) {
					makeEdgeRev(out_str, node, T, edgeWidthF, edgeWidthF/scaling, "green");
					makeEdgeRev(out_str, T, node, edgeWidthR, edgeWidthR/scaling, "red");
				}
				else if (clean == 1) {
					makeEdgeCleanRev(out_str, node, T, edgeWidthF, "green");
					makeEdgeCleanRev(out_str, T, node, edgeWidthR, "red");
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}

	//end reached - put a curly to end
	out_str << "}";	

	//free memory 
	out_str.close();
	delete []drawn;
}

void printClusterEq(std::ofstream& out_str, int index, double eq, double hit, double pw) {
	//print a node creation in graphviz - include hit and eq probs

	std::string sEq = std::to_string(eq*100);
	std::string sH = std::to_string(hit*100);
	sEq.erase(sEq.end()-4,sEq.end());
	sH.erase(sH.end()-4,sH.end());
	std::string P = "(E: " + sEq + ", H: " + sH + ")";

	out_str << "\"" + std::to_string(index) + "\" [label=\"\n\n\n\n\n\n" + P +" \", shape=circle, width = 1, regular = 1," +
	" penwidth=" + std::to_string(pw) + ", fontsize=50, style = filled, fillcolor=white," 
	+ "image=\"c" +std::to_string(index) + ".png\"]; \n";
} 

void printGraphEqHit(Graph* g, int source, double* eq, double* F, int reduce) {
	//construct a graphviz file to print the graph - compare the equilibrium probability 
	//with the hitting probability of each state, conditioned on the row (num bonds)

	//set usual input variables
	int clean = 1;

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphvizEqHit.txt";
	std::ofstream out_str(out);

	//probability tolerance for reduced graphs

	//write the graphviz header
	out_str << "graph p {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//print the source node
	double pw = getPW(F[source]);
	printClusterEq(out_str, source, 1, 1, pw);
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double rate = (*g)[source].getEdgeRate(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (reduce == 0 || prob > PTOL) {
			if (!drawn[T]) {
				double pw = getPW(F[T]);
				printClusterEq(out_str, T, eq[T], prob, pw);
				drawn[T] = 1; rankVec.push_back(T);
			}
			double edgeWidth = 40*node_prob*prob;
			if (clean == 0) {
				makeEdge(out_str, source, T, edgeWidth, rate);
			}
			else if (clean == 1) {
				makeEdgeClean(out_str, source, T, edgeWidth);
			}
			toVisit.push_back(T);
		}
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (rate > 0.01 && (reduce == 0 || prob > PTOL)) {
				if (!drawn[T]) {
					double pw = getPW(F[T]);
					printClusterEq(out_str, T, eq[T], (*g)[T].getProb(), pw);
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidth = 40*node_prob*prob;
				if (clean == 0) {
					makeEdge(out_str, node, T, edgeWidth, rate);
				}
				else if (clean == 1) {
					makeEdgeClean(out_str, node, T, edgeWidth);
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}

	//end reached - put a curly to end
	out_str << "}";	

	//free memory
	out_str.close();
	delete []drawn;
}

void printGraphPF(Graph* g, int source, double* Z, int draw, int clean, 
									 int reduce) {
	//construct a graphviz file to print the graph - irreversible with free enenrgy

	int scale = 3;
	int ew_scale = 90;

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphvizPF.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "digraph bd {\n nodesep = 1.5; ranksep = 4; ratio = 0.56\n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//print the source node
	printClusterRev(out_str, source, scale*getPW(Z[source]));
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double rate = (*g)[source].getEdgeRate(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (reduce == 0 || prob > PTOL) {
			if (!drawn[T]) {
				printClusterRev(out_str, T, scale*getPW(Z[T]));
				drawn[T] = 1; rankVec.push_back(T);
			}
			double edgeWidth = ew_scale*node_prob*prob;
			if (clean == 0) {
				makeEdge(out_str, source, T, edgeWidth, rate);
			}
			else if (clean == 1) {
				makeEdgeClean(out_str, source, T, edgeWidth);
			}
			toVisit.push_back(T);
		}
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (rate > 0.01 && (reduce == 0 || prob > PTOL)) {
				if (!drawn[T]) {
					printClusterRev(out_str, T, scale*getPW(Z[T]));
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidth = 3*ew_scale*node_prob*prob;
				if (clean == 0) {
					makeEdge(out_str, node, T, edgeWidth, rate);
				}
				else if (clean == 1) {
					makeEdgeClean(out_str, node, T, edgeWidth);
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}

	//end reached - put a curly to end
	out_str << "}";	

	//free memory
	out_str.close();
	delete []drawn;
}

void makeLegend(std::ofstream& out_str) {
	//make legend for colors

	out_str << "subgraph legend {\n";
	out_str << "label = \"Legend\"; fontsize = 72;\n";
	out_str << "\"AA\" [label=\"AA\", fontsize=72, shape=circle, width = 5, regular = 1, style = filled, fillcolor=blue];\n";
	out_str << "\"AB\" [label=\"AB\", fontsize=72, shape=circle, width = 5, regular = 1, style = filled, fillcolor=yellow];\n";
	out_str << "\"BB\" [label=\"BB\", fontsize=72, shape=circle, width = 5, regular = 1, style = filled, fillcolor=red];\n";
	out_str << "{rank = same; \"AA\";}\n";
	out_str << "\"AA\" -> \"AB\" [style=invis]\n";
	out_str << "{rank = same; \"AB\";}\n";
	out_str << "\"AB\" -> \"BB\" [style=invis]\n";
	out_str << "{rank = same; \"BB\";}\n";
	out_str << "}\n";

}

void printGraphQuenched(Graph* g, int source, int countQ, int layer, std::string color1, std::string color2,
											  std::string types, int draw, int clean, int reduce) {
	//construct a graphviz file to print the graph - for a quenching problem where 

	int scale = 3;
	int ew_scale = 50;

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphvizQ" + std::to_string(countQ) + ".txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "digraph bd {\n nodesep = 1.5; ranksep = 4; ratio = 0.56\n";
	out_str << "edge [ fontsize=48];\n";

	//make the legend and label of the graph
	makeLegend(out_str);
	out_str << " label = \"" + types + " \"; fontsize = 72;\n";



	//print the source node
	printCluster(out_str, source, draw);
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	std::string color = color1;
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double rate = (*g)[source].getEdgeRate(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (reduce == 0 || prob > PTOL) {
			if (!drawn[T]) {
				printCluster(out_str, T, draw);
				drawn[T] = 1; rankVec.push_back(T);
			}
			double edgeWidth = ew_scale*node_prob*prob;
			if (clean == 0) {
				makeEdgeRev(out_str, source, T, edgeWidth, rate, color);
			}
			else if (clean == 1) {
				makeEdgeCleanRev(out_str, source, T, edgeWidth, color);
			}
			toVisit.push_back(T);
		}
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0; int current_layer = 1;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (rate > 0.01 && (reduce == 0 || prob > PTOL)) {
				if (!drawn[T]) {
					printCluster(out_str, T, draw);
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidth = ew_scale*node_prob*prob;
				if (clean == 0) {
					makeEdgeRev(out_str, node, T, edgeWidth, rate, color);
				}
				else if (clean == 1) {
					makeEdgeCleanRev(out_str, node, T, edgeWidth, color);
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			current_layer++;
			if (current_layer == layer) {
				color = color2;
			}
			//std::cout << current_layer << "\n";
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}

	//end reached - put a curly to end
	out_str << "}";	

	//free memory
	out_str.close();
	delete []drawn;
}



}
