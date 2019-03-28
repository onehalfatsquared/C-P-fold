#pragma once
namespace bd { 
class Database;

//general functions
//perform tpt on initial and target states //output = ?
void performTPT(int N, int initial, int target, Database* db);

//transition matrix functions
//fill in rate matrix with data from mfpt estimates
void createTransitionMatrix(double* T, int num_states, Database* db);
//edit rate matrix to satisfy detailed balance
void satisfyDB(double* T, int num_states, Database* db);
//fill in diagonal of rate matrix
void fillDiag(double* T, int num_states);

//tpt functions
void computeCommittor(double* q, double* T, int num_states, int initial, int* target);























}