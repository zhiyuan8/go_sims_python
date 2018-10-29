#ifndef MAIN_H
#define MAIN_H

#include <chrono>  // for high_resolution_clock
#include <string>

#include "array_size.h"

// readout paths
string path_in;
string path_out;
string path_trials_out;

// all variables used before tests
vector<vector <int> > json;
vector <bool> self_contained; //it is a 0/1 vector, when it is 1, the corresponding vector has at least one rejection
vector <bool> competitive; //it is a 0/1 vector, when it is 1, the corresponding vector has at least one rejection
vector <double> pValue_ttest; // it has the same length as the total columns of pValue
vector <double> statistic_ttest; // it has the same length as the total columns of pValue

// all variables used in self-contained/Sims
vector <double> sims_test_statistic;
vector <bool> sims_result;
vector <bool> sims_result2;

// all variables used in competitive/geometric
vector <int> Bonf_reject;//a vector to temp rejections by mul_Bonf
vector <bool> geo_result;
vector <bool> geo_result2;
vector <double> geometric_p_value;

map<int, int> ID; //create a map
map<int, int>::const_iterator it; // it gives output of map

INPUT input = { 0 /*alpha*/, 0 /*ga_alpha*/,
0/*row*/,   0/*col*/,
0/*mu0*/,    1/*sig0*/,
0/*mu1*/,    1/*sig1*/,
0/*row_reject*/,  10/*seed*/,
1/*a*/, {}, {}, {}};

DECISION decision = {false, false, false, false, false, false};

READIN readin;
OUTPUT sims_bonf;
OUTPUT sims_BH;
OUTPUT hyper_bonf;
OUTPUT hyper_BH;

// definition for all dymanic arrays
double *data1;
double *data2;
double **data1_temp;
double **data2_temp;
double **sims_bonf_node;
double **sims_BH_node;
double **hyper_bonf_node;
double **hyper_BH_node;

// temps json.size() and json[i].size()
int json_vector = 0;
vector <int> json_length;

int position = 0;

// temp process time
double startTime = 0;
double endTime = 0;
double generate_time = 0;
double save_time = 0;
double hyper_time = 0;
double sims_time = 0;

auto start = std::chrono::high_resolution_clock::now();
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;

#endif
