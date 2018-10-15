#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <map>
#include <string>

#include "array_size.h" //refer to constant
using namespace std;

void Initialize_all(READIN &readin, INPUT &input, vector <double> &ttest, vector <double> &sims_test, vector <double> &geometric_p,
                    vector <bool> &sims_result, vector <bool> &sims_result2, vector <bool> &geo_result, vector <bool> &geo_result2, int col2);
int getKeyByValue(map <int, int> ID, int n);
double mean(double p[], int n);
double mean_vector(vector <double> data);
double deviation(double p[], double mean, int n);
double deviation_vector(vector <double> data, double average);

// functions to get set intersection, set difference and combination number
int set_difference_number(vector <int> A, vector <int> B);
int set_intersection_number(vector <int> A, vector <int> B);
void set_difference_vector(vector <int> &result, vector <int> A, vector <int> B);
double comb(int n, int k);

void show_final_result(OUTPUT output, int n_reps, int position, string str);
void write_out_node_specific_power(double **matrix, string path, int json_vector, READIN readin);
void write_out_json_bool (vector<bool> output, string path);
void write_out_json_double (vector<double> output, string path);
#endif
