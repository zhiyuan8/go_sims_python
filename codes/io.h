#ifndef IO_H
#define IO_H

#include <vector>
#include <map>
#include <string>

#include "array_size.h" //refer to constant
using namespace std;

void read_in_json(string filename, vector< vector<int> >&json, int &json_vector, vector <int> &json_length);
void read_in_json2(string filename, READIN &readin, DECISION & decision);
void create_ID(vector< vector<int> >&json, map<int, int> &ID, int &index, int json_vector, vector <int> json_length, bool dec);
void read_in_json3(string filename, vector<bool> &json_truth, vector<bool> &json_truth2, vector< vector<int>> json_data,
    map<int, int> ID, INPUT &input,int json_vector, vector <int> json_length, READIN readin);
void create_output(DECISION decision, string path_out, string path_out_trials);
void write_out_python(vector <double> pValue_ttest,vector <double> statistic_ttest, string path);
void write_out_matrix(double ** data1_temp, double ** data2_temp, int row, int col, int rep, string path);
void write_out_trails(vector <bool> result, vector <bool> result2, vector <double> p,
                      DECISION decision, string path, string alpha);
void write_out_results(OUTPUT sims_bonf, OUTPUT sims_BH, OUTPUT hyper_bonf, OUTPUT hyper_BH,
                       double **sims_bonf_node, double **sims_BH_node, double **hyper_bonf_node, double ** hyper_BH_node,
                       int json_vector, READIN readin, int position, string path_out, DECISION decision);
void finish_trial_save(string path);

#endif
