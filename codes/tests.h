#ifndef TESTS_H
#define TESTS_H

#include "array_size.h"

double t_test(double temp_1[], double temp_2[], int &Sn1, int &Sn2);
void mul_bonf(vector <double> pValue_ttest, vector <int> &Bonf_reject, INPUT input, map <int, int> ID, int dec);
void gl_geometric(vector<double> &result, vector<int> D, vector<vector <int> > json_data, int N, int dec);
void gl_sims(vector <double> pValue_ttest, vector <double> &sims_result, vector< vector<int>> json_data, map <int, int> ID, int dec);
void mul_BH(vector <bool> &result, vector <double> p, double alpha, double **node_power, int currentCol, int dec);
void mul_bonf2(vector <bool> &result, vector <double> p, double alpha, double **node_power, int currentCol, int dec);
void confusion(OUTPUT &output, vector <bool> result, vector <bool> truth, int regime, int experiment, int position, int dec);

#endif
