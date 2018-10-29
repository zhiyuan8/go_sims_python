#include <iostream>
#include <fstream>
#include <cmath> // use log
#include <vector>
#include <map>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>  //quantile of chi-square
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>    // std::sort , *std::min_element, std::set_difference , std::max
#include <string>

#include "array_size.h"

using namespace std;
namespace pt = boost::property_tree;

void Initialize_all(READIN &readin, INPUT &input, vector <double> &ttest, vector <double> &statistic_ttest, vector <double> &sims_test, vector <double> &geometric_p,
                    vector <bool> &sims_result, vector <bool> &sims_result2, vector <bool> &geo_result, vector <bool> &geo_result2, int col2)
{
    input.mu1 = readin.eff_size;
    input.alpha = readin.method_alpha;
    input.ga_alpha = readin.ga_alpha;
    if (input.mu0 <= input.mu1)
        input.a = -1; // use area in right tail as p - value
    else
        input.a = 1;  // use area in left tail as p - value
    ttest.resize(input.col);
    statistic_ttest.resize(input.col);
    sims_test.resize(col2);
    sims_result.resize(col2);
    sims_result2.resize(col2);
    geo_result.resize(col2);
    geo_result2.resize(col2);
    geometric_p.resize(col2);
}

double mean(double p[], int n)
{
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
        temp = temp + p[i];
    }
    return( (temp+0.0) / n);

}

double deviation(double p[], double mean, int n)
{
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
        temp = temp + (p[i] - mean)*(p[i] - mean);
    }
    return(sqrt((temp + 0.0) / (n - 1)));

}

double mean_vector(vector <double> data)
{
    int n = data.size();
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
        temp = temp + data[i];
    }
    return( (temp+0.0) / n);

}

double deviation_vector(vector <double> data, double mean)
{
    int n = data.size();
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
        temp = temp + (data[i] - mean)*(data[i] - mean);
    }
    return(sqrt((temp + 0.0) / (n - 1)));
}

int getKeyByValue(map <int, int> ID, int n)
{
    map<int, int>::iterator it;
    for (it = ID.begin(); it != ID.end(); it++)
    {
        if ( it->second == n )// value
        {
            return(it->first);   // key
        }
    }
    return(0);
}

int set_difference_number(vector <int> A, vector <int> B)
{
    int size_A = A.size();
    int size_B = B.size();
    vector<int> result( max(size_A, size_B) );//  must specify length
    vector<int>::iterator it;
    sort(A.begin(), A.end());     //  not address input here, so no worries for changing sequence of vector
    sort(B.begin(), B.end());   //  not address input here, so no worries for changing sequence of vector
    it = set_union(A.begin(), A.end(), B.begin(), B.end(), result.begin());
    return(it - result.begin());
}

int set_intersection_number(vector <int> A, vector <int> B)
{
    int size_A = A.size();
    int size_B = B.size();
    vector<int> result(max(size_A, size_B));//  must specify length
    vector<int>::iterator it;
    sort(A.begin(), A.end());     //  not address input here, so no worries for changing sequence of vector
    sort(B.begin(), B.end());   //  not address input here, so no worries for changing sequence of vector
    it = set_intersection(A.begin(), A.end(), B.begin(), B.end(), result.begin());
    return(it - result.begin());
}

void set_difference_vector(vector <int> &result, vector <int> A, vector <int> B)
{
    int size_A = A.size();
    int size_B = B.size();
    if (int( result.size() ) < max(size_A, size_B))
    {
        result.resize( max(size_A, size_B) );//  must specify length
    }
    vector<int>::iterator it;
    sort(A.begin(), A.end());     //  not address input here, so no worries for changing sequence of vector
    sort(B.begin(), B.end());   //  not address input here, so no worries for changing sequence of vector
    it = set_difference(A.begin(), A.end(), B.begin(), B.end(), result.begin());
    result.resize(it - result.begin());
}

double comb(int n, int k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    double result = n;
    for (int i = 2; i <= k; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

void show_final_result(OUTPUT output, int n_reps, int position, string str)
{
    double FDR = 0;
    double mean_beta = 0;
    for (int i = (position-n_reps); i < position; i++)
    {
        FDR = output.FDP[i] + FDR;
        mean_beta = output.power[i] + mean_beta;
    }
    cout << str << "  FDP= " << FDR/n_reps << " power= " << mean_beta/n_reps << endl;
}

void write_out_node_specific_power(double **matrix, string path, int json_vector, READIN readin)
{
    ofstream c_out;
    c_out.open(path.c_str());
    c_out << " ,";
    for (int i = 0; i<readin.n_regimes; i++)
    {
        c_out << readin.min_n + int((readin.max_n - readin.min_n + 0.0) / (readin.n_regimes - 1 + 0.0) *i) << ',';
    }
    c_out<<endl;
    for (int i = 0; i<json_vector; i++)
    {
        c_out << i <<',';
        for (int j = 0; j < readin.n_regimes; j++)
        {
            c_out << (matrix[i][j]+0.0)/(readin.n_reps+0.0) << ','; // make sure nominator and denominator are double
        }
        c_out << endl;
    }
    c_out.close(); //close the file
}

void write_out_json_bool (vector<bool> output, string path, string test, string adjust, string q)
{
    pt::ptree pt;
    pt::ptree temp;
    int length=output.size();
    int k=0; // temp # of rejections
    for (int i=0; i<length; i++)
    {
        if (output[i])
        {
            pt::ptree cell;
            cell.put_value(i);
            temp.push_back ( make_pair("", cell) );
            k++;
        }
    }

    pt.put ("test", test); // get "test"
    pt.put ("adjust", adjust); // get "num_pvals"
    pt.put ("q_threshold", q); // get "q_threshold"
    pt.put("num_rejections",k);//get "num_rejections"
    pt.add_child("rejections", temp);
    pt::write_json(path, pt);
}

void write_out_json_double (vector<double> output, string path, string name)
{
    pt::ptree pt;
    pt::ptree temp;
    int length=output.size();

    pt.put ("test", name); // get "test"
    //pt.put ("num_pvals", length); // get "num_pvals"
    for (int i=0; i<length; i++)
    {
        pt::ptree cell;
        cell.put_value(output[i]);
        temp.push_back ( make_pair("", cell));
    }
    pt.add_child("pvals",temp); // get "pvals" vector
    pt::write_json(path, pt);
}
