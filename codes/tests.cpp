#include <iostream>
#include <cmath> // use log
#include <vector>
#include <map>
#include <boost/math/distributions/chi_squared.hpp>  //quantile of chi-square
#include <boost/math/distributions/students_t.hpp>
#include <algorithm>    // std::sort , *std::min_element, , std::max

#include "array_size.h"
#include "utils.h"
using namespace std;

void t_test(double temp_1[], double temp_2[], int &Sn1, int &Sn2, int index,
              vector <double> &pValue_ttest, vector <double> &statistic_ttest)
{
    // declare parameters used in 2-sample t-test
    double Sm1;       // Sm1 = Sample 1 Mean.
    double Sd1;       // Sd1 = Sample 1 Standard Deviation.
    double Sm2;       // Sm2 = Sample 2 Mean.
    double Sd2;       // Sd2 = Sample 2 Standard Deviation.
    double v;         // Degrees of freedom.
    double t1;
    double t2;
    // All parameters used in 2-sample t-test
    Sm1 = mean(temp_1, Sn1);
    Sm2 = mean(temp_2, Sn2);
    Sd1 = deviation(temp_1, Sm1, Sn1);
    Sd2 = deviation(temp_2, Sm2, Sn2);
    v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;// Degrees of freedom
    v *= v;
    t1 = Sd1 * Sd1 / Sn1;
    t1 *= t1;
    t1 /= (Sn1 - 1);
    t2 = Sd2 * Sd2 / Sn2;
    t2 *= t2;
    t2 /= (Sn2 - 1);
    v /= (t1 + t2);
    statistic_ttest[index] = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);// t-statistic
    boost::math::students_t dist(v);
    pValue_ttest[index] =  2*cdf(complement(dist, fabs(statistic_ttest[index])));// return two-sided p value
}

void mul_bonf(vector <double> pValue_ttest, vector <int> &Bonf_reject, INPUT input, map <int, int> ID, int dec)
{
    Bonf_reject.clear();// Because I use push_back, clear the data first
    for (int l = 0; l < input.col; l++) // get the index of rejections after mul_Bonf
    {
        if (pValue_ttest[l] <= input.ga_alpha) //reject if p[i]<= alpha/n
        {
            Bonf_reject.push_back(getKeyByValue(ID, l)); //Bonf_reject temps the ID of rejected columns, not index
        }
    }
    if (dec == 1) // check Bonf_reject
    {
        cout << "Bonf_rejected" << "size=" << Bonf_reject.size() << " ID: ";
        int maxx = Bonf_reject.size();
        for (int x = 0; x < maxx; x++)
        {
            cout << Bonf_reject[x] << ' ';
        }
        cout << endl;
    }
}

void gl_geometric(vector<double> &result, vector<int> D, vector<vector <int> > json_data, int N, int dec)
{
    // scan through all vectors in json_data
    int maxi = json_data.size();
    int X0 = 0;
    int n = D.size();
    int K = 0;
    int X = 0;
    //cout << " n= " << n << " N= " << N << endl;
    for (int i = 0; i < maxi; i++)
    {
        double temp = 0; // temps p value[i]
                         // extract one line of json data
        vector <int> Gi = json_data[i];
        // get n, N, K, X for geometric sampling
        X0 = set_intersection_number(D, Gi);
        K = Gi.size();
        X = std::max(n + K - N, X0); // make sure max(0,n+K-N)=<X<=min(K,N) && X>=X0
                                     //cout << " X0= " <<X0 << " K= " << K << endl;
                                     // now get p value of hypergeometric term
                                     // p = Pr(X>=X0, n=|D|, N=|G|, K=|Gi|), X={X0,...,K}
                                     // Pr(X)=comb(K,X)*comb(N-K,n-X)/comb(N,n)
        while (X <= std::min(n, K))
        {
            temp = temp + (comb(K, X) + 0.0) *(comb(N - K, n - X) + 0.0) / (comb(N, n) + 0.0);
            X = X + 1;
        }
        result[i] = temp;
    }
    if (dec == 1)
    {
        cout << "geometric_p_value " << ' ';
        for (int x = 0; x < maxi; x++)
        {
            cout << result[x] << ' ';
        }
        cout << endl;
    }
}

void gl_sims(vector <double> pValue_ttest, vector <double> &sims_result, vector< vector<int>> json_data, map <int, int> ID, int dec)
{
    int maxj = json_data.size();
    vector <double> temp_p;
    for (int j = 0; j < maxj; j++)// loop over all vectors in json
    {
        int maxk = json_data[j].size();
        temp_p.resize(maxk); // resize temp_p and also initialize it
        for (int k = 0; k < maxk; k++)
        {
            temp_p[k] = pValue_ttest[ID[json_data[j][k]]];
        }
        // remeber to sort p-value, which is VERY important!
        sort(temp_p.begin(), temp_p.end()); //not sorted!
        //pass the columns correponding to pValue from 2-sample t-test to a new vector, then get p[i]*n/i
        for (int k = 0; k < maxk; k++)
        {
            temp_p[k] = temp_p[k] * (maxk + 0.0) / (k + 1.0); // p[i]*n/i
        }
        sims_result[j] = *min_element(temp_p.begin(), temp_p.end());
    }
    if (dec == 1)
    {
        cout << "sims_test_statistic  ";
        for (int x = 0; x < maxj; x++)
        {
            cout << sims_result[x] << ' ';
        }
        cout << endl;
    }
}

void mul_BH(vector <bool> &result, vector <double> p, double alpha, double **node_power, int currentCol, int dec)
{
    // remeber to return the unsorted 1/0 array!
    int threshold = 0;
    double criteria = 0;
    int col = p.size();
    vector <double> a = p;
    sort(a.begin(), a.end()); // sort a here
    for (int j = col; j > 0; j--) // the actual index is should -1
    {
        if (a[j - 1] <= j*alpha / col)
        {
            threshold = j;
            break;
        } //find the largest j, save it in threshold, and break out the loop
    }
    if (threshold > 0) // j exists, return criteria
    {
        criteria = a[threshold - 1];
    }
    else
    {
        criteria = 0;//if no such j exists, there is no rejections
    }
    // now get the 0/1 vector
    for (int i = 0; i < col; i++)
    {
        if (p[i] <= criteria)
        {
            result[i] = true; // rejected
            node_power[i][currentCol]=node_power[i][currentCol]+1;
        }
        else
        {
            result[i] = false; // not rejected
        }
    }
    if (dec == 1)
    {
        cout << "BH multiple test result " << ' ';
        for (int x = 0; x < col; x++)
        {
            cout << result[x] << ' ';
        }
        cout << endl;
    }
}

void mul_bonf2(vector <bool> &result, vector <double> p, double alpha, double **node_power, int currentCol, int dec)
{
    int col = p.size();
    // now get the 0/1 vector
    for (int i = 0; i < col; i++)
    {
        if (p[i] <= alpha ) // alpha = input.alpha / length(p)
        {
            result[i] = true; // rejected
            node_power[i][currentCol]=node_power[i][currentCol]+1;
        }
        else
        {
            result[i] = false; // not rejected
        }
    }
    if (dec == 1)
    {
        cout << "Bonferroni multiple test result " << ' ';
        for (int x = 0; x < col; x++)
        {
            cout << result[x] << ' ';
        }
        cout << endl;
    }
}


void confusion(OUTPUT & output, vector <bool> result, vector <bool> truth, int regime, int experiment, int position, int dec)
{
    int col = result.size();
    position = position-1; // the index starts from 0
    double temp[4] = { 0 }; //temp confusion matrix
    for (int j = 0; j < col; j++) // compare element in one row
    {
        if ((truth[j] == false) && (result[j] == false))
        {
            temp[0] = temp[0] + 1; // assume = 0 , result=0
        }
        if ((truth[j] == true) && (result[j] == false))
        {
            temp[1] = temp[1] + 1; // assume = 1 , result=0
        }
        if ((truth[j] == false) && (result[j] == true))
        {
            temp[2] = temp[2] + 1; // assume = 0 , result=1
        }
        if ((truth[j] == true) && (result[j] == true))
        {
            temp[3] = temp[3] + 1; // assume = 1 , result=1
        }
    }
    output.n_rejection[position] = int (temp[2] + temp[3]);
    if ((temp[2] + temp[3]) >= 1)// V+S>=1
    {
        output.FDP[position] = (temp[2] + 0.0) / (temp[2] + temp[3]);//FDR
    }
    else
    {
        output.FDP[position] = temp[2];//FDR
    }
    if ((temp[1] + temp[3]) >= 1) // T+S>=1
    {
        output.power[position] = (temp[3] + 0.0) / (temp[1] + temp[3]); // power
    }
    else // T=0, S=0
    {
        output.power[position] = 1; //power
    }
    if (dec == 1) //if 1, show result.
    {
        cout << " confusion matrix " << temp[0] << ' ' << temp[1] << ' ' << temp[2] << ' ' << temp[3] << endl;
        cout << "n_control " << regime << " experiment " << experiment+1
            << " FDP " << output.FDP[position] << " power " << output.power[position] << endl;
    }
}
