#include <iostream>
#include <random>  // default_random_engine generator
#include <boost/math/distributions/normal.hpp> // generate normal random variable
#include <boost/math/distributions/students_t.hpp> // get CDF of student t distribution
#include <vector>

#include "array_size.h"
#include "utils.h"
#include "tests.h"

using namespace std;

void generateN(double *data1, double *data2, INPUT input, vector <double> &pValue_ttest)
{
    seed_seq seq{ input.seed };
    mt19937 generator(seq);
    uniform_real_distribution<double> uniform_distribution(0, 1); // p value under H0 is uniform
    normal_distribution<double> distribution0(0, 1); // H0 distribution
    normal_distribution<double> distribution1(input.mu1, input.sig1); // H1 distribution
    for (int j = 0; j < input.col; j++) // generate p value column by column, then get 2-sample t-test value
    {
        if (input.if_reject[j]) // this column comes from non-null
        {
            for (int i = 0; i < input.row_reject; i++)
            {
                data1[i] = distribution0(generator); // top half from H0
                data2[i] = distribution1(generator); // bottom half from H1
            }
            pValue_ttest[j] = t_test(data1, data2, input.row_reject, input.row_reject);
        }
        else // this column comes from null
        {
            pValue_ttest[j] = uniform_distribution(generator); // t-test p value is uniform distributed
        }
    }
}

