#include <iostream>
#include <time.h> //show calculation time
#include <chrono>  // for high_resolution_clock
#include <map>
#include <vector>
#include <string>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "io.h"
#include "array_size.h"
#include "main.h"
#include "generate_data.h"
#include "utils.h"
#include "tests.h"

using namespace std;

int main(int argc, char* argv[])
//int main()
{
    cout << "The path to read in json files is:  " << argv[1] << "\n";
    path_in = argv[1]; // get path where we readin data
    cout << "The seed used in generating data is:  " << abs ( atoi( argv[2]) ) << "\n";
    input.seed = abs ( atoi( argv[2]) );
    if (atoi( argv[3])==1)
    {
        cout << "The trail folder will be generated and all intermediate results will be output" << "\n";
        decision.matrix=true;
    }
    path_out = path_in + "/summary";
    path_trials_out = path_in + "/trials";

    /*path_in = "/home/zhiyuan/Desktop/go_sim_for_python/simulation_output/case_heart-effect_0.5";
    input.seed = 100;*/

    // 1. read in data, intialize all vectors
    startTime = clock(); // begin documenting time
    read_in_json(path_in, json, json_vector, json_length); // get the ID data in json file
    read_in_json2(path_in, readin, decision); // read parameters for simulations
    if ( readin.n_reps  >100 )// if data is too large, stop program
    {
        cout<<"please do not output trail details, or make 'n_reps' <= 100" << endl;
        return 0;
    }
    create_ID(json, ID, input.col, json_vector, json_length, false);// build the map. If true, then show me all ID
    Initialize_all(readin, input, pValue_ttest, sims_test_statistic, geometric_p_value, sims_result, sims_result2, geo_result, geo_result2, json_vector);// initialize all vectors
    read_in_json3(path_in, self_contained, competitive, json, ID, input, json_vector, json_length, readin);// get self_contained and competitive non-null
    create_output(decision,path_out,path_trials_out); // create two folders, and create first line for files in summary folder

    // 2. Initialize all dynamic arrays
    sims_bonf.FDP = new double [readin.n_regimes*readin.n_reps];
    sims_bonf.power = new double [readin.n_regimes*readin.n_reps];
    sims_BH.FDP = new double [readin.n_regimes*readin.n_reps];
    sims_BH.power = new double [readin.n_regimes*readin.n_reps];
    hyper_bonf.FDP = new double [readin.n_regimes*readin.n_reps];
    hyper_bonf.power = new double [readin.n_regimes*readin.n_reps];
    hyper_BH.FDP = new double [readin.n_regimes*readin.n_reps];
    hyper_BH.power = new double [readin.n_regimes*readin.n_reps];
    sims_bonf.n_rejection = new int [readin.n_regimes*readin.n_reps];
    sims_BH.n_rejection = new int [readin.n_regimes*readin.n_reps];
    hyper_bonf.n_rejection = new int [readin.n_regimes*readin.n_reps];
    hyper_BH.n_rejection = new int [readin.n_regimes*readin.n_reps];
    data1 = new double [readin.max_n];
    data2 = new double [readin.max_n];// normal random variable 1D array for one row of control sample
    sims_bonf_node = new double *[json_vector];
    sims_BH_node = new double *[json_vector];// node specific power for all regimes
    hyper_bonf_node = new double *[json_vector];
    hyper_BH_node = new double *[json_vector];// node specific power for all regimes
    for (int i=0; i<json_vector; i++)
    {
        sims_bonf_node[i]=new double[readin.n_regimes];
        sims_BH_node[i]=new double[readin.n_regimes];//node specific power size = (json_vector * readin.n_regimes)
        hyper_bonf_node[i]=new double[readin.n_regimes];
        hyper_BH_node[i]=new double[readin.n_regimes];
        for (int j=0; j<readin.n_regimes; j++)
        {
            sims_bonf_node[i][j]=sims_BH_node[i][j]=hyper_bonf_node[i][j]=hyper_BH_node[i][j]=0; // node specific power needs to be 0 because it will be averaged
        }
    }

    // 3. do simulations and tests
    for (int i = 0; i < readin.n_regimes; i++)
    {
        input.row_reject = readin.min_n + int((readin.max_n - readin.min_n + 0.0) / (readin.n_regimes - 1 + 0.0) *i); // linspace
        cout << "\n" << "regime " << i + 1 << endl; // show current progress
        for (int j = 0; j < readin.n_reps; j++) // repeat tests
        {
            position++;
            input.seed = input.seed + (j + 1) * 100; // use different seed in each loop
            start = std::chrono::high_resolution_clock::now();
            generateN(data1,data2,input,pValue_ttest); // generate a column of p value and do t-test to get p value
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish - start;
            generate_time = generate_time + elapsed.count();

            // create "trialN" folder if required
            string command = "mkdir -p " + path_trials_out + "/trail_" +to_string(position) ; // path to "trialN" folder
            string path_trials_out_N = path_trials_out + "/trail_" +to_string(position);
            stringstream ss;
            ss << fixed << setprecision(2) << input.alpha; //print to string with specific precision
            string alpha_string = ss.str(); // get the '0.1' used in json file names
            if (decision.matrix)
            {
                system(command.c_str()); // generate "trialN" folder
            }

            // hypergeometric test for the competitive null
            if (decision.hyper) // do hypergeometric tests
            {
                start = std::chrono::high_resolution_clock::now();
                mul_bonf(pValue_ttest, Bonf_reject, input, ID, 0); // use  input.ga_alpha rather than alpha/n
                gl_geometric(geometric_p_value, Bonf_reject, json, input.col, 0); //fisher's exact test. if 1, show  result
                if (decision.BH) // choose bonferroni adjustment
                {
                     mul_BH(geo_result, geometric_p_value, input.alpha, hyper_BH_node, i, 0);// BH multiple testing. if 1, show result.
                     confusion(hyper_BH, geo_result, competitive, input.row_reject, j, position, 0);// get FDP & beta of hypergeometric test.  if 1, show result.
                }
                if (decision.bonf) // choose BH adjustment, the index always has '+readin.n_regimes*readin.n_reps'
                {
                     mul_bonf2(geo_result2, geometric_p_value, (input.alpha/json_vector), hyper_bonf_node, i, 0);// Bonferroni multiple testing. if 1, show result.
                     confusion(hyper_bonf, geo_result2, competitive, input.row_reject, j, position, 0);// get FDP & beta of hypergeometric test.  if 1, show result.
                }
                if (decision.matrix)
                {
                   write_out_trails(geo_result, geo_result2, geometric_p_value, decision, path_trials_out_N, alpha_string); //write out all details in 'trials' folder if required
                }

                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;
                hyper_time = hyper_time + elapsed.count();
            }

            // sims test for the self-contained null
            if (decision.sims) // do sims tests
            {
                start = std::chrono::high_resolution_clock::now();
                gl_sims(pValue_ttest, sims_test_statistic, json, ID, 0); // get the test statistic from sims tests
                if (decision.BH) // choose bonferroni adjustment
                {
                    mul_BH(sims_result, sims_test_statistic, input.alpha, sims_BH_node, i, 0);// BH multiple testing. if 1, show result.
                    confusion(sims_BH, sims_result, self_contained, input.row_reject, j, position, 0);// get FDP & beta of hypergeometric test.  if 1, show result.
                }
                if (decision.bonf) // choose bonferroni adjustment
                {
                    mul_bonf2(sims_result2, sims_test_statistic, (input.alpha/json_vector), sims_bonf_node, i, 0);// Bonferroni multiple testing. if 1, show result.
                    confusion(sims_bonf, sims_result2, self_contained, input.row_reject, j, position, 0);// get FDP & beta of hypergeometric test.  if 1, show result.
                }
                if (decision.matrix)
                {
                    write_out_trails(sims_result, sims_result2, sims_test_statistic, decision, path_trials_out_N,alpha_string); //write out all details in 'trials' folder if required
                }
                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;
                sims_time = sims_time + elapsed.count();
            }
        }
        //show averaged result for one regime
        show_final_result(sims_bonf, readin.n_reps, position,"sims test + bonferroni adjustment");
        show_final_result(sims_BH, readin.n_reps, position, "sims test + BH adjustment");
        show_final_result(hyper_bonf, readin.n_reps, position,"hypergeometric test + bonferroni adjustment");
        show_final_result(hyper_BH, readin.n_reps, position,"hypergeometric test + BH adjustment");
    }
    //write out power, FDP and node-specific power
    write_out_results(sims_bonf,sims_BH, hyper_bonf, hyper_BH, sims_bonf_node, sims_BH_node, hyper_bonf_node, hyper_BH_node,
                      json_vector,readin, position, path_out, decision);
    //4. delete dynamic array and show total time, I have to write separate delete [] expressions
    endTime = clock();
    cout << "generate_time + t_test_time : " << generate_time <<endl;
    cout << "save_time : " << save_time <<endl;
    cout << "hyper_time : " << hyper_time << endl;
    cout << "sims_time : " << sims_time << endl;
    cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    for (int i=0; i<json_vector; i++)
    {
        delete [] sims_BH_node[i];
        delete [] sims_bonf_node[i];
        delete [] hyper_BH_node[i];
        delete []hyper_bonf_node[i];
    }
    delete [] sims_BH_node;
    delete [] sims_bonf_node;
    delete [] hyper_BH_node;
    delete [] hyper_bonf_node;
    delete [] data1;
    delete [] data2;
    delete [] sims_BH.FDP;
    delete [] sims_BH.power;
    delete [] sims_BH.n_rejection;
    delete [] sims_bonf.FDP;
    delete [] sims_bonf.power;
    delete [] sims_bonf.n_rejection;
    delete [] hyper_BH.FDP;
    delete [] hyper_BH.power;
    delete [] hyper_BH.n_rejection;
    delete [] hyper_bonf.FDP;
    delete [] hyper_bonf.power;
    delete [] hyper_bonf.n_rejection;

    return 0;
}
