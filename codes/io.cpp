#include <iostream>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <sys/stat.h>
#include <vector>
#include <string>
#include <map>

#include "array_size.h" //refer to constant
#include "utils.h" // use set union

using namespace std;
namespace pt = boost::property_tree;

void read_in_json(string filename, vector< vector<int> >&json, int &json_vector, vector <int> &json_length)
{
    filename = filename + "/meta_gene_ids.json";
    pt::ptree pt;
    pt::read_json(filename, pt);
    for (pt::ptree::value_type &row : pt)// read in # of vector
    {
        int y = 0;// start from first ID in current gene
        for (pt::ptree::value_type &cell : row.second)
        {
            y++; // move to next ID
        }
        json_length.push_back(y); // temp number of ID in one gene
        json_vector++; // move to next vector
    }
    json.resize(json_vector);// read data into json_data 2D vector
    int x = 0;// start from first gene
    for (pt::ptree::value_type &row : pt)
    {
        json[x].resize(json_length[x]);
        int y = 0; // start from first ID in current gene
        for (pt::ptree::value_type &cell : row.second)
        {
            json[x][y] = cell.second.get_value<int>();
            y++;
        }
        x++; // move to next vector
    }
}

void read_in_json2(string filename, READIN &readin, DECISION & decision)
{
    filename = filename +"/meta_restore_params.json";
    string temp;
    pt::ptree pt;
    pt::read_json(filename, pt);
    readin.min_node_size = pt.get<int>("context_params.min_node_size", 0);
    readin.max_node_size = pt.get<int>("context_params.max_node_size", 0);
    readin.ga_alpha = pt.get<double>("test_params.nonnull_params.comp_nonnull.ga",0);
    for (pt::ptree::value_type &alphas : pt.get_child("test_params.method_alpha"))
    {
        readin.method_alpha = alphas.second.get_value<double>();
    }
    readin.n_regimes = pt.get<int>("oneway_params.n_regimes", 0);
    readin.n_reps = pt.get<int>("oneway_params.n_reps", 0);
    readin.min_n = pt.get<int>("oneway_params.min_n", 0);
    readin.max_n = pt.get<int>("oneway_params.max_n", 0);
    readin.eff_size = pt.get<double>("oneway_params.eff_size", 0);
    for ( pt::ptree::value_type &cell : pt.get_child("test_params.method_test") )
    {
        temp = cell.second.get_value<string>();
        if (temp == "simes")
        {
            decision.sims=true;
        }
        if (temp == "hypergeometric.ga")
        {
            decision.hyper=true;
        }
    }
    for ( pt::ptree::value_type &cell : pt.get_child("test_params.method_madj") )
    {
        temp = cell.second.get_value<string>();
        if (temp == "Bonferroni")
        {
            decision.bonf=true;
        }
        if (temp == "BH")
        {
            decision.BH=true;
        }
    }
    for ( pt::ptree::value_type &cell : pt.get_child("test_params.report_metrics") )
    {
        temp = cell.second.get_value<string>();
        if (temp == "NumRej")
        {
            decision.NumRej=true;
        }
    }
}

void create_ID(vector< vector<int> >&json, map<int, int> &ID, int &index, int json_vector, vector <int> json_length, bool dec)
{
    int ID_size = 0;
    for (int i = 0; i < json_vector; i++)
    {
        for (int j = 0; j < json_length[i]; j++)
        {
            ID_size = ID.size(); // temp the size of ID before insert
            ID.insert(pair<int, int>(json[i][j], index));
            if (int(ID.size()) == ID_size + 1) //check if size of ID is +1
            {
                index++;
            }
        }
    }
    if (dec)
    {
        map<int, int>::const_iterator it; // it gives output of map
        for (it = ID.begin(); it != ID.end(); ++it)
        {
            cout << it->first << "=" << it->second << endl;
        }
    }
}

void read_in_json3(string filename, vector<bool> &json_truth, vector<bool> &json_truth2, vector< vector<int>> json_data,
    map<int, int> ID, INPUT &input,int json_vector, vector <int> json_length, READIN readin)
{
    pt::ptree pt;
    filename = filename + "/meta_nonnull_gene_ids.json";
    pt::read_json(filename, pt);
    // To store key of rejection
    for (pt::ptree::value_type &keys : pt)
    {
        input.key_reject.push_back(keys.second.get_value<int>()); // must use push_back
    }
    // According to key, find the index of rejected columns
    int imax = input.key_reject.size();
    input.if_reject.resize(input.col); // initialize size of input.if_reject
    for (int i = 0; i < input.col; i++)
    {
        input.if_reject[i] = false; // initialize all columns as null
    }
    input.col_reject.resize(imax);// quicker to use "resize + loop" than to use push_back
    for (int i = 0; i < imax; i++)
    {
        input.col_reject[i] = ID[input.key_reject[i]];
        input.if_reject[input.col_reject[i]] = true; // some columns have non-null
    }
    // prepare for getting ground_truth vector
    int A1 = imax;
    json_truth.resize(json_vector); // define the length of self_contained ground_truth
    json_truth2.resize(json_vector); // define the length of competitive ground_truth
    double G = input.col; //the number of different elements in json file
    // now get the json_truth (self_contained), it is 1 if there is at least one element from rejection
    for (int j = 0; j < json_vector; j++)
    {
        json_truth[j] = false;
        for (int k = 0; k < json_length[j]; k++) // check each element one by one
        {
            for (int l = 0; l < A1; l++) // compare k th data with l th rejection element
            {
                if (ID[json_data[j][k]] == input.col_reject[l])
                {
                    json_truth[j] = true;
                    break;
                }
            }
            if (json_truth[j] == true)
            {
                break; // I have found a rejection, I don't need to scan the rest elements in this vector
            }
        }
    }
    // now get json_truth2, get 0 by comparing double(|A1^Gi| / |A1|) < double (|Gi|/|G|)
    for (int j = 0; j < json_vector; j++)
    {
        int Gi =json_length[j];
        vector <int> Gi_vector = json_data[j];
        int A1Gi = set_intersection_number(Gi_vector, input.key_reject);
        if (((A1Gi + 0.0) / (A1 + 0.0)) <= ((Gi + 0.0) / (G + 0.0)))
        {
            json_truth2[j] = false;
        }
        else
        {
            json_truth2[j] = true;
        }
    }    
    cout << "There are  " << input.col << " non-repeated index of genes in system" << endl;// check output
    cout << "There are " << json_vector << " genes in system" << endl;// check output
    cout << "n_regimes= " << readin.n_regimes << " n_reps= " << readin.n_reps << endl;// check output
}

void create_output(DECISION decision, string path_out, string path_out_trials)
{
    string command1 = "mkdir -p " + path_out;
    int dir1 = system(command1.c_str());
    if (dir1 == 0)
    {
         cout <<"Create " << path_out << endl;
    }
    if (decision.matrix) //trial folder is created
    {
        string command2 = "mkdir -p " + path_out_trials;
        int dir2 = system(command2.c_str());
        if (dir2 == 0)
        {
             cout <<"Create " << path_out_trials << endl;
        }
    }
}

void write_out_trails(vector <bool> result, vector <bool> result2, vector <double> p, DECISION decision, string path, string alpha)
{
    if (decision.sims)
    {
        string value_path = path +"/node_pvals_simes.json";
        write_out_json_double(p,value_path);
        if (decision.bonf)
        {
            string bonf_path = path +"/node_rej_"+alpha+"_simes_Bonferroni.json";
            write_out_json_bool(result,bonf_path);
        }
        if (decision.BH)
        {
            string BH_path = path +"/node_rej_"+alpha+"_simes_BH.json";
            write_out_json_bool(result2,BH_path);
        }
    }
    if (decision.hyper)
    {
        string value_path = path +"/node_pvals_hypergeometric.ga.json";
        write_out_json_double(p,value_path);
        if (decision.bonf)
        {
            string bonf_path = path +"/node_rej_"+alpha+"_hypergeometric.ga_Bonferroni.json";
            write_out_json_bool(result,bonf_path);
        }
        if (decision.BH)
        {
            string BH_path = path +"/node_rej_"+alpha+"_hypergeometric.ga_BH.json";
            write_out_json_bool(result2,BH_path);
        }
    }
} 

void write_out_results(OUTPUT sims_bonf, OUTPUT sims_BH, OUTPUT hyper_bonf, OUTPUT hyper_BH, 
                       double **sims_bonf_node, double **sims_BH_node, double **hyper_bonf_node, double ** hyper_BH_node,
                       int json_vector, READIN readin, int position, string path_out, DECISION decision)
{
    ofstream c_out;
    string summary = path_out + "/trial_summary.csv";
    c_out.open(summary.c_str());
    c_out<< " "<<','<<"adjustment_method,"<<"empirical_fdr,"<<"num_rejections,"<<"empirical_power,"
         <<"testing_method,"<<"nonnull_type,"<<"regime_id,"<<"repetition_id,"<<"trial_id"<<endl;
    int regime_id=0;
	int repetition_id=0;
    int k=0;
    for (int i = 0; i<position; i++)
    {
        regime_id = i / readin.n_reps;
        repetition_id = i % readin.n_reps;
        if (decision.sims && decision.BH)
        {
            c_out << k << ", BH , " << sims_BH.FDP[i]<<','<< sims_BH.n_rejection[i]<<','<<sims_BH.power[i]<<" , simes,self"<<','
                  <<regime_id<<','<<repetition_id<<','<<i<<endl;
            k++;
        }
        if (decision.sims && decision.bonf)
        {
            c_out << k << ",Bonferroni," << sims_bonf.FDP[i]<<','<< sims_bonf.n_rejection[i]<<','<<sims_bonf.power[i]<<" , simes,self"<<','
                  <<regime_id<<','<<repetition_id<<','<<i<<endl;
            k++;
        }
        if (decision.hyper && decision.BH)
        {
            c_out << k << ",BH," << hyper_BH.FDP[i]<<','<< hyper_BH.n_rejection[i]<<','<<hyper_BH.power[i]<<" , hypergeometric.ga,comp"<<','
                  <<regime_id<<','<<repetition_id<<','<<i<<endl;
            k++;
        }
        if (decision.hyper && decision.bonf)
        {
            c_out << k << ",Bonferroni," << hyper_bonf.FDP[i]<<','<< hyper_bonf.n_rejection[i]<<','<<hyper_bonf.power[i]<<" , hypergeometric.ga,comp"<<','
                  <<regime_id<<','<<repetition_id<<','<<i<<endl;
            k++;
        }
    }
    c_out.close(); //close the file

    if (decision.sims && decision.BH)
    {
        string sims_BH_node_path = path_out + "/node_simes_BH.csv";
        write_out_node_specific_power(sims_BH_node,sims_BH_node_path,json_vector,readin);
    }
    if (decision.sims && decision.bonf)
    {
        string sims_bonf_node_path = path_out + "/node_simes_Bonferroni.csv";
        write_out_node_specific_power(sims_bonf_node,sims_bonf_node_path,json_vector,readin);
    }
    if (decision.hyper && decision.bonf)
    {
        string hyper_bonf_node_path = path_out + "/node_hypergeometric.ga_Bonferroni.csv";
        write_out_node_specific_power(hyper_bonf_node,hyper_bonf_node_path,json_vector,readin);
    }
    if (decision.hyper && decision.BH)
    {
        string hyper_BH_node_path = path_out + "/node_hypergeometric.ga_BH.csv";
        write_out_node_specific_power(hyper_BH_node,hyper_BH_node_path,json_vector,readin);
    }
}
