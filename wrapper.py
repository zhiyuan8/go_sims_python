import os

# This function is a python wrapper for the C++ implementaion of GO-DAG
# Usage:
# import go_dag from wrapper
# go_dag(data_filepath, result_filepath, seed, boost_path)
#	data_filepath: the relative path to where you store the json data file
#	result_filepath: the relative path to output your result, which is a collection of csv files
#	seed: random variable seed
#	boost_path: boost library installation path, made visible to the system

def go_dag(filepath, seed=100, output =0, boost_path='/home/zhiyuan/Downloads/boost_1_68_0'):
	os.system("cd codes; g++ -I {} main.cpp generate_data.cpp io.cpp tests.cpp utils.cpp -o GO_DAG -std=c++11; cd -".format(boost_path))
	os.system("./codes/GO_DAG {} {} {}".format(filepath, seed, output))
go_dag('./simulation_output/case_heart-effect_0.5')
