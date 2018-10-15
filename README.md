# go_sim

Self-contained and competitive tests on simulated data. The two tests are used to detect non-null gene.

### System Requirement
g++ >= 4.2.1. (double check it! it is said g++ contains C++11 lib)
https://stackoverflow.com/questions/7482026/can-i-use-the-latest-features-of-c11-in-xcode-4-or-osx-lion

The following codes help to check version
```
g++ --version
```

### Before Installation

Boost library is used in normal random variable generation and chi-square quantile calculation. In terminal, firstly check if there exists boost library. Run the following codes:
```
locate /boost/version.hpp
```
If there exits boost library, then the path will be printed. For example, "/home/zhiyuan/Downloads/boost_1_68_0/boost/version.hpp". Define the path as "path to boost", then this line of code enables you to see your boost library version.
```
cat *path to boost* | grep "BOOST_LIB_VERSION"
```
If your boost library is installed, then you can skip this section. 
If boost library is not found, and you have the root access, you can easily install it via following commands in terminal:
```
sudo apt-get install libboost-all-dev
```
If boost library is not found, and you don't have the root access, you can find install instructions at http://masumhabib.com/blog/how-to-install-the-boost-library-with-mpi-and-without-root-access/

### Installation 

If you are running it on your PC, then use ```cd``` to locate at the place where you will store the codes. Download the codes on github.
```
git clone --https://github.com/zhiyuan8/go_sim.git
```
If you are running it on a remote server, download the codes on github.
```
git clone git@github.com:zhiyuan8/go_sim.git 
```
Next, create ```result``` folder and go to ```codes``` folder.
```
mkdir result
cd codes
```
complie all cpp files.
```
g++ main.cpp generate_data.cpp io.cpp tests.cpp utils.cpp -o GO_DAG -std=c++11
```
In this way, ```GO_DAG``` is created. 	This process takes a while.

## Customize your tests
The ```test_data``` folder contains datasets as a demo. You can replace them with your json files, but remember to keep the file name as the same. 
Open the file ```parameter.json```, in "report_metrics", there are four choices (1)"FDP"(2)"Power"(3)"NumRej"(4)"Data". They refer to (1) False Discovery Rate of two tests (2) Power of two tests (3) Node-specific power (4) original normal random variables in Self-contained and competitive tests. (1) and (2) are output mandatorily, and (3) and (4) are output if you type "NumRej" and "Data" in that json file.
Also in the file ```parameter.json```, you can change "method_alpha" to other criteria, such as 0.05, and you can change  "eff_size" which means the effective size of non-null random variable ( For example, eff_size=0.5, non-null data ~ N(0.5, 1) ).
In  "oneway_params", "n_rep" means number of repeated experiments in each regime. "n_regimes" means number of regimes. For example,  {"n_regimes": 12, "n_reps": 100, "min_n": 10, "max_n": 120} means "There are 12 regimes, each one will be repeated 100 times with different randomness, and in each regime, number of cases = number of control = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120}"

## Running the tests

Finally, input codes in terminal.
```
./codes/GO_DAG ./simulation_output/case_heart-effect_0.5 100 1
```
The first argument enables terminal to run the code. The second and third arguments are paths for reading in and writing out files. The fourth number is a seed used in generating data. If you keep the seed number unchanged, you will get same outputs.
noted, if you run my codes for a second time, the output in your first test will not be covered, and new output will be written at the end of previous results. So I suggest you copy the output .csv files in another path and delete them.

## Related publication


## Authors

* **Zhiyuan Jim Li** - *Initial work* - contact: zhiyuan.li1995@hotmail.com


