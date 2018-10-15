# go_sim
The codes is written by C++ but is wrapped by python.
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

## Running the tests
open ```wrapper.py``` file. 'seed' refers to random value generation and you can put any value. 'output' reflects your decison on whether all details of trials will be output. If 'output=1', then a ```trails``` folder will be generated and all details will be written in json files. 'boost_path' is path where you download boost library files and extract the compressed bag. The path within go_dag() refers to the relative path for readin data.
Finally, run this line and you will begin your test.
```
python wrapper.py
```

## Related publication


## Authors

* **Zhiyuan Jim Li** - *Initial work* - contact: zhiyuan.li1995@hotmail.com


