#include <iostream>
#include <exception>
#include <complex>
#include<string>
#include <cmath>
#include <algorithm>
#include<fstream>
#include<istream>
#include <string>
#include <Armadillo>
#include"ReadParam.h"

#include"sgsim.h"


using namespace std;
using namespace arma;

int main() {


	std::string fp = "sgsim.par";
	Param sgs_param;
	read_param(sgs_param, fp);


	sgsim SGS(sgs_param);
	SGS.simulate(sgs_param);

	arma::vec sim;
	sim.load("sim.bin");
	sim.save("simulation.asc", arma_ascii);
	




	
	return 0;
	system("pause");
	// define the data headers here

};



/*

sim -  starts from one
new_ind -  from 0
order  - from 0

all data starts from 0 - 
var, vr, vrgtr, lvm

ixnode, iynode,iznode,icnode starts from 1 

*/