#include <fstream>
#include "CPhys.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace CPhys;

int main(int argc, const char *argv[])
{	
	int 	walkers 	= 10000;
	double	walkProb 	= 0.5;
	double	T		= 0.01;
	double	dt		= 0.0001;
	double	minDist		= 0;
	double	maxDist		= 1;
	Walker walker = 
		MonteCarlo::getWalker(walkers,walkProb,T,dt,minDist,maxDist);
	walker.useGauss(true);
	walker.start1D();
	walker.output("1DT0_01Gauss");
}

