#include <fstream>
#include <iostream>
#include <cmath>
#include "CPhys.h"

inline double exact(double x, double t, int N){
	double pi = 3.14159;
	double sum = 0;
	for (int n = 1; n < N; n++) {
		sum += (1/n)*sin(n*pi*x)*exp(-n*n*pi*pi*t);
	}
	sum = 2*sum/pi;
	return 1 - x - sum;
}

using namespace std;
using namespace CPhys::Diffusion;

int main(int argc, const char *argv[])
{
	char	output[] = "T5x01t00005.txt";
	double 	dx 	= 0.01;
	double 	dt 	= 0.00005;
	double 	time 	= 0.50000;
	double 	alpha 	= dt/(dx*dx);
	int 	N 	= (1/dx)+1;
	Vector init1 	= Vector(N);
	Vector init2 	= Vector(N);
	Vector init3 	= Vector(N);
	Vector uEx	= Vector(N);
	// Inital values
	for (int i = 0; i < N; i++) {
		init1(i) = -1+dx*i;
		init2(i) = -1+dx*i;
		init3(i) = -1+dx*i;
		uEx(i)  = exact(dx*i,time,100);
	}
	init1(0) = 0;
	init2(0) = 0;
	init3(0) = 0;
	init1(N-1) = 0;
	init2(N-1) = 0;
	init3(N-1) = 0;

	Vector u1 = solveEuler(init1,time,dt,dx);
	Vector u2 = solveBackwardEuler(init2,time,dt,dx);
	Vector u3 = solveCrankNicolson(init3,time,dt,dx);

	for (int i = 0; i < N; i++) {
		u1(i) += 1 - dx*i;
		u2(i) += 1 - dx*i;
		u3(i) += 1 - dx*i;
	}

	// Write to file
	std::ofstream myfile;
	myfile.open(output);
	for (int i = 0; i < N; i++) {
		myfile << i*dx << " ";
		myfile << uEx(i) << " ";
		myfile << u1(i) << " ";
		myfile << u2(i) << " ";
		myfile << u3(i) << " ";
		myfile << endl;
	}
	myfile.close();
	return 0;
}

