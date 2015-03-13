#include <fstream>
#include "CPhys.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace CPhys::LinAlg;

inline double srcFunction(double x){
	return 100*exp(-10*x);
}

inline double solution(double x){
	return 1.-(1.-exp(-10))*x-exp(-10*x);
}

int main(int argc, const char *argv[])
{
	// We solve the equation Ay(x) = f(x)
	// We have the exact solution As(x) = f(x)
	int N = 11;
	Vector f = Vector(N);
	Vector x = Vector(N);
	Vector s = Vector(N);
	double h = x.linspace(0,1);
	for (int i = 0; i < N; i++) {
		f(i) = srcFunction(x(i)) * h * h;
		s(i) = solution(x(i));
	}
	f(0) = 0;
	f(N-1) = 0;
	Vector y = tridiagSolve(-1,2,-1,f);
	// Write to file
	std::ofstream myfile;
	myfile.open("test.txt");
	for (int i = 0; i < N; i++) {
		myfile << i*h << " ";
		myfile << s(i) << " ";
		myfile << y(i) << " ";
		myfile << endl;
	}
	myfile.close();
	return 0;
}

