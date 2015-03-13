#include <iostream>
#include <armadillo>
#include <time.h>
#include <fstream>
using namespace std;
using namespace arma;

inline vec srcFunction(vec x){
	return 100*exp(-10*x);
}

int main(int argc, const char *argv[])
{
	const int N = 1000;

	char 	fileName[] 	= "N1000LU.txt"; 
	double 	xStart		= 0;
	double 	xEnd  		= 1;
	double 	h		= (xEnd - xStart)/(N+1);
	mat A(N,N);
	mat L(N,N);
	mat U(N,N);
	// Init arrays
	A.diag( 0) += 2;
	A.diag( 1) +=-1;
	A.diag(-1) +=-1;
	vec x = linspace<vec>(xStart, xEnd, N);
	vec w = srcFunction(x)*h*h;
	// Init time
	clock_t start, finish;
	start = clock();
	// Solve with LU decomp
	lu(L,U,A);
	vec y = solve(L,w);
	vec v = solve(U,y);
	// Calculate time
	finish = clock();
	double time = (double)(finish - start)/CLOCKS_PER_SEC;
	cout << "Time : " << time << endl;
	// Write to file
	ofstream myfile;
	myfile.open(fileName);
	for (int i = 0; i < N; i++) {
		myfile << x[i] << " " << v[i] << endl;
	}
	myfile.close();
	return 0;

	// Execution times:
	// n = 10	: 0	
	// n = 100 	: 0
	// n = 1000 	: 0.04
	// n = 10000	: 4.31
	// n = 100000	: Too large
}


