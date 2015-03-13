#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

inline double srcFunction(double x){
	return 100*exp(-10*x);
}

inline double solution(double x){
	return 1.-(1.-exp(-10))*x-exp(-10*x);
}

int main(){
	const int N = 100000;
	char fileName[] = "N100000.txt"; 
	char fileNameError[] = "errorN100000.txt";

	double boundryStart = 0;
	double boundryEnd   = 0;
	double xStart = 0;
	double xEnd   = 1;
	double f[N];
	double x[N];
	double error[N];

	double a[N];
	double b[N];
	double c[N];
	double u[N];
	double v[N];
	double temp[N];
	double btemp;
	double diagA = -1;
	double diagB =  2;
	double diagC = -1;

	// Manual linspace
	double delta = xEnd - xStart;
	double step  = delta/(N-1);       // step defined here
	for (int i = 0; i < N; i++) {
		x[i] = xStart + i*step;
	}
	// Init arrays
	for (int i = 0; i < N; i++) {
		f[i] = srcFunction(x[i]) * step  * step;
		u[i] = solution(x[i]);
		a[i] = diagA;
		b[i] = diagB;
		c[i] = diagC;
	}
	// Init time
	clock_t start, finish;
	start = clock();
	// Forward substitution
	btemp = b[0];
	v[0]  = f[0]/btemp;
	for (int i = 1; i < N; i++) {
		temp[i] =  c[i-1]/btemp;
		btemp   =  b[i] - a[i]*temp[i];
		v[i]    = (f[i] - a[i]*v[i-1])/btemp;
	}
	// Backward substitution
	for (int i = N-2; i >= 0; i--) {
		v[i] -= temp[i+1]*v[i+1];
	}
	// Calculate time
	finish = clock();
	double time = (double)(finish - start)/CLOCKS_PER_SEC;
	cout << "Time : " << time << endl;
	// Compute error
	for (int i = 0; i < N; i++) { 
		error[i] = abs((v[i] - u[i])/u[i]);
		error[i] = log10(error[i]);
	}
	// Write solution to file
	ofstream myfile;
	myfile.open(fileName);
	for (int i = 0; i < N; i++) {
		myfile << x[i] << " " << v[i] << endl;
	}
	myfile.close();
	// Write error to file
	myfile.open(fileNameError);
	for (int i = 0; i < N; i++) {
		myfile << x[i] << " " << error[i] << endl;
	}
	myfile.close();
	return 0;

	// Execution times:
	// n = 10	: 0	
	// n = 100 	: 0
	// n = 1000 	: 0
	// n = 10000	: 0
	// n = 1e5	: 0
	// n = 1e6	: Segmentation fault
}
