#include "CPhys.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace CPhys;

main() {
	char	fileName[] = "2DT0_5Imp.dat";
	double 	alpha 	= 0.25;
	double 	time 	= 0.5;
	double 	h 	= 0.01;
	double 	dt 	= alpha*h*h;
	int 	N 	= int(1/h) + 1;

	Matrix init = Matrix(N,N);
	// Initialize the inital state
	for(int j = 0; j < N; j++){
		for(int i = 0; i < N; i++){
			init(i,j) = (1-j*h)*exp(i*h);
		}
	}
	// Find the solution
	Matrix solution = Diffusion::solveJacobi2D(init, time, dt, h);
	// Print to file
	ofstream outFile(fileName);
	for (int i=0; i < N; i++){
		for (int j=0; j < N; j++){
			outFile << setw(15) << setprecision(8) << solution(i,j);
		}
		outFile << endl;
	}
	outFile.close();
}	
