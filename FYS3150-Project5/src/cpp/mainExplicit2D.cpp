#include <fstream>
#include <iostream>
#include "CPhys.h"
#include <iomanip>
#include <cmath>

using namespace std;
using namespace CPhys;

main() {
	char	fileName[] = "2DT0_01Exp.dat";
	double 	alpha 	= 0.25;
	double	time	= 0.01;
	double 	h 	= 0.01;
	double 	dt 	= alpha*h*h;
	int 	N 	= int(1/h)+1;

	Matrix init = Matrix(N,N);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++){
			init(i,j) = (1-j*h)*exp(i*h);
		}
	}
	Matrix solution = Diffusion::solveEuler2D(init, time, dt, h);
	ofstream outFile(fileName);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			outFile << setw(15) << setprecision(8) << solution(i,j);
		}
		outFile << endl;
	}
	outFile.close();			
}
