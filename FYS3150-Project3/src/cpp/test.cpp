#include "CPhys.h"
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
using namespace CPhys::ODE;

class Deriv : public ODE_Interface{
public:
	//deriv(){}
	void	derivative(double* u, double t, double* v){
		v[0] =  u[1];
		v[1] = -u[0];
	}
};

int main(int argc, const char *argv[])
{
	int N = 100;
	int diff = 2;
	double h = 4.0*acos(-1.)/N;
	Matrix u = Matrix(N,diff);
	Vector uInit = Vector(diff);
	
	uInit(0) = 0.5;
	uInit(1) = 0;

	Deriv deriv = Deriv();
	solveRK4(u,uInit,h,&deriv);

	// Write to file
	std::ofstream myfile;
	myfile.open("test.txt");
	for (int i = 0; i < N; i++) {
		myfile << i*h << " ";
		myfile << u(i,0) << " ";
		myfile << u(i,1) << " ";
		myfile << endl;
	}
	myfile.close();
	return 0;
}

