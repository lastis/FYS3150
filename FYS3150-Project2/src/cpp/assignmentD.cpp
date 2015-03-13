#include <iostream>
#include "CPhys.h"
#include <vector>
#include <fstream>

using namespace std;
using namespace CPhys;

int main()
{
	char    fileName[] = "WaveFunc1.txt";
	int	N	 = 10;
	Matrix 	A 	 = Matrix(N,N);
	Matrix 	R	 = Matrix(N,N);
	Vector 	rho 	 = Vector(N);
	Vector 	tmp	 = Vector(N);
	double 	rho0	 = 0;
	double 	rhoN  	 = 5;
	double	omega	 = 0.01;
	double  h 	 = rho.linspace(rho0,rhoN,N);
	double 	d	 =  2/(h*h);
	double 	e	 = -1/(h*h);

	A.eye();
	// Make the diagonal elements and insert them
	tmp      = e;
	tmp(0)   = 0;
	tmp(N-1) = 0;
	A.diag(tmp, 1);
	A.diag(tmp,-1);
	tmp(0)   = 1;
	tmp(N-1) = 1;
	for (int i = 1; i < N; i++) {
		tmp(i) = d + omega*omega*rho(i)*rho(i) + 1.0/rho(i);
	}
	A.diag(tmp,0);
	// Solve the problem
	EigVal::jacobiMethod(A, R, N);
	// Sort the eigenvectors and values
	tmp = A.diag();
	MatOp::sortCol(R, tmp);
	// Put the sorted eigenvectors back into A
	A.diag(tmp);
	// Get the eigenvector we want
	Vector eig1 = R.getCol(1);
	VecOp::normalize(eig1,h);
	// Print
	cout << "Rho   :\t" << rhoN << endl;
	cout << "N     :\t" << N << endl;
	cout << "Omega :\t" << omega << endl;
	// Write to file
	std::ofstream myfile;
	myfile.open(fileName);
	for (int i = 0; i < N; i++) {
		myfile << rho(i) << " ";
		myfile << eig1(i) << " ";
		myfile << endl;
	}
	myfile.close();
	return 0;
}

