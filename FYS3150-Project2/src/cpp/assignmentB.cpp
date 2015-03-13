#include <cmath>
#include <iostream>
#include "CPhys.h"

using namespace std;
using namespace CPhys;

int main()
{
	int	N	= 200;
	Matrix 	A 	= Matrix(N,N);
	Matrix 	R	= Matrix(N,N);
	Vector 	rho 	= Vector(N);
	Vector 	tmp	= Vector(N);
	double 	rho0	= 0;
	double 	rhoN 	= 5;
	double  h 	= rho.linspace(rho0,rhoN,N);
	double 	d	=  2/(h*h);
	double 	e	= -1/(h*h);

	// Make the diagonal elements and insert them
	tmp = e;
	// Set the boundry conditions to 1 on either side
	// to get more accuray results
	tmp(0)   = 0;
	tmp(N-1) = 0;
	A.diag(tmp, 1);
	A.diag(tmp,-1);
	tmp(0)   = 1;
	tmp(N-1) = 1;
	for (int i = 1; i < N-1; i++) {
		tmp(i) = d + rho(i)*rho(i);
	}
	A.diag(tmp,0);
	// Solve the problem
	EigVal::jacobiMethod(A, R, N);
	// Print eigenvalues
	cout << "Rho:\t" << rhoN << endl;
	cout << "N  :\t" << N << endl;
	cout << "Eigenvalues :" << endl;
	Vector eigVal = A.diag();
	eigVal.sort();
	for (int i = 1; i < 6; i++) {
		cout << "R" << i << ":\t" << eigVal(i) << endl;
	}
	return 0;
}

