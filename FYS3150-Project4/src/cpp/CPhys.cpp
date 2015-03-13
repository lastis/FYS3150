#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "CPhys.h"


using namespace std;
using namespace CPhys;

Vector	LinAlg::tridiagSolve(double a, double b, double c, Vector y){
	int N = y.getLength();
	Vector x = Vector(N);
	double* xVec  = x.getArrayPointer();
	double* yVec  = y.getArrayPointer();
	pLinAlg::tridiagSolve(a,b,c,xVec,yVec,N);
	// return the vector
	return x;
}

Vector	Diffusion::solveEuler(Vector init, double time, double dt, double dx) {
	int N = init.getLength();
	Vector u1 = Vector(init);
	Vector u2 = Vector(N);
	double 	alpha = dt/(dx*dx);
	double* pU1 = u1.getArrayPointer();
	double* pU2 = u2.getArrayPointer();
	double  a = alpha;
	double  b = 1-2*alpha;
	for (int j = 1; j*dt < time; j++) {
		pDiffusion::stepEuler(pU1,pU2,a,b,N);
		// Copy values to pU1
		for (int i = 0; i < N; i++) {
			pU1[i] = pU2[i];
		}
	}
	return u1;
}

Vector	Diffusion::solveBackwardEuler(Vector init, double time, 
				      double   dt, double dx) {
	int N = init.getLength();
	Vector u1 = Vector(init);
	Vector u2 = Vector(N);
	double 	alpha = dt/(dx*dx);
	double* pU1 = u1.getArrayPointer();
	double* pU2 = u2.getArrayPointer();
	double  a = (-alpha);
	double  b = (1+2*alpha);
	for (int j = 1; j*dt < time; j++) {
		pDiffusion::stepBackwardEuler(pU1,pU2,a,b,N);
		// Copy values to pU1
		for (int i = 0; i < N; i++) {
			pU1[i] = pU2[i];
		}
	}
	return u1;
}

Vector	Diffusion::solveCrankNicolson(Vector init, double time, 
				      double   dt, double dx) {
	int N = init.getLength();
	Vector u1 = Vector(init);
	Vector u2 = Vector(N);
	double 	alpha = dt/(dx*dx);
	double* pU1 = u1.getArrayPointer();
	double* pU2 = u2.getArrayPointer();
	double  a1 = alpha;
	double  a2 = -alpha;
	double  b1 = 2-2*alpha;
	double  b2 = 2+2*alpha;
	for (int j = 1; j*dt < time; j++) {
		pDiffusion::stepEuler(pU1,pU2,a1,b1,N);
		// Copy values to pU1
		for (int i = 0; i < N; i++) {
			pU1[i] = pU2[i];
		}
		pDiffusion::stepBackwardEuler(pU1,pU2,a2,b2,N);
		// Copy values to pU1
		for (int i = 0; i < N; i++) {
			pU1[i] = pU2[i];
		}
	}
	return u1;
}

void	ODE::solveEuler(Matrix&  u, Vector& 	   uInit, 
		        double  dt, ODE_Interface*   ode){
	// Call a more general solve method
	using namespace pODE;
	solve(u,uInit,dt,ode,&stepEuler);
}

void	ODE::solveRK4(Matrix&  u, Vector& 	 uInit, 
		      double  dt, ODE_Interface*   ode){
	// Call more general solve method
	using namespace pODE;
	solve(u,uInit,dt,ode,&stepRK4);
}

void EigVal::jacobiMethod(Matrix& A, Matrix& R, int N){
	using namespace pEigVal;
	double** ppR = R.getArrayPointer();
	double** ppA = A.getArrayPointer();
	int 	 k, l;
	int 	 iterations  	= 0;
	double	 eps 		= 10e-8;
	double 	 maxIterations 	= (double) N * N * N;
	// Prepare eigenvector matrix
	R.eye();
	// Decide which elements to rotate
	double 	maxoff		= maxoffdiag(ppA, &k, &l, N);
	// Do the rotations until we converge at a solution
	while (		maxoff > eps 	
		 && iterations < maxIterations ) {
		
		int num = iterations;
		if(num%20000 == 0){
			cout << "Computing - "<< iterations << " iterations ";
			cout << " max off diag = " << maxoff << endl;
		}
		rotate(ppA, ppR, k, l, N );
		maxoff = maxoffdiag(ppA, &k, &l, N );
		iterations++;
	}
	cout << "Number of iterations: " << iterations << "\n";
}



double 	VecOp::normalize(Vector& v, double dx){
	int     N   = v.getLength();
	double* pV  = v.getArrayPointer();
	double  sum = 0;
	double  Z   = 0;
	for (int i = 0; i < N; i++) {
		sum += pV[i]*pV[i]*dx;
	}
	Z = 1/sqrt(sum);
	for (int i = 0; i < N; i++) {
		pV[i] *= Z;
	}
	return Z;

}

void	MatOp::sortCol(Matrix& A, Vector& v){
	// We want to sort the colums in A by the elements in v.
	int rows = A.getN();
	int cols = A.getM();
	// Make a new matrix wit dim N = A.N+1 and M = A.M
	Matrix tmp = Matrix(rows+1,cols);
	// Insert v in the first row
	tmp.setRow(0,v);
	// Copy A to tmp after the first row
	double** ppTmp = tmp.getArrayPointer();
	double** ppA   =    A.getArrayPointer();
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ppTmp[i+1][j] = ppA[i][j];
		}
	}
	// Our matrix is row major so we need to transpose
	// to be able to sort by columns
	// Now **arr[0,1,2, ... , i] points to whole arrays 
	// of the columns
	tmp.t();
	// Sort by first elements of A (sorting by vector v)
	ppTmp = tmp.getArrayPointer();
	qsort((void*)ppTmp,cols,sizeof(double),&pMatOp::compareTwoRows);
	// Transpose back
	tmp.t();
	ppTmp = tmp.getArrayPointer();
	// Now copy the sorted elements to A
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			 ppA[i][j] = ppTmp[i+1][j];
		}
	}
	// Then sort v
	v.sort();
}

void	pODE::solve(Matrix&  u, Vector&        uInit,  
		    double& dt, ODE_Interface*   ode,
	      void (*stepFunc) (double*,double,double*       ,
				double&,int&  ,ODE_Interface*)){
	// Number of equations
	int 	 eq    = u.getM();
	int	 N     = u.getN();
	double** ppU   = u.getArrayPointer();
	double*  pU    = new double[eq];
	double*  pUOut = new double[eq];

	// Insert the inital values
	for (int i = 0; i < eq; i++) {
		pU[i]   = uInit(i);
		ppU[0][i] = pU[i];
	}
	// Start calculating with the given step function
	for (int i = 0; i < N-1; i++) {
		// Paramaters are;
		// uIn, current time, uOut, dt, the number of equations and 
		// the object with the derivative function. 
		stepFunc(pU,i*dt,pUOut,dt,eq,ode);
		// We have now calculated step i+1
		for (int j = 0; j < eq; j++) {
			ppU[i+1][j] = pUOut[j];
			pU[j]       = pUOut[j];
		}
	}
}


void 	pODE::stepEuler(double* pU, double   t, double*        pUOut, 
			double& dt, int&    eq, ODE_Interface*   ode){
	double* pV  = new double[eq];
	// Calculate pV
	ode->derivative(pU,t,pV);
	for (int i = 0; i < eq; i++) {
		pUOut[i] = pU[i] + pV[i]*dt;
	}
}

void 	pODE::stepRK4(double* pU, double   t, double* 	     pUOut, 
		      double& dt, int&    eq, ODE_Interface*   ode){
	double* pV  = new double[eq];
	double* k1  = new double[eq];
	double* k2  = new double[eq];
	double* k3  = new double[eq];
	double* k4  = new double[eq];
	double* tmp = new double[eq];
	double dt2 = dt/2;
	// The derivative functions calculates pV
	// Calculate k2
	ode->derivative(pU,t,pV);
	for (int i = 0; i < eq; i++) {
		k1[i] = pV[i]*dt;
	}
	// Calculate k2
	for (int i = 0; i < eq; i++) {
		tmp[i] = pU[i] + 0.5*k1[i];
	}
	ode->derivative(tmp,t+dt2,pV);
	for (int i = 0; i < eq; i++) {
		k2[i] = pV[i]*dt;
	}
	// Calculate k3
	for (int i = 0; i < eq; i++) {
		tmp[i] = pU[i] + 0.5*k2[i];
	}
	ode->derivative(tmp,t+dt2,pV);
	for (int i = 0; i < eq; i++) {
		k3[i] = pV[i]*dt;
	}
	// Calculate k4
	for (int i = 0; i < eq; i++) {
		tmp[i] = pU[i] + k3[i];
	}
	ode->derivative(tmp,t+dt,pV);
	for (int i = 0; i < eq; i++) {
		k4[i] = pV[i]*dt;
	}
	// Calculate new value
	for (int i = 0; i < eq; i++) {
		pUOut[i] = pU[i] + (1/6.0)*
			(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	}
	
}

void	pLinAlg::tridiagSolve(double  a, double  b, double c, 
			      double* x, double* y, int    N){
	double* temp  = new double[N];
	double* aVec  = new double[N];
	double* bVec  = new double[N];
	double* cVec  = new double[N];
	double  bTemp = 0;
	// Init arrays
	for (int i = 0; i < N; i++) {
		temp[i] = 0;
		aVec[i] = a;
		bVec[i] = b;
		cVec[i] = c;
	}
	// Forward substitution
	bTemp = bVec[1];
	x[1]  = y[1]/bTemp;
	for (int i = 2; i < N; i++) {
		temp[i] =  cVec[i-1]/bTemp;
		bTemp   =  bVec[i] - aVec[i]*temp[i];
		x[i] = (y[i] - aVec[i]*x[i-1])/bTemp;
	}
	// Backward substitution
	for (int i = N-1; i >= 1; i--) {
		x[i] -= temp[i+1]*x[i+1];
	}
}

void	pDiffusion::stepEuler(double* u, double* un, 
			      double  a, double   b, int N){
	// Forward Euler scheme for the diffusion equation
	for (int i = 1; i < N-1; i++) {
		un[i] = a*u[i-1] + b*u[i] + a*u[i+1];
	}
	un[0] = u[0];
}

void	pDiffusion::stepBackwardEuler(double* u, double* un, 
				      double  a, double   b, int N){
	// Backward Euler scheme for the diffusion equation
	pLinAlg::tridiagSolve(a,b,a,un,u,N);
	un[N-1] = u[N-1];
}

int 	pMatOp::compareTwoRows(const void* rowA, const void* rowB){
	// Compare the first element in row A (this is what ** does) and 
	// compare it to the first element in row B
	return (**(double**) rowA - **(double**) rowB);
}

double 	pEigVal::maxoffdiag(double** A, int* k, int* l, int N){
	double max  = 0.0;
	double a_ij = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			a_ij = fabs(A[i][j]);
			if (a_ij > max){
				max	= a_ij;
				*l 	= i;
				*k 	= j;
			}
		}
	}
	return max;
}

void 	pEigVal::rotate(double** A, double** R,int k, int l, int N){
	double s, c;
	double tau, t;
	tau = (A[l][l] -A[k][k])/(2*A[k][l]);
	if(tau > 0) t = -tau + sqrt(1+tau*tau);
	else 	    t = -tau - sqrt(1+tau*tau);
	c = 1 / sqrt(1+t*t);
	s = t * c;
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A[k][k];
	a_ll = A[l][l];
	// changing the matrix elements with indices k and l
	A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
	A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
	A[k][l] = 0.0; // hard-coding of the zeros
	A[l][k] = 0.0;
	// and then we change the remaining elements
	for ( int i = 0; i < N; i++ ) {
		if ( i != k && i != l ) {
			a_ik 	= A[i][k];
			a_il 	= A[i][l];
			A[i][k]	= c*a_ik - s*a_il;
			A[k][i]	= A[i][k];
			A[i][l]	= c*a_il + s*a_ik;
			A[l][i]	= A[i][l];
		}
		// Finally, we compute the new eigenvectors
		r_ik 	= R[i][k];
		r_il 	= R[i][l];
		R[i][k]	= c*r_ik - s*r_il;
		R[i][l]	= c*r_il + s*r_ik;
	}
}
