#include <iostream>
#include <time.h>
#include <cmath>
#include "lib.h"

using namespace std;


int main()
{
	long seed  = 432423;
	int const N = 10000;
	double **A = (double **)matrix(N,N,sizeof(double));
	double **B = (double **)matrix(N,N,sizeof(double));
	double **C = (double **)matrix(N,N,sizeof(double));
	// Reset A
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 0;
		}
	}
	// Init matrices
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			B[i][j] = ran0(&seed);
			C[i][j] = ran0(&seed);
		}
	}
	// Init time
	clock_t start, finish = 0;;
	start = clock();
	// Multiplay with memory stride. 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				A[i][j] += B[i][k] * C[k][j];
			}
		}
	}
	// Stop time
	finish = clock();
	double time1 = (double)(finish - start)/CLOCKS_PER_SEC;
	start  = 0;
	finish = 0;
	// Reset A
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 0;
		}
	}
	// Start time
	start = clock();
	// Multiplay agains memory stride. 
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			for (int k = 0; k < N; k++) {
				A[i][j] += B[i][k] * C[k][j];
			}
		}
	}
	// Stop time
	finish = clock();
	double time2 = (double)(finish - start)/CLOCKS_PER_SEC;
	start  = 0;
	finish = 0;
	start = clock();
	// Print times
	cout << "Time 1 (with memory stride)     : " << time1 << endl;
	cout << "Time 2 (with out memory stride) : " << time2 << endl;
	cout << "Ratio Time 1 / Time 2           : " << time1/time2 << endl;
	return 0;
}

// Execution times:
// N =  500	Time1:	1.27	Time2:	0.86	Ratio:	1.47674
// N = 1000	Time1:	12.48	Time2:	9.18	Ratio:	1.35948
