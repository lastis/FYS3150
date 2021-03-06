#include "Matrix.h"
#include "ODE_Interface.cpp"
namespace CPhys{
	namespace LinAlg{
		Vector	tridiagSolve(double a, double b, double c, Vector y);
	}
	namespace Diffusion{
		Vector 	solveEuler(Vector init, double time, 
				   double   dt, double   dx);
		Vector 	solveBackwardEuler(Vector init, double time, 
				           double   dt, double   dx);
		Vector 	solveCrankNicolson(Vector init, double time, 
				           double   dt, double   dx);
	}
	namespace EigVal{
		void 	jacobiMethod(Matrix& A, Matrix& R, int N);
	}
	namespace MatOp{
		void	sortCol(Matrix& A, Vector& v);
	}
	namespace VecOp{
		double	normalize(Vector& v, double dx);
	}
	namespace ODE{
		void	solveEuler(Matrix&  u, Vector&        uInit,   
				   double  dt, ODE_Interface*   ode);

		void	solveRK4  (Matrix&  u, Vector& 	    uInit, 
				   double  dt, ODE_Interface*   ode);
	}


	// prefix p (abbr. for 'pointer') is used to avoid name collision
	namespace pMatOp{
		int 	compareTwoRows(const void* rowA, const void* rowB);
	}
	namespace pLinAlg{
		void	tridiagSolve(double  a, double  b, double c, 
				     double* x, double* y, int    N);
	}
	namespace pEigVal{
		double 	maxoffdiag        (double** A, int* k, int* l, int N);
		void	rotate(double** A, double** R, int  k, int  l, int N);
	}
	namespace pDiffusion{
		void	stepEuler(double* u, double* un, 
				  double  a, double   b, int N);
		void	stepBackwardEuler(double* u, double* un,
				          double  a, double   b, int N);
	}
	namespace pODE{
	      	// These function includes function pointer parameters
		void	solve(Matrix&  u, Vector&        uInit,  
			      double& dt, ODE_Interface*   ode,
		     	void (*stepFunc) (double*,double,double*       ,
			               	  double&,int&  ,ODE_Interface*));

		void 	stepEuler(double* pU, double   t, double*        pUOut, 
			          double& dt, int&    eq, ODE_Interface*   ode);

		void 	stepRK4(  double* pU, double  t, double* pUOut, 
			          double& dt, int&    eq, ODE_Interface*   ode);
	}
}

