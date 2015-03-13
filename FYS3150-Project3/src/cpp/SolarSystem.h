#include <iostream>
#include <cmath>
#include "Planet.h"
#include "CPhys.h"

class SolarSystem : public ODE_Interface{
public:
	 SolarSystem     (int steps, int capacity, int stepLength);
	~SolarSystem     ();
	void 	addPlanet(Planet& planet);
	void	simulate ();
	void	stabilize();
	Matrix	getArray();
	int	getPlanetCount();

private:
	void	appendPosVel(Planet& planet);
	void	derivative  (double*      u, double t, double* v);

	// Variables
	static const int	mDim = 2;
	double   G;
	Planet** mppPlanets;
	double*  mpMass;
	int	 mPlanetCnt;
	int	 mCapacity;
	int 	 mStepLength;
	// This variable contains the planet positions and velocities 
	Matrix	 mU;
	double** mppU; // "mpp" > my pointer pointer

	void inline force(   int& planet1,    int& planet2, double* u,
			  double&      Fx, double&      Fy){
		double pos[2];
		double r2Inv;
		// Calculate the indices, there are 4 coordinates
		// x, y, vx, vz
		int i1 = planet1*2*mDim;
		int i2 = planet2*2*mDim;
		// Get the relative coordinates
		pos[0] = u[i2+0] - u[i1+0];
		pos[1] = u[i2+1] - u[i1+1];
		double x = pos[0];
		double y = pos[1];
		// Calculate x^2, y^2, 1/r^2,
		pos[0] = pos[0]*pos[0];
		pos[1] = pos[1]*pos[1];
		r2Inv = 1/(pos[0]+pos[1]);
		//std::cout << "radius: " << sqrt(pos[0] + pos[1]) << std::endl;
		//std::cout << "radius2inv" << r2Inv << std::endl;
		// Calculate the total force
		double F = G*mpMass[planet1]*mpMass[planet2]*r2Inv;
		//std::cout << "F : " << F << std::endl;
		// Find the directions
		Fx = F*x*sqrt(r2Inv);
		Fy = F*y*sqrt(r2Inv);
		//std::cout << "F again : " << sqrt(Fx*Fx+Fy*Fy) << std::endl;
	}
	
	// First planet == 0
	void inline totForce(int& planet, double* u, double& Fx, double& Fy){
		// Add forces from all the planets
		double tmpFx, tmpFy;
		Fx = 0;
		Fy = 0;
		// Itterate over mPlanetCnt and not mCapacity
		// because ppU might not have been fully filled with planets
		// Calculate the force between planet p and planet i
		for (int i = 0; i < planet; i++) {
			// i cannot become p
			force(planet,i,u,tmpFx,tmpFy);
			Fx += tmpFx;
			Fy += tmpFy;
		}
		for (int i = planet+1; i < mPlanetCnt; i++) {
			// i cannot become p
			force(planet,i,u,tmpFx,tmpFy);
			Fx += tmpFx;
			Fy += tmpFy;
		}
	}
};
