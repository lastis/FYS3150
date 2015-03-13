#include "SolarSystem.h"
#include "constants.h"
#include <cmath>
using namespace std;

SolarSystem::SolarSystem(int points, int capacity, int stepLength){
	mStepLength = stepLength;
	G = 0.0002959;
	mPlanetCnt  = 0;
	mCapacity = capacity;
	// This needs the vector class
	mppPlanets  = new Planet*[capacity];
	mpMass	    = new  double[capacity];
	// u contains pos x,y and vel x,y for all the planets
	mU   = Matrix(points,2*mDim*mCapacity);
	mppU = mU.getArrayPointer();
}

SolarSystem::~SolarSystem(){
	delete mpMass;
}

Matrix SolarSystem::getArray(){
	return mU;
}

void	SolarSystem::simulate(){
	Vector uInit = mU.getRow(0);
	using namespace CPhys::ODE;
	solveRK4(mU,uInit,mStepLength,this);
}

void	SolarSystem::stabilize(){
	// We assume the sun is the first planet that is added
	double vx, vy;
	double energyX = 0;
	double energyY = 0;
	for (int i = 1; i < mPlanetCnt; i++) {
		vx = mppU[0][i*4+2];
		vy = mppU[0][i*4+3];
		energyX += mpMass[i]*vx;
		energyY += mpMass[i]*vy;
	}
	vx = -energyX/mpMass[0];
	vy = -energyY/mpMass[0];
	
	mppU[0][2] = vx;
	mppU[0][3] = vy;
}

void	SolarSystem::addPlanet(Planet& planet){
	// Pointer pointer array to the planets
	mppPlanets[mPlanetCnt] = &planet;
	mpMass[mPlanetCnt] = planet.getMass()/cMSun;
	mPlanetCnt++;
	appendPosVel(planet);
}

void	SolarSystem::appendPosVel(Planet& planet){
	// This method is only called by the add planet method
	int    index;
	double x,y;
	double vel = planet.getVelMin();
	double aph = planet.getAphelion();
	double ang = planet.getAphAngle();
	vel = vel * (60*60*24)/ cAU;
	aph = aph / cAU;
	ang = ang*2*cPI/360.0;

	x   = aph*cos(ang);
	y   = aph*sin(ang);
	// Figure out where we should append the values
	index = 2*mDim*(mPlanetCnt-1);
	mppU[0][index+0] = x;
	mppU[0][index+1] = y;
	// The angle to the velocity is 90 degrees of prior angle
	x = -vel*sin(ang);
	y =  vel*cos(ang);
	mppU[0][index+2] = x;
	mppU[0][index+3] = y;
}

int	SolarSystem::getPlanetCount(){
	return mPlanetCnt;
}

void	SolarSystem::derivative(double* u, double t, double* v){
	double Fx,Fy;
	for (int i = 0; i < mPlanetCnt; i++) {
		//cout << "Mass : " << mpMass[i] << endl;
		// This calculates Fx and Fy from the positions
		totForce(i,u,Fx,Fy);
		// x derivative is velocity x
		v[4*i+0] = u[4*i+2];
		// y derivative is velocity y
		v[4*i+1] = u[4*i+3];
		// vx derivative is Fx/m
		v[4*i+2] = Fx/mpMass[i];
		// vy derivative is Fy/m
		v[4*i+3] = Fy/mpMass[i];
		//cout << "Fx : " << Fx << endl;
		//cout << "Fy : " << Fy << endl;
		//cout << "Total force: " << sqrt(Fx*Fx+Fy*Fy) << endl;
	}
}
