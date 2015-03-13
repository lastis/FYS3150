#include <iostream>
#include "Planet.h"
using namespace std;

Planet::Planet(const char* name, double aphelion, double velMin,
	        double aphAngle, double mass		     ){

	mName 		= name;
	mAphelion 	= aphelion;
	mVelMin 	= velMin;
	mAphAngle 	= aphAngle;
	mMass 		= mass;

}

const char* Planet::getName(){
	return mName;
}

double	Planet::getMass(){
	return mMass;
}

double 	Planet::getAphelion(){
	return mAphelion;
}

double	Planet::getVelMin(){
	return mVelMin;
}

double	Planet::getAphAngle(){
	return mAphAngle;
}

