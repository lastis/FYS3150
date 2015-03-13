#include <fstream>
#include "SolarSystem.h"

using namespace std;
using namespace CPhys::ODE;


int main(int argc, const char *argv[])
{

	// Setting up the system
	// 	// Data: http://nssdc.gsfc.nasa.gov/planetary/factsheet
	//                             Label     Aphelion    
	//                             VMin      Aph Angle  Mass 
	Planet sun     = Planet(    "Sun",          0,          
				        0,          0,  1.9891e30);
	Planet mercury = Planet("Mercury",    69.82e9,
			          38.86e3,        350,  3.3022e23);
	Planet venus   = Planet(  "Venus",   108.94e9,   
			          34.79e3,         45,  4.8685e24);
	Planet earth   = Planet(  "Earth",   152.10e9,
			          29.29e3,         20,  5.9736e24);
	Planet moon    = Planet(   "Moon", 152.5055e9, 
				 30.254e3,         20,  7.349e22 );
	Planet mars    = Planet(   "Mars",   249.23e9, 
				  21.97e3,        250,  6.4185e23);
	Planet jupiter = Planet("Jupiter",   816.62e9,   
				  12.44e3,        285, 1.8986e27);
	Planet saturn  = Planet( "Saturn",  1514.50e9,  
				   9.09e3,          0,  5.6846e26);
	Planet uranus  = Planet( "Uranus",  3003.62e9,  
				   6.49e3,         75,  8.6810e25);
	Planet neptune = Planet("Neptune",  4545.67e9,  
				   5.37e3,        295,  1.0243e26);
	Planet pluto   = Planet(  "Pluto",  7375.93e9, 
				   3.71e3,        295,  1.305e22 );
//	Planet rocket  = new Planet( "Rocket",   152.10e9,   
//				      29.29e3,         20,  1e5);
//	Planet custom1 = new Planet("Custom1",   816.62e9,   
//				      12.44e3,        285,  1.8986e30);

	// Calculate the velocity of the sun to get the total momentum
	// of the system = 0
	
	int planets = 11;
	double year = 10;
	int stepLength = 5;
	int N = year*365/stepLength;

	// The paramters are points, capacity and steplength
	SolarSystem system = SolarSystem(N,planets,stepLength);
	system.addPlanet(sun);
	system.addPlanet(mercury);
	system.addPlanet(venus);
	system.addPlanet(earth);
	system.addPlanet(moon);
	system.addPlanet(mars);
	system.addPlanet(jupiter);
	system.addPlanet(saturn);
	system.addPlanet(uranus);
	system.addPlanet(neptune);
	system.addPlanet(pluto);
	
	system.stabilize();
	system.simulate();
	Matrix u = system.getArray();

	// Write to file
	std::ofstream myfile;
	myfile.open("SolarSystem.txt");
	N = u.getN();
	for (int i = 0; i < N; i++) {
		myfile << i << " ";
		for (int p = 0; p < planets; p++) {
			myfile << u(i,p*4+0) << " ";
			myfile << u(i,p*4+1) << " ";
		}
		myfile << endl;
	}
	myfile.close();
	return 0;
}

