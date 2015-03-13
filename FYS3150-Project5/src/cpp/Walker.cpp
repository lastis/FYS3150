#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CPhys.h"
#include "Walker.h"

using namespace std;
using namespace CPhys;

Walker::Walker(int walkers, double walkProb, double time, double dt, 
	       double minDist, double maxDist){
	mMinDist = minDist;
	mMaxDist = maxDist;
	mUseGauss = false;
	// This function initializes the vectors we need
	setTime(time, dt);

	mWalkProb	= walkProb;
	mWalkers 	= walkers;
}

void	Walker::useGauss(bool use){
	mUseGauss = use;
}

void	Walker::setTime(double time, double dt){
	mMaxSteps = int (time/dt+0.5);
	mdt = dt;
	ml = sqrt(2*mdt);
	mVecSize = 1/ml+1;

	mxPosComulative 	= Vector(mMaxSteps);
	mxPosComulative2	= Vector(mMaxSteps);
	mxProbabilities  	= Vector(mVecSize);
	mpxPosComulative  	= mxPosComulative.getArrayPointer();
	mpxPosComulative2 	= mxPosComulative2.getArrayPointer();
	mpxProbabilities  	= mxProbabilities.getArrayPointer();
	myPosComulative 	= Vector(mMaxSteps);
	myPosComulative2	= Vector(mMaxSteps);
	myProbabilities  	= Vector(mVecSize);
	mpyPosComulative  	= myPosComulative.getArrayPointer();
	mpyPosComulative2 	= myPosComulative2.getArrayPointer();
	mpyProbabilities  	= myProbabilities.getArrayPointer();
}

void	Walker::start1D(){
	mcSampling();
}

void	Walker::start2D(){
	mcSampling2();
}

void 	Walker::mcSampling(){
	long seed = -1;
	Vector posVec = Vector(mWalkers);
	double* pos = posVec.getArrayPointer();
	// Itterate over all the walkers
	for (int step = 0; step < mMaxSteps; step++) {
		// Itterate the steps
		for (int walker = 0; walker < mWalkers; walker++) {
			double 	ran 	= 0;
			double	gauss 	= 1;
			if (pos[walker] < mMinDist || pos[walker] > mMaxDist ) continue;
			if (mUseGauss) {
				gauss = Random::gauss(seed);
			}
			ran = Random::ran0(seed);
			if (ran < mWalkProb) {
				if (mUseGauss){
					pos[walker] += ml*gauss;
				}
				else {
					pos[walker] += ml;
				}
			}
			else {
				if (mUseGauss){
					pos[walker] -= ml*gauss;
				}
				else {
					pos[walker] -= ml;
				}
			}
			int index = (pos[walker] - mMinDist)/ml + 0.5;
			if (index < 0 || index >= mVecSize) continue;
			mpxPosComulative [step] += pos[walker];
			mpxPosComulative2[step] += pos[walker]*pos[walker];
			mpxProbabilities[index] += 1;
		}
	}
}

void 	Walker::mcSampling2(){
	long 	seed = -1;  // initialise random number generator
	Vector 	posxVec = Vector(mWalkers);
    	Vector 	posyVec = Vector(mWalkers);
    	double* posx = posxVec.getArrayPointer();
    	double* posy = posyVec.getArrayPointer();
	for (int step = 0; step < mMaxSteps; step++){
    		double x = 0; double y = 0;
    		for (int walker = 0; walker < mWalkers; walker++){  
			if (posy[walker] < mMinDist || posy[walker] > mMaxDist ) continue;
			if (posx[walker] < mMinDist || posx[walker] > mMaxDist ) continue;

      			double ran = Random::ran0(seed);
      			if (ran <= mWalkProb) {
    				posx[walker] += ml;
      			}
      			else if (ran <= 2*mWalkProb){
    				posx[walker] -= ml;
      			}
      			else if (ran <= 3*mWalkProb){
	        		posy[walker] += ml;
	      		}
	      		else {
	      			posy[walker] -= ml;
		      	}
			// Check if the walker is outside our boundries
			int indexX = (posx[walker] - mMinDist)/ml + 0.5;
			int indexY = (posy[walker] - mMinDist)/ml + 0.5;
			if (indexX < 0 || indexX >= mVecSize) continue;
			if (indexY < 0 || indexY >= mVecSize) continue;
	      		mpxPosComulative [step] += posx[walker];
	      		mpxPosComulative2[step] += posx[walker]*posx[walker];
	      		mpyPosComulative [step] += posy[walker];
	      		mpyPosComulative2[step] += posy[walker]*posy[walker];
			mpxProbabilities[indexX] += 1;
			mpyProbabilities[indexY] += 1;
		}
	}
}


void	Walker::output(string name){
	ofstream outFile1;
	string filename = name + "Walk.dat";
	outFile1.open(filename.c_str());
	for (int i = 0; i < mMaxSteps; i++) {
		double xAverage  = mpxPosComulative [i]/((double) mVecSize);
		double x2Average = mpxPosComulative2[i]/((double) mVecSize);
		double xVariance = x2Average - xAverage*xAverage;
		double yAverage  = mpyPosComulative [i]/((double) mVecSize);
		double y2Average = mpyPosComulative2[i]/((double) mVecSize);
		double yVariance = y2Average - yAverage*yAverage;
		outFile1 << setiosflags(ios::showpoint | ios::uppercase);
		outFile1 << setw(6) << i;
		outFile1 << setw(15) << setprecision(8) << xAverage;
		outFile1 << setw(15) << setprecision(8) << xVariance;
		outFile1 << setw(15) << setprecision(8) << yAverage;
		outFile1 << setw(15) << setprecision(8) << yVariance << endl;
	}
	outFile1.close();
	filename = name + "ProbWalk.dat";
	ofstream outFile2(filename.c_str());
	for (int i = 0; i < mVecSize; i++) {
		outFile2 << setw(6) << i*ml;
		outFile2 << setw(15) << setprecision(8) 
			     << mpxProbabilities[i];
		outFile2 << setw(15) << setprecision(8)
				 << mpyProbabilities[i] << endl;
	}
	outFile2.close();
}





