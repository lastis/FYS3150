#ifndef __WALKER_H_INCLUDED
#define __WALKER_H_INCLUDED 
using namespace std;

class Walker {
private:
	int	mWalkers;	// Number of walkers
	int	mMaxSteps;	// Steps the walkers should take
	int 	mVecSize;	// Vector size of mProb

	double 	mdt; 		// time step dt
	double	ml;		// step length
	double 	mMinDist;
	double 	mMaxDist;
	double	mWalkProb;
	Vector 	mxPosComulative;
	Vector 	mxPosComulative2;
	Vector	mxProbabilities;
	double*	mpxPosComulative;
	double*	mpxPosComulative2;
	double*	mpxProbabilities;
	Vector 	myPosComulative;
	Vector 	myPosComulative2;
	Vector	myProbabilities;
	double*	mpyPosComulative;
	double*	mpyPosComulative2;
	double*	mpyProbabilities;		
	bool	mUseGauss;
public:
	Walker(int walkers, double walkProb, double time, double dt, 
	       double minDist, double maxDist);
	void 	start1D();
	void	start2D();
	void 	mcSampling();
	void	mcSampling2();
	void	output(string name);
	void	setTime(double time, double dt);
	void	useGauss(bool use);
};
#endif
