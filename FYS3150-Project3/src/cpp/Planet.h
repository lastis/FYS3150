class Planet {
public:
	Planet(const char* name, double aphelion, double velMin,
	       double aphAngle, double mass                 );

	double 	getMass();
	double	getAphelion();
	double	getVelMin();
	double	getAphAngle();
	const char* getName();
private:
	// Variables
	const char*	mName;
	double		mAphelion;
	double		mVelMin;
	double 		mAphAngle;
	double 		mMass;
};
