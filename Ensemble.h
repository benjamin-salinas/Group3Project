#ifndef _Ensemble_
#define _Ensemble_
#include <vector>
#include <time.h>
#include "Flow.h"

using namespace std;

struct Particle
{
	double x, y, u, v;
};

class Ensemble
{
public:
	Flow* flow;
	double dt, dx, dy;
	vector<Particle> P;
	static double St; //Stokes number;
	static double m_in; //Particle inlet mass flow rate;
	int n_in; //Particle inlet number per time step;

	double R = 0.05;

public:
	Ensemble(Flow* flow_);
	~Ensemble();
	void updateEnsemble();
	void particleinsert(); //generate particle at the left boundary
	void Runge_Kutta_explicit();
	void particlecheck(); //remove particle leaving the domain

	double Runge_Kutta(double q);
	double getVx(Particle P_);
	double getVy(Particle P_);

	inline double getuleftBC(double y_)
	{
		return exp(-(y_ - 0.5) * (y_ - 0.5) / R / R);
	}
};

#endif