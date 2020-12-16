#ifndef _Flow_
#define _Flow_
#include "Mesh.h"

using namespace std;

class Flow
{
public:
	Mesh* mesh;
	double** uu, ** vv; //velocity
	double** pp; //pressure
	double dx, dy;
	double** ustar, ** vstar;
	double** ustar2, ** vstar2;

	static double Re; //Reynold number
	double** Hx; //n-1 time step convective term to update u;
	double** Hy; //n-1 time step convective term to update v;

	double R = 0.05; //jet radius
	static double T; //Simulation time
	static double dt; //Simulation timestep
	double alpha;

	int count = 0;
	double L2_D;

	static bool restartflow; //switch to use restart data
	double step1 = 0., step2 = 0., step3 = 0.;

public:
	Flow();
	~Flow();
	//Functions to update boundary conditions:
	double u(int i_, int j_);
	double v(int i_, int j_);
	double p(int i_, int j_);

	void updateFlow();
	double getHx(int i_, int j_); //get convection term
	double getHy(int i_, int j_); //get convection term

	void getustar(); //return delta ustar as a colume vector
	void getvstar(); //return delta vstar as a colume vector

	void getustar2(); //return delta ustar2 as a colume vector
	void getvstar2(); //return delta vstar2 as a colume vector

	void getuu();
	void getvv();
	void getpp();

	void writeRestart();
	void readRestart();

};

#endif