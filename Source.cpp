#include <iostream>
#include "Flow.h"
#include "Ensemble.h"

using namespace std;

/* Initialization of basic parameters */
//Scale of the domain [m]
double Mesh::Lx = 2.;
double Mesh::Ly = 1.;

//Number of cells at each dimension
int Mesh::Nx = 128;
int Mesh::Ny = 64;

//Simulation time
double Flow::T = 100.;

//Time step
double Flow::dt = 0.001;

//Reynold number
double Flow::Re = 1.;

//Stokes number
double Ensemble::St = 1.;

//Mass flow rate
double Ensemble::m_in = 1000;

int main()
{
	Flow* flow = new Flow;
	flow->updateFlow();
	/*
	double t = 0.;
	while (t < Flow::T)
	{
		flow->updateFlow();


		t += Flow::dt;
		if ((int)(t / Flow::dt) % 1 == 0) cout << "Time: " << t << endl;
	
	}
*/
	return 0;
}