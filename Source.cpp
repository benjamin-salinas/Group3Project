#include "Flow.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Ensemble.h"

using namespace std;

/* Initialization of basic parameters */
//Scale of the domain [m]
double Mesh::Lx = 4.;
double Mesh::Ly = 1.;

//Number of cells at each dimension
int Mesh::Nx = 256;
int Mesh::Ny = 64;

//Simulation time
double Flow::T = 60.000;

//Time step
double Flow::dt = 0.001;

//Reynold number
double Flow::Re = 10000.;

//Stokes number
double Ensemble::St = 1.;

//Mass flow rate
double Ensemble::m_in = 1000;

//switch to read restart data
bool Flow::restartflow = false;
bool Ensemble::restartensemble = false;

int main()
{
	Flow* flow = new Flow;
	Ensemble* ensemble = new Ensemble(flow);

	double t = 0.; //restart time
	int count = (int)(t/Flow::dt);
	cout << "count: " << count << endl;
	while (t < Flow::T)
	{
		flow->updateFlow();
		ensemble->updateEnsemble();
		count++;
		//output flow data
		ofstream fout;
			fout.open("./post/flow/niubility"+to_string(count)+".dat", ios::out | ios::trunc);
			fout << "VARIABLES = \"X\", \"Y\", \"VELOCITY_X\", \"VELOCITY_Y\", \"PRESSURE\"" << endl;
			fout << "ZONE I=" << Mesh::Nx << ", J=" << Mesh::Ny << ", F=POINT" << endl;
			Mesh* mesh1 = new Mesh;
			for (int j = 0; j < Mesh::Ny; j++)
				for (int i = 0; i < Mesh::Nx; i++)
				{
					fout << mesh1->xc[i] << " " << mesh1->yc[j] << " " << 0.5 * (flow->uu[i][j] +
						flow->uu[i + 1][j]) << " " << 0.5 * (flow->vv[i][j] + flow->vv[i][j + 1]) << " " <<
						flow->pp[i][j] << endl;
				}
			fout.close();
			
			//output particles data
			cout << "output particles number: " << ensemble->P.size() << endl;
			ofstream pout;
			pout.open("./post/ensemble/ensemble" + to_string(count) + ".dat", ios::out | ios::trunc);
			pout << "VARIABLES = \"X\", \"Y\", \"VELOCITY_X\", \"VELOCITY_Y\"" << endl;
			pout << "ZONE I=" << ensemble->P.size() << ", F=POINT" << endl;
			for (int i = 0; i < ensemble->P.size(); i++)
					pout << ensemble->P[i]->x << " " << ensemble->P[i]->y << " " << ensemble->P[i]->u << " "  <<
					ensemble->P[i]->v << endl;
				
			pout.close();
			t += Flow::dt;
			if ((int)(t / Flow::dt) % 1 == 0) cout << "Time: " << t << endl;
			if (((int)(t / Flow::dt)) % 200 == 0)
			{
				flow->writeRestart();
				ensemble->writeRestart();
			}

			//output L2
			ofstream lout;
			lout.open("./post/lnorm.dat", ios::app);
			//lout << "VARIABLES = \"L2\"" << endl;
			lout << sqrt(flow->L2_D) << endl;
			lout.close();
	}
	



	return 0;
}