#include <iostream>
#include <fstream>
#include<iomanip>
#include<sstream>
#include "Ensemble.h"

Ensemble::Ensemble(Flow* flow_)
{
	flow = flow_;
	dt = Flow::dt;
	dx = flow_->dx;
	dy = flow_->dy;
	n_in = 0;

	if (restartensemble == true) readRestart(); 
}

Ensemble::~Ensemble()
{
	delete flow;
}

void Ensemble::updateEnsemble()
{
	srand(time(0) + rand());
	//Generate number of particles inserted in the domain at each timestep
	double ntemp = m_in * dt;
	if (ntemp < 1)
	{
		double ran = 1. * rand() / RAND_MAX;
		if (ran <= ntemp) n_in = 1;
		else n_in = 0;
	}
	else n_in = m_in * dt;
	particleinsert();
	Runge_Kutta_explicit();
	particlecheck();

}

void Ensemble::particleinsert()
{
	for (int i = 0; i < n_in; i++)
	{
		Particle* temp = new Particle;
		temp->x = 0.;
		temp->y = 0.45 + 0.1 * rand() / RAND_MAX;
		temp->u = getuleftBC(temp->y);
		temp->v = 0.;

		P.push_back(temp);
	}
}

void Ensemble::Runge_Kutta_explicit()
{
	int j = 0;
	for (j = 0; j < P.size(); j++)
	{
		P[j]->x += Runge_Kutta(P[j]->u);
		P[j]->y += Runge_Kutta(P[j]->v);
		if (P[j]->x > Mesh::Lx || P[j]->y > Mesh::Ly || P[j]->y < 0. || P[j]->x < 0.) continue;
		P[j]->u += Runge_Kutta(1. / St * (getVx(P[j]) - P[j]->u));
		P[j]->v += Runge_Kutta(1. / St * (getVy(P[j]) - P[j]->v));
	}
}

double Ensemble::Runge_Kutta(double q)
{
	double K1 = q;
	double K2 = q + dt / 2 * K1;
	double K3 = q + dt / 2 * K2;
	double K4 = q + dt * K3;

	return (K1 + 2 * K2 + 2 * K3 + K4) * dt / 6.;
}

double Ensemble::getVx(Particle* P_)
{
	int ix = P_->x / dx;
	int iy = (P_->y - dy / 2.) / dy;

	//intepolate to get the velocity in the field
	double temp1 = flow->u(ix, iy) +
		(P_->x - ix * dx) / dx * (flow->u(ix, iy) + flow->u(ix + 1, iy));
	double temp2 = flow->u(ix, iy + 1) +
		(P_->x - ix * dx) / dx * (flow->u(ix, iy + 1) + flow->u(ix + 1, iy + 1));
	return temp1 + (temp2 - temp1) * (P_->y - iy * dy - dy / 2);
}

double Ensemble::getVy(Particle* P_)
{
	int ix = (P_->x - dx / 2.) / dx;
	int iy = P_->y / dy;

	//intepolate to get the velocity in the field
	double temp1 = flow->v(ix, iy) +
		(P_->y - iy * dy) / dy * (flow->v(ix, iy) + flow->v(ix, iy + 1));
	double temp2 = flow->v(ix + 1, iy) +
		(P_->y - iy * dy) / dy * (flow->v(ix + 1, iy) + flow->v(ix + 1, iy + 1));
	return temp1 + (temp2 - temp1) * (P_->x - ix * dx - dx / 2);
}

void Ensemble::particlecheck()
{
	int j = 0;
	for (j = 0; j < P.size(); j++)
	{
		if (P[j]->x > Mesh::Lx || P[j]->y > Mesh::Ly || P[j]->y < 0. || P[j]->x < 0.)
			P.erase(P.begin() + j);
	}
}

void Ensemble::writeRestart()
{
	ofstream restart;
	restart.open("restart_p.csv", ios::out);
	for (int i = 0; i < P.size(); i++)
		restart << P[i]->x << "," << P[i]->y << "," << P[i]->u << "," << P[i]->v << endl;
	restart.close();
}

void Ensemble::readRestart()
{
	//count number of rows in the file
	ifstream file;
	file.open("restart_p.csv", ios::in);

	int n = 0;
	string temp;
	while (getline(file, temp)) n++;	
	int LINES = n;
	cout << "n:" << LINES << endl;
	file.close();
	file.open("restart_p.csv", ios::in);
	
	for (int row = 0; row < LINES; row++) {

		string line;
		getline(file, line);
		stringstream iss(line);
		string val;

		double tempx, tempy, tempu, tempv;

		for (int col = 0; col < 4; col++) {

			getline(iss, val, ',');

			stringstream convertor(val);

			double inv;
			convertor >> inv;
			if (col == 0) tempx = inv;
			else if (col == 1) tempy = inv;
			else if (col == 2) tempu = inv;
			else if (col == 3) tempv = inv;
		}
		Particle* temp = new Particle;
		temp->x = tempx;
		temp->y = tempy;
		temp->u = tempu;
		temp->v = tempv;

		P.push_back(temp);
	}
	file.close();
}