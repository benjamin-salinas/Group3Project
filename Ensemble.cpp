#include "Ensemble.h"

Ensemble::Ensemble(Flow* flow_)
{
	flow = flow_;
	dt = Flow::dt;
	dx = flow_->dx;
	dy = flow_->dy;
}

Ensemble::~Ensemble()
{
	delete flow;
}

void Ensemble::updateEnsemble()
{
	srand(time(0) + rand());
	n_in = m_in * dt;
	particleinsert();
	Runge_Kutta_explicit();
	particlecheck();
}

void Ensemble::particleinsert()
{
	for (int i = 0; i < n_in; i++)
	{
		Particle* temp;
		temp->x = 0.;
		temp->y = 0.45 + 0.1 * rand() / RAND_MAX;
		temp->u = getuleftBC(temp->y);
		temp->v = 0.;

		P.push_back(*temp);
	}
}

void Ensemble::Runge_Kutta_explicit()
{
	int j = 0;
	for (j = 0; j < P.size; j++)
	{
		P[j].x += Runge_Kutta(P[j].u);
		P[j].y += Runge_Kutta(P[j].v);
		P[j].u += Runge_Kutta(1. / St * (getVx(P[j]) - P[j].u));
		P[j].v += Runge_Kutta(1. / St * (getVy(P[j]) - P[j].v));
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

double Ensemble::getVx(Particle P_)
{
	int ix = P_.x / dx;
	int iy = (P_.y - dy / 2.) / dy;

	//intepolate to get the velocity in the field
	double temp1 = flow->u(ix, iy) +
		(P_.x - ix * dx) / dx * (flow->u(ix, iy) + flow->u(ix + 1, iy));
	double temp2 = flow->u(ix, iy + 1) +
		(P_.x - ix * dx) / dx * (flow->u(ix, iy + 1) + flow->u(ix + 1, iy + 1));
	return temp1 + (temp2 - temp1) * (P_.y - iy * dy - dy / 2);
}

double Ensemble::getVy(Particle P_)
{
	int ix = (P_.x - dx / 2.) / dx;
	int iy = P_.y / dy;

	//intepolate to get the velocity in the field
	double temp1 = flow->v(ix, iy) +
		(P_.y - iy * dy) / dy * (flow->v(ix, iy) + flow->v(ix, iy + 1));
	double temp2 = flow->v(ix + 1, iy) +
		(P_.y - iy * dy) / dy * (flow->v(ix + 1, iy) + flow->v(ix + 1, iy + 1));
	return temp1 + (temp2 - temp1) * (P_.x - ix * dx - dx / 2);
}

void Ensemble::particlecheck()
{
	int j = 0;
	for (j = 0; j < P.size; j++)
	{
		if (P[j].x > Mesh::Lx || P[j].y > Mesh::Ly || P[j].y < 0. || P[j].x < 0.)
			P.erase(P.begin + j);
	}
}