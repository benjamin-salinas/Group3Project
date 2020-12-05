#include <math.h>
#include "Flow.h"
#include "Solver.h"
#include <iostream>

using namespace std;

Flow::Flow()
{
	dx = Mesh::Lx / Mesh::Nx;
	dy = Mesh::Ly / Mesh::Ny;
	mesh = new Mesh;
	//Set up for (Nx+1)*Ny storage
	uu = new double* [Mesh::Nx + 1];
	ustar = new double* [Mesh::Nx + 1];
	for (int i = 0; i <= Mesh::Nx; i++)
	{
		uu[i] = new double[Mesh::Ny];
		ustar[i] = new double[Mesh::Ny];
		for (int j = 0; j < Mesh::Ny; j++)
		{
			uu[i][j] = 0.;
			ustar[i][j] = 0.;
		}
	}

	//Set up for Nx*(Ny+1) storage
	vv = new double* [Mesh::Nx];    
	vstar = new double* [Mesh::Nx];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		vv[i] = new double[Mesh::Ny + 1];
		vstar[i] = new double[Mesh::Ny + 1];

		for (int j = 0; j <= Mesh::Ny; j++)
		{
			vv[i][j] = 0.;
			vstar[i][j] = 0.;
		}
	}

	//Set up for Nx*Ny storage
	pp = new double* [Mesh::Nx];
	Hx = new double* [Mesh::Nx];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		pp[i] = new double[Mesh::Ny];
		Hx[i] = new double[Mesh::Ny];
		for (int j = 0; j < Mesh::Ny; j++)
		{
			pp[i][j] = 0.;
			Hx[i][j] = 0.;
		}
	}
	
	//Set up for Nx*(Ny-1) storage
	Hy = new double* [Mesh::Nx];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		Hy[i] = new double[Mesh::Ny - 1];
		for (int j = 0; j < Mesh::Ny - 1; j++)
		{
			Hy[i][j] = 0.;
		}
	}


	//Set up left boundary condition
	for (int j = 0; j < Mesh::Ny; j++) uu[0][j] = exp(-(mesh->yc[j] - 0.5) * (mesh->yc[j] - 0.5) / R / R); 
	
	alpha = dt / (2. * Re * dx * dx);

}

Flow::~Flow() //Ning delete?
{
	for (int i = 0; i <= Mesh::Nx; i++)
	{
		delete[] uu[i];
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		delete[] vv[i];
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		delete[] pp[i];
	}
	delete[] uu;
	delete[] vv;
	delete[] pp;
}

double Flow::u(int i_, int j_)
{
	if (i_ == 0) return exp(-(mesh->yf[j_] - 0.5) * (mesh->yf[j_] - 0.5) / R / R);
	//else if (i_ == -1) return 0.;
	else if (j_ == -1) return uu[i_][j_ + 1];
	else if (j_ == Mesh::Ny) return uu[i_][j_ - 1];
	else if (i_ == Mesh::Nx + 1) return uu[i_ - 2][j_];
	else return uu[i_][j_];
}

double Flow::v(int i_, int j_)
{
	if (i_ == -1) return -vv[i_ + 1][j_];
	else if (j_ == 0 || j_ == Mesh::Ny) return 0.;
	else if (i_ == Mesh::Nx) return vv[i_ - 1][j_];
	else return vv[i_][j_];
}

double Flow::p(int i_, int j_)
{
	if (i_ == -1) return pp[i_ + 1][j_];
	else if (j_ == -1) return pp[i_][j_ + 1];
	else if (j_ == Mesh::Ny) return pp[i_][j_ - 1];
	else if (i_ == Mesh::Nx) return pp[i_ - 1][j_];
	else return pp[i_][j_];
}

void Flow::updateFlow()
{
	double t = 0.;
	int count = 0;
	while (t < Flow::T)
	{
		


	
		//Step 1 to get ustar and vstar
		getustar();
		getvstar();
		//Correct for mass conservation Ning
		double sum1 = 0., sum2 = 0., kk;
		for (int j = 0; j < Mesh::Ny; j++)
		{
			sum1 += ustar[0][j];
			sum2 += ustar[Mesh::Nx][j];
		}
		kk = sum1 / sum2;
		for (int j = 0; j < Mesh::Ny; j++)
			ustar[Mesh::Nx][j] = kk * ustar[Mesh::Nx][j];
		//Step 2
		getpp();

		//Step 3
		getuu();
		getvv();
			t += Flow::dt;
			count++;
		//if (count % 1 == 0) 
			cout << "Time: " << t << endl;

	}
}

double Flow::getHx(int i_, int j_)
{
	return
		- (u(i_ + 1, j_) - u(i_ - 1, j_)) * (u(i_ + 1, j_) + u(i_ - 1, j_) + 2. * u(i_, j_)) / (4. * dx)
			+ ((u(i_, j_) + u(i_, j_ - 1)) * (v(i_ - 1, j_) + v(i_, j_))
				- (u(i_, j_ + 1) + u(i_, j_)) * (v(i_ - 1, j_ + 1) + v(i_, j_ + 1))) / (4. * dy); //Ning
}

double Flow::getHy(int i_, int j_)
{
	return
		- ((u(i_ + 1, j_ - 1) + u(i_ + 1, j_)) * (v(i_, j_) + v(i_ + 1, j_))
			- (u(i_, j_ - 1) + u(i_, j_)) * (v(i_ , j_) + v(i_ - 1, j_))) / (4. * dx)
				-(v(i_, j_ + 1) - v(i_, j_ - 1)) * (v(i_, j_ + 1) + v(i_, j_ - 1) + 2. * v(i_, j_)) / (4. * dy); //Ning
				
}

void Flow::getustar()
{
	double** ustar2 = new double*[Mesh::Ny];
	for (int j = 0; j < Mesh::Ny; j++)
	{
		ustar2[j] = new double[Mesh::Nx];
		ustar2[j] = getustar2(j);
	}
	
	double* lower = new double[Mesh::Ny];
	double* diag = new double[Mesh::Ny];
	double* upper = new double[Mesh::Ny];

	for (int j = 0; j < Mesh::Ny; j++) 
	{
		if (j == 0) diag[j] = 1 - 3. * alpha;	// Add vector value at the boundary
		else if (j == Mesh::Ny - 1) diag[j] = 1 - 3. * alpha;
		else
		{
			lower[j] = -alpha;
			diag[j] = 1 - 2. * alpha;
			upper[j] = -alpha;
		}
	}

	for (int i = 1; i <= Mesh::Nx; i++) //
	{
		double* rhs = new double[Mesh::Ny];
		for (int j = 0; j < Mesh::Ny; j++) rhs[j] = ustar2[j][i - 1];	
		Solver* solver2 = new Solver;
		ustar[i] = solver2->GaussElimination(diag, upper, lower, rhs, Mesh::Ny);
	}
	/*for (int j = 0; j < Mesh::Ny; j++)
	{
		ustar[0] = new double[Mesh::Ny];
		ustar[0][j] = 0.;
	}*/
	//get ustar from delta_ustar
	for (int i = 0; i < Mesh::Nx + 1; i++)
	{
		for (int j = 0; j < Mesh::Ny; j++) ustar[i][j] += u(i, j);
	}
}

double* Flow::getustar2(int j)
{
	double* rhsu2 = new double[Mesh::Nx];
	double* lower = new double[Mesh::Nx];
	double* diag = new double[Mesh::Nx];
	double* upper = new double[Mesh::Nx];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		if (i == Mesh::Nx - 1) lower[i] = -2. * alpha;// Ning
		else
		{
			lower[i] = -alpha;
			diag[i] = 1 - 2. * alpha;
			upper[i] = -alpha;
		}
	}

	for (int i = 1; i <= Mesh::Nx; i++)
	{
		double newH = getHx(i, j);
		rhsu2[i - 1] = dt / 2. * (3 * newH - Hx[i - 1][j])
			+ dt / Re * (u(i + 1, j) + u(i - 1, j) - 2. * u(i, j)) / dx / dx
			+ dt / Re * (u(i, j + 1) + u(i, j - 1) - 2. * u(i, j)) / dy / dy;
		Hx[i - 1][j] = newH;
	}
	
	Solver* solver1 = new Solver;
	return solver1->GaussElimination(diag, upper, lower, rhsu2, Mesh::Nx); //Ning
}


//Ning
void Flow::getvstar()
{
	double** vstar2 = new double*[Mesh::Ny - 1];
	for (int j = 0; j < Mesh::Ny - 1; j++)
	{
		vstar2[j] = new double[Mesh::Nx];
		vstar2[j] = getustar2(j);
	}
	
	double* lower = new double[Mesh::Ny - 1];
	double* diag = new double[Mesh::Ny - 1];
	double* upper = new double[Mesh::Ny - 1];

	for (int j = 0; j < Mesh::Ny - 1; j++)
	{
			lower[j] = -alpha;
			diag[j] = 1 - 2. * alpha;
			upper[j] = -alpha;

	}

	for (int i = 0; i < Mesh::Nx; i++) 
	{
		double* rhs = new double[Mesh::Ny - 1];
		for (int j = 0; j < Mesh::Ny - 1; j++) rhs[j] = vstar2[j][i];	
		Solver* solver2 = new Solver;
		vstar[i] = vstar[i] + 1;
		vstar[i] = solver2->GaussElimination(diag, upper, lower, rhs, Mesh::Ny - 1);
	}
	/*for (int i = 0; i < Mesh::Nx; j++)
	{
		vstar[i][0] = 0.;
		vstar[i][Ny] = 0.;
	}*/
	//get vstar from delta_vstar
	for (int i = 0; i < Mesh::Nx; i++)
	{
		for (int j = 0; j < Mesh::Ny + 1; j++) vstar[i][j] += v(i, j);
	}
}

double* Flow::getvstar2(int j)  //j from 1 to Ny - 1
{
	double* rhsv2 = new double[Mesh::Nx];
	double* lower = new double[Mesh::Nx];
	double* diag = new double[Mesh::Nx];
	double* upper = new double[Mesh::Nx];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		if (i == Mesh::Nx - 1) diag[i] = 1 - 3. * alpha;
		else if (i == 0) diag[i] = 1 - alpha;
		else
		{
			lower[i] = -alpha;
			diag[i] = 1 - 2. * alpha;
			upper[i] = -alpha;
		}
	}

	for (int i = 0; i < Mesh::Nx; i++)
	{
		double newH = getHy(i, j);
		rhsv2[i] = dt / 2. * (3 * newH - Hy[i][j - 1])
			+ dt / Re * (v(i + 1, j) + v(i - 1, j) - 2. * v(i, j)) / dx / dx
			+ dt / Re * (v(i, j + 1) + v(i, j - 1) - 2. * v(i, j)) / dy / dy;
		Hy[i][j - 1] = newH;
	}
	
	Solver* solver3 = new Solver;
	return solver3->GaussElimination(diag, upper, lower, rhsv2, Mesh::Nx);
}


void Flow::getpp()
{

	Solver* solver5 = new Solver;
	double** Ae = new double* [Mesh::Nx];
	double**  Aw = new double* [Mesh::Nx];
	double**  As = new double* [Mesh::Nx];
	double**  An = new double* [Mesh::Nx];
	double**  Ap = new double* [Mesh::Nx];
	double**  rhs = new double* [Mesh::Nx];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		Ae[i] = new double[Mesh::Ny];
		Aw[i] = new double[Mesh::Ny];
		As[i] = new double[Mesh::Ny];
		An[i] = new double[Mesh::Ny];
		Ap[i] = new double[Mesh::Ny];
		rhs[i] = new double[Mesh::Ny];
		for (int j = 0; j < Mesh::Ny; j++)
		{
			Ae[i][j] = 1. / dx / dx;
			Aw[i][j] = 1. / dx / dx;
			As[i][j] = 1. / dy / dy;
			An[i][j] = 1. / dy / dy;
			Ap[i][j] = -2. / dx / dx - 2. / dy / dy;
			rhs[i][j] = dt * ((ustar[i + 1][j] - ustar[i][j]) / dx + (vstar[i][j + 1] - vstar[i][j]) / dy);
		}
	}
	pp = solver5->SOR(Ae, Aw, An, As, Ap, rhs, Mesh::Nx, Mesh::Ny);
	//Reshape for pressure to a 2D matrix
	

}

void Flow::getuu()
{

	for (int i = 1; i <= Mesh::Nx; i++)
	{
		for (int j = 0; j < Mesh::Ny; j++)
		{
			uu[i][j] = ustar[i][j] - dt * (p(i, j) - p(i - 1, j)) / dx;
		}
	}
}

void Flow::getvv()
{

	for (int i = 0; i < Mesh::Nx; i++)
	{
		for (int j = 1; j < Mesh::Ny; j++)
		{
			vv[i][j] = vstar[i][j] - dt * (p(i, j) - p(i, j - 1)) / dy;
		}
	}
}