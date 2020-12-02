#include <math.h>
#include "Flow.h"
#include "Solver.h"

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
	vv = new double* [Mesh::Ny];
	vstar = new double* [Mesh::Ny];
	for (int i = 0; i < Mesh::Nx; i++)
	{
		vv[i] = new double[Mesh::Ny + 1];
		Hy[i] = new double[Mesh::Ny + 1];
		vstar[i] = new double[Mesh::Ny + 1];

		for (int j = 0; j <= Mesh::Ny; j++)
		{
			vv[i][j] = 0.;
			Hy[i][j] = 0.;
			vstar[i][j] = 0.;
		}
	}

	//Set up for Nx*Ny storage
	pp = new double* [Mesh::Nx];
	Hx = new double* [Mesh::Nx + 1];
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
	for (int j = 0; j < Mesh::Ny; j++) uu[0][j] = exp(-(mesh->yf[j] - 0.5) * (mesh->yf[j] - 0.5) / R / R);
	
	alpha = dt / (2. * Re * dx * dx);

}

Flow::~Flow()
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
	//Step 1 to get ustar and vstar
	getustar();
	getvstar();
	//Correct for mass conservation
	
	//Step 2
	getpp();

	//Step 3
	getuu();
	getvv();
}

double Flow::getHx(int i_, int j_)
{
	return
		-u(i_ + 1, j_) - u(i_ - 1, j_) * (u(i_ + 1, j_) + u(i_ - 1, j_) + 2. * u(i_, j_)) / (4. * dx)
			+ (u(i_, j_) + u(i_, j_ - 1) * (v(i_ - 1, j_) + v(i_, j_))
				- (u(i_, j_ + 1) + u(i_, j_)) * (v(i_ - 1, j_ + 1) + v(i_, j_ + 1))) / (4. * dy);
}

double Flow::getHy(int i_, int j_)
{
		
}

void Flow::getustar()
{
	double** ustar2 = new double*[Mesh::Ny];
	for (int j = 0; j < Mesh::Ny; j++)
	{
		ustar2[j] = new double[Mesh::Nx];
		ustar2[j] = getustar2(j);
	}
	
	double* rhsu = new double[Mesh::Ny];
	double* lower = new double[Mesh::Ny];
	double* diag = new double[Mesh::Ny];
	double* upper = new double[Mesh::Ny];

	for (int j = 0; j < Mesh::Nx; j++)
	{
		if (j == 0);	// Add vector value at the boundary
		else if (j == Mesh::Ny - 1);
		else
		{
			lower[j] = -alpha;
			diag[j] = 1 - alpha;
			upper[j] = -alpha;
		}
	}

	for (int i = 1; i < Mesh::Nx; i++)
	{
		ustar[i] = new double[Mesh::Ny];
		double* rhs = new double[Mesh::Ny];
		for (int j = 0; j < Mesh::Ny; j++) rhs[j] = ustar2[j][i-1];	
		Solver* solver2 = new Solver;
		ustar[i] = solver2->GaussElimination(diag, upper, lower, rhs);
	}
	for (int j = 0; j < Mesh::Ny; j++)
	{
		ustar[0] = new double[Mesh::Ny];
		ustar[0][j] = 0.;
	}
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
		if (i == 0);	// Add vector value at the boundary
		else if (i == Mesh::Nx - 1);
		else
		{
			lower[i] = -alpha;
			diag[i] = 1 - alpha;
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
	return solver1->GaussElimination(diag, upper, lower, rhsu2);
}

void Flow::getpp()
{

	Solver* solver3 = new Solver;
	double* temp = solver3->SOR();
	//Reshape for pressure to a 2D matrix
	

}

void Flow::getuu()
{

	for (int i = 0; i <= Mesh::Nx; i++)
	{
		for (int j = 0; j < Mesh::Ny; j++)
		{
			uu[i][j] = ustar[i][j] - dt * (p(i, j) - p(i - 1, j)) / dx;
		}
	}
}