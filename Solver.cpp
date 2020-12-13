#include "Solver.h"

Solver::Solver(int Nx)
{
	ss = new double*[1];
	ss[0] = new double[Nx];
	for (int i = 1; i <= Nx - 1; i++) ss[0][i] = 0.;
	ix = 1;
	iy = Nx;
}

Solver::Solver(int Nx, int Ny)
{
	ss = new double* [Nx];
	for (int i = 0; i < Nx; i++)
	{
		ss[i] = new double[Ny];
		for (int j = 0; j < Ny; j++) ss[i][j] = 0.;
	}
	ix = Nx; iy = Ny;
}

Solver::~Solver()
{
	if (ix == 1)
	{

		//for (int i = 0; i < ix; i++) delete[] ss[i];
		//delete ss;
	}
	else
	{
	//for (int i = 0; i < ix; i++) delete[] ss[i];
	//delete[] ss;
	}
}

double* Solver::GaussElimination(double* diag, double* upper, double* lower, double* rhs, int N)
{
	for (int i = 1; i <= N - 1; i++)
	{
		diag[i] = diag[i] - lower[i] * upper[i - 1] / diag[i - 1];
		rhs[i] = rhs[i] - rhs[i - 1] * lower[i] / diag[i - 1];
	}
	ss[0][N - 1] = rhs[N - 1] / diag[N - 1];
	for (int i = N - 2; i >= 0; i--)
		ss[0][i] = (rhs[i] - upper[i] * ss[0][i + 1]) / diag[i];
	return ss[0];
}



double** Solver::SOR(double** Ae, double** Aw, double** An, double** As, double** Ap, double** rhs, int Nxx, int Nyy )
{
	double** ss0 = new double* [Nxx];
	double tol = 1e-7, w = 1.8;
	int iw, jw, ie, je, is, js, in, jn;
	int tt = 0;
	for(int i = 0; i < Nxx; i++) ss0[i] = new double [Nyy];
	double sum, res = 1.0;
	while (res > tol && tt <= 100000)
	{
		sum = 0.0;
		for (int i = 0; i < Nxx; i++)
			for (int j = 0; j < Nyy; j++)
				ss0[i][j] = ss[i][j];
		for(int i = 0; i < Nxx; i++)
		{
			for(int j = 0; j < Nyy; j++)
			{
				//ss0[i][j] = ss[i][j];
				iw = i - 1; jw = j;
				if (iw == -1) iw = 0;
				ie = i + 1; je = j;
				if (ie == Nxx) ie = Nxx - 1;
				is = i; js = j - 1;
				if (js == -1) js = 0;
				in = i; jn = j + 1;
				if (jn == Nyy) jn = Nyy -1;
				ss[i][j] = (1 - w) * ss[i][j] + w / Ap[i][j] * (rhs[i][j] -Ae[i][j] * ss[ie][je] - Aw[i][j] * ss[iw][jw]
					-An[i][j] * ss[in][jn] - As[i][j] * ss[is][js]);
				//ss[i][j] =1 / Ap[i][j] * (rhs[i][j] - Ae[i][j] * ss0[ie][je] - Aw[i][j] * ss0[iw][jw]
					//- An[i][j] * ss0[in][jn] - As[i][j] * ss0[is][js]);
			}
		}
		for (int i = 0; i < Nxx; i++)
		{
			for (int j = 0; j < Nyy; j++)
				sum += (ss[i][j] - ss0[i][j]) * (ss[i][j] - ss0[i][j]);
			//cout << ss[i][0] << " " << ss[i][63] << endl;
		}
		res = sqrt(sum);
		//cout << sum << endl;
		tt = tt + 1;
		//if (tt >= 30000)
		//cout << sum << endl;
	}	
	for (int i = 0; i < Nxx; i++) delete[] ss0[i];
	delete[] ss0;
	cout << tt << endl;
	return ss;
}

double** Solver::CG(double** Ae, double** Aw, double** An, double** As, double** Ap, double** rhs, int Nxx, int Nyy)
{
	double** dd = new double* [Nxx];
	double** epsilon = new double* [Nxx];
	double tol = 1e-7;
	int k = 0, i1, i2, j1, j2;
	int tt = 0;
	for (int i = 0; i < Nxx; i++)
	{
		ss[i] = new double[Nyy];
		dd[i] = new double[Nyy];
		epsilon[i] = new double[Nyy];
		for (int j = 0; j < Nyy; j++)
		{
			epsilon[i][j] = 0;
			ss[i][j] = 0;
		}
	}
	double rho = 0.0, rho_old, alpha;
	for (int i = 0; i < Nxx; i++)
		for (int j = 0; j < Nyy; j++)
			rho += rhs[i][j] * rhs[i][j];
	while (sqrt(rho) > tol&& tt <= 100000)
	{
		k++;
		for (int i = 0; i < Nxx; i++)
			for (int j = 0; j < Nyy; j++)
			{
				if (k == 1)
					dd[i][j] = rhs[i][j];
				else
					dd[i][j] = rhs[i][j] + rho / rho_old * dd[i][j];
			}
		alpha = 0;
		for (int i = 0; i < Nxx; i++)
			for (int j = 0; j < Nyy; j++)
			{
				i1 = i - 1;
				if (i1 == -1) i1 = 0;
				i2 = i + 1;
				if (i2 == Nxx) i2 = Nxx - 1;
				j1 = j - 1;
				if (j1 == -1) j1 = 0;
				j2 = j + 1;
				if (j2 == Nyy) j2 = Nyy - 1;
				epsilon[i][j] = Ap[i][j] * dd[i][j] + Aw[i][j] * dd[i1][j] + Ae[i][j] * dd[i2][j]
					+ As[i][j] * dd[i][j1] + An[i][j] * dd[i][j2];
				alpha += epsilon[i][j] * dd[i][j];
			}
		alpha = rho / alpha;
		for (int i = 0; i < Nxx; i++)
			for (int j = 0; j < Nyy; j++)
			{
				ss[i][j] += alpha * dd[i][j];
				rhs[i][j] -= alpha * epsilon[i][j];
			}
		rho_old = rho;
		rho = 0.0;
		for (int i = 0; i < Nxx; i++)
			for (int j = 0; j < Nyy; j++)
				rho += rhs[i][j] * rhs[i][j];
		tt++;
		//if (tt >= 30000)
		//cout << sum << endl;
	}
	for (int i = 0; i < Nxx; i++) delete[] dd[i];
	for (int i = 0; i < Nxx; i++) delete[] epsilon[i];
	delete[] dd;
	delete[] epsilon;
	cout << tt << endl;
	return ss;
	return 0;
}