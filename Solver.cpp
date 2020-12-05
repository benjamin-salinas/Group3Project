#include "Solver.h"
double* Solver::GaussElimination(double* diag, double* upper, double* lower, double* rhs, int N)
{
	double* ss = new double[N];
	for (int i = 1; i <= N - 1; i++)
	{
		diag[i] = diag[i] - lower[i] * upper[i - 1] / diag[i - 1];
		rhs[i] = rhs[i] - rhs[i - 1] * lower[i] / diag[i - 1];
	}
	ss[N - 1] = rhs[N - 1] / diag[N - 1];
	for (int i = N - 2; i >= 0; i--)
		ss[i] = (rhs[i] - upper[i] * ss[i + 1]) / diag[i];
	return ss;
}



double** Solver::SOR(double** Ae, double** Aw, double** An, double** As, double** Ap, double** rhs, int Nxx, int Nyy )
{
	double** ss = new double* [Nxx];
	double** ss0 = new double* [Nxx];
	double tol = 1e-7, w = 1.71;
	int iw, jw, ie, je, is, js, in, jn;
	for(int i = 0; i < Nxx; i++)
	{
		ss[i] = new double [Nyy];
		ss0[i] = new double [Nyy];
		for(int j = 0; j < Nyy; j++)
			ss[i][j] = 0;
	}
	double sum, res = 1.0;
	while (res > tol)
	{
		sum = 0.0;
		for(int i = 0; i < Nxx; i++)
		{
			for(int j = 0; j < Nyy; j++)
			{
				ss0[i][j] = ss[i][j];
				iw = i - 1; jw = j;
				if (iw == -1) iw = 0;
				ie = i + 1; je = j;
				if (ie == Nxx) ie = Nxx - 1;
				is = i; js = Nyy - 1;
				if (js == -1) js = 0;
				in = i; jn = 0;
				if (jn == Nyy) jn = Nyy -1;
				ss[i][j] = (1 - w) * ss[i][j] + w / Ap[i][j] * (rhs[i][j] -Ae[i][j] * ss[ie][je] - Aw[i][j] * ss[iw][jw]
					-An[i][j] * ss[in][jn] - As[i][j] * ss[is][js]);
			}
		}
		for(int i = 0; i < Nxx; i++)
			for(int j = 0; j < Nyy; j++)
				sum += (ss[i][j] - ss0[i][j]) * (ss[i][j] - ss0[i][j]);
		res = sqrt(sum);
	}	
	return ss;
}