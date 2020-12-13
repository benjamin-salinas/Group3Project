#ifndef _Solver_
#define _Solver_
#include <math.h>
#include <iostream>
using namespace std;
class Solver
{
public:
	//double tol;
	double** ss;
	int ix, iy;
public:
	Solver(int Nx);
	Solver(int Nx, int Ny);
	~Solver();
	double* GaussElimination(double* diag, double* upper, double* lower, double* rhs, int N);
	double** SOR(double** Ae, double** Aw, double** An, double** As, double** Ap, double** rhs, int Nxx, int Nyy );
	double* CG();
	double* multigrid();
};

#endif