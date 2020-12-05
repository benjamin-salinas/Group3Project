#ifndef _Solver_
#define _Solver_
#include <math.h>
class Solver
{
public:
	//double tol;

public:
	double* GaussElimination(double* diag, double* upper, double* lower, double* rhs, int N);
	double** SOR(double** Ae, double** Aw, double** An, double** As, double** Ap, double** rhs, int Nxx, int Nyy );
	double* CG();
	double* multigrid();
};

#endif