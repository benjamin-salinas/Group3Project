#ifndef _Solver_
#define _Solver_

class Solver
{
public:

public:
	double* GaussElimination(double* diag, double* upper, double* lower, double* rhs);
	double* SOR();
	double* CG();
	double* multigrid();
};

#endif