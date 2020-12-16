#include <math.h>
#include "Flow.h"
#include "Solver.h"
#include <iostream>
#include <fstream>
#include<sstream>
#include<iomanip>

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
	for (int j = 0; j < Mesh::Ny; j++)
	{
		uu[0][j] = exp(-(mesh->yc[j] - 0.5) * (mesh->yc[j] - 0.5) / R / R); 
		//uu[Mesh::Nx][j] = exp(-(mesh->yc[j] - 0.5) * (mesh->yc[j] - 0.5) / R / R);
	}
		
	
	alpha = dt / (2. * Re * dx * dx);
	//alpha = 0.;

	ustar2 = new double* [Mesh::Ny];
	for (int j = 0; j < Mesh::Ny; j++)
	{
		ustar2[j] = new double[Mesh::Nx];
	}
	vstar2 = new double* [Mesh::Ny - 1];
	for (int j = 0; j < Mesh::Ny - 1; j++)
	{
		vstar2[j] = new double[Mesh::Nx];
	}

	if (restartflow == true) readRestart(); 
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
	if (j_ == -1) j_ = 0;
	if (j_ == Mesh::Ny) j_ = Mesh::Ny - 1;
	if (i_ == Mesh::Nx + 1) i_ = Mesh::Nx - 1;
	if (uu[i_][j_] < 0) uu[i_][j_] = 0.;
	if (i_ == 0) 
		return exp(-(mesh->yc[j_] - 0.5) * (mesh->yc[j_] - 0.5) / R / R);
	else
		return uu[i_][j_];
}

double Flow::v(int i_, int j_)
{
	if (i_ == -1) i_ = 0;
	if (i_ == Mesh::Nx) i_ = Mesh::Nx - 1;
	if (vv[i_][j_] < -0.001) vv[i_][j_] = -0.001; //Only apply for Re = inf
	if (vv[i_][j_] > 0.001) vv[i_][j_] = 0.001;
	if (j_ == 0 || j_ == Mesh::Ny) return 0.;
	else
		return vv[i_][j_];
}

double Flow::p(int i_, int j_)
{
	if (i_ == -1) i_ = 0;
	if (i_ == Mesh::Nx) i_ = Mesh::Nx - 1;
	if (j_ == -1) j_ = 0;
	if (j_ == Mesh::Ny) j_ = Mesh::Ny - 1;
	return pp[i_][j_];

}

void Flow::updateFlow()
{
	clock_t time1, time2, time3, time4;
	//Step 1 to get ustar and vstar
	time1 = clock();
	getustar();
	getvstar();
	count++;

	//Correct for mass conservation Ning
	double sum1 = 0., sum2 = 0., kk, mm;
	L2_D = 0.;
	for (int j = 0; j < Mesh::Ny; j++)
	{
		sum1 += ustar[0][j];
		sum2 += ustar[Mesh::Nx][j];
	}
	mm = (sum1 + 1e-100) / (sum2 + 1e-100);
	for (int j = 0; j < Mesh::Ny; j++) ustar[Mesh::Nx][j] = mm * ustar[Mesh::Nx][j];

	time2 = clock();
	step1 = (double)(time2 - time1) / CLOCKS_PER_SEC;
	//Step 2
	getpp();
	time3 = clock();
	step2 = (double)(time3 - time2) / CLOCKS_PER_SEC;

	//Step 3
	getuu();
	getvv();
	time4 = clock();
	step3 = (double)(time4 - time3) / CLOCKS_PER_SEC;

	for (int i = 0; i < Mesh::Nx; i++)
		for (int j = 0; j < Mesh::Ny; j++)
			L2_D = L2_D + ((uu[i + 1][j] - uu[i][j]) / dx + (vv[i][j + 1] - vv[i][j]) / dy) *
			((uu[i + 1][j] - uu[i][j]) / dx + (vv[i][j + 1] - vv[i][j]) / dy);
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
	getustar2();
	
	double* lower = new double[Mesh::Ny];
	double* diag = new double[Mesh::Ny];
	double* upper = new double[Mesh::Ny];
	double* rhs = new double[Mesh::Ny];

	for (int j = 0; j < Mesh::Ny; j++) 
	{
		lower[j] = -alpha;
		diag[j] = 1 - 2. * alpha;
		upper[j] = -alpha;
		if (j == 0) diag[j] = 1 - 3. * alpha;	// Add vector value at the boundary
		else if (j == Mesh::Ny - 1) diag[j] = 1 - 3. * alpha; 
	}

	for (int i = 1; i <= Mesh::Nx; i++) 
	{
		for (int j = 0; j < Mesh::Ny; j++) rhs[j] = ustar2[j][i - 1];
		Solver* solver2 = new Solver(Mesh::Ny);
		ustar[i] = solver2->GaussElimination(diag, upper, lower, rhs, Mesh::Ny);
	}
	// delete temporal arrays
	delete[] lower; lower = NULL;
	delete[] upper; upper = NULL;
	delete[] diag; diag = NULL;
	delete[] rhs; rhs = NULL;
	for (int j = 0; j < Mesh::Ny; j++) ustar[0][j] = 0.;

	//get ustar from delta_ustar
	for (int i = 0; i < Mesh::Nx + 1; i++)
		for (int j = 0; j < Mesh::Ny; j++) ustar[i][j] += u(i, j);

}

void Flow::getustar2()
{
	
	for (int j = 0; j < Mesh::Ny; j++)
	{
		double* rhsu2 = new double[Mesh::Nx];
		double* lower = new double[Mesh::Nx];
		double* diag = new double[Mesh::Nx];
		double* upper = new double[Mesh::Nx];
		for (int i = 0; i < Mesh::Nx; i++)
		{
			lower[i] = -alpha;
			diag[i] = 1 - 2. * alpha;
			upper[i] = -alpha;
			if (i == Mesh::Nx - 1) lower[i] = -2. * alpha;
		}

		for (int i = 1; i <= Mesh::Nx; i++)
		{
			double newH = getHx(i, j);
			rhsu2[i - 1] = dt / 2. * (3 * newH - Hx[i - 1][j])
				+ dt / Re * (u(i + 1, j) + u(i - 1, j) - 2. * u(i, j)) / dx / dx
				+ dt / Re * (u(i, j + 1) + u(i, j - 1) - 2. * u(i, j)) / dy / dy;
			Hx[i - 1][j] = newH;
		}
	
		Solver* solver1 = new Solver(Mesh::Nx);
	
		ustar2[j] = solver1->GaussElimination(diag, upper, lower, rhsu2, Mesh::Nx); //Ning
		delete[] rhsu2; rhsu2 = NULL;
		delete[] lower; lower = NULL;
		delete[] diag; diag = NULL;
		delete[] upper; upper = NULL;
	}
}

void Flow::getvstar()
{
	getvstar2();
	
	double* lowerv = new double[Mesh::Ny - 1];
	double* diagv = new double[Mesh::Ny - 1];
	double* upperv = new double[Mesh::Ny - 1];
	double* rhsv = new double[Mesh::Ny - 1];
	for (int j = 0; j < Mesh::Ny - 1; j++)
	{
			lowerv[j] = -alpha;
			diagv[j] = 1 - 2. * alpha;
			upperv[j] = -alpha;

	}

	for (int i = 0; i < Mesh::Nx; i++) 
	{
		for (int j = 0; j < Mesh::Ny - 1; j++) rhsv[j] = vstar2[j][i];	
		Solver* solver2 = new Solver(Mesh::Ny - 1);
		vstar[i] = solver2->GaussElimination(diagv, upperv, lowerv, rhsv, Mesh::Ny - 1);
	}

	
	for (int i = 0; i < Mesh::Nx; i++)
	{
		for (int j = Mesh::Ny - 1; j >= 1; j--)
			vstar[i][j] = vstar[i][j - 1];
		vstar[i][0] = 0.;
		vstar[i][Mesh::Ny] = 0.;
 	}
	//get vstar from delta_vstar
	for (int i = 0; i < Mesh::Nx; i++)
		for (int j = 0; j < Mesh::Ny + 1; j++) vstar[i][j] += v(i, j);
	
}

void Flow::getvstar2()  //j from 1 to Ny - 1
{
	for (int j = 0; j < Mesh::Ny - 1; j++)
	{
		int jj = j + 1;
		double* rhsv2 = new double[Mesh::Nx];
		double* lower = new double[Mesh::Nx];
		double* diag = new double[Mesh::Nx];
		double* upper = new double[Mesh::Nx];
		for (int i = 0; i < Mesh::Nx; i++)
		{
			lower[i] = -alpha;
			diag[i] = 1 - 2. * alpha;
			upper[i] = -alpha;
			if (i == Mesh::Nx - 1) diag[i] = 1 - 3. * alpha;
			else if (i == 0) diag[i] = 1 - alpha;  //Ning **
		}

		for (int i = 0; i < Mesh::Nx; i++)
		{
			double newH = getHy(i, jj);
			rhsv2[i] = dt / 2. * (3 * newH - Hy[i][jj - 1])
				+ dt / Re * (v(i + 1, jj) + v(i - 1, jj) - 2. * v(i, jj)) / dx / dx
				+ dt / Re * (v(i, jj + 1) + v(i, jj - 1) - 2. * v(i, jj)) / dy / dy;
			Hy[i][jj - 1] = newH;
		}

		Solver* solver3 = new Solver(Mesh::Nx);
		vstar2[j] = solver3->GaussElimination(diag, upper, lower, rhsv2, Mesh::Nx);
		delete[] rhsv2; rhsv2 = NULL;
		delete[] lower; lower = NULL;
		delete[] diag; diag = NULL;
		delete[] upper; upper = NULL;
	}
}


void Flow::getpp()
{

	Solver* solver5 = new Solver(Mesh::Nx, Mesh::Ny);
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
	for (int i = 0; i < Mesh::Nx; i++)
	{
		delete[] Ae[i];
		delete[] Aw[i];
		delete[] As[i];
		delete[] An[i];
		delete[] Ap[i];
		delete[] rhs[i];
	}
	delete[] Ae;
	delete[] Aw;
	delete[] As;
	delete[] An;
	delete[] Ap;
	delete[] rhs;
}

void Flow::getuu()
{
	for (int i = 1; i <= Mesh::Nx; i++)
		for (int j = 0; j < Mesh::Ny; j++) 
			uu[i][j] = ustar[i][j] - dt * (p(i, j) - p(i - 1, j)) / dx;
}

void Flow::getvv()
{
	for (int i = 0; i < Mesh::Nx; i++)
		for (int j = 1; j < Mesh::Ny; j++)
			vv[i][j] = vstar[i][j] - dt * (p(i, j) - p(i, j - 1)) / dy;
}

void Flow::writeRestart()
{
	ofstream restartuu;
	restartuu.open("./restart/restartuu.csv", ios::out);
	for (int j = 0; j < Mesh::Ny; j++)
		for (int i = 0; i <= Mesh::Nx; i++) restartuu << uu[i][j] << endl;
	restartuu.close();

	ofstream restartvv;
	restartvv.open("./restart/restartvv.csv", ios::out);
	for (int j = 0; j <= Mesh::Ny; j++)
		for (int i = 0; i < Mesh::Nx; i++) restartvv << vv[i][j] << endl;
	restartvv.close();

	ofstream restarthx;
	restarthx.open("./restart/restarthx.csv", ios::out);
	for (int j = 0; j < Mesh::Ny; j++)
		for (int i = 0; i < Mesh::Nx; i++) restarthx << Hx[i][j] << endl;
	restarthx.close();

	ofstream restarthy;
	restarthy.open("./restart/restarthy.csv", ios::out);
	for (int j = 0; j < Mesh::Ny - 1; j++)
		for (int i = 0; i < Mesh::Nx; i++) restarthy << Hy[i][j] << endl;
	restarthy.close();
}

void Flow::readRestart()
{
	//read velocity u
	ifstream file;
	file.open("./restart/restartuu.csv", ios::in);
	int n = 0;
	string temp;
	while (getline(file, temp)) n++;
	int LINES = n;
	cout << "n:" << LINES << endl;
	file.close();
	file.open("./restart/restartuu.csv", ios::in);
	for (int row = 0; row < LINES; row++) {

		string line;
		getline(file, line);
		stringstream iss(line);
		string val;

		int ii = row % (Mesh::Nx + 1);
		int jj = row / (Mesh::Nx + 1);

		getline(iss, val, ',');

		stringstream convertor(val);

		double inv;
		convertor >> inv;
		uu[ii][jj] = inv;

	}
	file.close();


	//read velocity v
	ifstream file2;
	file2.open("./restart/restartvv.csv", ios::in);
	n = 0;
	string temp2;
	while (getline(file2, temp2)) n++;
	LINES = n;
	file2.close();
	file2.open("./restart/restartvv.csv", ios::in);
	for (int row = 0; row < LINES; row++) {

		string line;
		getline(file2, line);
		stringstream iss(line);
		string val;
		int ii = row % Mesh::Nx;
		int jj = row / Mesh::Nx;

		getline(iss, val, ',');

		stringstream convertor(val);

		double inv;
		convertor >> inv;
		vv[ii][jj] = inv;
	}
	file2.close();

	//read Hx data
	ifstream filehx;
	filehx.open("./restart/restarthx.csv", ios::in);
	n = 0;
	string temp5;
	while (getline(filehx, temp5)) n++;
	LINES = n;
	filehx.close();
	filehx.open("./restart/restarthx.csv", ios::in);
	for (int row = 0; row < LINES; row++) {

		string line;
		getline(filehx, line);
		stringstream iss(line);
		string val;
		int ii = row % Mesh::Nx;
		int jj = row / Mesh::Nx;

		getline(iss, val, ',');

		stringstream convertor(val);

		double inv;
		convertor >> inv;
		Hx[ii][jj] = inv;
		//cout << "v: " << vv[ii][jj] << endl;
	}
	filehx.close();

	//read Hy data
	ifstream filehy;
	filehy.open("./restart/restarthy.csv", ios::in);
	n = 0;
	string temp6;
	while (getline(filehy, temp6)) n++;
	LINES = n;
	filehy.close();
	filehy.open("./restart/restarthy.csv", ios::in);
	for (int row = 0; row < LINES; row++) {

		string line;
		getline(filehy, line);
		stringstream iss(line);
		string val;
		int ii = row % Mesh::Nx;
		int jj = row / Mesh::Nx;

		getline(iss, val, ',');

		stringstream convertor(val);

		double inv;
		convertor >> inv;
		Hy[ii][jj] = inv;
		//cout << "v: " << vv[ii][jj] << endl;
	}
	filehy.close();
}