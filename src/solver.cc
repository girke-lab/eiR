#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
using namespace std;
#include <f2c.h>
#include "solver.h"

extern "C" {
	typedef long int integer;
	typedef double doublereal;
	typedef long int ftnlen;


	int s_copy(char *, char *, ftnlen, ftnlen);
  int setulb_(integer *n, integer *m, doublereal *x, 
	doublereal *l, doublereal *u, integer *nbd, doublereal *f, doublereal 
	*g, doublereal *factr, doublereal *pgtol, doublereal *wa, integer *
	iwa, char *task, integer *iprint, char *csave, logical *lsave, 
	integer *isave, doublereal *dsave, ftnlen task_len, ftnlen csave_len);
}

Solver::Solver(int dim, int size, int param_m, double *coords)
{
	k = dim;
	m = size;
	nmax = k + 1;					// nmax is the max dimensionality
	mmax = param_m + 1;		// mmax is the max value on an internal parameter m
	p = new doublereal[m * k];
	for (int i = 0; i < m*k; i ++)
		p[i] = coords[i];
	// prepare workspace buffers
	nbd = new integer[nmax];
	iwa = new integer[3*nmax];
	l = new doublereal[nmax];
	u = new doublereal[nmax];
	g = new doublereal[nmax];
	wa = new doublereal[2*mmax*nmax + 4*nmax + 12*mmax*mmax + 12*mmax];
}

Solver::~Solver()
{
	delete p, nbd, iwa, l, u, g, wa;
}

void Solver::init_values(doublereal *x, doublereal *l, doublereal *u,
	integer *nbd)
{
	for (int i = 0; i < k; i ++) {
		x[i] = 0;
		l[i] = -.5; u[i] = .5; nbd[i] = 2;
	}
}
//euclidean distance
doublereal Solver::udist(doublereal* x, doublereal* y)
{
	doublereal dist = 0;
	doublereal t;
	for (int i = 0; i < k; i ++) {
		t = (*(x+i) - *(y+i));
		dist += (t * t);
	}
	return sqrt(dist);
}

doublereal Solver::calc_f(doublereal *x, doublereal *d)
{
	doublereal e = 0;
	doublereal dd = 0;
	doublereal t = 0;
	for (int i = 0; i < m; i ++) {
		dd = udist(x, p + k * i);
		t = dd - d[i];
		e += (t * t);
	}
	return e;
}

void Solver::calc_g(doublereal *x, doublereal *d, doublereal *g)
{
	for (int i = 0; i < k; i ++) {
		doublereal s = 0, tt;
		for (int j = 0; j < m; j ++) {
			doublereal dd = udist(x, p + k * j);
			tt = (x[i] - p[k * j + i]) * (1 - d[j] / dd);
			s += tt;
		}
		g[i] = 2 * s;
	}
}

int Solver::optim(doublereal *x, doublereal *d, doublereal factr,
	doublereal pgtol)
{
	char task[60] = "";
	char buf[60] = "";
	char csave[60] = "";
   char start[]="START";
	ftnlen task_len = 60;;
	ftnlen csave_len = 60;
	logical lsave[4];
	integer n, _m, iprint, isave[44];

	_m = mmax - 1;
	doublereal f, dsave[29];

	iprint = -1;

	n = k;

	init_values(x, l, u, nbd);
 
	s_copy(task, start, (ftnlen)60, (ftnlen)5);
	while(1) {
		setulb_(&n, &_m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task,
				&iprint, csave, lsave, isave, dsave, task_len, csave_len);
 
		if (!strncmp(task, "FG", 2)) {
			f = calc_f(x, d);
			calc_g(x, d, g);
		}
		else if (!strncmp(task, "NEW_X", 5)) continue;
    else return 0;
	}
}

// vim:tabstop=2:shiftwidth=2:smartindent

