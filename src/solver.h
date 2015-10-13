#include "config.h"
#include <f2c.h>

using namespace std;

class Solver
{
	private:
		int k;			// dimensionality
		int m;			// trait size
		int nmax, mmax;
		doublereal *p;	// trait
		// workspace buffers
		integer *nbd, *iwa; 
		doublereal *l, *u, *g, *wa, *x, *d; 
		// internal methods
		void init_values(doublereal *x, doublereal *l, doublereal *u,
				integer *nbd);
		doublereal udist(doublereal *x, doublereal *y);
		doublereal calc_f(doublereal *x, doublereal *d);
		void calc_g(doublereal *x, doublereal *d, doublereal *g);
	public:
		Solver(int dim, int size, int param_m, double *coords);
		~Solver();
		int optim(doublereal *x, doublereal *d,
			doublereal factr=FACTR, doublereal pgtol=PGTOL);
};
