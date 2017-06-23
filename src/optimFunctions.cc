#include <R.h>
#include <Rinternals.h>

extern "C" {
	SEXP calcF(SEXP xS,SEXP dS,SEXP pS);
	SEXP calcG(SEXP xS,SEXP dS,SEXP pS);
}

typedef long int integer;
typedef double doublereal;
typedef long int ftnlen;
typedef long int logical;


double udist(int k,doublereal* x, doublereal* y)
{
//	printf("----- udist -----------------\n ");
//	printf("x: ");
//	for(int i=0;i<k;i++) printf("%f,",*(x+i));
//	printf("\n");
//
//	printf("y: ");
//	for(int i=0;i<k;i++) printf("%f,",*(y+i));
//	printf("\n");
//	printf("-----------------------------\n ");


	doublereal dist = 0;
	doublereal t;
	for (int i = 0; i < k; i ++) {
		t = (*(x+i) - *(y+i));
		dist += (t * t);
	}
	return sqrt(dist);
}

SEXP calcF(SEXP xS,SEXP dS,SEXP pS)
{

	int k = length(xS); //length of result vector
	int m = length(dS); //length of distance vector
	double *x=REAL(xS);
	double *d=REAL(dS);
	double *p=REAL(pS);

//	printf("x: ");
//	for(int i=0;i<m;i++) printf("%f,",x[i]);
//	printf("\n");
	//printf("p: "); //see what major order is
	//for(int i=0;i<length(pS);i++) printf("%f,",p[i]);
	
	/*
	 * d=2, 
	 * N1=4, N2=10
	 *
	 * k=1: n1
	 * k=2: N1*n2
	 * 4*n2+n1
	 */


//	for(int i=0;i<4;i++){
//		for(int j=0;j<10;j++)
//			printf("%f,",p[10*i + j]); //row-major
//			//printf("%f,",p[4*j + i]); //col-major
//		printf("\n");
//	}
//	printf("\n");





	double e = 0;
	double t = 0;
	//double *row;
	for (int i = 0; i < m; i ++) {

	//	row = p+k*i;
	//	printf("coords[%d]: ",i);
	//	for(int j=0;j<k;j++) printf("%f,",*(row+j));
	//	printf("\n");


		t = udist(k,x, p + k * i) - d[i];
		e += (t * t);
	}

	//printf("f(x): %f\n",e);

	SEXP result = PROTECT(allocVector(REALSXP, 1));
	REAL(result)[0] = e;
	UNPROTECT(1);

	return result;
}
SEXP calcG(SEXP xS,SEXP dS,SEXP pS)
{
	int k = length(xS); //length of result vector
	int m = length(dS); //length distance vector

	SEXP result = PROTECT(allocVector(REALSXP, k));
	double *x=REAL(xS);
	double *d=REAL(dS);
	double *p=REAL(pS);
	double *g=REAL(result);

	for (int i = 0; i < k; i ++) {
		doublereal s = 0, tt;
		for (int j = 0; j < m; j ++) {
			doublereal dd = udist(k,x, p + k * j);
			tt = (x[i] - p[k * j + i])   *    (1 - d[j] / dd);
			s += tt;
		}
		g[i] = 2 * s;
	}
	//calc_g(REAL(x),REAL(d),REAL(result));

	UNPROTECT(1);
	return result;
}

