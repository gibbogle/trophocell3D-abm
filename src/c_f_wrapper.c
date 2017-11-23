#include <stdio.h>
//#include "slu_mt_ddefs.h"

void solve(int nprocs, int m, int n, int nnz, double *a, int *asub, int *xa, double *rhs);

/* Wrapper code to call solve() in dlinsol.dll from Fortran
  m    = nrow
  n    = ncol
  nnz  = nonz
  a    = nzval
  asub = rowind	(0-based)
  xa   = colptr	(0-based)
*/
void slu_solve_(int *nprocs, int *n, int *nnz, double *a, int *asub, int *xa, double *rhs)
{
	int np, mm, nn, nz;
//	printf("call solve: nprocs,n,nnz: %d %d %d\n",*nprocs,*n,*nnz);
	np = *nprocs;
	mm = *n;
	nn = *n;
	nz = *nnz;
	solve(*nprocs, *n, *n, *nnz, a, asub, xa, rhs);
}
