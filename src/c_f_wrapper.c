#include <stdio.h>

void solve(int nprocs, int m, int n, int nnz, double *a, int *asub, int *xa, double *rhs);

/* Wrapper code to call solve() in dlinsol.dll from Fortran
  m    = nrow
  n    = ncol
  nnz  = nonz
  a    = nzval
  asub = rowind	(0-based)
  xa   = colptr	(0-based)
*/

void slu_solve(int *nprocs, int *n, int *nnz, double a[], int asub[], int xa[], double rhs[])
{
//	printf("call solve: nprocs,n,nnz: %d %d %d\n",*nprocs,*n,*nnz);
	solve(*nprocs, *n, *n, *nnz, a, asub, xa, rhs);
}
