#ifndef SCALAPACK_H
#define SCALAPACK_H

#include <complex>
typedef complex<float> complex_s;
typedef complex<double> complex_d;

extern "C" {
	void Cblacs_pinfo( int* mypnum, int* nprocs);
	void Cblacs_get(int context, int request, int* value);
	int Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
	void Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
	void Cblacs_gridexit( int context);
	void Cblacs_exit( int error_code);
	void Cblacs_gridmap( int* context, int* map, int ld_usermap, int np_row, int np_col);

	int npreroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
	int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

	void descinit_(int *idescal, int *m, int *n, int *mb, int *nb, int *dummy1 , int *dummy2 , int *icon, int *procRows, int *info);

	void psgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *ia, int *ja, int *desca, float *s, float *u, int *iu, int *ju, int *descu, float *vt, int *ivt, int *jvt, int *descvt, float *work, int *lwork, int *info);
	void pdgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *s, double *u, int *iu, int *ju, int *descu, double *vt, int *ivt, int *jvt, int *descvt, double *work, int *lwork, int *info);
	void pcgesvd_(char *jobu, char *jobvt, int *m, int *n, complex_s *a, int *ia, int *ja, int *desca, float *s, complex_s *u, int *iu, int *ju, int *descu, complex_s *vt, int *ivt, int *jvt, int *descvt, complex_s *work, int *lwork, float *rwork, int *info);
	void pzgesvd_(char *jobu, char *jobvt, int *m, int *n, complex_d *a, int *ia, int *ja, int *desca, double *s, complex_d *u, int *iu, int *ju, int *descu, complex_d *vt, int *ivt, int *jvt, int *descvt, complex_d *work, int *lwork, double *rwork, int *info);
	void pzheevd_(char *jobz, char *uplo, int *n, complex_d *a, int *ia, int *ja, int *desca, double *w, complex_d *z, int *iz, int *jz, int *descz, complex_d *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
	void pztranc_ (int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc,
int *descc);

/*
 * sub( C ) := alpha*op( sub( A ) )*op( sub( B ) ) + beta*sub( C )
 */
	void pdgemm_ (char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
	void pzgemm_ (	char *transa, char *transb, int *m, int *n, int *k, 
					double *alpha, 	
					complex_d *a, int *ia, int *ja, int *desca, 
					complex_d *b, int *ib, int *jb, int *descb, 
					double *beta, 
					complex_d *c, int *ic, int *jc, int *descc);

}

#endif
