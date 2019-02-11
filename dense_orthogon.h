#ifndef _DENSE_ORTHOGON_H_
#define _DENSE_ORTHOGON_H_

#define ERR 1.0E-08

void gs(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R);
void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R);
void get_z(double *z, double *x, int k, int n);
void vec_scal_prod(double *x, double *xhat, double y, int n, int div);
void vec_vec_add(double *x, double *xhat, double *y, int n, int sub);
double norm_2(double *x, int n);
double inner_prod(double *x, double *y, int n);
double fwd_err(DenseMatrix* A);
DenseMatrix *mat_mul(DenseMatrix *A, DenseMatrix *B);
DenseMatrix *mat_mat_add(DenseMatrix *A, DenseMatrix *B, int sub);

#endif
