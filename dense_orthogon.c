#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "dense_orthogon.h"

void gs(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R){
	int i,j;
	int n = A->I, m = A->J;
	double *qhat = (double*)calloc(n,sizeof(double));
	double *tmp_sum = (double*)calloc(n,sizeof(double));

	R->col_ptr[0][0] = norm_2(A->col_ptr[0], n);
	if(R->col_ptr[0][0] < ERR)
		return;
	vec_scal_prod(Q->col_ptr[0], A->col_ptr[0], R->col_ptr[0][0], n, 1);

	for(j=1;j<m;j++){
		for(i=0;i<n;i++)
			qhat[i] = A->col_ptr[j][i];
		for(i=0;i<j;i++){
			R->col_ptr[j][i] = inner_prod(qhat, Q->col_ptr[i], n);
			vec_scal_prod(tmp_sum, Q->col_ptr[i], R->col_ptr[j][i], n, 0);
			vec_vec_add(qhat, qhat, tmp_sum, n, 1);
		}
		R->col_ptr[j][j] = norm_2(qhat,n);
		if(R->col_ptr[j][j] < ERR)
			return;
		vec_scal_prod(Q->col_ptr[j], qhat, R->col_ptr[j][j], n, 1);
	}

	free(tmp_sum);
	free(qhat);
}

void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R){
	int i,j;
	int n = A->I, m = A->J;
    double *w_vals = (double*)calloc(n*m,sizeof(double));
	double **w = (double**)malloc(n*sizeof(double*));
	double *z = (double*)calloc(n,sizeof(double));
    double *wv = (double*)calloc(n,sizeof(double));
    double v;

    for(j=0;j<m;j++){
        w[j] = &w_vals[j*n];
    }

	for(j=0;j<m;j++){
		for(i=0;i<n;i++)
            R->col_ptr[j][i] = A->col_ptr[j][i];
        if(j > 0){
            for(i=0;i<j;i++){
                v = 2.0*inner_prod(R->col_ptr[j], w[i], n);
                vec_scal_prod(wv, w[i], v, n, 0); 
                vec_vec_add(R->col_ptr[j], R->col_ptr[j], wv, n, 1);
            }
		}
		get_z(z, R->col_ptr[j], j, n);
		vec_scal_prod(w[j], z, norm_2(z, n), n, 1);
        v = 2.0*inner_prod(R->col_ptr[j], w[j], n);
        vec_scal_prod(wv, w[j], v, n, 0);
        vec_vec_add(R->col_ptr[j], R->col_ptr[j], wv, n, 1);
        Q->col_ptr[j][j] = 1.0;
        for(i=j;i>=0;i--){
            v = 2.0*inner_prod(Q->col_ptr[j], w[i], n);
            vec_scal_prod(wv, w[i], v, n, 0); 
            vec_vec_add(Q->col_ptr[j], Q->col_ptr[j], wv, n, 1);
        }
    }
    free(w_vals);
    free(w);
    free(z);
    free(wv);
}

void get_z(double *z, double *x,  int k, int n){
	int i,j;
	double beta, tmp = 0.0;

	for(i=0;i<n;i++){
		if(i<k){
			z[i] = 0;
		} else if(i==k){
			for(j=k;j<n;j++){
				tmp += x[j]*x[j];
			}
            beta = x[k] > 0.0 ? sqrt(tmp) : -sqrt(tmp);
			z[i] = beta + x[i];
		} else {
			z[i] = x[i];
		}
	}
}

void vec_scal_prod(double *xhat, double *x, double y, int n, int div){
	int i;
	for(i=0;i<n;i++)
		if(div)
			xhat[i] = x[i]/y;
		else
			xhat[i] = x[i]*y;
}

void vec_vec_add(double *xhat, double *x, double *y, int n, int sub){
	int i;
	for(i=0;i<n;i++){
		if(sub)
			xhat[i] = x[i] - y[i];
		else 
			xhat[i] = x[i] + y[i];
	}
}

double norm_2(double *x, int n){
	double norm=0;
	int i;
	for(i=0;i<n;i++)
		norm += x[i]*x[i];
	norm = sqrt(norm);

	return norm;
}

double inner_prod(double *x, double *y, int n){
	double prod;
	int i;
	for(i=0;i<n;i++)
		prod += x[i]*y[i];
	
	return prod;
}

double fwd_err(DenseMatrix *A){
	int i,j;
	double fwd_err = 0.0;
	for(j=0;j<A->J;j++)
		for(i=0;i<A->I;i++)
			fwd_err = A->col_ptr[j][i] > fwd_err ? A->col_ptr[j][i] : fwd_err;
	return fwd_err;
}

DenseMatrix *mat_mul(DenseMatrix *A, DenseMatrix *B){
	int i,j,k;
	DenseMatrix *AB;
	if(A->J != B->I)
		return NULL;

	AB = gen_mat(A->I,B->J,0);

	for(j=0;j<AB->J;j++)
		for(i=0;i<AB->I;i++)
			for(k=0;k<A->J;k++)
				AB->col_ptr[j][i] += A->col_ptr[k][i]*B->col_ptr[j][k];

	return AB;
}

/*
DenseMatrix *mat_mat_diff(DenseMatrix *A, DenseMatrix *B){

}
*/
