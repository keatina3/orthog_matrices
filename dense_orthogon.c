#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "dense_orthogon.h"

void gs(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R){
	int i,j;
	int n = A->I, m = A->J;
	double *qhat = (double*)calloc(n,sizeof(double));
	double *tmp_sum = (double*)calloc(n,sizeof(double));	// to store sum(Rij*Qi)

	R->col_ptr[0][0] = norm_2(A->col_ptr[0], n);			// getting initial R11 = ||x1||
	vec_scal_prod(Q->col_ptr[0], A->col_ptr[0], R->col_ptr[0][0], n, 1);	// Q1 = x1/R11

	for(j=1;j<m;j++){
		for(i=0;i<n;i++)
			qhat[i] = A->col_ptr[j][i];			// since qhat = Xj - tmpsum, storing Xj first
		for(i=0;i<j;i++){
			R->col_ptr[j][i] = inner_prod(qhat, Q->col_ptr[i], n);		// Rij = (Xj,Qi)
			vec_scal_prod(tmp_sum, Q->col_ptr[i], R->col_ptr[j][i], n, 0);	// see above
			vec_vec_add(qhat, qhat, tmp_sum, n, 1);		// qhat = Xj - sum(Rij*Qi)
		}
		R->col_ptr[j][j] = norm_2(qhat,n);			// normalising Rij
		vec_scal_prod(Q->col_ptr[j], qhat, R->col_ptr[j][j], n, 1);		// Qi = qhat/Rij
	}
	// free allocated memory //
	free(tmp_sum);
	free(qhat);
}

void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R){
	int i,j;
	int n = A->I, m = A->J;
	double *w_vals = (double*)calloc(n*m,sizeof(double));
	double **w = (double**)malloc(n*sizeof(double*));		// to store matrix of w's
	double *z = (double*)calloc(n,sizeof(double));			// stores z vector
	double *wv = (double*)calloc(n,sizeof(double));			// temp storage for Wj*V)
	double v;
	double *Xk = (double*)calloc(n,sizeof(double));		// temp placeholder of X
	
	// initialisign pointers for W //
	for(j=0;j<m;j++){
    	w[j] = &w_vals[j*n];
	}
	for(j=0;j<m;j++){
		for(i=0;i<n;i++)
			Xk[i] = A->col_ptr[j][i];			// setting Xj = Aj
        if(j > 0){
            for(i=0;i<j;i++){			// iterating through (Pj-1)(Pj-2)...(P1)(Xj)
                v = 2.0*inner_prod(Xk, w[i], n);	// getting v=2*X(^T)w
                vec_scal_prod(wv, w[i], v, n, 0); 	// getting wv(^T)
                vec_vec_add(Xk, Xk, wv, n, 1);		// PXk = Xk - wv(^T)
            }
		}
		get_z(z, Xk, j, n);					// getting z to be used in getting Wj
        vec_scal_prod(w[j], z, norm_2(z, n), n, 1);		// w = z / ||z||
        v = 2.0*inner_prod(Xk, w[j], n);		// see above
        vec_scal_prod(wv, w[j], v, n, 0);
        vec_vec_add(Xk, Xk, wv, n, 1);
		for(i=0;i<=j;i++)
			R->col_ptr[j][i] = Xk[i];			//	rk = Pk(Pk-1)...Xk
        Q->col_ptr[j][j] = 1.0;					//	Qej
        for(i=j;i>=0;i--){						// Qk = P1..(Pk)(ek)
            v = 2.0*inner_prod(Q->col_ptr[j], w[i], n);
            vec_scal_prod(wv, w[i], v, n, 0); 
            vec_vec_add(Q->col_ptr[j], Q->col_ptr[j], wv, n, 1);
        }
    }
	// freeing allocated memory //
	free(Xk);
    free(w_vals);
    free(w);
    free(z);
    free(wv);
}

// applies piecewise function of z //
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

// scalar product //
// division if div=TRUE //
void vec_scal_prod(double *xhat, double *x, double y, int n, int div){
	int i;
    if(fabs(y) < ERR)
		return;
	
	for(i=0;i<n;i++){
		if(div)
			xhat[i] = x[i]/y;
		else
			xhat[i] = x[i]*y;
	}
}

// vector-vector add //
// subtract if sub=TRUE //
void vec_vec_add(double *xhat, double *x, double *y, int n, int sub){
	int i;
	for(i=0;i<n;i++){
		if(sub)
			xhat[i] = x[i] - y[i];
		else 
			xhat[i] = x[i] + y[i];
	}
}

// matrix-matrix multiply //
void mat_mul(DenseMatrix *AB, DenseMatrix *A, DenseMatrix *B){
	int i,j,k;
	if(A->J != B->I)
		return;
	for(j=0;j<AB->J;j++)
		for(i=0;i<AB->I;i++)
			for(k=0;k<A->J;k++)
				AB->col_ptr[j][i] += A->col_ptr[k][i]*B->col_ptr[j][k];

}

// matrix-matrix addition //
void mat_mat_add(DenseMatrix *A_B, DenseMatrix *A, DenseMatrix* B, int sub){
    int j;
    if(A->J != B->J || A->I != B->I)
        return;
    for(j=0;j<A_B->J;j++)
        vec_vec_add(A_B->col_ptr[j], A->col_ptr[j], B->col_ptr[j], A_B->I, sub);
}

// matrix transpose //
void transpose(DenseMatrix *A, DenseMatrix *At){
	int i,j;
	
	for(j=0;j<At->J;j++)
		for(i=0;i<At->I;i++)
			At->col_ptr[j][i] = A->col_ptr[i][j];
}

// euclidian norm //
double norm_2(double *x, int n){
	double norm=0;
	int i;
	for(i=0;i<n;i++)
		norm += x[i]*x[i];
	norm = sqrt(norm);

	return norm;
}

// inner-product of 2 vectors //
double inner_prod(double *x, double *y, int n){
	double prod;
	int i;
	for(i=0;i<n;i++)
		prod += x[i]*y[i];
	
	return prod;
}

// calculates max( |Aij| )
double fwd_err(DenseMatrix *A){
	int i,j;
	double fwd_err = 0.0;
	for(j=0;j<A->J;j++)
		for(i=0;i<A->I;i++)
			fwd_err = A->col_ptr[j][i] > fwd_err ? A->col_ptr[j][i] : fwd_err;
	return fwd_err;
}
