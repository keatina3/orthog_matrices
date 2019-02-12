#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"
#include "dense_orthogon.h"

int main(int argc, char *argv[]){
	int row = 100, vecs = 90;
	int i,j,k;
	double fwd_err1, fwd_err2, fwd_err3, fwd_err4;
	DenseMatrix A, Q, R, QR, A_QR, Qt, QtQ, I, I_QtQ;
	srand48(time(NULL));

	A = gen_mat(row,vecs,1);
	Q = gen_mat(row,vecs,0);
	R = gen_mat(vecs,vecs,0);
	Qt = gen_mat(vecs,row,0);
	QR = gen_mat(row,vecs,0);
	A_QR = gen_mat(row,vecs,0);
	I = gen_mat(vecs,vecs,0);
	QtQ = gen_mat(vecs,vecs,0);
	I_QtQ = gen_mat(vecs,vecs,0);
	
	for(i=0;i<vecs;i++)
		I.col_ptr[i][i] = 1.0;
	printf("GS: A-QR\tGS: I-QtQ\tHH: A-QR\tHH: I-QtQ\n");
	for(i=0;i<10;i++){
		gs(&A, &Q, &R);
    	mat_mul(&QR, &Q, &R);
   		mat_mat_add(&A_QR, &A, &QR, 1);
		fwd_err1 = fwd_err(&A_QR);
		
		transpose(&Q, &Qt);
		mat_mul(&QtQ, &Qt, &Q);
		mat_mat_add(&I_QtQ, &I, &QtQ, 1);
		fwd_err2 = fwd_err(&I_QtQ);

		for(j=0;j<vecs;j++){
			for(k=0;k<vecs;k++){
				R.col_ptr[j][k] = 0.0;
				QtQ.col_ptr[j][k] = 0.0;
				I_QtQ.col_ptr[j][k] = 0.0;
			}
		}
		
		for(j=0;j<vecs;j++){
			for(k=0;k<row;k++){
				A_QR.col_ptr[j][k] = 0.0;
				Q.col_ptr[j][k] = 0.0;
				QR.col_ptr[j][k] = 0.0;
			}
		}
		
		hhorth(&A, &Q, &R);
    	mat_mul(&QR, &Q, &R);
   		mat_mat_add(&A_QR, &A, &QR, 1);
		fwd_err3 = fwd_err(&A_QR);
		
		transpose(&Q, &Qt);
		mat_mul(&QtQ, &Qt, &Q);
		mat_mat_add(&I_QtQ, &I, &QtQ, 1);
		fwd_err4 = fwd_err(&I_QtQ);
		
		if(row <= 10 && vecs <= 10){
			printf("Matrix A\n");
			print_mat(&A);
			printf("Matrix Q\n");
			print_mat(&Q);
			printf("Matrix R\n");
			print_mat(&R);
			printf("Matrix QR\n");
			print_mat(&QR);
			printf("Matrix I\n");
			print_mat(&I);
			printf("Matrix QtQ\n");
			print_mat(&QtQ);
		} else { 
		printf("%0.15lf\t%0.15lf\t%0.15lf\t%0.15lf\n",fwd_err1,fwd_err2,fwd_err3,fwd_err4);
		}
		assign_rand(&A);
    }

	free_mat(&A);
	free_mat(&Q);
	free_mat(&R);
	free_mat(&QR);
	free_mat(&A_QR);
	free_mat(&I_QtQ);
	free_mat(&QtQ);
	free_mat(&I);
	free_mat(&Qt);
	
	return 0;
}
