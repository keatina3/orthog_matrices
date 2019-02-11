#include <stdio.h>
#include "utils.h"
#include "dense_orthogon.h"

int main(int argc, char *argv[]){
	int row = 3, vecs = 3;
	DenseMatrix *A, *Q, *R, *QR;

	A = gen_mat(row,vecs,1);
	Q = gen_mat(row,vecs,0);
	R = gen_mat(vecs,vecs,0);
	
	// test matrix //
	A->col_ptr[0][0] = 12.0;
	A->col_ptr[0][1] = 6.0;
	A->col_ptr[0][2] = -4.0;
	A->col_ptr[1][0] = -51.0;
	A->col_ptr[1][1] = 167.0;
	A->col_ptr[1][2] = 24.0;
	A->col_ptr[2][0] = 4.0;
	A->col_ptr[2][1] = -68.0;
	A->col_ptr[2][2] = -41.0;
	///////////////
	
	printf("Matrix A:\n");
	print_mat(A);
	
	//gs(A,Q,R);
	//printf("Matrix Q:\n");
	//print_mat(Q);
	//printf("Matrix R:\n");
	//print_mat(R);

    hhorth(A,Q,R);
	printf("Matrix Q:\n");
	print_mat(Q);
	printf("Matrix R:\n");
	print_mat(R);
    
    QR = mat_mul(Q,R);
	printf("Matrix QR:\n");
	print_mat(QR);
    
	// mat_diff //
	free_mat(A);
	free_mat(Q);
	free_mat(R);
	free_mat(QR);

	return 0;
}
