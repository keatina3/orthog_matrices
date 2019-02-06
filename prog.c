#include <stdio.h>
#include "utils.h"
#include "dense_orthogon.h"

int main(int argc, char *argv[]){
	int row = 100, col = 90;
	DenseMatrix *A, *Q, *R;

	A = gen_mat(row,col,1);
	Q = gen_mat(row,col,0);
	R = gen_mat(row,col,0);

	printf("test1\n");
	print_mat(A);
	printf("test2\n");
	
	free_mat(A);
	free_mat(Q);
	free_mat(R);

	return 0;
}
