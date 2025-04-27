

#define ARMA_USE_SUPERLU
#define ARMA_USE_ARPACK
#include <armadillo>

using namespace arma;
int main() {


	sp_cx_mat A =spdiags(cx_vec(100 , fill::value(cx_double(1,0.11))),
		ivec{ 0 } , 100 , 100 );

	cx_vec eigval;
	cx_mat eigvec;

	eigs_gen(eigval, eigvec, A, 5);  // find 5 eigenvalues/eigenvectors

	eigs_opts opts;
	opts.maxiter = 10000;            // increase max iterations to 10000

	eigs_gen(eigval, eigvec, A, 5, "lm", opts);

	eigval.print();
	return 0;
}