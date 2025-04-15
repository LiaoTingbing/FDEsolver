#pragma once


#include "sparseEigs.h"

void sparseEigs(VectorXcd& eigval, MatrixXcd& eigvec,
	SparseMatrix<double>& M, int  searchNum, double sigma) {

	SparseGenRealShiftSolve<double> op(M); // Matrix-vector product operation for A
	op.set_shift(sigma);

	GenEigsRealShiftSolver<SparseGenRealShiftSolve<double>> eigs(
		op, searchNum, 2 * searchNum + 1, sigma);

	eigs.init();
	int nconv = eigs.compute(SortRule::LargestMagn, 1000,
		1e-10, SortRule::LargestMagn);

	if (eigs.info() == CompInfo::Successful)
	{
		eigval = eigs.eigenvalues();
		eigvec = eigs.eigenvectors();
	}
}