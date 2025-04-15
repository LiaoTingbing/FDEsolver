#pragma once

#include "common.h"

//(eigval, eigvec, P4, searchNum, sigma);

void sparseEigs(VectorXcd& eigval, MatrixXcd& eigvec,
	SparseMatrix<double>& M, int  searchNum, double sigma);