#pragma once


#include "common.h"

using namespace Eigen;

SparseMatrix<double> horizontalConcatenation(const SparseMatrix<double>& A, const SparseMatrix<double>& B);

SparseMatrix<double> verticalConcatenation(const SparseMatrix<double>& A, const SparseMatrix<double>& B);