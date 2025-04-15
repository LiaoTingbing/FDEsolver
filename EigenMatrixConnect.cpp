

#pragma once
#include "EigenMatrixConnect.h"

SparseMatrix<double> horizontalConcatenation(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    // 1. ��������Ƿ�һ��
    if (A.rows() != B.rows()) {
        std::cerr << "Error: Matrices have incompatible row sizes for horizontal concatenation." << std::endl;
        return SparseMatrix<double>(0, 0); // ����һ���յ�ϡ�����
    }
    // 2. ���������� C �Ĵ�С
    Eigen::Index rows = A.rows();
    Eigen::Index cols = A.cols() + B.cols();
    // 3. ����������� C
    SparseMatrix<double> C(rows, cols);
    // 4. Ԥ�ȷ��� C �ķ���Ԫ�ؿռ� (��ѡ����ǿ���Ƽ�)
    //    ����ȷ�Ĺ��ƿ��������������
    //    ������һ�ָ���ȷ�Ĺ��Ʒ�����
    C.reserve(A.nonZeros() + B.nonZeros());
    // 5. ƴ�Ӿ��� A �� B �� C
    std::vector<Triplet<double>> triplets; // ���ڹ���ϡ������ triplets
    triplets.reserve(A.nonZeros() + B.nonZeros()); // Ԥ���� triplets �ռ�
    // ���� A ��Ԫ��
    for (Eigen::Index k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    // ���� B ��Ԫ�أ���������ƫ��
    for (Eigen::Index k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col() + A.cols(), it.value()));
        }
    }
    // 6. ʹ�� triplets ���� C
    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}
SparseMatrix<double> verticalConcatenation(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    // 1. ��������Ƿ�һ��
    if (A.cols() != B.cols()) {
        std::cerr << "Error: Matrices have incompatible column sizes for vertical concatenation." << std::endl;
        return SparseMatrix<double>(0, 0); // ����һ���յ�ϡ�����
    }
    // 2. ���������� C �Ĵ�С
    Eigen::Index rows = A.rows() + B.rows();
    Eigen::Index cols = A.cols();
    // 3. ����������� C
    SparseMatrix<double> C(rows, cols);
    // 4. Ԥ�ȷ��� C �ķ���Ԫ�ؿռ� (��ѡ����ǿ���Ƽ�)
    C.reserve(A.nonZeros() + B.nonZeros());
    // 5. ƴ�Ӿ��� A �� B �� C (ʹ�� Triplet)
    std::vector<Triplet<double>> triplets;
    triplets.reserve(A.nonZeros() + B.nonZeros());
    // ���� A ��Ԫ��
    for (Eigen::Index k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    // ���� B ��Ԫ�أ���������ƫ��
    for (Eigen::Index k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row() + A.rows(), it.col(), it.value()));
        }
    }
    // 6. ʹ�� triplets ���� C
    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}