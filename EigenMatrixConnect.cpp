

#pragma once
#include "EigenMatrixConnect.h"

SparseMatrix<double> horizontalConcatenation(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    // 1. 检查行数是否一致
    if (A.rows() != B.rows()) {
        std::cerr << "Error: Matrices have incompatible row sizes for horizontal concatenation." << std::endl;
        return SparseMatrix<double>(0, 0); // 返回一个空的稀疏矩阵
    }
    // 2. 计算结果矩阵 C 的大小
    Eigen::Index rows = A.rows();
    Eigen::Index cols = A.cols() + B.cols();
    // 3. 创建结果矩阵 C
    SparseMatrix<double> C(rows, cols);
    // 4. 预先分配 C 的非零元素空间 (可选，但强烈推荐)
    //    更精确的估计可以显著提高性能
    //    以下是一种更精确的估计方法：
    C.reserve(A.nonZeros() + B.nonZeros());
    // 5. 拼接矩阵 A 和 B 到 C
    std::vector<Triplet<double>> triplets; // 用于构建稀疏矩阵的 triplets
    triplets.reserve(A.nonZeros() + B.nonZeros()); // 预分配 triplets 空间
    // 复制 A 的元素
    for (Eigen::Index k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    // 复制 B 的元素，并进行列偏移
    for (Eigen::Index k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col() + A.cols(), it.value()));
        }
    }
    // 6. 使用 triplets 构造 C
    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}
SparseMatrix<double> verticalConcatenation(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    // 1. 检查列数是否一致
    if (A.cols() != B.cols()) {
        std::cerr << "Error: Matrices have incompatible column sizes for vertical concatenation." << std::endl;
        return SparseMatrix<double>(0, 0); // 返回一个空的稀疏矩阵
    }
    // 2. 计算结果矩阵 C 的大小
    Eigen::Index rows = A.rows() + B.rows();
    Eigen::Index cols = A.cols();
    // 3. 创建结果矩阵 C
    SparseMatrix<double> C(rows, cols);
    // 4. 预先分配 C 的非零元素空间 (可选，但强烈推荐)
    C.reserve(A.nonZeros() + B.nonZeros());
    // 5. 拼接矩阵 A 和 B 到 C (使用 Triplet)
    std::vector<Triplet<double>> triplets;
    triplets.reserve(A.nonZeros() + B.nonZeros());
    // 复制 A 的元素
    for (Eigen::Index k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    // 复制 B 的元素，并进行行偏移
    for (Eigen::Index k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row() + A.rows(), it.col(), it.value()));
        }
    }
    // 6. 使用 triplets 构造 C
    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}