#include "armaToEigen.h"

Eigen::SparseMatrix<double> armaToEigenSparseManual(const arma::sp_mat& arma_sparse) {
    Eigen::SparseMatrix<double> eigen_sparse(arma_sparse.n_rows, arma_sparse.n_cols);

    // 预分配空间（提高性能）
    eigen_sparse.reserve(arma_sparse.n_nonzero);

    // 遍历所有非零元素并插入
    int c = 0;
    for (auto it = arma_sparse.begin(); it != arma_sparse.end(); ++it) {
        eigen_sparse.insert(it.row(), it.col()) = *it;
        c++;
    }
    //cout << c;
    eigen_sparse.makeCompressed(); // 压缩矩阵
    return eigen_sparse;
}

