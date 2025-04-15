#include "armaToEigen.h"

Eigen::SparseMatrix<double> armaToEigenSparseManual(const arma::sp_mat& arma_sparse) {
    Eigen::SparseMatrix<double> eigen_sparse(arma_sparse.n_rows, arma_sparse.n_cols);

    // Ԥ����ռ䣨������ܣ�
    eigen_sparse.reserve(arma_sparse.n_nonzero);

    // �������з���Ԫ�ز�����
    int c = 0;
    for (auto it = arma_sparse.begin(); it != arma_sparse.end(); ++it) {
        eigen_sparse.insert(it.row(), it.col()) = *it;
        c++;
    }
    //cout << c;
    eigen_sparse.makeCompressed(); // ѹ������
    return eigen_sparse;
}

