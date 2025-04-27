#include "fde_solve.h"

FdeSolve::FdeSolve() {

}

FdeSolve::FdeSolve(map<string, cube>* dev)
{
	dev_ = dev;
}

FdeSolve::~FdeSolve()
{
}

void FdeSolve::initialize()
{
	lambda_ = (*dev_)["lambda"](0);
	cout << "\tlambda\t" <<lambda_;
	dx_ = (*dev_)["x"](1) - (*dev_)["x"](0);
	dy_ = (*dev_)["y"](1) - (*dev_)["y"](0);
	nx_ = (*dev_)["x"].n_elem;
	ny_ = (*dev_)["y"].n_elem;
	nt_ = nx_ * ny_;
	k0_ = 2 * PI / lambda_;

}

void FdeSolve::calculateIsotropicPMatrix( )
{
	vec eps = vectorise( (*dev_)["indexZ"] ) ;
	vec full1ColumnVector(nt_, fill::ones);

	Pxx_ = dxdxFuncReal(full1ColumnVector,
		eps, eps, nx_, ny_, dx_, dy_)
		+ dydyFuncReal(full1ColumnVector, full1ColumnVector,
			full1ColumnVector, nx_, ny_, dx_, dy_);
	Pxx_.diag(0) += eps * k0_ * k0_;

	Pyy_ = dxdxFuncReal(full1ColumnVector, full1ColumnVector, full1ColumnVector,
		nx_, ny_, dx_, dy_) + dydyFuncReal(full1ColumnVector, eps, eps,
			nx_, ny_, dx_, dy_);
	Pyy_.diag(0) += eps * k0_ * k0_;

	Pxy_ = dxdyFuncReal(full1ColumnVector, eps, eps, nx_, ny_, dx_, dy_) -
		dxdyFuncReal(full1ColumnVector, full1ColumnVector, full1ColumnVector,
			nx_, ny_, dx_, dy_);

	Pyx_ = dydxFuncReal(full1ColumnVector, eps, eps, nx_, ny_, dx_, dy_) -
		dydxFuncReal(full1ColumnVector, full1ColumnVector, full1ColumnVector,
			nx_, ny_, dx_, dy_);
}

void FdeSolve::calculateCharacteristicValues()
{
	cout << "\n\t计算特征值\n";
	sp_mat pMatrix(2 * nt_, 2 * nt_);
	pMatrix(0, 0, size(nt_, nt_)) = Pxx_;
	pMatrix(0, nt_, size(nt_, nt_)) = Pxy_;
	pMatrix(nt_, 0, size(nt_, nt_)) = Pyx_;
	pMatrix(nt_, nt_, size(nt_, nt_)) = Pyy_;


	cx_vec betaValue;
	cx_mat fieldValue;
 
	double sigma = 3.4*3.4*k0_*k0_;
	arma::eigs_gen(betaValue, fieldValue,   (pMatrix), 20, sigma);
	betaValue = sqrt(betaValue) / k0_;
	betaValue.print();

	fieldValue.save("text.txt", csv_ascii);

}
