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
	cout << "\tlambda\t" << lambda_;
	x_ = vectorise((*dev_)["x"]);
	y_ = vectorise((*dev_)["y"]);

	dx_ = (*dev_)["x"](1) - (*dev_)["x"](0);
	dy_ = (*dev_)["y"](1) - (*dev_)["y"](0);
	nx_ = (*dev_)["x"].n_elem;
	ny_ = (*dev_)["y"].n_elem;
	nt_ = nx_ * ny_;
	k0_ = 2 * PI / lambda_;

}

void FdeSolve::calculateIsotropicPMatrix()
{
	cx_vec epsZ = vectorise((*dev_)["indexZ"]) + 0.0 * IU;
	//cx_vec epsX = vectorise((*dev_)["indexX"]) + 0.0 * IU;
	//cx_vec epsY = vectorise((*dev_)["indexY"]) + 0.0 * IU;
	cx_vec epsX = epsZ;
	cx_vec epsY = epsZ;
	cx_vec epsXY = vectorise((*dev_)["indexXY"]) + 0.0 * IU;
	cx_vec epsYX = vectorise((*dev_)["indexYX"]) + 0.0 * IU;




	cx_vec full1ColumnVector(nt_, fill::ones);

	Pxx_ = dxdxFunc(full1ColumnVector,
		epsZ, epsX, nx_, ny_, dx_, dy_)
		+ dydyFunc(full1ColumnVector, full1ColumnVector,
			full1ColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsX, ivec{ 0 }, nt_, nt_)
		+ dxdyFunc(full1ColumnVector, epsZ, epsYX, nx_, ny_, dx_, dy_);

	Pxy_ = dxdyFunc(full1ColumnVector, epsZ, epsY, nx_, ny_, dx_, dy_)
		- dxdyFunc(full1ColumnVector, full1ColumnVector, full1ColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsXY, ivec{ 0 }, nt_, nt_)
		+ dxdxFunc(full1ColumnVector, epsZ, epsXY, nx_, ny_, dx_, dy_);

	Pyx_ = dydxFunc(full1ColumnVector, epsZ, epsX, nx_, ny_, dx_, dy_)
		- dydxFunc(full1ColumnVector, full1ColumnVector, full1ColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsYX, ivec{ 0 }, nt_, nt_)
		+ dydyFunc(full1ColumnVector, epsZ, epsYX, nx_, ny_, dx_, dy_);

	Pyy_ = dydyFunc(full1ColumnVector, epsZ, epsY, nx_, ny_, dx_, dy_)
		+ dxdxFunc(full1ColumnVector, full1ColumnVector, full1ColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsY, ivec{ 0 }, nt_, nt_)
		+ dydxFunc(full1ColumnVector, epsZ, epsXY, nx_, ny_, dx_, dy_);
}

void FdeSolve::calculateCharacteristicValues()
{
	cout << "\n\t计算特征值\n";
	sp_cx_mat pMatrix(2 * nt_, 2 * nt_);
	pMatrix(0, 0, size(nt_, nt_)) = Pxx_;
	pMatrix(0, nt_, size(nt_, nt_)) = Pxy_;
	pMatrix(nt_, 0, size(nt_, nt_)) = Pyx_;
	pMatrix(nt_, nt_, size(nt_, nt_)) = Pyy_;


	cx_vec betaValue;
	cx_mat fieldValue;
 
	cx_double sigma = 3.4*3.4*k0_*k0_ + 0.0*IU;
	arma::eigs_gen(betaValue, fieldValue, real (pMatrix),20, sigma);
	betaValue = sqrt(betaValue) / k0_;
	betaValue.print();

	string outFilePath = "matlab/out.h5";

	mat fieldAbs = abs(fieldValue );
	fieldAbs.save(hdf5_name(outFilePath,"/fieldAbs"));
	x_.save(hdf5_name(outFilePath, "/x",hdf5_opts::append));
	y_.save(hdf5_name(outFilePath, "/y", hdf5_opts::append));




}
