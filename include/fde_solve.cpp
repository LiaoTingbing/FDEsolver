#include "fde_solve.h"

FdeSolve::FdeSolve() {

}

 
FdeSolve::FdeSolve(map<string, cube>* dev, double lambda, int nmodes, double search)
{
	dev_ = dev;
	lambda_ = lambda;
	nmodes_ = nmodes;
	search_ = search;
	k0_ = 2 * PI / lambda_;

}

FdeSolve::~FdeSolve()
{
}

void FdeSolve::initialize()
{
	
	cout << "\tlambda\t" << lambda_ << endl;;
	x_ = vectorise((*dev_)["x"]);
	y_ = vectorise((*dev_)["y"]);

	dx_ = (*dev_)["x"](1) - (*dev_)["x"](0);
	dy_ = (*dev_)["y"](1) - (*dev_)["y"](0);

	nx_ = (*dev_)["x"].n_elem;
	ny_ = (*dev_)["y"].n_elem;
	nt_ = nx_ * ny_;

}

//void FdeSolve::computePML(int layersPML)
//{
//	//int layersPML = 10 ;
//	cout << "\t初始化PML参数\tPML层数：" << layersPML << "\n";
//	cx_mat sx(nx_, ny_, fill::ones);
//	cx_mat sy(nx_, ny_, fill::ones);
//	if (layersPML != 0) {
//		int m = 3;
//		double omega = C0 * k0_;
//		const double R0 = 1e-6;
//		double sigmaXMax = -(m + 1) * log10(R0) / (2 * 373 * dx_ * layersPML);
//		double sigmaYMax = -(m + 1) * log10(R0) / (2 * 373 * dy_ * layersPML);
//
//		cx_double value;
//		for (int i = 0;i < layersPML + 1;i++) {
//
//			value = 1.0 + sigmaXMax * pow((layersPML - i + 0.0) / layersPML, m) / (IU * omega * EPS0);
//			sx.row(i) = cx_rowvec(ny_, fill::value(value));
//			sx.row(nx_ - 1 - i) = sx.row(i);
//
//			value = 1.0 + sigmaYMax * pow((layersPML - i + 0.0) / layersPML, m) / (IU * omega * EPS0);
//			sy.col(i) = cx_vec(nx_, fill::value(value));
//			sy.col(ny_ - 1 - i) = sy.col(i);
//
//		}
//	}
//
//	sx_ = sx.as_col();
//	sy_ = sy.as_col();
// 
//}
 
//void FdeSolve::calculateIsotropicPMatrix()
//{
//	cx_vec epsZ = vectorise((*dev_)["indexZ"]) + 0.0 * IU;
//	cx_vec epsX = vectorise((*dev_)["indexX"]) + 0.0 * IU;
//	cx_vec epsY = vectorise((*dev_)["indexY"]) + 0.0 * IU;
//	cx_vec epsXY = vectorise((*dev_)["indexXY"]) + 0.0 * IU;
//	cx_vec epsYX = vectorise((*dev_)["indexYX"]) + 0.0 * IU;
//
//	cx_vec full1ColumnVector(nt_, fill::ones);
//
//	Pxx_ = dxdxFunc(sx_,sx_%epsZ, epsX, nx_, ny_, dx_, dy_)
//		+ dydyFunc(sy_, sy_,full1ColumnVector, nx_, ny_, dx_, dy_)
//		+ spdiags(k0_ * k0_ * epsX, ivec{ 0 }, nt_, nt_)
//		+ dxdyFunc(sx_, sy_%epsZ, epsYX, nx_, ny_, dx_, dy_);
//
//	Pxy_ = dxdyFunc(sx_, sy_%epsZ, epsY, nx_, ny_, dx_, dy_)
//		- dxdyFunc(sx_, sy_, full1ColumnVector, nx_, ny_, dx_, dy_)
//		+ spdiags(k0_ * k0_ * epsXY, ivec{ 0 }, nt_, nt_)
//		+ dxdxFunc(sx_, sx_%epsZ, epsXY, nx_, ny_, dx_, dy_);
//
//	Pyx_ = dydxFunc(sy_, sx_%epsZ, epsX, nx_, ny_, dx_, dy_)
//		- dydxFunc(sy_, sx_, full1ColumnVector, nx_, ny_, dx_, dy_)
//		+ spdiags(k0_ * k0_ * epsYX, ivec{ 0 }, nt_, nt_)
//		+ dydyFunc(sy_, sy_%epsZ, epsYX, nx_, ny_, dx_, dy_);
//
//	Pyy_ = dydyFunc(sy_, sy_%epsZ, epsY, nx_, ny_, dx_, dy_)
//		+ dxdxFunc(sx_, sx_, full1ColumnVector, nx_, ny_, dx_, dy_)
//		+ spdiags(k0_ * k0_ * epsY, ivec{ 0 }, nt_, nt_)
//		+ dydxFunc(sy_, sx_%epsZ, epsXY, nx_, ny_, dx_, dy_);
//}

//void FdeSolve::calculateCharacteristicValues()
//{
//	cout << "\n\t计算特征值\n";
//	sp_cx_mat pMatrix(2 * nt_, 2 * nt_);
//	pMatrix(0, 0, size(nt_, nt_)) = Pxx_;
//	pMatrix(0, nt_, size(nt_, nt_)) = Pxy_;
//	pMatrix(nt_, 0, size(nt_, nt_)) = Pyx_;
//	pMatrix(nt_, nt_, size(nt_, nt_)) = Pyy_;
//
//
//	cx_vec betaValue;
//	cx_mat fieldValue;
// 
//	double sigma = 3.4*3.4*k0_*k0_ ;
//	arma::eigs_gen(betaValue, fieldValue, pMatrix,20, sigma);
//	cx_vec neff = sqrt(betaValue) / k0_;
//	neff.print();
//
//	string outFilePath = "matlab/out.h5";
//
//	mat fieldAbs = abs(fieldValue );
//	fieldAbs.save(hdf5_name(outFilePath,"/fieldAbs"));
//	x_.save(hdf5_name(outFilePath, "/x",hdf5_opts::append));
//	y_.save(hdf5_name(outFilePath, "/y", hdf5_opts::append));
//
//}

void FdeSolve::calculatePECBoundary()
{
	cout << "\t计算PEC边界\n";
	clock_t t1 = clock();
	vec epsZ = vectorise((*dev_)["indexZ"]) ;
	vec epsX = vectorise((*dev_)["indexX"]) ;
	vec epsY = vectorise((*dev_)["indexY"])  ;
	vec epsXY = vectorise((*dev_)["indexXY"]) ;
	vec epsYX = vectorise((*dev_)["indexYX"])  ;

	vec unitColumnVector(nt_, fill::ones);
	vec sx = ones(nt_);
	vec sy = ones(nt_);

	//计算P矩阵
	cout << "\t计算P矩阵\n";
	sp_mat Pxx = dxdxFunc(sx, sx % epsZ, epsX, nx_, ny_, dx_, dy_)
		+ dydyFunc(sy, sy, unitColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsX, ivec{ 0 }, nt_, nt_)
		+ dxdyFunc(sx, sy % epsZ, epsYX, nx_, ny_, dx_, dy_);

	sp_mat Pxy = dxdyFunc(sx, sy % epsZ, epsY, nx_, ny_, dx_, dy_)
		- dxdyFunc(sx, sy, unitColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsXY, ivec{ 0 }, nt_, nt_)
		+ dxdxFunc(sx, sx % epsZ, epsXY, nx_, ny_, dx_, dy_);

	sp_mat Pyx = dydxFunc(sy, sx % epsZ, epsX, nx_, ny_, dx_, dy_)
		- dydxFunc(sy, sx, unitColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsYX, ivec{ 0 }, nt_, nt_)
		+ dydyFunc(sy, sy % epsZ, epsYX, nx_, ny_, dx_, dy_);

	sp_mat Pyy = dydyFunc(sy, sy % epsZ, epsY, nx_, ny_, dx_, dy_)
		+ dxdxFunc(sx, sx, unitColumnVector, nx_, ny_, dx_, dy_)
		+ spdiags(k0_ * k0_ * epsY, ivec{ 0 }, nt_, nt_)
		+ dydxFunc(sy, sx % epsZ, epsXY, nx_, ny_, dx_, dy_);

	cout << "\t计算特征值\n";
	sp_mat pMatrix(2 * nt_, 2 * nt_);
	pMatrix(0, 0, size(nt_, nt_)) = Pxx;
	pMatrix(0, nt_, size(nt_, nt_)) = Pxy;
	pMatrix(nt_, 0, size(nt_, nt_)) = Pyx;
	pMatrix(nt_, nt_, size(nt_, nt_)) = Pyy;

	cx_vec betaValue;
	cx_mat fieldValue;

	double sigma = search_ * search_ * k0_ * k0_;
	arma::eigs_gen(betaValue, fieldValue, pMatrix, nmodes_, sigma);
	cx_vec neff = sqrt(betaValue) / k0_;


	// 写入文件
	cout << "\t写入文件\n";
	string outFilePath = "matlab/out.h5";

	fieldValue.save(arma::hdf5_name(outFilePath, "/field"));
	neff.save(arma::hdf5_name(outFilePath, "/neff", arma::hdf5_opts::append));

	x_.save(arma::hdf5_name(outFilePath, "/x", arma::hdf5_opts::append));
	y_.save(arma::hdf5_name(outFilePath, "/y", arma::hdf5_opts::append));



	clock_t t2 = clock();
	cout << "\t耗时\t" << double(t2 - t1) / CLOCKS_PER_SEC << "s";

}

sp_mat FdeSolve::DX()
{
	vec a = ones(nt_);
	return sp_mat();
}
