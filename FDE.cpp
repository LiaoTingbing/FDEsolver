#pragma once

#include "FDE.h"

FDE::FDE()
{
}


FDE::FDE(Device device_, double lambda_, double searchIndex_, int searchNum_)
{
	device = device_;
	lambda = lambda_;
	guess = searchIndex_;
	nmodes = searchNum_;
	//this->device.lambda = lambda;
	//this->device.erx = pow(this->device.index_x, 2);
	//this->device.ery = pow(this->device.index_y, 2);
	//this->device.erz = pow(this->device.index_z, 2);
}



FDE::~FDE()
{
}

void FDE::solve()
{
	size_t ng = device.Nx * device.Ny;

	sp_mat I = speye<sp_mat>(ng, ng);
	sp_mat Ux(ng, ng), Uy(ng, ng), Vx(ng, ng), Vy(ng, ng);

	sp_mat  Erx(ng, ng), Erx_inv(ng, ng), Ery(ng, ng);
	sp_mat Ery_inv(ng, ng), Erz(ng, ng), Erz_inv(ng, ng);

	//Erx = spdiags(device.erx.as_col());

	for (size_t i = 0; i < ng; i++)
	{
		Erx(i, i) = device.erx(i);
		Ery(i, i) = device.ery(i);
		Erz(i, i) = device.erz(i);

		Erx_inv(i, i) = 1.0 / device.erx(i);
		Ery_inv(i, i) = 1.0 / device.ery(i);
		Erz_inv(i, i) = 1.0 / device.erz(i);

		Ux(i, i) = -1.0;
		if (i + 1 < ng and i % device.Nx != device.Nx - 1)
		{
			Ux(i, i + 1) = 1.0;
		}

		Uy(i, i) = -1.0;
		if (i + device.Nx < ng)
		{
			Uy(i, i + device.Nx) = 1.0;
		}
	}

	Ux = Ux / device.dx;
	Uy = Uy / device.dy;

	Vx = -Ux.st();
	Vy = -Uy.st();

	double k0 = 2.0 * pi / lambda;
	//% Calculate Pxx
	sp_mat	Pxx = -pow(k0, -2) * Ux * Erz_inv * Vy * Vx * Uy + (pow(k0, 2) * I + Ux * Erz_inv * Vx) * (Erx + pow(k0, -2) * Vy * Uy);
	//% Calculate Pyy
	sp_mat	Pyy = -pow(k0, -2) * Uy * Erz_inv * Vx * Vy * Ux + (pow(k0, 2) * I + Uy * Erz_inv * Vy) * (Ery + pow(k0, -2) * Vx * Ux);
	//% Calculate Pxy
	sp_mat	Pxy = Ux * Erz_inv * Vy * (Ery + pow(k0, -2) * Vx * Ux) - pow(k0, -2) * (pow(k0, 2) * I + Ux * Erz_inv * Vx) * Vy * Ux;
	//% Calculate Pyx
	sp_mat	Pyx = Uy * Erz_inv * Vx * (Erx + pow(k0, -2) * Vy * Uy) - pow(k0, -2) * (pow(k0, 2) * I + Uy * Erz_inv * Vy) * Vx * Uy;

	sp_mat P4 = join_cols(join_rows(Pxx, Pxy), join_rows(Pyx, Pyy));

	double sigma = pow(guess, 2.0);

	SparseMatrix<double> P4_ = armaToEigenSparseManual(P4 / k0 / k0);


	VectorXcd eigval;
	MatrixXcd eigvec;
	sparseEigs(eigval, eigvec,
		P4_, nmodes, sigma);

	VectorXcd NEFF = eigval.array().sqrt().matrix();

	cout << NEFF << endl;


	//cx_vec neff = sqrt(eigval) / k0;
	//neff.print();

}

void FDE::solve_eigen()
{
	clock_t t1, t2;
	t1 = clock();

	const size_t  Nx = device.Nx;
	const size_t  Ny = device.Ny;

	size_t ng = Nx * Ny;

	SparseMatrix<double> I(ng, ng);
	SparseMatrix<double> Erx(ng, ng);
	SparseMatrix<double> Ery(ng, ng);
	SparseMatrix<double> Erz(ng, ng);
	SparseMatrix<double> Erx_inv(ng, ng);
	SparseMatrix<double> Ery_inv(ng, ng);
	SparseMatrix<double> Erz_inv(ng, ng);

	SparseMatrix<double> Ux(ng, ng);
	SparseMatrix<double> Uy(ng, ng);
	SparseMatrix<double> Vx(ng, ng);
	SparseMatrix<double> Vy(ng, ng);

	for (size_t i = 0; i < ng; i++)
	{
		I.insert(i, i) = 1.0;

		Erx.insert(i, i) = device.erx(i);
		Ery.insert(i, i) = device.ery(i);
		Erz.insert(i, i) = device.erz(i);
		Erx_inv.insert(i, i) = 1.0 / device.erx(i);
		Ery_inv.insert(i, i) = 1.0 / device.ery(i);
		Erz_inv.insert(i, i) = 1.0 / device.erz(i);

		Ux.insert(i, i) = -1.0;
		if (i + 1 < ng and i % Nx != Nx - 1)
		{
			Ux.insert(i, i + 1) = 1.0;
		}
		Uy.insert(i, i) = -1.0;
		if (i + device.Nx < ng)
		{
			Uy.insert(i, i + device.Nx) = 1.0;
		}
	}

	I.makeCompressed();
	Erx.makeCompressed();
	Ery.makeCompressed();
	Erz.makeCompressed();
	Erx_inv.makeCompressed();
	Ery_inv.makeCompressed();
	Erz_inv.makeCompressed();
	Ux.makeCompressed();
	Uy.makeCompressed();


	Ux = Ux / device.dx;
	Uy = Uy / device.dy;

	Vx = -Ux.transpose();
	Vy = -Uy.transpose();

	//double lambda = device.lambda(0);
	double k0 = 2 * pi / lambda;


	SparseMatrix<double> Pxx = -pow(k0, -2) * Ux * Erz_inv * Vy * Vx * Uy + (pow(k0, 2) * I + Ux * Erz_inv * Vx) * (Erx + pow(k0, -2) * Vy * Uy);
	SparseMatrix<double> Pyy = -pow(k0, -2) * Uy * Erz_inv * Vx * Vy * Ux + (pow(k0, 2) * I + Uy * Erz_inv * Vy) * (Ery + pow(k0, -2) * Vx * Ux);
	SparseMatrix<double> Pxy = Ux * Erz_inv * Vy * (Ery + pow(k0, -2) * Vx * Ux) - pow(k0, -2) * (pow(k0, 2) * I + Ux * Erz_inv * Vx) * Vy * Ux;
	SparseMatrix<double> Pyx = Uy * Erz_inv * Vx * (Erx + pow(k0, -2) * Vy * Uy) - pow(k0, -2) * (pow(k0, 2) * I + Uy * Erz_inv * Vy) * Vx * Uy;

	SparseMatrix<double> P12 = horizontalConcatenation(Pxx, Pxy);
	SparseMatrix<double> P34 = horizontalConcatenation(Pyx, Pyy);

	SparseMatrix<double> P4 = verticalConcatenation(P12, P34);

	//cout << P4.coeff(0, 0);
	SparseMatrix<double> P4_ = P4 / k0 / k0;

	double sigma = guess * guess;    // Shift 值 (要查找接近 10.0 的特征值)

	VectorXcd eigval;
	MatrixXcd eigvec;
	sparseEigs(eigval, eigvec,
		P4_, nmodes, sigma);

	VectorXcd NEFF = eigval.array().sqrt().matrix();

	cout << NEFF << endl;

	t2 = clock();
	cout << double(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;
}

void FDE::solve_eigen_meta()
{

	size_t Nx = device.Nx;
	size_t Ny = device.Ny;
	size_t ng = Nx * Ny;

	double dx = device.dx;
	double dy = device.dy;

	vec Iv(ng, fill::ones);
	sp_mat I = speye<sp_mat>(ng, ng);

	sp_mat Ux(ng, ng), Uy(ng, ng), Vx(ng, ng), Vy(ng, ng);
	sp_mat  ERX(ng, ng), ERY(ng, ng), ERZ(ng, ng);
	sp_mat  ERX_inv(ng, ng), ERY_inv(ng, ng), ERZ_inv(ng, ng);

	ERX.diag() = vectorise(device.erx);
	ERY.diag() = vectorise(device.ery);
	ERZ.diag() = vectorise(device.erx);
	ERX_inv.diag() = vectorise(1.0 / device.erx);
	ERY_inv.diag() = vectorise(1.0 / device.ery);
	ERZ_inv.diag() = vectorise(1.0 / device.erz);

	sp_mat URX = I;
	sp_mat URY = I;
	sp_mat URZ = I;
	sp_mat URX_inv = I;
	sp_mat URY_inv = I;
	sp_mat URZ_inv = I;

	Ux.diag(0) = -Iv;
	Ux.diag(1) = vec(ng-1,fill::ones);
	Uy.diag(0) = -Iv;
	Uy.diag(Nx) = vec(ng - Nx, fill::ones);

	for (size_t i = Nx-1; i +1< ng; i = i + Nx)
	{
		Ux(i,  i + 1) = 0;
	}

	Ux = Ux / dx;
	Uy = Uy / dy;

	Vx = -Ux.st();
	Vy = -Uy.st();

	double k0 = 2 * pi / lambda;
 
	sp_mat P(2 * ng, 2 * ng);
	P.submat(0, 0, size(ng, ng)) = Ux * ERZ_inv * Vy;
	P.submat(0, ng, size(ng, ng)) = -(k0 * k0 * URY + Ux * ERZ_inv * Vx);
	P.submat(ng, 0, size(ng, ng)) = k0 * k0 * URX + Uy * ERZ_inv * Vy;
	P.submat(ng, ng, size(ng, ng)) = -Uy * ERZ_inv * Vx;
	P = P / k0;

	sp_mat Q(2 * ng, 2 * ng);
	Q.submat(0, 0, size(ng, ng)) = -Vx * URZ_inv * Uy;
	Q.submat(0, ng, size(ng, ng)) = (k0 * k0 * ERY + Vx * URZ_inv * Ux);
	Q.submat(ng, 0, size(ng, ng)) = -(k0 * k0 * ERX + Vy * URZ_inv * Uy);
	Q.submat(ng, ng, size(ng, ng)) = Vy * URZ_inv * Ux;
	Q = Q / k0;

	sp_mat PQ = P * Q / k0 / k0   ; 

	SparseMatrix<double> P4_ = armaToEigenSparseManual(PQ  );

	double sigma = guess * guess  ;    // Shift 值 (要查找接近 10.0 的特征值)
	VectorXcd eigval;
	MatrixXcd eigvec;
	sparseEigs(eigval, eigvec,
		P4_, nmodes, sigma);

	//VectorXcd NEFF = eigval.array().sqrt().matrix()  ;

	cx_vec D(eigval.size());
	cx_mat Et(eigvec.rows(), eigvec.cols());

	for (size_t i = 0; i < eigval.size(); i++)
	{
		D(i) = eigval(i);
		for (size_t j = 0; j < eigvec.rows(); j++)
		{
			Et(j,i) = eigvec(j, i);
		}
	}

	cx_vec NEFF = sqrt(D)   ;
	cx_vec BETA = k0 * NEFF; 




	cout << NEFF << endl;

 
	cx_mat Ht(2 * ng, nmodes);
	cx_mat Ez(ng, nmodes);
	cx_mat Hz(ng, nmodes);

	cx_vec ex,ey,hx,hy;

	for (size_t i = 0; i < nmodes; i++)
	{
		Ht.col(i) = 1.0 / BETA(i) * Q * Et.col(i);
		ex = Et.submat(0, i, size(ng, 1));
		ey = Et.submat(ng, i, size(ng, 1));
		hx = Ht.submat(0, i, size(ng, 1));
		hy = Ht.submat(ng, i, size(ng, 1));
		Hz.col(i) = -iu * URZ_inv * (-Ux * ey + Uy * ex) / k0;
		Ez.col(i) = -iu * ERZ_inv * ( Vx * hy - Vy * hx) / k0;
	}


	cx_mat Ex = Et.submat(0, 0, size(ng, Et.n_cols));
	cx_mat Ey = Et.submat(ng, 0, size(ng, Et.n_cols));
	cx_mat Hx = Ht.submat(0, 0, size(ng, Ht.n_cols));
	cx_mat Hy = Ht.submat(ng, 0, size(ng, Ht.n_cols));

 
	string Opfile = "Output/";
	system("mkdir Output");

	Ex = Ex * sqrt(mu0 / eps0);
	Ey = Ey * sqrt(mu0 / eps0);
	Ez = Ez * sqrt(mu0 / eps0);

	mat Ex_real = real(Ex);
	mat Ex_imag = imag(Ex);
	mat Ey_real = real(Ey);
	mat Ey_imag = imag(Ey);
	mat Ez_real = real(Ez);
	mat Ez_imag = imag(Ez);
 
	Ex_real.save(Opfile + "/Ex_real.txt", raw_ascii);
	Ex_imag.save(Opfile + "/Ex_imag.txt", raw_ascii);
	Ey_real.save(Opfile + "/Ey_real.txt", raw_ascii);
	Ey_imag.save(Opfile + "/Ey_imag.txt", raw_ascii);
	Ez_real.save(Opfile + "/Ez_real.txt", raw_ascii);
	Ez_imag.save(Opfile + "/Ez_imag.txt", raw_ascii);


	
	mat Hx_real = real(Hx);
	mat Hx_imag = imag(Hx);
	mat Hy_real = real(Hy);
	mat Hy_imag = imag(Hy);
	mat Hz_real = real(Hz);
	mat Hz_imag = imag(Hz);

	Hx_real.save(Opfile + "/Hx_real.txt", raw_ascii);
	Hx_imag.save(Opfile + "/Hx_imag.txt", raw_ascii);
	Hy_real.save(Opfile + "/Hy_real.txt", raw_ascii);
	Hy_imag.save(Opfile + "/Hy_imag.txt", raw_ascii);
	Hz_real.save(Opfile + "/Hz_real.txt", raw_ascii);
	Hz_imag.save(Opfile + "/Hz_imag.txt", raw_ascii);

	device.x.save(Opfile + "/x.txt", raw_ascii);
	device.y.save(Opfile + "/y.txt", raw_ascii);
	device.lambda.save(Opfile + "/lambda.txt", raw_ascii);

	
	vec neff_real = real(NEFF);
	vec neff_imag = imag(NEFF);
	
	neff_real.save(Opfile + "/neff_real.txt", raw_ascii);
	neff_imag.save(Opfile + "/neff_imag.txt", raw_ascii);
















}
