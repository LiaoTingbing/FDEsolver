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

void FDE::init()
{
	Nx = device.Nx;
	Ny = device.Ny;
	ng = Nx * Ny;
	dx = device.dx;
	dy = device.dy;
	k0 = 2 * pi / lambda;
}

void FDE::solve_eigen_meta()
{

	// Construct Difference Matrix
	sp_mat I = speye<sp_mat>(ng, ng);

	sp_mat  ERX(ng, ng), ERY(ng, ng), ERZ(ng, ng);
	sp_mat  ERX_inv(ng, ng), ERY_inv(ng, ng), ERZ_inv(ng, ng);
	sp_mat Ux(ng, ng), Uy(ng, ng), Vx(ng, ng), Vy(ng, ng);

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

	Ux.diag(0) = -vec(ng, fill::ones);
	Ux.diag(1) = vec(ng - 1, fill::ones);
	Uy.diag(0) = -vec(ng, fill::ones);;
	Uy.diag(Nx) = vec(ng - Nx, fill::ones);

	for (size_t i = Nx - 1; i + 1 < ng; i = i + Nx)
	{
		Ux(i, i + 1) = 0;
	}

	Ux = Ux / dx;
	Uy = Uy / dy;

	Vx = -Ux.st();
	Vy = -Uy.st();


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

	sp_mat PQ = P * Q / k0 / k0;

	// Solve Eigs
	SparseMatrix<double> P4_ = armaToEigenSparseManual(PQ);
	double sigma = guess * guess;    // Shift 值 (要查找接近 10.0 的特征值)
	VectorXcd eigval;
	MatrixXcd eigvec;
	sparseEigs(eigval, eigvec,
		P4_, nmodes, sigma);

	// ARMA eigs

	cx_vec eigvala;
	cx_mat eigveca;

	eigs_gen(eigvala, eigveca, PQ, nmodes, sigma);  // find 5 eigenvalues/eigenvectors
	eigvala = sqrt(eigvala);
	eigvala.print();
	//eigs_opts opts;
	//opts.maxiter = 10000;            // increase max iterations to 10000

	//eigs_gen(eigval, eigvec, A, 5, "lm", opts);
	//

	cx_vec D(eigval.size());
	cx_mat Et(eigvec.rows(), eigvec.cols());

	for (size_t i = 0; i < eigval.size(); i++)
	{
		D(i) = eigval(i);
		for (size_t j = 0; j < eigvec.rows(); j++)
		{
			Et(j, i) = eigvec(j, i);
		}
	}

	// Solve Fields
	NEFF = sqrt(D);
	cx_vec BETA = k0 * NEFF;
	cx_mat Ht(2 * ng, nmodes);

	Ez = cx_mat(ng, nmodes);
	Hz = cx_mat(ng, nmodes);

	cx_vec ex, ey, hx, hy;

	for (size_t i = 0; i < nmodes; i++)
	{
		Ht.col(i) = 1.0 / BETA(i) * Q * Et.col(i);
		ex = Et.submat(0, i, size(ng, 1));
		ey = Et.submat(ng, i, size(ng, 1));
		hx = Ht.submat(0, i, size(ng, 1));
		hy = Ht.submat(ng, i, size(ng, 1));
		Hz.col(i) = -iu * URZ_inv * (-Ux * ey + Uy * ex) / k0;
		Ez.col(i) = -iu * ERZ_inv * (Vx * hy - Vy * hx) / k0;
	}

	Ex = Et.submat(0, 0, size(ng, Et.n_cols));
	Ey = Et.submat(ng, 0, size(ng, Et.n_cols));
	Hx = Ht.submat(0, 0, size(ng, Ht.n_cols));
	Hy = Ht.submat(ng, 0, size(ng, Ht.n_cols));

	Ex = Ex * sqrt(mu0 / eps0);
	Ey = Ey * sqrt(mu0 / eps0);
	Ez = Ez * sqrt(mu0 / eps0);

	cout << NEFF << endl;
}

void FDE::saveTXT(string OutputPath)
{
	string Opfile = OutputPath;
	string command = "mkdir ";
	string s =  command +  OutputPath;
	system(s.c_str());

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

	system("pause");
}

 


