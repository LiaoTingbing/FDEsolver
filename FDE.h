#pragma once

#include "Device.h"
#include "EigenMatrixConnect.h"
#include "sparseEigs.h"
#include "armaToEigen.h"

class FDE 
{
	Device  device;
	double lambda;
	double guess;
	int nmodes;

	size_t Nx, Ny, ng;
	double dx, dy, k0;

	cx_vec NEFF;
	cx_mat Ex, Ey, Ez;
	cx_mat Hx, Hy, Hz;

	//string OutputPath;

	//vec neff_real;
	//vec neff_imag;

	//mat Ex_real;
	//mat Ex_imag ;
	//mat Ey_real;
	//mat Ey_imag ;
	//mat Ez_real ;
	//mat Ez_imag ;

	//mat Hx_real ;
	//mat Hx_imag ;
	//mat Hy_real ;
	//mat Hy_imag ;
	//mat Hz_real ;
	//mat Hz_imag ;

public:
	FDE();

	FDE(Device device_,double lambda_  = 1.55*1e-6, double searchIndex_ = 3.4, int searchNum_ = 4);

 
	~FDE();

	void init();
	void solve_eigen_meta();
	void saveTXT(string OutputPath = "Output/");

	void run();



};

