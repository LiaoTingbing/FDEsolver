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

public:
	FDE();

	FDE(Device device_,double lambda_  = 1.55*1e-6, double searchIndex_ = 3.4, int searchNum_ = 4);

 
	~FDE();

	void solve();
	void solve_eigen();


	void solve_eigen_meta();



};

