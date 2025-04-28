#pragma once

#include "common.h"
#include "dxdx_function.h"
#include "dxdy_function.h"
#include "dydx_function.h"
#include "dydy_function.h"



class FdeSolve
{

public:
	FdeSolve();
	FdeSolve(map<string, cube>* dev);

	~FdeSolve();

	void initialize();
 

	// PEC边界
	void calculatePECBoundary();


private:
	map<string, cube>* dev_;

	double lambda_;
	double k0_;
	int nx_;
	int ny_;
	int nt_;
	double dx_;
	double dy_;
	vec x_; 
	vec y_;

};
 