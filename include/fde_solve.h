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
	FdeSolve(map<string, cube>* dev,
		double lambda = 1.55e-6,	//波长
		int nmodes = 10,			//模式数量
		double search = 3.4);		//搜索附近值);



	~FdeSolve();

	void initialize();

	// PML边界
	//void calculatePMLBoundary();
	// PEC边界
	void calculatePECBoundary();

	// 计算DX
	sp_mat DX();

private:
	map<string, cube>* dev_;
	double lambda_ = 1.55e-6;	//波长
	int nmodes_ = 10;			//模式数量
	double search_ = 3.4;			//搜索附近值
	double k0_;
	int nx_;
	int ny_;
	int nt_;
	double dx_;
	double dy_;
	vec x_;
	vec y_;

};
