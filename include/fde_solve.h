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

	void computePML(int layersPML = 10);

	//传入各向同性材料
	void calculateIsotropicPMatrix( );

	void calculateCharacteristicValues();


private:
	map<string, cube>* dev_;

	sp_cx_mat Pxx_;
	sp_cx_mat Pyy_;
	sp_cx_mat Pxy_;
	sp_cx_mat Pyx_;

	double lambda_;
	double k0_;
	int nx_;
	int ny_;
	int nt_;
	double dx_;
	double dy_;
	vec x_; 
	vec y_;

	cx_vec sx_;
	cx_vec sy_;
	cx_vec isx_;
	cx_vec isy_;


};
 