#pragma once

#include "common.h"



sp_cx_mat dxdxFunc(const cx_vec& p, const cx_vec& q, const cx_vec& r,
	int nx, int ny,
	double dx, double dy);


sp_mat dxdxFuncReal(const vec& p,  vec& q,  vec& r,
	int nx, int ny,
	double dx, double dy);

void dxdxFuncBoth( sp_mat &result,
	const vec& p, vec& q, vec& r,
	int nx, int ny,
	double dx, double dy);