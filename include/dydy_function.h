﻿#pragma once

#include "common.h"


sp_cx_mat dydyFunc(const cx_vec& p, const  cx_vec& q, const cx_vec& r,
	int nx, int ny, double dx, double dy);

sp_mat dydyFunc(const  vec& p, const   vec& q, const  vec& r,
	int nx, int ny, double dx, double dy);