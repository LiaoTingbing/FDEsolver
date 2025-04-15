#pragma once
#include "common.h"
class Device
{
public:
	// 设备参数

	double Lx, Ly, Lz; // 设备尺寸
	size_t Nx, Ny, Nz; // 网格数目
	double dx, dy, dz; // 网格大小
	
	vec lambda;

	vec x, y, z;
	mat index_x, index_y, index_z;
	mat erx, ery, erz;
	mat MURX, MURY, MURZ;

	// 边界条件
	int BCtype; // 边界条件类型

public:
	Device();
	Device(string InputFilePath);
	~Device();


	
};

