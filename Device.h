#pragma once
#include "common.h"
class Device
{
public:
	// �豸����

	double Lx, Ly, Lz; // �豸�ߴ�
	size_t Nx, Ny, Nz; // ������Ŀ
	double dx, dy, dz; // �����С
	
	vec lambda;

	vec x, y, z;
	mat index_x, index_y, index_z;
	mat erx, ery, erz;
	mat MURX, MURY, MURZ;

	// �߽�����
	int BCtype; // �߽���������

public:
	Device();
	Device(string InputFilePath);
	~Device();


	
};

