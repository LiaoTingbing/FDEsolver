

#pragma once

#include "Device.h"

Device::Device()
{
}

Device::Device( string InputFilePath)
{
	index_x.load(InputFilePath+"/index_x.txt");
	index_y.load(InputFilePath+"/index_y.txt");
	index_z.load(InputFilePath+"/index_z.txt");

	erx = pow(index_x, 2);
	ery = pow(index_y, 2);
	erz = pow(index_z, 2);

	//erx.print();

	lambda.load(InputFilePath+"/lambda.txt");
	//lambda.print();

	x.load(InputFilePath+"/x.txt");
	y.load(InputFilePath+"/y.txt");
	z.load(InputFilePath+"/z.txt");

	dx = x(1) - x(0);
	dy = y(1) - y(0);
	//dx = z(1) - z(0);

	Nx = x.size();
	Ny = y.size();
	Nz = z.size();

	Lx = x(x.size() - 1) - x(0);
	Ly = y(y.size() - 1) - y(0);

}

Device::~Device()
{
}
