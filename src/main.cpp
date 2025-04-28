


#include"../include/load_hdf5.h"
#include "../include/fde_solve.h"

int main() {

	string filePath = "lumerical/lumerical.h5";

	field<string> dataSetName{
		"indexX","indexY","indexZ","indexXY","indexYX",
		"x","y" 
	};
	
	double lambda = 1.55e-6;	//波长
	int nmodes = 10;			//模式数量
	double search = 1.45;			//搜索附近值
	map<string, cube> dev;
	loadHdf5Data(dev , filePath , dataSetName);

	FdeSolve fde(&dev  ,lambda,nmodes,search);
	fde.initialize();
	fde.calculatePECBoundary();
 

	return 0;
}