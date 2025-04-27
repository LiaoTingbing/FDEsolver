


#include"../include/load_hdf5.h"
#include "../include/fde_solve.h"

int main() {

	string filePath = "lumerical/lumerical.h5";

	field<string> dataSetName{
		"indexX","indexY","indexZ",
		"x","y","lambda"
	};
	
	map<string, cube> dev;
	loadHdf5Data(dev , filePath , dataSetName);

	FdeSolve fde(&dev);
	fde.initialize();
	fde.calculateIsotropicPMatrix();
	fde.calculateCharacteristicValues();

	return 0;
}