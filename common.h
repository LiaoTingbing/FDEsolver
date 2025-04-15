#pragma once





#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>



//#define ARMA_USE_ARPACK 
//#define ARMA_USE_SUPERLU
//#define ARMA_USE_SUPERLU    
////#define ARMA_USE_ARPACK   
//#define ARMA_SUPERLU_INCLUDE_DIR  E:/VisualStudioFiles/VCPKG/vcpkg-master/vcpkg-master/packages/superlu_x64-windows/include/

//#include <arpack/arpack.h>
//#include <slu_cdefs.h>
#include <armadillo>
#include <iostream>
#include <cstring>
#include <ctime>

using namespace std;
using namespace arma;
using namespace Eigen;
using namespace Spectra;



const std::complex<double> iu(0.0, 1.0);
const double eps0 = 8.85419e-12;
const double mu0 = 1.25663706e-6;
const double c0 = 2.99792458e8;
const double pi = 3.14159265358979323846;




