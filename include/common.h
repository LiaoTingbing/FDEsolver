
#pragma once

#define ARMA_USE_HDF5
#define ARMA_USE_SUPERLU
#define ARMA_USE_ARPACK
#include <armadillo>
#include<cstring>
#include<ctime>

using namespace std;
using namespace arma;

const double EPS0 = 8.85419e-12;
const double PI = 3.14159;
const double MU0 = PI * 4e-7;
const double C0 = 1 / sqrt(EPS0 * MU0);
const cx_double IU(0, 1);

 