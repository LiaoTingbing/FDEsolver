// FDEsolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//


#pragma once
#include "Device.h"
#include "FDE.h"

int main()
{
    string model;
     //model  = "waveguide";
     model = "fiber";
     double lambda = 1.5e-6;
     double guess = 1.45; 
     int nmodes = 20;


    Device device("Input/" + model);
    FDE fde(device,lambda,guess,nmodes);
 
    fde.init();
    fde.solve_eigen_meta();
    fde.saveTXT("Output\\" + model);  // windows system 创建目录反斜杠

    return 0;
}


 
