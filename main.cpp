// FDEsolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//


#pragma once
#include "Device.h"
#include "FDE.h"

int main()
{

    Device device("E:\\研究生\\FDE\\20250412_3D\\lumerical");
    FDE fde(device,1.5e-6,3.4,20);
 
    fde.run();

    return 0;
}


 
