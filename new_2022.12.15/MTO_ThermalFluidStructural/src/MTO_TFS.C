//Author: Yu Minghao    Updated: 2021.8.28 

static char help[] = "topology optimization of thermal-fluid-structural problem\n";
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"//
#include "MMA/MMA.h"
#include <diff.c>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"
    #include "createFields.H"//一些物理场的声明
    #include "readTransportProperties.H" //一些物理场的声明
    #include "readThermalProperties.H" //一些物理场的声明
    #include "readMechanicalProperties.H" //一些物理场的声明
     #include "initContinuityErrs.H"
    #include "opt_initialization.H"//一些参数的声明和初始化
    while (simple.loop(runTime))
    {
        #include "update.H"//根据设计变量更新材料属性
        #include "NS.H"//NS方程求解
        #include "HeatTransfer.H"//传热方程求解
        #include "LinearElasticity.H"//平衡方程求解
        #include "AdjHeatTransfer.H"//伴随能量方程(of 平均温度)
        #include "AdjNS_HT.H"//伴随NS方程(of 平均温度)
        #include "AdjNS_FF.H"//伴随能量方程(of 能量耗散)
        #include "costfunction.H"//计算目标和约束函数值             
        #include "sensitivity.H"//灵敏度分析及MMA
    }
    #include "finalize.H"//程序终止，delete变量
    return 0;
}
