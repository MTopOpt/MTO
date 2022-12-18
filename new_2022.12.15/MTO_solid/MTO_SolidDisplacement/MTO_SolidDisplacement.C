//Author: Yu Minghao    Updated: May 2020     Email:yuminghao_dlut@163.com

static char help[] = "topology optimization of linear elasticity problem\n";
#include "fvCFD.H"
#include "simpleControl.H"
#include "MMA/MMA.h"
#include <diff.c>
#include <vector>
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "readMechanicalProperties.H" 
    #include "opt_initialization.H"

    while (simple.loop(runTime))
    {
        #include "update.H"
        #include "LinearElasticity.H"
        #include "costfunction.H"              
        #include "sensitivity.H"
    }
    #include "finalize.H"
    return 0;
}
