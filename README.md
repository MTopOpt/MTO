 MTO——Multiphysics Topology Optimization
=========================================
Description
-----------
The code presented is a parallel solver for multi-physics topology optimization on structured or unstructured  grids. Density-based method is adopted for solving 2D or 3D thermal-fluid-structural topology optimization problems based on OpenFOAM. 

Please cite: Yu, M., Ruan, S., Gu, J. et al. Three-dimensional topology optimization of thermal-fluid-structural problems for cooling system design. Struct Multidisc Optim (2020). https://doi.org/10.1007/s00158-020-02731-z 

My PHD thesis: 《液冷散热器流道与热源分布拓扑优化研究》


Installation
------------
**The old version of MTO is difficult to install, so I made some changes in the code.**

Before running this solver, following softwares are needed.  

(a) **OpenFOAM 6.0**  (www.openfoam.org)

(b) **swak4Foam** (if you use thermal-fluid solver) (https://openfoamwiki.net/index.php/Contrib/swak4Foam)

**Here are some installation instructions:**

1. OpenFOAM 6.0 is only valid on Ubuntu16.04 - 18.04, so you should use the correct Ubuntu version.

2. Install OpenFOAM 6.0 (see appendix documents OpenFOAM_Ubuntu.md for details)

3. cd /opt/openfoam6/wmake/rules/linux64Gcc then open the file c and relpace the fifth line gcc with mpicc

4. cd /opt/openfoam6/wmake/rules/linux64Gcc then open the file c++ and relpace the fifth line g++ with mpic++

Run the solver
--------------
 After finishing above works, users should run **wmake** in the src folder, and then run **blockMesh** and **decomposePar** in the app folder. Finally, the multi-physics topology optimization solver can run in the app folder by **mpirun –n 20 MTO_ThermalFluid –parallel**, where 20 is the default number of cores for parallel computing.
 
Postplot
--------
After the optimization, users should run **reconstructPar** in the app folder and then run **paraFoam** to view the results with Paraview.  
