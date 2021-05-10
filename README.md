 MTO——Multiphysics Topology Optimization
=========================================
Description
-----------
The code presented is a parallel solver for multi-physics topology optimization on structured grids. Density-based method is adopted for solving 2D or 3D thermal-fluid-structural topology optimization problems based on OpenFOAM. 

Yu, M., Ruan, S., Gu, J. et al. Three-dimensional topology optimization of thermal-fluid-structural problems for cooling system design. Struct Multidisc Optim (2020). https://doi.org/10.1007/s00158-020-02731-z 

Installation
------------
Before running this solver, following softwares are needed.  
(a) **OpenFOAM 5.0**  
(b) **PETSc 3.12**  
(c) **swak4Foam**  development version（my version is 0.4.3）
(d) **openMPI**   (my version is 1.10）

Parallel version of MMA written by Aage is also needed for updating the design variables (https://github.com/topopt/TopOpt_in_PETSc/MMA.cc). We have changed this code to apply it in our solver.Users should compile the C++ script (MMA.C MMA.H) into a dynamic library (**libMMA_yu.so**) and put this file in **FOAM_USER_LIBBIN** after installing OpenFOAM.  

Run the solver
--------------
 After finishing above works, users should run **wmake** in the src folder, and then run **blockMesh** and **decomposePar** in the app folder. Finally, the multi-physics topology optimization solver can run in the app folder by **mpirun –n 20 MTO_solid –parallel**, where 20 is the default number of cores for parallel computing.
 
Postplot
--------
After the optimization, users should run **reconstructPar** in the app folder and then run **paraFoam** to view the results with Paraview.  

Solid displacement problem-2D  
-----------------------------
![image](https://github.com/MTopOpt/MTO/blob/master/MTO/beam_2D.gif)  

Solid displacement problem-3D  
-----------------------------
![image](https://github.com/MTopOpt/MTO/blob/master/MTO/beam_3D.gif)  
Thermal-fluid problem-2D  
------------------------
![image](https://github.com/MTopOpt/MTO/blob/master/MTO/heatsink_2D.gif)  
Thermal-fluid problem-3D  
------------------------
![image](https://github.com/MTopOpt/MTO/blob/master/MTO/heatsink_3D.gif)  

![image](https://github.com/MTopOpt/MTO/blob/master/MTO/12.gif)  

![image](https://github.com/MTopOpt/MTO/blob/master/MTO/13.gif)  


