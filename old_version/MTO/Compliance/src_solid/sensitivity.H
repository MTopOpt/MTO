
      Info<<"sensitivity analysis"<<endl;
     
     volScalarField fsens(-gradD&&(rho*(3*gamma*gamma*E*(1-eps)/(2.0*(1.0 + Po))*twoSymm(gradD) + (3*gamma*gamma*Po*E*(1-eps)/((1.0 + Po)*(1.0 - 2.0*Po))*I)*tr(gradD)))); 
     if(solid_area==1)
     {
       setCells(fsens, cells, 0);
     }    

      #include "Heaviside_sens.H"

      //constraint value
     gx[0]=(PetscScalar)(V);

     VecGetArray(dfdx,&dfdx_array);
     VecGetArrays(dgdx,m,&dgdx_array);
     for(i=0;i<allcells;i++)
     {
       dfdx_array[i]=fsens[i]/scalef/nallcells;//sensitivity of objective function
       dgdx_array[0][i]=gsens[i];//sensitivity of constraint function
     }
     VecRestoreArray(dfdx,&dfdx_array);
     VecRestoreArrays(dgdx,m,&dgdx_array);

     Info<< "\n run MMA \n" << endl;
     mma->SetOuterMovelimit(Xmin,Xmax,movlim,xmma,xmin,xmax);
     mma->Update(xmma,dfdx,gx,dgdx,xmin,xmax,raa0,mma_dec,mma_inc);

     VecGetArray(xmma,&xmma_array);
     for(i=0;i<allcells;i++)
     {
       gamma[i]=xmma_array[i];
     }
     VecRestoreArray(xmma,&xmma_array);

     #include "Heaviside_rho.H"

     //update mu and lambda
     opt++;
     mu=(E*eps+gamma*gamma*gamma*(1-eps)*E)/(2.0*(1.0 + Po));
     lambda=(E*eps+gamma*gamma*gamma*(1-eps)*E)*Po/((1.0 + Po)*(1.0 - 2.0*Po));
     if (planeStress)
     {
        lambda = (E*eps+gamma*gamma*gamma*(1-eps)*E)*Po/((1.0 + Po)*(1.0 - Po));
     }




