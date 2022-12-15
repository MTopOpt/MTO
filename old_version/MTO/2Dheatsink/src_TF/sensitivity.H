
Info<<"sensitivity analysis"<<endl;

fsens=-alphaMax*(1+qu)*qu/((qu+gamma)*(qu+gamma))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+gamma)*(qu+gamma))*(fvc::grad(T)&fvc::grad(Tb))/rhoc;

gsens_U=-alphaMax*(1+qu)*qu/((qu+gamma)*(qu+gamma))*(U&Uc);

#include "Heaviside_sens.H"
     
gx[0]=(PetscScalar)(DissPower/Foam::max(D0-opt*0.2,D1)-1.0);
gx[1]=(PetscScalar)(V);

VecGetArray(dfdx,&dfdx_array);
VecGetArrays(dgdx,m,&dgdx_array);
for(i=0;i<allcells;i++)
{
  dfdx_array[i]=fsens[i]/scalef/nallcells;
  dgdx_array[0][i]=gsens_U[i]/scaleg/nallcells;
  dgdx_array[1][i]=gsens[i];    
}
VecRestoreArray(dfdx,&dfdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

Info<< "\nrun mmasub\n" << endl;
mma->SetOuterMovelimit(Xmin,Xmax,movlim,xmma,xmin,xmax);
mma->Update(xmma,dfdx,gx,dgdx,xmin,xmax,raa0,mma_dec,mma_inc);

VecGetArray(xmma,&xmma_array);
for(i=0;i<allcells;i++)
{
  gamma[i]=xmma_array[i];
}
VecRestoreArray(xmma,&xmma_array);

#include "Heaviside_rho.H"

opt++;

if(opt<=100)
{
  alphaMax=alphamax*(opt/7.0+1.0);
}
if(opt>100)
{
  alphaMax=alphaMax*1.05;
}
alphaMax.value()=Foam::min(alphaMax.value(),10000000.0);

if(opt>=60)
{
  qu=0.005+(opt-60)*1e-4;
}
qu=Foam::min(qu,0.01);
qu.correctBoundaryConditions();
alpha=alphaMax*qu*(1-gamma)/(qu+gamma);
DT=(ks+(kf-ks)*gamma*(1+qu)/(qu+gamma))/rhoc;





