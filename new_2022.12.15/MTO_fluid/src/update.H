//***********************************//update mu and lambda
opt++;
if(solid_area)
{
   forAll(cells_solid, i)
   {
      xh[cells_solid[i]]=0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      xh[cells_fluid[i]]=1.0;
   }
}
xh.correctBoundaryConditions();
if(runTime.writeTime())
{
   xh.write();  
}

if(opt<=27)
{
  alphaMax=alphamax*(opt/3.0+1.0);
}

alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
// alphaMax.value()=alphaMax.value();

alpha=alphaMax*qu*(1-xh)/(qu+xh);
