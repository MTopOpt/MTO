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

if(opt<=63)
{
  alphaMax=alphamax*(opt/7.0+1.0);
}
if(opt>100)
{
  alphaMax=alphaMax*1.05;
}
alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
if(opt>=80)
{
  qu=0.005+(opt-80)*1e-4;
}
qu=Foam::min(qu,0.01);
alpha=alphaMax*qu*(1-xh)/(qu+xh);
DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoc;
