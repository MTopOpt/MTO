
opt++;

//不可设计域
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
   T.write();
}

//流体问题一些参数的放松
rex=63;
if(opt<=rex)//45,95,145
{
  alphaMax=alphamax*(opt/7.0+1.0);
}
if(opt>rex)
{
  alphaMax=alphaMax*1.05;
}

alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
if(opt>=rex)
{
  qu=0.01+(opt-rex)*1e-4;
}
qu=Foam::min(qu,0.02);

//材料参数更新
alpha=alphaMax*qu*(1-xh)/(qu+xh);
DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoc;
mu=qd*(1-xh)/(qd+xh)*(E-Emin)/(2.0*(1.0+Po))+Emin/(2.0*(1.0+Po));
lambda=qd*(1-xh)/(qd+xh)*Po*(E-Emin)/((1.0+Po)*(1.0-2.0*Po))+Po*Emin/((1.0+Po)*(1.0-2.0*Po));
if (planeStress)
{
   lambda=qd*(1-xh)/(qd+xh)*Po*(E-Emin)/((1.0+Po)*(1.0-Po))+Po*Emin/((1.0+Po)*(1.0-Po));
}