 Info<< "\n Heaviside_rho\n" << endl;

localr=gamma.internalField();
AllrList[Pstream::myProcNo()] = localr;
Pstream::gatherList(AllrList);
Pstream::scatterList(AllrList);
Allgamma=ListListOps::combine<scalarField>(AllrList,accessOp<scalarField>());

for(i=0;i<nallcells;i++)
{
   add1=0;
   add2=0;
   for(j=0;j<Nnum[i];j++)
   {
      add1=add1+Nweight[i][j]*Allgamma[Nlabel[i][j]];
      add2=add2+Nweight[i][j];
   }
   gamma_densityfilter[i]=add1/add2;
}

if(opt>1)
{
  del=del+0.2;
}
del=Foam::min(del,200);

eta0=0.0001;
eta1=1.0;
y0=fun(gamma_densityfilter,del,eta0,nallcells);
y1=fun(gamma_densityfilter,del,eta1,nallcells);
for(i=1;i<17;i++)//solve eta5
{
    eta5=(eta0+eta1)/2.0;
    y5=fun(gamma_densityfilter,del,eta5,nallcells);
    if(y0*y5<0)
    { 
       eta1=eta5;
       y1=y5;
    }
    else
    { 
       eta0=eta5;
       y0=y5;
    }
}

for(i=0;i<nallcells;i++)
{
   if(gamma_densityfilter[i]<=eta5)
   {
     Allgamma[i]=eta5*(Foam::exp(-del*(1.0-gamma_densityfilter[i]/eta5))-(1.0-gamma_densityfilter[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     Allgamma[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(gamma_densityfilter[i]-eta5)/(1.0-eta5))+(gamma_densityfilter[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}

for(i=0;i<allcells;i++)
{
   gamma[i]=Allgamma[i+mycells];
}
if(solid_area==1)
{
  setCells(gamma, cells_solid, 0);
}
if(fluid_area==1)
{
  setCells(gamma, cells_fluid, 1);
}
gamma.correctBoundaryConditions();
