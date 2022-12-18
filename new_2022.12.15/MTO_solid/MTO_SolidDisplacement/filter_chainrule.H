
fsensh=fsensh/gMax(mag(fsensh.primitiveField()));
gsensh.primitiveFieldRef() = 1;
//***********************************//Heaviside_filter
if(opt>1)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-xp[i]/eta5))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+Foam::exp(-del); 
     } 
  }
  fsensh=fsensh*drho;
  gsensh=gsensh*drho;
}

//***********************************//PDE filter
solve(fvm::laplacian(fsens)-fvm::Sp(b,fsens)+fsensh*b);
solve(fvm::laplacian(gsens)-fvm::Sp(b,gsens)+gsensh*b);
fsens.primitiveFieldRef()=fsens.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsens.primitiveFieldRef()=gsens.primitiveFieldRef()*mesh.V()/gMax(mesh.V());