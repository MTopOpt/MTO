
Info<<"sensitivity analysis"<<endl;

gsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);

#include "filter_chainrule.H"

gx[0]=V;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=gsensPowerDiss[i]/N;//sensitivities of objective function
 dgdx[0][i]=gsensVol[i]/N;//sensitivities of constraint function
}


Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}

#include "filter_x.H"
