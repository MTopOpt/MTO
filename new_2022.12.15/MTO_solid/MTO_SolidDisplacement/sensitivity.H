
Info<<"sensitivity analysis"<<endl;
fsensh=-gradD&&(rho*(3*xh*xh*(E-Emin)/(2.0*(1.0 + Po))*twoSymm(gradD) + (3*xh*xh*Po*(E-Emin)/((1.0 + Po)*(1.0 - 2.0*Po))*I)*tr(gradD)));  

#include "filter_chainrule.H"

gx[0]=V;//constraint function


for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsens[i]/N;//sensitivities of objective function
 dgdx[0][i]=gsens[i]/N;//sensitivities of constraint function
}


Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}

#include "filter_x.H"
