
Info<<"sensitivity analysis"<<endl;
//计算灵敏度
fsenshMeanT=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoc;
gsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);
mu=-(1+qd)*qd/((qd+xh)*(qd+xh))*(E-Emin)/(2.0*(1.0+Po));
lambda=-(1+qd)*qd/((qd+xh)*(qd+xh))*Po*(E-Emin)/((1.0+Po)*(1.0-2.0*Po));
gsenshCompliance=-gradD && (rho*(mu*twoSymm(gradD)+ (lambda*I)*tr(gradD)));

//链式求导计算对物理密度的灵敏度
#include "filter_chainrule.H"

//往MMA中传入数据
gx[0]=(PowerDiss/Foam::max(PowerDissRelax-opt*0.2,PowerDissMax)-1.0);
gx[1]=(C/Foam::max(CRelax-opt*0.2,CMax)-1.0);
gx[2]=V;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
 dgdx[0][i]=gsensPowerDiss[i]/N;//sensitivities of constraint function
 dgdx[1][i]=gsensCompliance[i]/N;//sensitivities of constraint function
 dgdx[2][i]=gsensVol[i]/N;//sensitivities of constraint function
}

Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}
//PDE过滤和Heaviside过滤计算物理密度
#include "filter_x.H"
