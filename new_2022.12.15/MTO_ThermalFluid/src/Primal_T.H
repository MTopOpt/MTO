T=T2;
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(DT, T)
);
if(DDHS)
{
 TEqn-=100000000.0/4180000*Q*(1-xh);
}
else
{
 TEqn-=100000000.0/4180000*Q;
}
TEqn.solve();

