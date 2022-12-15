//Tb=T2;
sumT=0;
maxT=gMax(T);
for(i=0;i<n;i++)
{
   sumT=sumT+std::pow(T[i]/maxT,Pnorm);
}      
reduce(sumT, sumOp<scalar>());
sumT=std::pow(sumT,1.0/Pnorm-1);

fvScalarMatrix TEqnb
(
   -fvm::div(-phi,Tb) 
   +fvm::laplacian(DT, Tb)
);
if (objFunction==1)
{
   TEqnb -= Q;
}
if (objFunction==2)
{
   TEqnb -= b2*sumT*Foam::pow(T/maxT,Pnorm-1);
}
TEqnb.solve();
