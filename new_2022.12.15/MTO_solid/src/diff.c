
scalar diff(volScalarField gamma,scalarField V,double del,double eta,int n)
{
     int i;
     scalar z=0;
     double *x =new double[n];
     
     for(i=0;i<n;i++)
     {
        if(gamma[i]<=eta)
        {
          x[i]=eta*(Foam::exp(-del*(1-gamma[i]/eta))-(1-gamma[i]/eta)*Foam::exp(-del));
        }
        else
        {
          x[i]=eta+(1-eta)*(1-Foam::exp(-del*(gamma[i]-eta)/(1-eta))+(gamma[i]-eta)*Foam::exp(-del)/(1-eta));
        }
     }
     for(i=0;i<n;i++)
     {
        z=z+(gamma[i]-x[i])*V[i];
     }
     delete x;
     return {z};
}

