// define scalar field Allfsens
List<scalarField> AllfList(Pstream::nProcs());
scalarField localfsens(fsens.internalField());
AllfList[Pstream::myProcNo()] = localfsens;
Pstream::gatherList(AllfList);
Pstream::scatterList(AllfList);
scalarField Allfsens 
(
    ListListOps::combine<scalarField>
    (
        AllfList,
        accessOp<scalarField>()
    )
);

scalef=0;      
for(i=0;i<nallcells;i++)
{
   gsens[i]=1.0/nallcells;
   scalef=std::max(scalef,mag(Allfsens[i]));
}

if(opt>1)
{
  for(i=0;i<nallcells;i++)
  {
     if(gamma_densityfilter[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-gamma_densityfilter[i]/eta5))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(gamma_densityfilter[i]-eta5)/(1-eta5))+Foam::exp(-del); 
     } 
     Allfsens[i]=Allfsens[i]*drho[i];
     gsens[i]=gsens[i]*drho[i];
  }
}

for(i=0;i<nallcells;i++)
{
   add1=0;
   Faa[i]=0;
   Gbb[i]=0;
   for(j=0;j<Nnum[i];j++)
   {
      add2=0;
      for(k=0;k<Nnum[Nlabel[i][j]];k++)
      {
         add2=add2+Nweight[Nlabel[i][j]][k];
      }
      add1=Nweight[i][j];
      Faa[i]=Faa[i]+Allfsens[Nlabel[i][j]]*add1/add2;
      Gbb[i]=Gbb[i]+gsens[Nlabel[i][j]]*add1/add2;
    }
 }

if(opt>1)
{
   for(i=0;i<allcells;i++)
   {
      fsens[i]=Faa[i+mycells];
      gsens[i]=Gbb[i+mycells];
   }          
}