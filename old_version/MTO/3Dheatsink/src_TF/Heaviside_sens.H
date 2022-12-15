
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
// define scalar field GSENS_U
List<scalarField> AllguList(Pstream::nProcs());
scalarField localgusens(gsens_U.internalField());
AllguList[Pstream::myProcNo()] = localgusens;
Pstream::gatherList(AllguList);
Pstream::scatterList(AllguList);
scalarField Allgusens 
(
    ListListOps::combine<scalarField>
    (
        AllguList,
        accessOp<scalarField>()
    )
);  

scalef=0;
scaleg=0;    
for(i=0;i<allcells;i++)
{
   scalef=std::max(scalef,mag(fsens[i]));
   scaleg=std::max(scaleg,mag(gsens_U[i]));
}
reduce(scalef, maxOp<scalar>()); 
reduce(scaleg, maxOp<scalar>()); 

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
     Allgusens[i]=Allgusens[i]*drho[i];
     gsens[i]=1.0/nallcells*drho[i];
  }
}

for(i=0;i<nallcells;i++)//ith node
     {
        add1=0;
        Faa[i]=0;
        Gbb[i]=0;
        Gdd[i]=0;
        for(j=0;j<Nnum[i];j++)//j nodes around the ith node
        {
           add2=0;
           for(k=0;k<Nnum[Nlabel[i][j]];k++)//k nodes around the jth node
           {
              add2=add2+Nweight[Nlabel[i][j]][k];
           }
           add1=Nweight[i][j];
           Faa[i]=Faa[i]+Allfsens[Nlabel[i][j]]*add1/add2;
           Gbb[i]=Gbb[i]+Allgusens[Nlabel[i][j]]*add1/add2;
           Gdd[i]=Gdd[i]+gsens[Nlabel[i][j]]*add1/add2;
         }
      }
     if(opt>1)
     {
        for(i=0;i<allcells;i++)
        {
           fsens[i]=Faa[i+mycells];
           gsens_U[i]=Gbb[i+mycells];
           gsens[i]=Gdd[i+mycells];
        }          
     }