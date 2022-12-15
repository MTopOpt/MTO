
 V=0;
 meanT=0;
 for (i=0;i<allcells;i++)
 {
     V=gamma[i]+V;
     meanT=meanT+T[i];
 }
 reduce(meanT, sumOp<scalar>());
 reduce(V, sumOp<scalar>());     
 V=V/nallcells-voluse;

 //DissPower
 DissPower=0;
 for (i=0; i<nObjPatch; i++)
 {
      DissPower =DissPower - sum(
      phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
      + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
 }
 DissPower=DissPower/D_normalization;
 if(geo_dim==2) 
 {
   DissPower=DissPower/len;
 }
 reduce(DissPower, sumOp<scalar>()); 
 Info<<"constraint function (DissPower): "<<DissPower<<endl; 

 if (Pstream::master())
{
ofstream outfile1("Voluse.txt",std::ios::app);
outfile1<<V<<"\n";
outfile1.close();

ofstream outfile2("meanT.txt",std::ios::app);
outfile2<<meanT/nallcells<<"\n";
outfile2.close(); 
   
ofstream outfile3("Disspower.txt",std::ios::app);
outfile3<<DissPower<<"\n";
outfile3.close();

ofstream outfile4("Time.txt",std::ios::app);
outfile4<<runTime.elapsedClockTime()<<"\n";
outfile4.close(); 

}