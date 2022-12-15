Info<<"cost function"<<endl;
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));  

PowerDiss=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDiss=PowerDiss/PowerDiss0;
if(GeoDim==2)
{
 PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>()); 

if (Pstream::master())
{ 
   ofstream outfile1("Voluse.txt",std::ios::app);
   outfile1<<V<<"\n";
   outfile1.close();

   ofstream outfile2("PowerDiss.txt",std::ios::app);
   outfile2<<PowerDiss<<"\n";
   outfile2.close(); 

   ofstream outfile4("Time.txt",std::ios::app);
   outfile4<<runTime.elapsedClockTime()<<"\n";
   outfile4.close();    

}
if(runTime.writeTime())
{
   U.write();
}