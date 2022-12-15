Info<<"cost function"<<endl;

gradD=fvc::grad(D);  
sigmaD=mu*twoSymm(gradD) + (lambda*I)*tr(gradD);  
volScalarField Energy(0.25*(gradD+gradD.T())&&sigmaD); 

C=fvc::domainIntegrate(Energy).value();
V=fvc::domainIntegrate(xh).value()/area-voluse; 

if (Pstream::master())
{ 
   ofstream outfile1("Voluse.txt",std::ios::app);
   outfile1<<V<<"\n";
   outfile1.close();

   ofstream outfile2("Compliance.txt",std::ios::app);
   outfile2<<C<<"\n";
   outfile2.close(); 

   ofstream outfile3("Time.txt",std::ios::app);
   outfile3<<runTime.elapsedClockTime()<<"\n";
   outfile3.close();    
}
if(runTime.writeTime())
{
  // D.write();     
}