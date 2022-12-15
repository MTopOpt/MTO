  gradD=fvc::grad(D);  
  volSymmTensorField sigma(rho*(mu*twoSymm(gradD) + (lambda*I)*tr(gradD)));  
  volScalarField Energy((gradD+gradD.T())&&sigma); 
  //volume use     
  V=0;     
  for(i=0;i<allcells;i++)
  {
     V=gamma[i]+V;
  }
  reduce(V, sumOp<scalar>()); 
  V=V/nallcells-voluse;
  //compliance    
  C=0;  
  for (i=0;i<allcells;i++)
  {
      C=C+Energy[i];
  }
  reduce(C, sumOp<scalar>());     
  C=C/2.0*len*len*len; 

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