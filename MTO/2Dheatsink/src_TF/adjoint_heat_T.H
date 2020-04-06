Info<< "\nslove adj-heat-b\n" << endl;
     fvScalarMatrix TEqnb
     (
        -fvm::div(-phi,Tb) 
        +fvm::laplacian(DT, Tb)
      ==
         fvOptions(Tb)+Q
     );
     TEqnb.relax();
     fvOptions.constrain(TEqnb);
     TEqnb.solve();
     fvOptions.correct(Tb);
