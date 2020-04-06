Info<< "\nslove pri-heat\n" << endl;
      fvScalarMatrix TEqn
      (
          fvm::div(phi,T)
         -fvm::laplacian(DT, T)
       ==
          fvOptions(T)+Q
      );
      TEqn.relax();
      fvOptions.constrain(TEqn);
      TEqn.solve();
      fvOptions.correct(T);

