
        Info<< "\nslove primal-equation\n" << endl;
        nCorr =1;
        convergenceTolerance = 1e-6;
        int iCorr = 0;
        scalar initialResidual = 0;
        
        do
        {
            {
                fvVectorMatrix DEqn
                (
                    fvm::d2dt2(D)
                 ==
                    fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
                  + divSigmaExp
                );

                initialResidual = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                volTensorField gradD(fvc::grad(D));
                sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

                if (compactNormalStress)
                {
                    divSigmaExp = fvc::div
                    (
                        sigmaD - (2*mu + lambda)*gradD,
                        "div(sigmaD)"
                    );
                }
                else
                {
                    divSigmaExp += fvc::div(sigmaD);
                }
            }

        } while (initialResidual > convergenceTolerance&& ++iCorr < nCorr );

