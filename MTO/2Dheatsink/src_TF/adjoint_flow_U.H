            Info<< "\nslove adj-flow-U\n" << endl;  
            volVectorField adjointTransposeConvectionc((fvc::grad(U) & Uc));

            tmp<fvVectorMatrix> tUcEqn
            (
                fvm::div(-phi, Uc)
              + adjointTransposeConvectionc
              - fvm::laplacian(nu,Uc)
              + fvm::Sp(alpha, Uc)
             ==
                fvOptions(Uc)
            );
            fvVectorMatrix& UcEqn = tUcEqn.ref();

            UcEqn.relax();

            fvOptions.constrain(UcEqn);

            solve(UcEqn == -fvc::grad(pc));

            fvOptions.correct(Uc);
//****************************************
            volScalarField rAUc(1.0/UcEqn.A());
            volVectorField HbyAc("HbyAc", Uc);
            HbyAc = rAUc*UcEqn.H();
            tUcEqn.clear();          
            surfaceScalarField phiHbyAc("phiHbyAc", fvc::flux(HbyAc));
            adjustPhi(phiHbyAc, Uc, pc);
            
            tmp<volScalarField> rAtUc(rAUc);
            
             if (simple.consistent())
            {
                rAtUc = 1.0/(1.0/rAUc - UcEqn.H1());
                phiHbyAc +=
                    fvc::interpolate(rAtUc() - rAUc)*fvc::snGrad(pc)*mesh.magSf();
                HbyAc -= (rAUc - rAtUc())*fvc::grad(pc);
            }

            tUcEqn.clear();

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(pc, Uc, phiHbyAc, rAtUc(), MRF);           
            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pcEqn
                (
                    fvm::laplacian(rAUc, pc) == fvc::div(phiHbyAc)
                );

                pcEqn.setReference(pcRefCell, pcRefValue);
                pcEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phic = phiHbyAc - pcEqn.flux();
                }
            }

            //#include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pc.relax();

            // Adjoint momentum corrector
            Uc = HbyAc - rAUc*fvc::grad(pc);
            Uc.correctBoundaryConditions();
            fvOptions.correct(Uc);
            
            Uc.storePrevIter();
            pc.storePrevIter();
            phic.storePrevIter();

