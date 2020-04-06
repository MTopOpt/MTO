    Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dictionary& rhoDict(mechanicalProperties.subDict("rho"));
    word rhoType(rhoDict.lookup("type"));

    autoPtr<volScalarField> rhoPtr;

    IOobject rhoIO
    (
        "rho",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (rhoType == "uniform")
    {
        scalar rhoValue(readScalar(rhoDict.lookup("value")));

        rhoPtr.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh,
                dimensionedScalar
                (
                    "rho",
                    dimMass/dimVolume,
                    rhoValue
                )
            )
        );
    }
    else if (rhoType == "field")
    {
        rhoIO.readOpt() = IOobject::MUST_READ;

        rhoPtr.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for rho"
            << abort(FatalError);
    }

    volScalarField& rho = rhoPtr();

    const dictionary& EDict(mechanicalProperties.subDict("E"));
    word EType(EDict.lookup("type"));

    autoPtr<volScalarField> EPtr;

    IOobject EHeader
    (
        "E",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (EType == "uniform")
    {
        scalar rhoEValue(readScalar(EDict.lookup("value")));

        EPtr.reset
        (
            new volScalarField
            (
                EHeader,
                mesh,
                dimensionedScalar
                (
                    "Erho",
                    dimMass/dimLength/sqr(dimTime),
                    rhoEValue
                )
            )
        );
    }
    else if (EType == "field")
    {
        EHeader.readOpt() = IOobject::MUST_READ;

        EPtr.reset
        (
            new volScalarField
            (
                EHeader,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for E"
            << abort(FatalError);
    }

    volScalarField& rhoE = EPtr();

    autoPtr<volScalarField> PPtr;

    IOobject PIO
    (
        "Po",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& PDict(mechanicalProperties.subDict("Po"));
    word PType(PDict.lookup("type"));

    if (PType == "uniform")
    {
        scalar PValue(readScalar(PDict.lookup("value")));
        PPtr.reset
        (
            new volScalarField
            (
                PIO,
                mesh,
                dimensionedScalar
                (
                    "Po",
                    dimless,
                    PValue
                )
            )
        );
    }
    else if (PType == "field")
    {
        PIO.readOpt() = IOobject::MUST_READ;
        PPtr.reset
        (
            new volScalarField
            (
                PIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for P"
            << abort(FatalError);
    }

    volScalarField& Po = PPtr();

    Info<< "Normalising E : E/rho\n" << endl;
    volScalarField E(rhoE/rho);

    Info<< "Calculating Lame's coefficients\n" << endl;

    volScalarField mu(gamma*gamma*gamma*E/(2.0*(1.0 + Po)));
    volScalarField lambda(gamma*gamma*gamma*Po*E/((1.0 + Po)*(1.0 - 2.0*Po)));

    Switch planeStress(mechanicalProperties.lookup("planeStress"));


//************************************************************************//pri-solid   
Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volVectorField Dorg(D);

volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(D)
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(fvc::grad(D)) + lambda*(I*tr(fvc::grad(D)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

Switch compactNormalStress("yes");
int nCorr = 1;
scalar convergenceTolerance(1e-6);

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());
