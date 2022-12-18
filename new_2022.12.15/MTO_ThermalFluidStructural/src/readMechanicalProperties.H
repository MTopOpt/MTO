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
int iCorr=0,nCorr=150;
scalar convergenceTolerance(1e-5),initialResidual(0);

Switch planeStress(mechanicalProperties.lookup("planeStress"));
dimensionedScalar rhoE("rhoE", dimMass/dimLength/sqr(dimTime),readScalar(mechanicalProperties.lookup("rhoE")));
dimensionedScalar rho("rho",  dimMass/dimVolume,readScalar(mechanicalProperties.lookup("rho")));
scalar Po(readScalar(mechanicalProperties.lookup("Po")));
dimensionedScalar E(rhoE/rho);
dimensionedScalar Emin(E*1e-9);
double qd=0.1;
volScalarField mu(qd*(1-xh)/(qd+xh)*(E-Emin)/(2.0*(1.0+Po))+Emin/(2.0*(1.0+Po)));
volScalarField lambda(qd*(1-xh)/(qd+xh)*Po*(E-Emin)/((1.0+Po)*(1.0-2.0*Po))+Po*Emin/((1.0+Po)*(1.0-2.0*Po)));

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
    mu*twoSymm(gradD) + (lambda*I)*tr(gradD)
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

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());

volScalarField gsensCompliance
(
    IOobject
    (
        "gsensCompliance",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -gradD&&(rho*(3*xh*xh*(E-Emin)/(2.0*(1.0 + Po))*twoSymm(gradD))),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshCompliance(gsensCompliance);


