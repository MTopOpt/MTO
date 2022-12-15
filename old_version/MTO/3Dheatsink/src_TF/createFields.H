//******************************************************************//pri-fluid
#include "createMRF.H"
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "createPhi.H"
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
//************************************************************************//
Info<< "Reading field pb\n" << endl;
volScalarField pb
(
    IOobject
    (
        "pb",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field Ub\n" << endl;
volVectorField Ub
(
    IOobject
    (
        "Ub",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading/calculating face flux field phib\n" << endl;

surfaceScalarField phib
(
    IOobject
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ub) & mesh.Sf()
);

label pbRefCell = 0;
scalar pbRefValue = 0.0;
setRefCell
(
    pb,
    simple.dict(),
    pbRefCell,
    pbRefValue
);
mesh.setFluxRequired(pb.name());
//************************************************************************//
Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field Uc\n" << endl;
volVectorField Uc
(
    IOobject
    (
        "Uc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading/calculating face flux field phic\n" << endl;

surfaceScalarField phic
(
    IOobject
    (
        "phic",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);
label pcRefCell = 0;
scalar pcRefValue = 0.0;
setRefCell
(
    pc,
    simple.dict(),
    pcRefCell,
    pcRefValue
);
mesh.setFluxRequired(pc.name());
//************************************************************************//property
Info<< "Reading Material Properties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


//fluid
dimensionedScalar nu(transportProperties.lookup("nu"));//Viscosity of fluid
dimensionedScalar alphaMax(transportProperties.lookup("alphaMax"));//Flow resistance
dimensionedScalar alphamax(transportProperties.lookup("alphaMax"));//Flow resistance

//optimization
scalar raa0(readScalar(transportProperties.lookup("raa0")));
scalar mma_init(readScalar(transportProperties.lookup("mma_init")));
scalar mma_dec(readScalar(transportProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(transportProperties.lookup("mma_inc")));
scalar movlim(readScalar(transportProperties.lookup("movlim")));

scalar voluse(readScalar(transportProperties.lookup("voluse")));
scalar filter_Radius(readScalar(transportProperties.lookup("filter_Radius")));
scalar solid_area(readScalar(transportProperties.lookup("solid_area")));
scalar fluid_area(readScalar(transportProperties.lookup("fluid_area")));
scalar test_PD(readScalar(transportProperties.lookup("test_PD")));

//  constraint=DissPower/D_normalization/Foam::max(D0-opt*delta,D1)-1.0
scalar D_normalization(readScalar(transportProperties.lookup("D_normalization")));
scalar D0(readScalar(transportProperties.lookup("D0")));
scalar D1(readScalar(transportProperties.lookup("D1")));
scalar geo_dim(readScalar(transportProperties.lookup("geo_dim")));

dimensionedScalar b1(transportProperties.lookup("b1"));

volScalarField qu
(
    IOobject
    (
        "qu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    readScalar(transportProperties.lookup("qu")),
    zeroGradientFvPatchScalarField::typeName
);
//gamma
Info<< "Reading field gamma\n" << endl;
volScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
labelList cells1,cells_solid,cells_fluid;
if(test_PD==1)
{
  word zoneName1="zone_test";
  label zoneID1=mesh.cellZones().findZoneID(zoneName1);
  cells1=mesh.cellZones()[zoneID1];
  setCells(gamma, cells1, 1);
}
if(solid_area==1)
{
  word zoneName2="zone_solid";
  label zoneID2=mesh.cellZones().findZoneID(zoneName2);
  cells_solid=mesh.cellZones()[zoneID2];
  setCells(gamma, cells_solid, 0);
}
if(fluid_area==1)
{
  word zoneName3="zone_fluid";
  label zoneID3=mesh.cellZones().findZoneID(zoneName3);
  cells_fluid=mesh.cellZones()[zoneID3];
  setCells(gamma, cells_fluid, 1);
}
volScalarField alpha(alphaMax*qu*(1-gamma)/(qu+gamma)); 

//Allgamma
List<scalarField> AllrList(Pstream::nProcs());
scalarField localr(gamma.internalField());
AllrList[Pstream::myProcNo()] = localr;
Pstream::gatherList(AllrList);
Pstream::scatterList(AllrList);
scalarField Allgamma 
(
    ListListOps::combine<scalarField>
    (
        AllrList,
        accessOp<scalarField>()
    )
);

volScalarField fsens
    (
        IOobject
        (
            "fsens",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -alphaMax*(1+qu)*qu/((qu+gamma)*(qu+gamma))*(U&Ub)
    );
volScalarField gsens_U
    (
        IOobject
        (
            "gsens_U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        -alphaMax*(1+qu)*qu/((qu+gamma)*(qu+gamma))*(U&Uc)
    );  
       
//
dictionary Fun1=mesh.solutionDict().subDict("power_dissipation");
int nObjPatch=Fun1.lookupOrDefault<scalar>("numberConstraintPatches",2);//number of patch 
wordList conPatchNames = Fun1.lookup("constraintPatchesNames");//function name
label conPatchList[nObjPatch];//patch list
int iLoop;
for (iLoop=0; iLoop<nObjPatch; iLoop++)
{
    conPatchList[iLoop] = mesh.boundaryMesh().findPatchID(conPatchNames[iLoop]);
}   

