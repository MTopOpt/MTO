Info<< "Reading optimization parameters\n" << endl;
IOdictionary optProperties
(
   IOobject
   (
       "optProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
//optimization parameters
scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));
scalar voluse(readScalar(optProperties.lookup("voluse")));
scalar filterR(readScalar(optProperties.lookup("filterR")));
scalar Hdel(readScalar(optProperties.lookup("Hdel")));
scalar volinit(readScalar(optProperties.lookup("volinit")));
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),0);
Switch solid_area(optProperties.lookup("solid_area"));

Info<< "Reading pseudo-density field\n" << endl;
volScalarField xp
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(volinit),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(volinit),
    zeroGradientFvPatchScalarField::typeName
);


volScalarField x(xh);
volScalarField drho(xh);

// find the cells in the non-design area
labelList cells;
if(solid_area)
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells=mesh.cellZones()[zoneID];
  forAll(cells, i)
  {
     x[cells[i]]=1.0;
  }
}
