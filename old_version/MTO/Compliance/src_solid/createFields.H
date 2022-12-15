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

scalar filter_Radius(readScalar(optProperties.lookup("filter_Radius")));
scalar solid_area(readScalar(optProperties.lookup("solid_area")));
scalar voluse(readScalar(optProperties.lookup("voluse")));

//pseudo-density field gamma
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

// find the cells in the non-design area
labelList cells;
if(solid_area==1)
{
  word zoneName="zone_solid";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells=mesh.cellZones()[zoneID];
}
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

