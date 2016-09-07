    Info<< "Reading field C\n" << endl;

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Cs
    (
        IOobject
        (
            "Cs",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
       mesh
    );

    volScalarField T
    (
        IOobject
        (
            "T",
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

    Info<< "Reading field j\n" << endl;

    Info<< "Reading DCdict\n" << endl;

    IOdictionary DCdict
    (
        IOobject
        (
            "DCdict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar DC
    (
        DCdict.lookup("DC")
    );

    Info<< "Reading A B d Ce\n" << endl;

    dimensionedScalar A
    (
        DCdict.lookup("A")
    );

    dimensionedScalar B
    (
        DCdict.lookup("B")
    );

    dimensionedScalar d
    (
        DCdict.lookup("d")
    );

    dimensionedScalar kb
    (
        DCdict.lookup("kb")
    );
    dimensionedScalar rho
    (
        DCdict.lookup("rho")
    );
    dimensionedScalar nu
    (
        DCdict.lookup("nu")
    );

#   include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );