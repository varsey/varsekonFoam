   Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar DC
    (
        DCdict.lookup("DC")
    );
    dimensionedScalar Sc
    (
        DCdict.lookup("Sc")
    );
    dimensionedScalar Sct
    (
        DCdict.lookup("Sct")
    );

    Info<< "Reading constants: A B d\n" << endl;

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

    Info<< "Reading fields: rho nu\n" << endl;

    dimensionedScalar rho
    (
        DCdict.lookup("rho")
    );
    dimensionedScalar nu
    (
        DCdict.lookup("nu")
    );

