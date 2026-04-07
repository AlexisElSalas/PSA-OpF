// 1. INICIALIZACIÓN TERMODINÁMICA Y DE FLUIDO BASE
Info<< "Reading thermophysical properties\n" << endl;
//autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
//fluidThermo& thermo = pThermo();
autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
fluidMulticomponentThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

// Referencias a las variables termodinámicas de OpenFOAM
volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
// const volScalarField& mu = thermo.mu();

Info<< "Reading field rho (density)\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    thermo.rho()
);

Info<< "Reading field U (velocity)\n" << endl;
volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// Crea el campo de flujos en las caras (phi) vital para el operador fvm::div
#include "createPhi.H"


// 2. FRACCIONES MÁSICAS (ESPECIES)
// Info<< "Reading species fields YCO2 and YN2\n" << endl;
// volScalarField YCO2
// (
//     IOobject("YCO2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
//     mesh
// );

// volScalarField YN2
// (
//     IOobject("YN2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
//     mesh
// );

Info<< "Extracting species fields CO2 and N2 from thermo\n" << endl;

//basicMultiComponentMixture& composition = thermo.composition();
//volScalarField& YCO2 = thermo.composition().Y("CO2");
//volScalarField& YN2 = thermo.composition().Y("N2");

//basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = thermo.composition().Y();

// Ahora Y[0] será CO2 y Y[1] será N2 (dependiendo del orden en thermophysicalProperties).
// Puedes obtener una referencia directa a la especie por nombre:
volScalarField& YCO2 = Y[0];
volScalarField& YN2 = Y[1];


// 3. DICCIONARIO DE PROPIEDADES DE ADSORCIÓN (HENRY 2024)
Info<< "Reading Adsorption Properties\n" << endl;
IOdictionary adsorptionProperties
(
    IOobject
    (
        "adsorptionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Parámetros de la cama porosa
dimensionedScalar epsil("epsil", dimless, adsorptionProperties);
dimensionedScalar dp("dp", dimLength, adsorptionProperties);
dimensionedScalar rhoSolid("rhoSolid", dimDensity, adsorptionProperties);
dimensionedScalar kap("kappaSolid", dimPower/dimLength/dimTemperature, adsorptionProperties);
dimensionedScalar CpSolid("CpSolid", dimEnergy/dimMass/dimTemperature, adsorptionProperties);
dimensionedScalar D_L("D_L", dimArea/dimTime, adsorptionProperties);

// Parámetros térmicos y cinéticos
dimensionedScalar deltaH_CO2("deltaH_CO2", dimEnergy/dimMass, adsorptionProperties);
dimensionedScalar deltaH_N2("deltaH_N2", dimEnergy/dimMass, adsorptionProperties);
dimensionedScalar LDF_CO2("LDF_CO2", dimless/dimTime, adsorptionProperties); // [1/s]
dimensionedScalar LDF_N2("LDF_N2", dimless/dimTime, adsorptionProperties); // [1/s]

dimensionedScalar MW_CO2("MW_CO2", dimMass/dimMoles, adsorptionProperties);
dimensionedScalar MW_N2("MW_N2", dimMass/dimMoles, adsorptionProperties);

//dimensionedScalar qm_CO2("qm_CO2", dimless, adsorptionProperties); // Capacidad max
//dimensionedScalar b_CO2("b_CO2", dimless/dimPressure, adsorptionProperties); // Afinidad
//dimensionedScalar qm_N2("qm_N2", dimless, adsorptionProperties);
//dimensionedScalar b_N2("b_N2", dimless/dimPressure, adsorptionProperties);

// dimSet(Mass, Length, Time, Temperature, Moles, Current, LuminousIntensity)
dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), adsorptionProperties);

// ---- Parámetros de Isoterma DSL para CO2 ----
dimensionedScalar qsb_CO2("qsb_CO2", dimensionSet(-1, 0, 0, 0, 1, 0, 0), adsorptionProperties);
dimensionedScalar b0_CO2("b0_CO2", dimensionSet(0, 3, 0, 0, -1, 0, 0), adsorptionProperties);
dimensionedScalar dU_b_CO2("dU_b_CO2", dimensionSet(1, 2, -2, 0, -1, 0, 0), adsorptionProperties);

dimensionedScalar qsd_CO2("qsd_CO2", dimensionSet(-1, 0, 0, 0, 1, 0, 0), adsorptionProperties);
dimensionedScalar d0_CO2("d0_CO2", dimensionSet(0, 3, 0, 0, -1, 0, 0), adsorptionProperties);
dimensionedScalar dU_d_CO2("dU_d_CO2", dimensionSet(1, 2, -2, 0, -1, 0, 0), adsorptionProperties);

// ---- Parámetros de Isoterma DSL para N2 ----
dimensionedScalar qsb_N2("qsb_N2", dimensionSet(-1, 0, 0, 0, 1, 0, 0), adsorptionProperties);
dimensionedScalar b0_N2("b0_N2", dimensionSet(0, 3, 0, 0, -1, 0, 0), adsorptionProperties);
dimensionedScalar dU_b_N2("dU_b_N2", dimensionSet(1, 2, -2, 0, -1, 0, 0), adsorptionProperties);

dimensionedScalar qsd_N2("qsd_N2", dimensionSet(-1, 0, 0, 0, 1, 0, 0), adsorptionProperties);
dimensionedScalar d0_N2("d0_N2", dimensionSet(0, 3, 0, 0, -1, 0, 0), adsorptionProperties);
dimensionedScalar dU_d_N2("dU_d_N2", dimensionSet(1, 2, -2, 0, -1, 0, 0), adsorptionProperties);
// 1. Capacidad máxima [mol/kg_solido]
//dimensionedScalar qm_CO2("qm_CO2", dimensionSet(-1, 0, 0, 0, 1, 0, 0), adsorptionProperties);
//dimensionedScalar qm_N2("qm_N2", dimensionSet(-1, 0, 0, 0, 1, 0, 0), adsorptionProperties);

// 2. Afinidad [1/Pa] 
//dimensionedScalar b_CO2("b_CO2", dimless/dimPressure, adsorptionProperties);
//dimensionedScalar b_N2("b_N2", dimless/dimPressure, adsorptionProperties);

// 4. CAMPOS DE FASE ADSORBIDA (Masa en el sólido)
Info<< "Reading adsorbed phase fields qCO2 and qN2\n" << endl;
volScalarField qCO2
(
    IOobject("qCO2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField qN2
(
    IOobject("qN2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// 5. VARIABLES INTERMEDIAS (Fuentes y Equilibrio) - No se leen de carpeta, se calculan
Info<< "Creating source term and equilibrium fields\n" << endl;

volScalarField qStar_CO2
(
    IOobject("qStar_CO2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, 0, 0, 0, 1, 0, 0), 0.0) // [mol/kg]
);

volScalarField qStar_N2
(
    IOobject("qStar_N2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, 0, 0, 0, 1, 0, 0), 0.0) // [mol/kg]
);

volScalarField dqdt_CO2
(
    IOobject("dqdt_CO2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, 0, -1, 0, 1, 0, 0), 0.0) // [mol/kg/s]
);

volScalarField dqdt_N2
(
    IOobject("dqdt_N2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, 0, -1, 0, 1, 0, 0), 0.0) // [mol/kg/s]
);

volScalarField S_Y_CO2
(
    IOobject("S_Y_CO2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0) // [kg/(m3 s)]
);

volScalarField S_Y_N2
(
    IOobject("S_Y_N2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
);

volScalarField S_rho
(
    IOobject("S_rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
);

// 6. CAMPOS AUXILIARES PARA CÁLCULOS (Solución a errores de memoria/constructores)
Info<< "Creating auxiliary fields for calculations\n" << endl;

volScalarField MW_mix(IOobject("MW_mix", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("MW_mix", dimMass/dimMoles, 0.0));
volScalarField y_CO2(IOobject("y_CO2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimless, 0.0));
volScalarField y_N2(IOobject("y_N2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimless, 0.0));
volScalarField P_CO2(IOobject("P_CO2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimPressure, 0.0));
volScalarField P_N2(IOobject("P_N2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimPressure, 0.0));
volScalarField denLangmuir(IOobject("denLangmuir", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("one", dimless, 1.0));
volScalarField kappaSolid(IOobject("kappaSolid", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimPower/dimLength/dimTemperature, 0.0));

volScalarField Cp_CO2(IOobject("Cp_CO2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimEnergy/dimMass/dimTemperature, 0.0));
volScalarField Cp_N2(IOobject("Cp_N2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimEnergy/dimMass/dimTemperature, 0.0));
volScalarField rhoCpEff(IOobject("rhoCpEff", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimEnergy/dimVolume/dimTemperature, 0.0));
volScalarField kappaEff(IOobject("kappaEff", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimPower/dimLength/dimTemperature, 0.0));
volScalarField Sh_T(IOobject("Sh_T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0));