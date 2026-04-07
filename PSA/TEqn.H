// TEqn.H

// Obtenemos los índices internos de las especies de la mezcla térmica
label idCO2 = thermo.composition().species()["CO2"];
label idN2  = thermo.composition().species()["N2"];

surfaceScalarField phiCp
(
    "phiCp",
    fvc::interpolate(thermo.rho() * thermo.Cp()) * phi
);

// 1. Poblar los campos de Cp celda por celda
forAll(T, cellI)
{
    Cp_CO2[cellI] = thermo.composition().Cp(idCO2, p[cellI], T[cellI]);
    Cp_N2[cellI]  = thermo.composition().Cp(idN2,  p[cellI], T[cellI]);
}

// 2. Cálculo de propiedades efectivas (Inercia térmica)
// Matriz base: fase gas + fase sólida pura
rhoCpEff = epsil * thermo.rho() * thermo.Cp() + (1.0 - epsil) * rhoSolid * CpSolid;

// Sumar la contribución de masa térmica de la fase adsorbida (cambia en el tiempo)
rhoCpEff += (1.0 - epsil) * rhoSolid * (qCO2 * MW_CO2 * Cp_CO2 + qN2 * MW_N2 * Cp_N2);

// Conductividad térmica efectiva (Eq. 23)
kappaEff = epsil * thermo.kappa() + (1.0 - epsil) * kappaSolid;

// 3. TÉRMINO FUENTE DE CALOR (Entalpía de adsorción)
// Como deltaH está como valor positivo en el diccionario, y la adsorción 
// (dqdt > 0) genera calor, el término se suma directamente.
Sh_T = (1.0 - epsil) * rhoSolid * (dqdt_CO2 * MW_CO2 * deltaH_CO2 + dqdt_N2 * MW_N2 * deltaH_N2);

// 4. Preparar el flujo convectivo de calor
// Se interpola el Cp a las caras de los volúmenes finitos y se multiplica por el flujo másico (phi)
//surfaceScalarField phiCp("phiCp", fvc::interpolate(thermo.Cp()) * phi);

// 5. Matriz de la ecuación de energía (Resuelve Eq. 22)
fvScalarMatrix TEqn
(
    fvm::ddt(rhoCpEff, T)
  + fvm::div(phiCp, T)
  - fvm::laplacian(kappaEff, T)
 ==
    Sh_T
);

// Relajar la matriz para mayor estabilidad y resolver
TEqn.relax();
TEqn.solve();

// 6. Actualizar termodinámica de la mezcla basándose en el nuevo T
thermo.correct();