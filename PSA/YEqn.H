// YEqn.H
{
    // 0. Crear el flujo másico advectivo interpolando la densidad a las caras
    // Esto asegura dimensiones [1 0 -1 0 0 0 0] (kg/s)
    surfaceScalarField rhoPhi = fvc::interpolate(rho) * phi;

    // 1. Ecuación de conservación de especies SÓLO para CO2
    fvScalarMatrix Y_CO2Eqn
    (
        fvm::ddt(epsil*rho, YCO2)
      + fvm::div(rhoPhi, YCO2) // <--- USAR rhoPhi AQUÍ
      - fvm::laplacian(epsil * rho * D_L, YCO2)
     ==
        S_Y_CO2
    );
    
    Y_CO2Eqn.relax();
    Y_CO2Eqn.solve();

    // 2. Imponer límites físicos (Sintaxis universal en reemplazo de .clamp())
    YCO2 = max(0.0, min(1.0, YCO2));

    // 3. Conservación Estricta: Calcular el gas portador algebraicamente
    YN2 = 1.0 - YCO2;
    YN2 = max(0.0, min(1.0, YN2));
}