// ========================================================================= //
// adsorptionEqn.H
// Resolución del Equilibrio Dual-Site Langmuir y Cinética LDF
// ========================================================================= //

dimensionedScalar dte = runTime.deltaT();


// 1. CÁLCULO DE FRACCIONES Y PRESIONES PARCIALES
MW_mix = 1.0 / ((YCO2 / MW_CO2) + (YN2 / MW_N2));
y_CO2 = YCO2 * (MW_mix / MW_CO2);
y_N2  = YN2  * (MW_mix / MW_N2);

P_CO2 = y_CO2 * p;
P_N2  = y_N2  * p;

// 2. CONVERSIÓN A CONCENTRACIÓN MOLAR (c_i en mol/m^3)
// Usando la ley de los gases ideales: c = P / (R * T)
volScalarField c_CO2 = P_CO2 / (R * T);
volScalarField c_N2  = P_N2  / (R * T);

// 3. DEPENDENCIAS TÉRMICAS DE ARRHENIUS PARA AFINIDADES
// Se calculan celda por celda basándose en la Temperatura (T)
volScalarField b_CO2 = b0_CO2 * exp(-dU_b_CO2 / (R * T));
volScalarField d_CO2 = d0_CO2 * exp(-dU_d_CO2 / (R * T));

volScalarField b_N2  = b0_N2  * exp(-dU_b_N2 / (R * T));
volScalarField d_N2  = d0_N2  * exp(-dU_d_N2 / (R * T));

// 4. MODELO DUAL-SITE LANGMUIR
// Denominadores competitivos para los sitios 'b' (fuertes) y 'd' (débiles)
volScalarField den_b = 1.0 + (b_CO2 * c_CO2) + (b_N2 * c_N2);
volScalarField den_d = 1.0 + (d_CO2 * c_CO2) + (d_N2 * c_N2);

// Capacidad de equilibrio (q*) usando Ecuación 27 del paper
qStar_CO2 = (qsb_CO2 * b_CO2 * c_CO2) / den_b 
          + (qsd_CO2 * d_CO2 * c_CO2) / den_d;

qStar_N2  = (qsb_N2 * b_N2 * c_N2) / den_b 
          + (qsd_N2 * d_N2 * c_N2) / den_d;
	
// 5. RESOLUCIÓN DE CINÉTICA LDF (Linear Driving Force) mediante Euler Implícito
qCO2 = (qCO2.oldTime() + dte * LDF_CO2 * qStar_CO2) / (1.0 + dte * LDF_CO2);
qN2  = (qN2.oldTime()  + dte * LDF_N2  * qStar_N2)  / (1.0 + dte * LDF_N2);

// 6. DERIVADAS TEMPORALES DE ADSORCIÓN [mol/(kg*s)]
dqdt_CO2 = (qCO2 - qCO2.oldTime()) / dte;
dqdt_N2  = (qN2  - qN2.oldTime())  / dte;

// 7. TÉRMINOS FUENTE PARA TRANSPORTE DE ESPECIES [kg/(m3*s)]
// Sumidero negativo porque la fase fluida pierde masa
S_Y_CO2 = - (1.0 - epsil) * rhoSolid * dqdt_CO2 * MW_CO2;
S_Y_N2  = - (1.0 - epsil) * rhoSolid * dqdt_N2  * MW_N2;

// 8. TÉRMINO FUENTE PARA ECUACIÓN DE CONTINUIDAD (Masa total)
S_rho = S_Y_CO2 + S_Y_N2;