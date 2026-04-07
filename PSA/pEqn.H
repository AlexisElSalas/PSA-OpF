volScalarField rAU = 1.0 / UEqn.A();
volVectorField HbyA("HbyA", U);
HbyA = rAU * UEqn.H();

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho) * fvc::flux(HbyA));

while (pimple.correctNonOrthogonal())
{
 fvScalarMatrix pEqn
(
    fvm::ddt(epsil * psi, p) // Transitorio con compresibilidad
  + fvc::div(phiHbyA)
  - fvm::laplacian(rho * rAU, p)
 ==
    S_rho // El sumidero total extrae fluido y, por ende, reduce la presión localmente
);
    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        // phi = phiHbyA - pEqn.flux();
		phi = (phiHbyA - pEqn.flux()) / fvc::interpolate(rho);
    }
}

U = HbyA - rAU * fvc::grad(p);
U.correctBoundaryConditions();

rho = thermo.rho();
rho.relax();