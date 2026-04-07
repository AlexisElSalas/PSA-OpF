volScalarField mu = thermo.mu();


surfaceScalarField rhoPhi = fvc::interpolate(rho) * phi;


volScalarField D_Ergun = 150.0 * mu * pow(1.0 - epsil, 2) / (pow(dp, 2) * pow(epsil, 3));
volScalarField F_Ergun = 1.75 * rho * (1.0 - epsil) / (dp * pow(epsil, 3));

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mu, U)
  + fvm::Sp(D_Ergun + F_Ergun * mag(U), U)
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
