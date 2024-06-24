

fvVectorMatrix UEqn
(
    fvm::ddt(alpha, U) + fvm::div(alphaPhi, U)
  - fvm::Sp(fvc::ddt(alpha) + fvc::div(alphaPhi), U)
  + turbulence->divDevReff()
 ==
    (dimSet/rho)*cloudSU
);

UEqn.relax();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

surfaceScalarField phiForces
(

    (fvc::interpolate(rAU*cloudVolSUSu/rho)& mesh.Sf())*dimSet + rAUf*(g & mesh.Sf())
);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            phiForces/rAUf - fvc::snGrad(p)*mesh.magSf()
        )
    );
}
