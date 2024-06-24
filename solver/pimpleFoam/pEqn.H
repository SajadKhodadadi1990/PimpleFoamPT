{
    volVectorField HbyA(rAU*UEqn.H());
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           (fvc::interpolate(U) & mesh.Sf())
         + alphaf*fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    if (p.needReference())
    {
        adjustPhi(phiHbyA, U, p);
    }

    phiHbyA += phiForces;

    // Update the pressure BCs to ensure flux consistency
   // constrainPressure(p, U, phiHbyA, rAUf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphaf*rAUf, p)
         ==
            fvc::ddt(alpha) + fvc::div(alphaf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
        pEqn.solve(
            		mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
		  );
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux()/alphaf;

            p.relax();

            U = HbyA
              + rAU*fvc::reconstruct((phiForces - pEqn.flux()/alphaf)/rAUf);
            U.correctBoundaryConditions();
        }
    }
}

#include "continuityErrs.H"
