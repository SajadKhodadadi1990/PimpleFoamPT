/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::lambda
(
    const vector& from,
    const vector& to,
    const label facei,
    const scalar stepFraction
) const
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    bool movingMesh = mesh.moving();

    if (movingMesh)
    {
        vector Sf = mesh.faceAreas()[facei];
        Sf /= mag(Sf);
        vector Cf = mesh.faceCentres()[facei];

        // move reference point for wall
        if
        (
            !cloud_.internalFace(facei)
         && isA<wallPolyPatch>(patches[cloud_.facePatch(facei)])
        )
        {
            const vector& C = mesh.cellCentres()[celli_];
            scalar CCf = mag((C - Cf) & Sf);
            // check if distance between cell centre and face centre
            // is larger than wallImpactDistance
            if
            (
                CCf
              > static_cast<const ParticleType&>(*this).wallImpactDistance(Sf)
            )
            {
                Cf -= static_cast<const ParticleType&>(*this)
                    .wallImpactDistance(Sf)*Sf;
            }
        }

        // for a moving mesh we need to reconstruct the old
        // Sf and Cf from oldPoints (they aren't stored)

        const vectorField& oldPoints = mesh.oldPoints();

        vector Cf00 = mesh.faces()[facei].centre(oldPoints);
        vector Cf0 = Cf00 + stepFraction*(Cf - Cf00);

        vector Sf00 = mesh.faces()[facei].normal(oldPoints);

        // for layer addition Sf00 = vector::zero and we use Sf
        if (mag(Sf00) > SMALL)
        {
            Sf00 /= mag(Sf00);
        }
        else
        {
            Sf00 = Sf;
        }

        scalar magSfDiff = mag(Sf - Sf00);

        // check if the face is rotating
        if (magSfDiff > SMALL)
        {
            vector Sf0 = Sf00 + stepFraction*(Sf - Sf00);

            // find center of rotation
            vector omega = Sf0 ^ Sf;
            scalar magOmega = mag(omega);
            omega /= magOmega + SMALL;
            vector n0 = omega ^ Sf0;
            scalar lam = ((Cf - Cf0) & Sf)/(n0 & Sf);
            vector r0 = Cf0 + lam*n0;

            // solve '(p - r0) & Sfp = 0', where
            // p = from + lambda*(to - from)
            // Sfp = Sf0 + lambda*(Sf - Sf0)
            // which results in the quadratic eq.
            // a*lambda^2 + b*lambda + c = 0
            vector alpha = from - r0;
            vector beta = to - from;
            scalar a = beta & (Sf - Sf0);
            scalar b = (alpha & (Sf - Sf0)) + (beta & Sf0);
            scalar c = alpha & Sf0;

            if (mag(a) > SMALL)
            {
                // solve the second order polynomial
                scalar ap = b/a;
                scalar bp = c/a;
                scalar cp = ap*ap - 4.0*bp;

                if (cp < 0.0)
                {
                    // imaginary roots only
                    return GREAT;
                }
                else
                {
                    scalar l1 = -0.5*(ap - ::sqrt(cp));
                    scalar l2 = -0.5*(ap + ::sqrt(cp));

                    // one root is around 0-1, while
                    // the other is very large in mag
                    if (mag(l1) < mag(l2))
                    {
                        return l1;
                    }
                    else
                    {
                        return l2;
                    }
                }
            }
            else
            {
                // when a==0, solve the first order polynomial
                return (-c/b);
            }
        }
        else // no rotation
        {
            vector alpha = from - Cf0;
            vector beta = to - from - (Cf - Cf0);
            scalar lambdaNominator = alpha & Sf;
            scalar lambdaDenominator = beta & Sf;

            // check if trajectory is parallel to face
            if (mag(lambdaDenominator) < SMALL)
            {
                if (lambdaDenominator < 0.0)
                {
                    lambdaDenominator = -SMALL;
                }
                else
                {
                    lambdaDenominator = SMALL;
                }
            }

            return (-lambdaNominator/lambdaDenominator);
        }
    }
    else
    {
        // mesh is static and stepFraction is not needed
        return lambda(from, to, facei);
    }
}


template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::lambda
(
    const vector& from,//C
    const vector& to,//Position
    const label facei //faceI
) const
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //Info<<"patches:"<<patches<<endl;
    
    vector Sf = mesh.faceAreas()[facei];
//Info<<"Sf faceI["<<facei<<"]:"<<Sf<<endl;
    Sf /= mag(Sf);
//Info<<"Sf faceI["<<facei<<"]:"<<Sf<<endl;
    vector Cf = mesh.faceCentres()[facei];
//Info<<"Cf faceI["<<facei<<"]:"<<Cf<<endl;
//Info<<"cloud_.internalFace["<<facei<<"]:"<<cloud_.internalFace(facei)<<endl;
//Info<<"cloud_.write["<<facei<<"]:"<<isA<wallPolyPatch>(patches[cloud_.facePatch(facei)])<<endl;
    // move reference point for wall
    if
    (
        !cloud_.internalFace(facei)
     && isA<wallPolyPatch>(patches[cloud_.facePatch(facei)])
    )
    {
	//Info<<"hello["<<facei<<"]:"<<endl;
        const vector& C = mesh.cellCentres()[celli_];
        scalar CCf = mag((C - Cf) & Sf);
        // check if distance between cell centre and face centre
        // is larger than wallImpactDistance
        if
        (
            CCf
          > static_cast<const ParticleType&>(*this).wallImpactDistance(Sf)
        )
        {
            Cf -= static_cast<const ParticleType&>(*this)
                .wallImpactDistance(Sf)*Sf;
        }
    }
	//Info<<"from["<<facei<<"]:"<<from<<endl;
    scalar lambdaNominator = (Cf - from) & Sf;
	//Info<<"lambdaNominator["<<facei<<"]:"<<lambdaNominator<<endl;
    scalar lambdaDenominator = (to - from) & Sf;
	//Info<<"lambdaDenominator["<<facei<<"]:"<<lambdaDenominator<<endl;
    // check if trajectory is parallel to face
    if (mag(lambdaDenominator) < SMALL)
    {
        if (lambdaDenominator < 0.0)
        {
            lambdaDenominator = -SMALL;
        }
        else
        {
            lambdaDenominator = SMALL;
        }
    }

    return lambdaNominator/lambdaDenominator;
}


template<class ParticleType>
inline bool Foam::Particle<ParticleType>::inCell() const
{
    dynamicLabelList& faces = cloud_.labels_;
    findFaces(position_, faces);

    return (!faces.size());
}


template<class ParticleType>
inline bool Foam::Particle<ParticleType>::inCell
(
    const vector& position,
    const label celli,
    const scalar stepFraction
) const
{
    dynamicLabelList& faces = cloud_.labels_;
    findFaces(position, celli, stepFraction, faces);

    return (!faces.size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
inline Foam::Particle<ParticleType>::trackData::trackData
(
    Cloud<ParticleType>& cloud
)
:
    cloud_(cloud)
{}


template<class ParticleType>
inline Foam::Cloud<ParticleType>&
Foam::Particle<ParticleType>::trackData::cloud()
{
    return cloud_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
inline const Foam::Cloud<ParticleType>&
Foam::Particle<ParticleType>::cloud() const
{
    return cloud_;
}


template<class ParticleType>
inline const Foam::vector& Foam::Particle<ParticleType>::position() const
{
    return position_;
}


template<class ParticleType>
inline Foam::vector& Foam::Particle<ParticleType>::position()
{
    return position_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::cell() const
{
    return celli_;
}


template<class ParticleType>
inline Foam::label& Foam::Particle<ParticleType>::cell()
{
    return celli_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::face() const
{
    return facei_;
}


template<class ParticleType>
inline bool Foam::Particle<ParticleType>::onBoundary() const
{
    return facei_ != -1 && facei_ >= cloud_.pMesh().nInternalFaces();
}


template<class ParticleType>
inline Foam::scalar& Foam::Particle<ParticleType>::stepFraction()
{
    return stepFraction_;
}


template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::stepFraction() const
{
    return stepFraction_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::origProc() const
{
    return origProc_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::origId() const
{
    return origId_;
}


template<class ParticleType>
inline bool Foam::Particle<ParticleType>::softImpact() const
{
    return false;
}


template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::currentTime() const
{
    return
        cloud_.pMesh().time().value()
      + stepFraction_*cloud_.pMesh().time().deltaT().value();
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::patch(const label facei) const
{
    return cloud_.facePatch(facei);
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::patchFace
(
    const label patchi,
    const label facei
) const
{
    return cloud_.patchFace(patchi, facei);
}


template<class ParticleType>
inline Foam::scalar
Foam::Particle<ParticleType>::wallImpactDistance(const vector&) const
{
    return 0.0;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::faceInterpolation() const
{
    return facei_;
}


// ************************************************************************* //
