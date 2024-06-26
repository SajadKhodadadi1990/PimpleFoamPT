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

#include "fvmSup.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::label Foam::KinematicCloud<ParcelType>::parcelTypeId() const
{
    return parcelTypeId_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::KinematicCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}

template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<CloudType>::outputProperties() const //add by sajad
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::KinematicCloud<CloudType>::outputProperties() //add by sajad
{
    return outputProperties_;
}

template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::KinematicCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::Switch Foam::KinematicCloud<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline const Foam::Switch Foam::KinematicCloud<ParcelType>::coupled() const
{
    return coupled_;
}


template <class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cellValueSourceCorrection() const
{
    return cellValueSourceCorrection_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::KinematicCloud<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::KinematicCloud<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::KinematicCloud<ParcelType>::mu() const
{
    return mu_;
}


template<class ParcelType>
inline const Foam::dimensionedVector&
Foam::KinematicCloud<ParcelType>::g() const
{
    return g_;
}


template<class ParcelType>
inline const Foam::particleForces&
Foam::KinematicCloud<ParcelType>::forces() const
{
    return forces_;
}


template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicCloud<ParcelType>::interpolationSchemes() const
{
    return interpolationSchemes_;
}


template<class ParcelType>
inline const Foam::DispersionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::dispersion() const
{
    return dispersionModel_;
}


template<class ParcelType>
inline Foam::DispersionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::dispersion()
{
    return dispersionModel_();
}


template<class ParcelType>
inline const Foam::DragModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::drag() const
{
    return dragModel_;
}


template<class ParcelType>
inline const Foam::InjectionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::injection() const
{
    return injectionModel_;
}

template<class ParcelType>
inline const Foam::PatchInteractionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::patchInteraction() const
{
    return patchInteractionModel_;
}

template<class ParcelType>
inline Foam::InjectionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::injection()
{
    return injectionModel_();
}

template<class ParcelType>
inline Foam::PatchInteractionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::patchInteraction() 
{
    return patchInteractionModel_();
}


template<class ParcelType>
inline Foam::PostProcessingModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::postProcessing()
{
    return postProcessingModel_();
}


template<class ParcelType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicCloud<ParcelType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
         const ParcelType& p = iter();
         sysMass += p.mass()*p.nParticle();
    }

    return sysMass;
}


template<class ParcelType>
inline Foam::Random& Foam::KinematicCloud<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UTrans()
{
    return UTrans_;
}


template<class ParcelType>
inline Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh> >
Foam::KinematicCloud<ParcelType>::SU() const
{
    tmp<DimensionedField<vector, volMesh> > tSU
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + "SU",
                this->db().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedVector
            (
                 "zero",
                 dimDensity*dimVelocity/dimTime,
                 vector::zero
            )
        )
    );

    vectorField& SU = tSU().field();
    SU = UTrans_/(mesh_.V()*this->db().time().deltaT());

    return tSU;
}

template<class CloudType>
inline Foam::fvVectorMatrix
Foam::KinematicCloud<CloudType>::SU(volVectorField& U) const
{

            //tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix fvm(U, dimForce);

            fvm.source() = -UTrans_/(this->db().time().deltaT());

            return fvm;
}

template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<ParcelType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& theta = ttheta().internalField();
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta /= mesh().V();

    return ttheta;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<ParcelType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha().internalField();
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh().V()*rho_);

    return talpha;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<ParcelType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "RhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff().internalField();
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        rhoEff[cellI] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh().V();

    return trhoEff;
}


// ************************************************************************* //
