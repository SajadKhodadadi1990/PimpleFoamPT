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

#include "KinematicCloudTemplate.H"
#include "IntegrationScheme.H"
#include "interpolation.H"

#include "DispersionModel.H"
#include "DragModel.H"
#include "InjectionModel.H"
#include "PatchInteractionModel.H"
#include "PostProcessingModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::preEvolve()
{
    this->dispersion().cacheFields(true);
    forces_.cacheFields(true);
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolveCloud()
{
    autoPtr<interpolation<scalar> > rhoInterpolator =
        interpolation<scalar>::New
        (
            interpolationSchemes_,
            rho_
        );

    autoPtr<interpolation<vector> > UInterpolator =
        interpolation<vector>::New
        (
            interpolationSchemes_,
            U_
        );

    autoPtr<interpolation<scalar> > muInterpolator =
        interpolation<scalar>::New
        (
            interpolationSchemes_,
            mu_
        );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterpolator(),
        UInterpolator(),
        muInterpolator(),
        g_.value()
    );

    this->injection().inject(td);

    if (coupled_)
    {
        resetSourceTerms();
    }

    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::postEvolve()
{
    if (debug)
    {
        this->writePositions();
    }

    this->dispersion().cacheFields(false);
    forces_.cacheFields(false);

    this->postProcessing().post();
    if (this->db().time().outputTime())
    {
        outputProperties_.writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            this->db().time().writeCompression()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::KinematicCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    Cloud<ParcelType>(rho.mesh(), cloudName, false),
    kinematicCloud(),
    mesh_(rho.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    outputProperties_
    (
        IOobject
        (
            cloudName + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/cloudName,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    constProps_(particleProperties_),
    active_(particleProperties_.lookup("active")),
    parcelTypeId_(readLabel(particleProperties_.lookup("parcelTypeId"))),
    coupled_(particleProperties_.lookup("coupled")),
    cellValueSourceCorrection_
    (
        particleProperties_.lookup("cellValueSourceCorrection")
    ),
    rndGen_(label(0)),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    forces_(mesh_, particleProperties_, g_.value()),
    interpolationSchemes_(particleProperties_.subDict("interpolationSchemes")),
    dispersionModel_
    (
        DispersionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    dragModel_
    (
        DragModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    injectionModel_
    (
        InjectionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    patchInteractionModel_
    (
        PatchInteractionModel<KinematicCloud<ParcelType> >::New
        (
            this->particleProperties_,
            *this
        )
    ),
    postProcessingModel_
    (
        PostProcessingModel<KinematicCloud<ParcelType> >::New
        (
            this->particleProperties_,
            *this
        )
    ),
    UIntegrator_
    (
        vectorIntegrationScheme::New
        (
            "U",
            particleProperties_.subDict("integrationSchemes")
        )
    ),
    UTrans_
    (
        IOobject
        (
            this->name() + "UTrans",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimMass*dimVelocity, vector::zero)
    )
{
    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::~KinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    if (!fullyDescribed)
    {
        parcel.rho() = constProps_.rho0();
    }

    scalar carrierDt = this->db().time().deltaTValue();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::resetSourceTerms()
{
    UTrans_.field() = vector::zero;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolve()
{
    if (active_)
    {
        preEvolve();

        evolveCloud();

        postEvolve();

        info();
        Info<< endl;
    }
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::info() 
{
    Info<< "Cloud: " << this->name() << nl
        << "    Total number of parcels added   = "
        << this->injection().parcelsAddedTotal() << nl
        << "    Total mass introduced           = "
        << this->injection().massInjected() << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl;
        
        this->patchInteraction().info(Info); //Add by sajad

}


// ************************************************************************* //
