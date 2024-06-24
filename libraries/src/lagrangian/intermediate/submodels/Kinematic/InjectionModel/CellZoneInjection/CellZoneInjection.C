/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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

#include "CellZoneInjection.H"
#include "DataEntry.H"
#include "pdf.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::CellZoneInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0) && (time0 < duration_))
    {
        return round(fraction_*(time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::CellZoneInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return fraction_*volumeFlowRate_().integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellZoneInjection<CloudType>::CellZoneInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    cellZoneId_
    (
        this->coeffDict().lookup("cellZoneName"),
        owner.mesh().cellZones()
    ),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    U0_(this->coeffDict().lookup("U0")),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            this->coeffDict()
        )
    ),
    parcelPDF_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    fraction_(1.0)
{
    if (!cellZoneId_.active())
    {
        FatalErrorInFunction
            << "Requested cellZone " << cellZoneId_.name() << " not found" << nl
            << "Available zones are: " << owner.mesh().cellZones().names()
            << nl << abort(FatalError);
    }

    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellZoneInjection<CloudType>::~CellZoneInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::CellZoneInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::CellZoneInjection<CloudType>::updateMesh()
{
    const cellZone& cz = this->owner().mesh().cellZones()[cellZoneId_.index()];

    label czSize = cz.size();
    label totalCzSize = czSize;
    reduce(totalCzSize, sumOp<label>());
    fraction_ = scalar(czSize)/totalCzSize;

    // Set total volume/mass to inject
    this->volumeTotal_ = fraction_*volumeFlowRate_().integrate(0.0, duration_);
    this->massTotal_ *= fraction_;
}


template<class CloudType>
Foam::scalar Foam::CellZoneInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::CellZoneInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    const cellZone& cz = this->owner().mesh().cellZones()[cellZoneId_.index()];

    if (!cz.empty())
    {
        label cellI = this->owner().rndGen().integer(0, cz.size() - 1);

        cellOwner = cz[cellI];
        position = this->owner().mesh().C()[cellOwner];
    }
    else
    {
        cellOwner = -1;
        // dummy position
        position = pTraits<vector>::max;
    }
}


template<class CloudType>
void Foam::CellZoneInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = parcelPDF_->sample();
}


template<class CloudType>
bool Foam::CellZoneInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::CellZoneInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
