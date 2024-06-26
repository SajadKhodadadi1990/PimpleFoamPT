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

Class
    Foam::InjectionModel

Description
    Templated injection model class.

    The injection model nominally describes the parcel:
    - position
    - diameter
    - velocity
    In this case, the fullyDescribed() flag should be set to 0 (false). When
    the parcel is then added to the cloud, the remaining properties are
    populated using values supplied in the constant properties.

    If, however, all of a parcel's properties are described in the model, the
    fullDescribed() flag should be set to 1 (true).


SourceFiles
    InjectionModel.C
    NewInjectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef InjectionModel_H
#define InjectionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class InjectionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectionModel
{
public:

    // Enumerations

        //- Parcel basis representation options
        //  i.e constant number of particles OR constant mass per parcel
        enum parcelBasis
        {
            pbNumber,
            pbMass,
            pbFixed
        };


private:

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


    // Private member functions

        //- Read injector properties from previous run (if applicable)
        void readProps();

        //- Write injector properties
        void writeProps();


protected:

    // Protected data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        // Global injection properties

            //- Start of injection [s]
            const scalar SOI_;

            //- Total volume of particles introduced by this injector [m^3]
            //  - scaled to ensure massTotal is achieved
            scalar volumeTotal_;

            //- Total mass to inject [kg]
            scalar massTotal_;

            //- Total mass injected to date [kg]
            scalar massInjected_;


        // Counters

            //- Number of injections counter
            label nInjections_;

            //- Running counter of total number of parcels added
            label parcelsAddedTotal_;


        // Injection properties per Lagrangian time step

            //- Parcel basis enumeration
            parcelBasis parcelBasis_;

            //- Continuous phase time at start of injection time step [s]
            scalar time0_;

            //- Time at start of injection time step [s]
            scalar timeStep0_;

            scalar nParticleFixed_;//sajad added 18 Feb 2023
    // Protected member functions

        //- Number of parcels to introduce over the time step relative to SOI
        virtual label parcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const = 0;

        //- Volume of parcels to introduce over the time step relative to SOI
        virtual scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const = 0;

        //- Additional flag to identify whether or not injection of parcelI is
        //  permitted
        virtual bool validInjection(const label parcelI) = 0;

        //- Determine properties for next time step/injection interval
        virtual void prepareForNextTimeStep
        (
            const scalar time,
            label& newParcels,
            scalar& newVolume
        );

        //- Find the cell that contains the supplied position
        //  Will modify position slightly towards the owner cell centroid to
        //  ensure that it lies in a cell and not edge/face
        virtual void findCellAtPosition(label& cellI, vector& position);

        //- Set number of particles to inject given parcel properties
        virtual scalar setNumberOfParticles
        (
            const label parcels,
            const scalar volume,
            const scalar diameter,
            const scalar rho
        );

        //- Post injection checks
        virtual void postInjectCheck
        (
            const label parcelsAdded,
            const scalar massAdded
        );


public:

    //- Runtime type information
    TypeName("InjectionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InjectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        InjectionModel(CloudType& owner);

        //- Construct from dictionary
        InjectionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~InjectionModel();


    //- Selector
    static autoPtr<InjectionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return the owner cloud dictionary
        inline const dictionary& dict() const;

        //- Return const access the owner cloud object
        inline const CloudType& owner() const;

        //- Return non-const access the owner cloud object for manipulation
        inline CloudType& owner();

        //- Return the coefficients dictionary
        inline const dictionary& coeffDict() const;


    // Member Functions

        //- Flag to indicate whether model activates injection model
        virtual bool active() const = 0;


        // Mapping

            //- Update mesh
            virtual void updateMesh();


        // Global information

            //- Return the start-of-injection time
            inline scalar timeStart() const;

            //- Return the total volume to be injected across the event
            inline scalar volumeTotal() const;

            //- Return mass of particles to introduce
            inline scalar massTotal() const;

            //- Return mass of particles injected (cumulative)
            inline scalar massInjected() const;

            //- Return the end-of-injection time
            virtual scalar timeEnd() const = 0;

            // Counters

                //- Return the number of injections
                inline label nInjections() const;

                //- Return the total number parcels added
                inline label parcelsAddedTotal() const;


        // Per-injection event functions

            //- Main injection loop
            template<class TrackData>
            void inject(TrackData& td);


        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner
            ) = 0;

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            ) = 0;

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInjectionModel(CloudType)                                         \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(InjectionModel<CloudType>, 0);        \
                                                                              \
    defineTemplateRunTimeSelectionTable(InjectionModel<CloudType>, dictionary);


#define makeInjectionModelType(SS, CloudType, ParcelType)                     \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    InjectionModel<CloudType<ParcelType> >::                                  \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


#define makeInjectionModelThermoType(SS, CloudType, ParcelType, ThermoType)   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    InjectionModel<CloudType<ParcelType<ThermoType> > >::                     \
        adddictionaryConstructorToTable                                       \
        <SS<CloudType<ParcelType<ThermoType> > > >                            \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InjectionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InjectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
