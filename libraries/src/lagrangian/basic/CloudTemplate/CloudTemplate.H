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

Class
    Foam::Cloud

Description

SourceFiles
    CloudTemplate.C
    CloudTemplateIO.C

\*---------------------------------------------------------------------------*/

#ifndef CloudTemplate_H
#define CloudTemplate_H

#include "cloud.H"
#include "IDLList.H"
#include "IOField.H"
#include "polyMesh.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of functions
template<class ParticleType>
class Cloud;

template<class ParticleType>
class IOPosition;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const Cloud<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                           Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class Cloud
:
    public cloud,
    public IDLList<ParticleType>
{
    // Private data

      //  const fvMesh& fvmesh_;
        const polyMesh& polyMesh_;
        //- Overall count of particles ever created. Never decreases.
        mutable label particleCount_;

        //- Temporary storage for addressing. Used in findFaces.
        mutable dynamicLabelList labels_;


    // Private member functions

        //- Initialise cloud on IO constructor
        void initCloud(const bool checkClass);

        //- Read cloud properties dictionary
        void readCloudUniformProperties();

        //- Write cloud properties dictionary
        void writeCloudUniformProperties() const;


public:

    template<class ParticleT>
    friend class Particle;
    template<class ParticleT>
    friend class IOPosition;

    typedef ParticleType particleType;

    typedef typename IDLList<ParticleType>::iterator iterator;
    typedef typename IDLList<ParticleType>::const_iterator const_iterator;

    //-Runtime type information
    TypeName("Cloud");


    // Static data

        //- Name of cloud properties dictionary
        static word cloudPropertiesName;


    // Constructors

        //- Construct from mesh and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh, cloud name, and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh by reading from file
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& mesh,
            const bool checkClass = true
        );


        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        // Access

            //- Return the polyMesh reference
            /*const fvMesh& fvmesh() const
            {
                return fvmesh_;
            }*/

            const polyMesh& pMesh() const
            {
                return polyMesh_;
            }

            //- Is this global face an internal face?
            bool internalFace(const label facei) const
            {
                return polyMesh_.isInternalFace(facei);
            }

            //- Is this global face a boundary face?
            bool boundaryFace(const label facei) const
            {
                return !internalFace(facei);
            }

            //- Which patch is this global face on
            label facePatch(const label facei) const
            {
                return polyMesh_.boundaryMesh().whichPatch(facei);
            }

            //- Which face of this patch is this global face
            label patchFace(const label patchi, const label facei) const
            {
                return polyMesh_.boundaryMesh()[patchi].whichFace(facei);
            }

            label size() const
            {
                return IDLList<ParticleType>::size();
            };


            // Iterators

                const const_iterator begin() const
                {
                    return IDLList<ParticleType>::begin();
                };

                const const_iterator cbegin() const
                {
                    return IDLList<ParticleType>::cbegin();
                };

                const const_iterator end() const
                {
                    return IDLList<ParticleType>::end();
                };

                const const_iterator cend() const
                {
                    return IDLList<ParticleType>::cend();
                };

                iterator begin()
                {
                    return IDLList<ParticleType>::begin();
                };

                iterator end()
                {
                    return IDLList<ParticleType>::end();
                };


        // Edit

            void clear()
            {
                return IDLList<ParticleType>::clear();
            };

            //- Get unique particle creation id
            label getNewParticleID() const;

            //- Transfer particle to cloud
            void addParticle(ParticleType* pPtr);

            //- Remove particle from cloud and delete
            void deleteParticle(ParticleType&);

            //- Move the particles
            //  passing the TrackingData to the track function
            template<class TrackingData>
            void move(TrackingData& td);

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            virtual void autoMap(const mapPolyMesh&);


        // Read

            //- Helper to construct IOobject for field and current time.
            IOobject fieldIOobject
            (
                const word& fieldName,
                const IOobject::readOption r
            ) const;

            //- Check lagrangian data field
            template<class DataType>
            void checkFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const IOField<DataType>& data
            ) const;

            //- Read the field data for the cloud of particles. Dummy at
            //  this level.
            virtual void readFields();


        // Write

            //- Write the field data for the cloud of particles Dummy at
            //  this level.
            virtual void writeFields() const;

            //- Write using given format, version and compression.
            //  Only writes the cloud file if the Cloud isn't empty
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;

            //- Write positions to <cloudName>_positions.obj file
            void writePositions() const;


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const Cloud<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
