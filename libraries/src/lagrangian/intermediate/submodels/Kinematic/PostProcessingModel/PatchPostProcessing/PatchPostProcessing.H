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
    Foam::PatchPostProcessing

Description
    Standard post-processing

SourceFiles
    PatchPostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef PatchPostProcessing_H
#define PatchPostProcessing_H

#include "PostProcessingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PatchPostProcessing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchPostProcessing
:
    public PostProcessingModel<CloudType>
{
    // Private data

        typedef typename CloudType::parcelType parcelType;

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Maximum number of parcels to store
        label maxStoredParcels_;


        //- List of patch names
        wordList patchNames_;

        //- List of output data per patch
        List<DynamicList<string> > patchData_;

        //- Mapping from local to global patch ids
        labelList patchIds_;


    // Private member functions

        //- Returns local patchI if patch is in patchIds_ list
        label applyToPatch(const label globalPatchI) const;


protected:

    // Protected member functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("PatchPostProcessing");


    // Constructors

        //- Construct from dictionary
        PatchPostProcessing(const dictionary& dict, CloudType& owner);


    //- Destructor
    virtual ~PatchPostProcessing();


    // Member Functions

        // Access

            //- Return const access to the mesh
            inline const polyMesh& mesh() const;

            //- Return maximum number of parcels to store per patch
            inline label maxStoredParcels() const;

            //- Return const access to the list of patch names
            inline const wordList& patchNames() const;

            //- Return const mapping from local to global patch ids
            inline const labelList& patchIds() const;


        // Evaluation

            //- Flag to indicate whether model activates post-processing model
            virtual bool active() const;

            //- Gather post-processing info on patch
            virtual void postPatch
            (
                const parcelType& p,
                const label patchI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchPostProcessingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PatchPostProcessing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
