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
    Foam::PatchInteractionModel

Description
    Templated patch interaction model class

SourceFiles
    PatchInteractionModel.C
    NewPatchInteractionModel.C

\*---------------------------------------------------------------------------*/

#ifndef PatchInteractionModel_H
#define PatchInteractionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class PatchInteractionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchInteractionModel
{
public:

    // Public enumerations

        // Interaction types
        enum interactionType
        {
            itRebound,
            itStick,
            itEscape,
            itOther
        };

        static wordList interactionTypeNames_;


private:

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("PatchInteractionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PatchInteractionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from components
        PatchInteractionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~PatchInteractionModel();


    //- Selector
    static autoPtr<PatchInteractionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        //- Convert interaction result to word
        static word interactionTypeToWord(const interactionType& itEnum);

        //- Convert word to interaction result
        static interactionType wordToInteractionType(const word& itWord);

        //- Flag to indicate whether model activates patch interaction model
        virtual bool active() const = 0;

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            const polyPatch& pp,
            const label faceId,
            bool& keepParticle,
            vector& U
        ) const = 0;
        virtual bool correct1
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            const label faceId,
            bool& keepParticle,
            vector& U
        )  = 0;
        virtual void info(Ostream& os);//Add by sajad
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePatchInteractionModel(CloudType)                                  \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(PatchInteractionModel<CloudType>, 0); \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        PatchInteractionModel<CloudType>,                                     \
        dictionary                                                            \
    );


#define makePatchInteractionModelType(SS, CloudType, ParcelType)              \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    PatchInteractionModel<CloudType<ParcelType> >::                           \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


#define makePatchInteractionModelThermoType(SS, CloudType, ParcelType, ThermoType)\
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    PatchInteractionModel<CloudType<ParcelType<ThermoType> > >::              \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PatchInteractionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
