/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::CrossPowerLawTempLog

Description
    An incompressible Cross-Power law non-Newtonian viscosity model.

SourceFiles
    CrossPowerLawTempLog.C

\*---------------------------------------------------------------------------*/

#ifndef CrossPowerLawTempLog_H
#define CrossPowerLawTempLog_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class CrossPowerLawTempLog Declaration
\*---------------------------------------------------------------------------*/

const Foam::dimensionSet dimInvTemperature(0, 0, 0, -1, 0, 0, 0);
const Foam::dimensionSet dimTemperature2(0, 0, 0, 2, 0, 0, 0);

class CrossPowerLawTempLog
:
    public viscosityModel
{
    // Private data

        dictionary CrossPowerLawTempLogCoeffs_;

        dimensionedScalar srate_cutoff_;
        dimensionedScalar nu0_;
        dimensionedScalar nuInf_;
        dimensionedScalar m_;
        dimensionedScalar n_;
        dimensionedScalar a1_;
        dimensionedScalar a2_;
        dimensionedScalar b1_;
        dimensionedScalar T0_;
        dimensionedScalar k_;
        dimensionedScalar timeFactor_;

        volScalarField nu_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("CrossPowerLawTempLog");


    // Constructors

        //- Construct from components
        CrossPowerLawTempLog
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~CrossPowerLawTempLog()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
