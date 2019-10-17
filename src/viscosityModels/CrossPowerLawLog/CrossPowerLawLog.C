/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CrossPowerLawLog.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(CrossPowerLawLog, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        CrossPowerLawLog,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::CrossPowerLawLog::calcNu() const
{
    
    // Added 1.0 to avoid log10(0) undefined behaviour
    // Info << Foam::log10(timeFactor_*strainRate()+scalar(1)) << " " << timeFactor_*strainRate() << " " << nl << endl;
    Foam::tmp<Foam::volScalarField> nu_out =\
                                   (nu0_ - nuInf_)/(scalar(1) + pow(m_*Foam::log10(timeFactor_*strainRate()+scalar(1)), n_)) + nuInf_;
    Foam::tmp<Foam::volScalarField> srate = strainRate();
    forAll(srate(), celli) {
        // cout << "cell: " << celli << ", " << nu_out()[celli] << std::endl;
        if (srate()[celli] > srate_cutoff_.value())
            nu_out()[celli] = 0.0;
    }
    // Info << nu_out() << nl << endl;
    return nu_out;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::CrossPowerLawLog::CrossPowerLawLog
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CrossPowerLawLogCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    srate_cutoff_("srateCutoff", dimTime, CrossPowerLawLogCoeffs_),
    nu0_("nu0", dimViscosity, CrossPowerLawLogCoeffs_),
    nuInf_("nuInf", dimViscosity, CrossPowerLawLogCoeffs_),
    timeFactor_("timeFactor", dimTime, CrossPowerLawLogCoeffs_),
    m_("m", dimless, CrossPowerLawLogCoeffs_),
    n_("n", dimless, CrossPowerLawLogCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::CrossPowerLawLog::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    CrossPowerLawLogCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    CrossPowerLawLogCoeffs_.lookup("nu0") >> nu0_;
    CrossPowerLawLogCoeffs_.lookup("nuInf") >> nuInf_;
    CrossPowerLawLogCoeffs_.lookup("m") >> m_;
    CrossPowerLawLogCoeffs_.lookup("n") >> n_;

    return true;
}


// ************************************************************************* //
