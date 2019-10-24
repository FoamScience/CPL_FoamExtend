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

#include "CrossPowerLawTempLog.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(CrossPowerLawTempLog, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        CrossPowerLawTempLog,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::CrossPowerLawTempLog::calcNu() const
{
    const volScalarField& T= U_.mesh().lookupObject<volScalarField>("T"); 
    volScalarField dT = 1.0/T - 1.0/T0_;
    volScalarField prefactorT = Foam::exp(a1_*dT + a2_*dT*dT);
    volScalarField sratePrefactorT = pow(b1_ * T, k_);
    
    // Added 1.0 to avoid log10(0) undefined behaviour
    // Info << Foam::log10(timeFactor_*strainRate()+scalar(1)) << " " << timeFactor_*strainRate() << " " << nl << endl;
    Foam::tmp<Foam::volScalarField> nu_out =\
                                   prefactorT * ((nu0_ - nuInf_)/(scalar(1) + pow(m_*sratePrefactorT*Foam::log10(timeFactor_*strainRate()+scalar(1)), n_)) + nuInf_);
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

Foam::viscosityModels::CrossPowerLawTempLog::CrossPowerLawTempLog
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CrossPowerLawTempLogCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    srate_cutoff_("srateCutoff", dimTime, CrossPowerLawTempLogCoeffs_),
    nu0_("nu0", dimViscosity, CrossPowerLawTempLogCoeffs_),
    nuInf_("nuInf", dimViscosity, CrossPowerLawTempLogCoeffs_),
    timeFactor_("timeFactor", dimTime, CrossPowerLawTempLogCoeffs_),
    m_("m", dimless, CrossPowerLawTempLogCoeffs_),
    n_("n", dimless, CrossPowerLawTempLogCoeffs_),
    a1_("a1", dimTemperature, CrossPowerLawTempLogCoeffs_),
    a2_("a2", dimTemperature2, CrossPowerLawTempLogCoeffs_),
    b1_("b1", dimInvTemperature, CrossPowerLawTempLogCoeffs_),
    T0_("T0", dimTemperature, CrossPowerLawTempLogCoeffs_),
    k_("k", dimless, CrossPowerLawTempLogCoeffs_),
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

bool Foam::viscosityModels::CrossPowerLawTempLog::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    CrossPowerLawTempLogCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    CrossPowerLawTempLogCoeffs_.lookup("nu0") >> nu0_;
    CrossPowerLawTempLogCoeffs_.lookup("nuInf") >> nuInf_;
    CrossPowerLawTempLogCoeffs_.lookup("m") >> m_;
    CrossPowerLawTempLogCoeffs_.lookup("n") >> n_;
    CrossPowerLawTempLogCoeffs_.lookup("a1") >> a1_;
    CrossPowerLawTempLogCoeffs_.lookup("a2") >> a2_;
    CrossPowerLawTempLogCoeffs_.lookup("b1") >> b1_;
    CrossPowerLawTempLogCoeffs_.lookup("T0") >> T0_;
    CrossPowerLawTempLogCoeffs_.lookup("k") >> k_;

    return true;
}


// ************************************************************************* //
