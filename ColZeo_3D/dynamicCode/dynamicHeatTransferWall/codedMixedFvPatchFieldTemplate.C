/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedMixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
#line 43 "/home/alexis/OpenFOAM/Tesis/run/ColZeo/0/T/boundaryField/wall"
#include "physicoChemicalConstants.H" // Necesario para sigma
            #include <cmath>
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 7460da6a7de59036b3a94aaa81860346263d4b05
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void dynamicHeatTransferWall_7460da6a7de59036b3a94aaa81860346263d4b05(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    dynamicHeatTransferWallMixedValueFvPatchScalarField
);


const char* const dynamicHeatTransferWallMixedValueFvPatchScalarField::SHA1sum =
    "7460da6a7de59036b3a94aaa81860346263d4b05";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicHeatTransferWallMixedValueFvPatchScalarField::
dynamicHeatTransferWallMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct dynamicHeatTransferWall sha1: 7460da6a7de59036b3a94aaa81860346263d4b05"
            " from patch/dictionary\n";
    }
}


dynamicHeatTransferWallMixedValueFvPatchScalarField::
dynamicHeatTransferWallMixedValueFvPatchScalarField
(
    const dynamicHeatTransferWallMixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct dynamicHeatTransferWall sha1: 7460da6a7de59036b3a94aaa81860346263d4b05"
            " from patch/DimensionedField/mapper\n";
    }
}


dynamicHeatTransferWallMixedValueFvPatchScalarField::
dynamicHeatTransferWallMixedValueFvPatchScalarField
(
    const dynamicHeatTransferWallMixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct dynamicHeatTransferWall sha1: 7460da6a7de59036b3a94aaa81860346263d4b05 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicHeatTransferWallMixedValueFvPatchScalarField::
~dynamicHeatTransferWallMixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy dynamicHeatTransferWall sha1: 7460da6a7de59036b3a94aaa81860346263d4b05\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicHeatTransferWallMixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs dynamicHeatTransferWall sha1: 7460da6a7de59036b3a94aaa81860346263d4b05\n";
    }

//{{{ begin code
    #line 48 "/home/alexis/OpenFOAM/Tesis/run/ColZeo/0/T/boundaryField/wall"
// 1. Parámetros
            const scalar T_inf = 298.15;
            const scalar epsilon_rad = 0.85;
            const scalar sigma = Foam::constant::physicoChemical::sigma.value();
            const scalar g = 9.81;
            const scalar L_char = 1.2;
            const scalar nu_air = 1.568e-5;
            const scalar Pr_air = 0.71;
            const scalar k_air  = 0.0262;
            const scalar beta   = 1.0 / T_inf;

            // 2. Recuperación de campos
            const scalarField Tw = this->patchInternalField();
            const scalarField& delta = this->patch().deltaCoeffs();

            // PROTECCIÓN: Comprobar si kappaEff ya fue creado en createFields.H
            bool hasKappa = this->db().foundObject<volScalarField>("kappaEff");
            
            // 3. Cálculo celda por celda
            forAll(Tw, i)
            {
                scalar deltaT = mag(Tw[i] - T_inf) + 1e-6;
                scalar Ra = (g * beta * deltaT * pow(L_char, 3.0) / pow(nu_air, 2.0)) * Pr_air;

                // Churchill-Chu
                scalar Nu = pow(0.825 + (0.387 * pow(Ra, 1.0/6.0)) / 
                            pow(1.0 + pow(0.492 / Pr_air, 9.0/16.0), 8.0/27.0), 2.0);

                scalar h_conv = (Nu * k_air) / L_char;
                scalar h_rad = epsilon_rad * sigma * (pow(Tw[i], 2) + pow(T_inf, 2)) * (Tw[i] + T_inf);
                scalar h_total = h_conv + h_rad;

                // 4. Acoplamiento Protegido
                scalar k_over_d = 0.0;
                
                if (hasKappa)
                {
                    const volScalarField& kappa = this->db().lookupObject<volScalarField>("kappaEff");
                    const scalarField& kappa_w = kappa.boundaryField()[this->patch().index()];
                    k_over_d = kappa_w[i] * delta[i];
                }

                this->refValue()[i] = T_inf;
                // Agregamos + 1e-10 para asegurar que nunca haya división por cero
                this->valueFraction()[i] = h_total / (h_total + k_over_d + 1e-10);
            }
//}}} end code

    this->mixedFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

