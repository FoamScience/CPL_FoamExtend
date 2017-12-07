#ifndef ICOFOAM_INCOMING_FIELD_H
#define ICOFOAM_INCOMING_FIELD_H

#include "cpl/TransmittingField.h"
#include "fvCFD.H"
#include "meshSearch.H"
#include <valarray>

class IcoFOAMIncomingField : public CPL::IncomingField {
    public:
        IcoFOAMIncomingField() : IncomingField() {}
        IcoFOAMIncomingField(std::string name, const PortionField& portion_field, 
                             const PortionField& field, 
                             volVectorField& velocity,
                             dimensionedScalar k_viscosity, 
                             fvMesh& mesh, double density): 
                             CPL::IncomingField(name, portion_field, field),
                             velocityField(&velocity), kViscosity(k_viscosity),
                             foamMesh(&mesh), foamDensity(density) {}
        
        volVectorField* velocityField;
        dimensionedScalar kViscosity;
        fvMesh* foamMesh; 
        Foam::meshSearch* meshSearcher;
        double foamDensity;

        virtual ~IcoFOAMIncomingField(){delete meshSearcher;}


};

//TODO: Maybe use dependencies to compute stress once
class VelIncomingField:  public IcoFOAMIncomingField {
    public:
        using IcoFOAMIncomingField::IcoFOAMIncomingField;
        void unpack_(const std::vector<int>& glob_cell,
                     const std::vector<int>& loc_cell,
                     const std::valarray<double>& coord);
        void setup();
        void update();
        int mode;
        int applyBCx, applyBCy, applyBCz;
        Foam::vectorField bcFaceCenters;
        Foam::fvPatchVectorField* bcPatch;
        virtual ~VelIncomingField() {}

};


#endif //ICOFOAM_INCOMING_FIELD_H
