#ifndef ICOFOAM_OUTGOING_FIELD_H
#define ICOFOAM_OUTGOING_FIELD_H

#include "cpl/TransmittingField.h"
#include "fvCFD.H"
#include "meshSearch.H"
#include <valarray>

class IcoFOAMOutgoingField : public CPL::OutgoingField {
    public:
        IcoFOAMOutgoingField() : OutgoingField() {}
        IcoFOAMOutgoingField(std::string name, const PortionField& portion_field,
                             const PortionField& field,
                             volVectorField& velocity,
                             dimensionedScalar k_viscosity, 
                             fvMesh& mesh, double density): 
                             CPL::OutgoingField(name, portion_field, field),
                             velocityField(&velocity), kViscosity(k_viscosity),
                             foamMesh(&mesh), foamDensity(density) {}
        
        volVectorField* velocityField;
        dimensionedScalar kViscosity;
        fvMesh* foamMesh; 
        Foam::meshSearch* meshSearcher;
        double foamDensity;

        virtual ~IcoFOAMOutgoingField(){delete meshSearcher;}


};

//TODO: Maybe use dependencies to compute stress once
class StressOutgoingField:  public IcoFOAMOutgoingField {
    public:
        using IcoFOAMOutgoingField::IcoFOAMOutgoingField;
        void pack_(const std::vector<int>& glob_cell,
                   const std::vector<int>& loc_cell,
                   const std::valarray<double>& coord);
        void setup();
        void update();
        Foam::symmTensorField stressField;
        int mode;
        virtual ~StressOutgoingField() {}
};


#endif //ICOFOAM_OUTGOING_FIELD_H
