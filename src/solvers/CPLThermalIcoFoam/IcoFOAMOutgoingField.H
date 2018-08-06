#ifndef ICOFOAM_OUTGOING_FIELD_H
#define ICOFOAM_OUTGOING_FIELD_H

#include "cpl/TransmittingField.h"
#include "cpl/cpl.h"
#include "fvCFD.H"
#include "meshSearch.H"
#include <valarray>

class IcoFOAMOutgoingField : public CPL::OutgoingField {
    public:
        IcoFOAMOutgoingField() : OutgoingField() {}
        IcoFOAMOutgoingField(std::string name, const PortionField& portion_field,
                             const PortionField& field, const volSymmTensorField& sigma,
                             const fvMesh& mesh): 
                             CPL::OutgoingField(name, portion_field, field),
                             sigma(&sigma), foamMesh(&mesh) {}
        
        const fvMesh* foamMesh; 
        Foam::meshSearch* meshSearcher;
        const Foam::volSymmTensorField* sigma;
        virtual ~IcoFOAMOutgoingField(){delete meshSearcher;}
        // virtual ~IcoFOAMOutgoingField(){}


};

class StressOutgoingField:  public IcoFOAMOutgoingField {
    public:
        using IcoFOAMOutgoingField::IcoFOAMOutgoingField;
        void pack_(const std::vector<int>& glob_cell,
                   const std::vector<int>& loc_cell,
                   const std::valarray<double>& coord);
        void setup();
        void update();
        Foam::symmTensorField stressField;
        std::string compute_mode;
        virtual ~StressOutgoingField() {}
};


#endif //ICOFOAM_OUTGOING_FIELD_H
