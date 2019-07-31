#ifndef ICOFOAM_INCOMING_FIELD_H
#define ICOFOAM_INCOMING_FIELD_H

#include "cpl/TransmittingField.h"
#include "fvCFD.H"
#include "meshSearch.H"
#include "cpl/cpl.h"
#include <valarray>

template <class FieldType>
class IncompressIncomingField : public CPL::IncomingField {
    public:
        IcoFOAMIncomingField() : IncomingField() {}
        IcoFOAMIncomingField(std::string name, const PortionField& portion_field, 
                             const PortionField& field, 
                             FieldType physical_field,
                             fvMesh& mesh, double density): 
                             CPL::IncomingField(name, portion_field, field),
                             physicalField(&physical_field),
                             temperatureField(&temperature),
                             foamMesh(&mesh), foamDensity(density) {}
        fvMesh* foamMesh; 
        Foam::meshSearch* meshSearcher;
        double foamDensity;
        std::vector<int> applyBc;
        Foam::vectorField bcFaceCenters;
        Foam::fvPatchVectorField* bcPatch;
        FieldType* physicalField;
        virtual ~IcoFOAMIncomingField(){delete meshSearcher;}



};

//TODO: Maybe use dependencies to compute stress once
class TempIncomingField:  public IncompressIncomingField<volScalarField> {
    public:
        using IcoFOAMIncomingField::IcoFOAMIncomingField;
        void unpack_(const std::vector<int>& glob_cell,
                     const std::vector<int>& loc_cell,
                     const std::valarray<double>& coord);
        void setup();
        void update();
        virtual ~VelIncomingField() {}

};


#endif //ICOFOAM_INCOMING_FIELD_H
