// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "FieldFilterProject.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::FieldFilterProject::FieldFilterProject(void) :
    _basisOrder(1) {
    PyreComponent::setName("fieldfilterproject");
    _passThruFns[0] = passThruSoln;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::FieldFilterProject::~FieldFilterProject(void) {}


// ------------------------------------------------------------------------------------------------
// Set basis order for projected field.
void
pylith::meshio::FieldFilterProject::setBasisOrder(const int value) {
    PYLITH_METHOD_BEGIN;

    if (value < 0) {
        std::ostringstream msg;
        msg << "Basis order (" << value << ") for field filter must be nonnegative.";
        throw std::runtime_error(msg.str());
    } // if

    _basisOrder = value;

    PYLITH_METHOD_END;
} // basisOrder


// ------------------------------------------------------------------------------------------------
// Create DM associated with filtered field.
PetscDM
pylith::meshio::FieldFilterProject::createDM(PetscDM dm,
                                             const pylith::topology::FieldBase::Description& description,
                                             const pylith::topology::FieldBase::Discretization& discretization) const {
    PYLITH_METHOD_BEGIN;

    if (discretization.quadOrder < _basisOrder) {
        PYLITH_COMPONENT_WARNING(
            "Projecting subfield '"
                << description.label << "' in field '" << description.label << " from basis order "
                << discretization.basisOrder << " to basis order " << _basisOrder
                << " with quadrature order " << discretization.quadOrder << " will result in under integration of the "
                << "subfield. Accurate projection requires a quadrature order of at least " << _basisOrder << "."
            );
    } // if

    PetscErrorCode err;

    PetscDM dmOut = NULL;
    err = DMClone(dm, &dmOut);PYLITH_CHECK_ERROR(err);

    pylith::topology::FieldBase::Discretization discretizationFiltered(discretization);
    discretizationFiltered.basisOrder = _basisOrder;

    PetscFE fe = pylith::topology::FieldOps::createFE(discretizationFiltered, dmOut, description.numComponents);assert(fe);
    err = PetscFESetName(fe, description.label.c_str());PYLITH_CHECK_ERROR(err);
    err = DMSetField(dmOut, 0, NULL, (PetscObject)fe);PYLITH_CHECK_ERROR(err);
    err = DMSetFieldAvoidTensor(dmOut, 0, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    err = PetscFEDestroy(&fe);PYLITH_CHECK_ERROR(err);

    err = DMCreateDS(dmOut);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(dmOut);
}


// ------------------------------------------------------------------------------------------------
// Apply filter to global PETSc vector.
void
pylith::meshio::FieldFilterProject::apply(PetscVec* vectorOut,
                                          PetscDM dmOut,
                                          PetscVec vectorIn) const {
    PYLITH_METHOD_BEGIN;
    assert(vectorOut);
    assert(*vectorOut);

    PetscErrorCode err = 0;
    const PylithReal t = 0.0;
    err = DMProjectField(dmOut, t, vectorIn, const_cast<PetscPointFunc*>(_passThruFns),
                         INSERT_VALUES, *vectorOut);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ---------------------------------------------------------------------------------------------------------------------
// Identify function kernel.
void
pylith::meshio::FieldFilterProject::passThruSoln(const PylithInt dim,
                                                 const PylithInt numS,
                                                 const PylithInt numA,
                                                 const PylithInt sOff[],
                                                 const PylithInt sOff_x[],
                                                 const PylithScalar s[],
                                                 const PylithScalar s_t[],
                                                 const PylithScalar s_x[],
                                                 const PylithInt aOff[],
                                                 const PylithInt aOff_x[],
                                                 const PylithScalar a[],
                                                 const PylithScalar a_t[],
                                                 const PylithScalar a_x[],
                                                 const PylithReal t,
                                                 const PylithScalar x[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar field[]) {
    assert(s);
    assert(field);
    assert(1 == numS);

    const PylithInt sEnd = sOff[1];
    for (PylithInt i = sOff[0]; i < sEnd; ++i) {
        field[i] = s[i];
    } // for
} // passThruSoln


// End of file
