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

#include "FieldFilterNone.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::FieldFilterNone::FieldFilterNone(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::FieldFilterNone::~FieldFilterNone(void) {}


// ------------------------------------------------------------------------------------------------
// Create DM associated with filtered field.
PetscDM
pylith::meshio::FieldFilterNone::createDM(PetscDM dm,
                                          const pylith::topology::FieldBase::Description& description,
                                          const pylith::topology::FieldBase::Discretization& discretization) const {
    PetscErrorCode err = PetscObjectReference((PetscObject)dm);PYLITH_CHECK_ERROR(err);
    return dm;
}


// ------------------------------------------------------------------------------------------------
// Apply filter to global PETSc vector.
void
pylith::meshio::FieldFilterNone::apply(PetscVec* vectorOut,
                                       PetscDM dmOut,
                                       PetscVec vectorIn) const {
    assert(vectorOut);
    assert(*vectorOut == vectorIn);
}


// End of file
