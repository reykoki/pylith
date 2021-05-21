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

#include "FieldFilter.hh" // Implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::FieldFilter::FieldFilter(void) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::FieldFilter::~FieldFilter(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::FieldFilter::deallocate(void) {}


// End of file
