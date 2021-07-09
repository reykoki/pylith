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

#include "MeshIOPETSc.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "spatialdata/utils/LineParser.hh" // USES LineParser

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iomanip> // USES setw(), setiosflags(), resetiosflags()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <fstream> // USES std::ifstream, std::ofstream
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES std::typeid

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
    } // meshio
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOPETSc::MeshIOPETSc(void) :
    _useIndexZero(true) { // constructor
    PyreComponent::setName("meshiopetsc");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOPETSc::~MeshIOPETSc(void) { // destructor
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOPETSc::deallocate(void) { // deallocate
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// label mesh.
void
pylith::meshio::MeshIOPETSc::_label(void) { // _label
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_label()");
    
    assert(_mesh);
    MPI_Comm comm = _mesh->comm();
    PetscErrorCode err;
} // label


// ---------------------------------------------------------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOPETSc::_read(void) { // _read
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");
    
    assert(_mesh);
    MPI_Comm comm = _mesh->comm();
    PetscErrorCode err;

    PetscDM dm = NULL;
    err = DMCreate(comm, &dm);PYLITH_CHECK_ERROR(err);
    err = DMSetType(dm, DMPLEX);PYLITH_CHECK_ERROR(err);
    err = DMSetFromOptions(dm);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dm, NULL, "-dm_view");PYLITH_CHECK_ERROR(err);

    char labelName[] = "material-id";
    DMLabel label;
    err = DMCreateLabel(dm, labelName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, labelName, &label);PYLITH_CHECK_ERROR(err);
    
    PetscInt p, pStart, pEnd;
    err = DMPlexGetHeightStratum(dm, 0, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    for (p=pStart; p < pEnd; p++) {
        PetscInt val;
        err = DMLabelSetValue(label, p, 1);PYLITH_CHECK_ERROR(err);
    }
    
    DMLabel face_sets_label;
    err = DMGetLabel(dm, "Face Sets", &face_sets_label);PYLITH_CHECK_ERROR(err);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelComplete(dm, face_sets_label);PYLITH_CHECK_ERROR(err);

    IS is;
    DMLabel xnegLabel;
    char xnegName[] = "boundary_xneg";
    err = DMCreateLabel(dm, xnegName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, xnegName, &xnegLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(face_sets_label, 1, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(xnegLabel, 1, is);PYLITH_CHECK_ERROR(err);
    
    DMLabel xposLabel;
    char xposName[] = "boundary_xpos";
    err = DMCreateLabel(dm, xposName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, xposName, &xposLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(face_sets_label, 3, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(xposLabel, 1, is);PYLITH_CHECK_ERROR(err);
    
    DMLabel domainLabel;
    char domainName[] = "domain_all";
    err = DMCreateLabel(dm, domainName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, domainName, &domainLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(face_sets_label, 1, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(face_sets_label, 2, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(face_sets_label, 3, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(face_sets_label, 4, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);

    _mesh->dmMesh(dm);

    PYLITH_METHOD_END;
} // read


// ---------------------------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOPETSc::_write(void) const {  
    PYLITH_JOURNAL_LOGICERROR("Writing meshes via MeshIOPetsc not implemented.");
} // write


