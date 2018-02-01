// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "NeumannNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

extern "C" PetscErrorCode DMPlexComputeBdResidual_Internal(PetscDM dm,
                                                           PetscInt cStart,
                                                           PetscInt cEnd,
                                                           PetscReal time,
                                                           PetscVec locX,
                                                           PetscVec locX_t,
                                                           PetscVec locF,
                                                           void *user);


// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannNew::NeumannNew(void) :
    _boundaryMesh(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannNew::~NeumannNew(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannNew::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    IntegratorPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate

// Set up direction to discriminate among shear directions in 3-D.
// ----------------------------------------------------------------------
void
pylith::bc::NeumannNew::upDir(const double vec[3]) {
    // Set up direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    for (int i = 0; i < 3; ++i) {
        _upDir[i] = vec[i] / mag;
    } // for
} // upDir

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::NeumannNew::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_field.c_str())) {
        std::ostringstream msg;
        msg << "Cannot apply Neumann boundary condition to field '"<< _field
            << "; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::NeumannNew::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    _setFEKernelsRHSResidual(solution);

    _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxField; _auxField = new pylith::topology::Field(*_boundaryMesh); assert(_auxField);
    _auxField->label("auxiliary fields");
    _auxFieldSetup(solution);
    _auxField->subfieldsSetup();
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->initializeSubfields();

#if 0 // Use low-level function to set kernels
    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // :TODO: @brad @matt Expect this to need updating once we associate point functions with a domain.
    void* context = NULL;
    const int labelId = 1;
    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    err = PetscDSAddBoundary(prob, DM_BC_NATURAL, label(), label(), info.index, 0, NULL,
                             NULL, 1, &labelId, context); PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::bc::NeumannNew::computeRHSResidual(pylith::topology::Field* residual,
                                           const PylithReal t,
                                           const PylithReal dt,
                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    _setFEKernelsRHSResidual(solution);
    _setFEConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");

    assert(residual);
    assert(_auxField);

    PetscDS prob = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();
    PetscDMLabel dmLabel;

    // Pointwise function have been set in DS
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    // Compute the local residual
#if 0 // :TODO: @brad @matt Update this for Matt's interface for computing the residual on a boundary.
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, "material-id", &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd); PYLITH_CHECK_ERROR(err);

    PYLITH_COMPONENT_DEBUG("DMPlexComputeResidual_Internal() with material-id '"<<id()<<"' and cells ["<<cStart<<","<<cEnd<<")");
    err = DMPlexComputeResidual_Internal(dmSoln, cStart, cEnd, PETSC_MIN_REAL, solution.localVector(), solutionDot.localVector(), residual->localVector(), NULL); PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian and preconditioner for G(t,s).
void
pylith::bc::NeumannNew::computeRHSJacobian(PetscMat jacobianMat,
                                           PetscMat preconMat,
                                           const PylithReal t,
                                           const PylithReal dt,
                                           const pylith::topology::Field& solution) {}

// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::bc::NeumannNew::computeLHSResidual(pylith::topology::Field* residual,
                                           const PylithReal t,
                                           const PylithReal dt,
                                           const pylith::topology::Field& solution,
                                           const pylith::topology::Field& solutionDot) {}

// ----------------------------------------------------------------------
// Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
void
pylith::bc::NeumannNew::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                   PetscMat precondMat,
                                                   const PylithReal t,
                                                   const PylithReal dt,
                                                   const PylithReal tshift,
                                                   const pylith::topology::Field& solution,
                                                   const pylith::topology::Field& solutionDot) {}


// ----------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
void
pylith::bc::NeumannNew::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                    const PylithReal t,
                                                    const PylithReal dt,
                                                    const pylith::topology::Field& solution) {}


// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::bc::NeumannNew::_setFEConstants(const pylith::topology::Field& solution,
                                        const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(solution="<<solution.label()<<", dt="<<dt<<")");

    const PetscInt numConstants = 3;
    PetscScalar constants[3];
    for (int i = 0; i < 3; ++i) {
        constants[i] = _upDir[i];
    } // for

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err); assert(prob);
    err = PetscDSSetConstants(prob, numConstants, constants); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;

    // Put up direction into constants.
} // _setFEConstants


// End of file