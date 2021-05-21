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

/**
 * @file libsrc/meshio/FieldFilter.hh
 *
 * @brief C++ object for filtering field upon output, etc.
 */

#if !defined(pylith_meshio_fieldfilter_hh)
#define pylith_meshio_fieldfilter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase
#include "pylith/utils/petscfwd.h" // USES PetscDM, PetscVec

// FieldFilter ---------------------------------------------------------
/** @brief C++ object for filtering field upon output, etc.
 */
class pylith::meshio::FieldFilter : public pylith::utils::PyreComponent {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    FieldFilter(void);

    /// Destructor
    virtual ~FieldFilter(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Create DM associated with filtered field.
     *
     * @param[in] dm PETSc DM of unfiltered field.
     * @param[in] description Description of input field.
     * @param[in] discretization Discretization of input field.
     */
    virtual
    PetscDM createDM(PetscDM dm,
                     const pylith::topology::FieldBase::Description& description,
                     const pylith::topology::FieldBase::Discretization& discretization) const = 0;

    /** Apply filter to global PETSc vector.
     *
     * @param[out] vectorOut Filtered global PETSc vector.
     * @param[in] dmOut PETSc DM associated with filtered global PETSc vector.
     * @param[in] vectorIn PETSc global vector to filter.
     */
    virtual
    void apply(PetscVec* vectorOut,
               PetscDM dmOut,
               PetscVec vectorIn) const = 0;

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

private:

    FieldFilter(const FieldFilter&); ///< Not implemented
    const FieldFilter& operator=(const FieldFilter&); ///< Not implemented.

}; // FieldFilter

#endif // pylith_meshio_fieldfilter_hh

// End of file
