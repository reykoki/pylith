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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSubfield.hh
 *
 * @brief Manager for extracting subfield for output from a PyLith field.
 */

#if !defined(pylith_meshio_outputsubfield_hh)
#define pylith_meshio_outputsubfield_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA PyreComponent

#include "pylith/topology/FieldBase.hh" // HASA Description, Discretization

#include "pylith/topology/topologyfwd.hh" // USES Field, FieldFilter
#include "pylith/utils/petscfwd.h" // HASA PetscVec

class pylith::meshio::OutputSubfield : public pylith::utils::GenericComponent {
    friend class TestOutputSubfield; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Create OutputSubfield from Field.
     *
     * @param[in] field Field with subfields.
     * @param[in] name Name of subfield that will be extracted.
     * @param[in] filter Filter to apply to subfield.
     */
    static
    OutputSubfield* create(const pylith::topology::Field& field,
                           const char* name,
                           const pylith::meshio::FieldFilter* filter);

    /// Destructor
    ~OutputSubfield(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get description of subfield.
     *
     * @returns Description of subfield.
     */
    const pylith::topology::FieldBase::Description& getDescription(void) const;

    /** Get filtered PETSc global vector.
     *
     * @returns PETSc global vector.
     */
    PetscVec getVector(void) const;

    /** Get PETSc DM for filtered vector.
     *
     * @returns PETSc DM.
     */
    PetscDM getDM(void) const;

    /** Extract subfield data from global PETSc vector with subfields.
     *
     * @param[in] fieldVector PETSc vector with subfields.
     */
    void extract(const PetscVec& fieldVector);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /** Constructor.
     *
     * @param[in] description Description of subfield.
     * @param[in] discretization Discretization of subfield.
     * @param[in] dm PETSc DM for subfield.
     * @param[in] is PETSc IS for subfield.
     * @param[in] filter Filter to apply to subfield.
     */
    OutputSubfield(const pylith::topology::FieldBase::Description& description,
                   const pylith::topology::FieldBase::Discretization& discretization,
                   PetscDM dm,
                   PetscIS is,
                   const pylith::meshio::FieldFilter* filter);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::FieldBase::Description _description; ///< Description of subfield.
    pylith::topology::FieldBase::Discretization _discretization; ///< Discretization of subfield.
    PetscDM _dm; ///< PETSc DM for subfield.
    PetscIS _is; ///< PETSc IS for subfield.
    PetscVec _vector; ///< PETSc global vector for subfield.

    const pylith::meshio::FieldFilter* _filter; ///< Filter to apply to field.
    PetscDM _filteredDM; ///< PETSc DM for filtered vector.
    PetscVec _filteredVector; ///< PETSc global vector after filtering.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    OutputSubfield(const OutputSubfield&); ///< Not implemented.
    const OutputSubfield& operator=(const OutputSubfield&); ///< Not implemented

}; // OutputSubfield

#endif // pylith_meshio_outputsubfield_hh

// End of file
