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
 * @file libsrc/meshio/FieldFilterProject.hh
 *
 * @brief C++ object for projecting a field to another basis.
 */

#if !defined(pylith_meshio_fieldfilterproject_hh)
#define pylith_meshio_fieldfilterproject_hh

#include "FieldFilter.hh" // ISA FieldFilter

class pylith::meshio::FieldFilterProject : public FieldFilter {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    FieldFilterProject(void);

    /// Destructor
    ~FieldFilterProject(void);

    /** Set basis order for projected field.
     *
     * @param[in] value Basis order.
     */
    void setBasisOrder(const int value);

    /** Create DM associated with filtered field.
     *
     * @param[in] dm PETSc DM of unfiltered field.
     * @param[in] description Description of input field.
     * @param[in] discretization Discretization of input field.
     */
    PetscDM createDM(PetscDM dm,
                     const pylith::topology::FieldBase::Description& description,
                     const pylith::topology::FieldBase::Discretization& discretization) const;

    /** Apply filter to global PETSc vector.
     *
     * @param[out] vectorOut Filtered global PETSc vector.
     * @param[in] dmOut PETSc DM associated with filtered global PETSc vector.
     * @param[in] vectorIn PETSc global vector to filter.
     */
    void apply(PetscVec* vectorOut,
               PetscDM dmOut,
               PetscVec vectorIn) const;

    /** Pass thru solution in pointwise function.
     *
     * We pass through the solution to the resulting field. The auxiliary field is ignored.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] field [dim].
     */
    static
    void passThruSoln(const PylithInt dim,
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
                      PylithScalar field[]);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PetscPointFunc _passThruFns[1]; ///< Pass through point functions.
    int _basisOrder; ///< Basis order for projected field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    FieldFilterProject(const FieldFilterProject& f); ///< Not implemented.
    const FieldFilterProject& operator=(const FieldFilterProject&); ///< Not implemented.

}; // FieldFilterProject

#endif // pylith_meshio_fieldfilterproject_hh

// End of file
