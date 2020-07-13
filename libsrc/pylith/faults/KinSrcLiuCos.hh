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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/KinSrcLiuCos.hh
 *
 * @brief C++ implementation of a Liu cosine-sine slip time function.
 */

#if !defined(pylith_faults_kinsrcliucos_hh)
#define pylith_faults_kinsrcliucos_hh

// Include directives ---------------------------------------------------
#include "KinSrc.hh"

// KinSrcLiuCos ------------------------------------------------------
/** @brief Slip function time history from Liu, Archuleta, and Hartzell,
 * BSSA, 2006 (doi:10.1785/0120060036) which has a rapid rise and then
 * a gradual falloff with a finite duration.
 */
class pylith::faults::KinSrcLiuCos : public KinSrc {
    friend class TestKinSrcLiuCos; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcLiuCos(void);

    /// Destructor.
    ~KinSrcLiuCos(void);

    /** Slip time function kernel.
     *
     * The "solution" field s is ignored.
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
     * @param[out] slip [dim].
     */
    static
    void slipFn(const PylithInt dim,
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
                PylithScalar slip[]);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                        const spatialdata::geocoords::CoordSys* cs);

    /** Set kernel for slip time function.
     *
     * @param[in] auxField Auxiliary field for fault with prescribed slip.
     */
    void _setSlipFnKernel(const pylith::topology::Field& auxField) const;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    KinSrcLiuCos(const KinSrcLiuCos&); ///< Not implemented
    const KinSrcLiuCos& operator=(const KinSrcLiuCos&); ///< Not implemented

}; // class KinSrcLiuCos

#endif // pylith_faults_kinsrcliucos_hh

// End of file
