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
 * @file libsrc/meshio/MeshIOPETSc.hh
 *
 * @brief C++ input/output manager for PyLith PETSc gmsh mesh files.
 */

#if !defined(pylith_meshio_meshiopetsc_hh)
#define pylith_meshio_meshiopetsc_hh

#include "MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

class pylith::meshio::MeshIOPETSc : public MeshIO {
    friend class TestMeshIOPETSc; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOPETSc(void);

    /// Destructor
    ~MeshIOPETSc(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);


    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    bool _useIndexZero; ///< Flag indicating if indicates start at 0 (T) or 1 (F)

}; // MeshIOPETSc


#endif // pylith_meshio_meshiopetsc_hh

// End of file
