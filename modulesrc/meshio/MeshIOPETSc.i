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
 * @file modulesrc/meshio/MeshIOPETSc.i
 *
 * @brief Python interface to C++ MeshIOPETSc object.
 */

namespace pylith {
  namespace meshio {

    class MeshIOPETSc : public MeshIO
    { // MeshIOPETSc

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      MeshIOPETSc(void);

      /// Destructor
      ~MeshIOPETSc(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for ASCII file.
       *
       * @param filename Name of file
       */
      void filename(const char* name);
      
      /** Get filename of ASCII file.
       *
       * @returns Name of file
       */
      const char* filename(void) const;

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /// Write mesh
      void _write(void) const;
      
      /// Read mesh
      void _read(void);

    }; // MeshIOPETSc

  } // meshio
} // pylith


// End of file 
