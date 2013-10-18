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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestReverseCuthillMcKee.hh" // Implementation of class methods

#include "pylith/topology/ReverseCuthillMcKee.hh" // USES ReverseCuthillMcKee

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Jacobian.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestReverseCuthillMcKee );

// ----------------------------------------------------------------------
// Test reorder() with tri3 cells and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTri3(void)
{ // testReorderTri3
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_tri3.mesh");

  PYLITH_METHOD_END;
} // testReorderTri3

// ----------------------------------------------------------------------
// Test reorder() with level 2, tri3 cells, and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTri3Fault(void)
{ // testReorderTri3Fault
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_tri3.mesh", "fault");

  PYLITH_METHOD_END;
} // testReorderTri3Fault

// ----------------------------------------------------------------------
// Test reorder() with level 2, quad4 cells, and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderQuad4(void)
{ // testReorderQuad4
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_quad4.mesh");

  PYLITH_METHOD_END;
} // testReorderQuad4

// ----------------------------------------------------------------------
// Test reorder() with level 2, quad4 cells, and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderQuad4Fault(void)
{ // testReorderQuad4Fault
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_quad4.mesh", "fault");

  PYLITH_METHOD_END;
} // testReorderQuad4Fault

// ----------------------------------------------------------------------
// Test reorder() with level 2, tet4 cells, and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTet4(void)
{ // testReorderTet4
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_tet4.mesh");

  PYLITH_METHOD_END;
} // testReorderTet4

// ----------------------------------------------------------------------
// Test reorder() with level 2, tet4 cells, and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTet4Fault(void)
{ // testReorderTet4Fault
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_tet4.mesh", "fault");

  PYLITH_METHOD_END;
} // testReorderTet4Fault

// ----------------------------------------------------------------------
// Test reorder() with hex8 cells and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderHex8(void)
{ // testReorderHex8
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_hex8.mesh");

  PYLITH_METHOD_END;
} // testReorderHex8

// ----------------------------------------------------------------------
// Test reorder() with hex8 cells and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderHex8Fault(void)
{ // testReorderHex8Fault
  PYLITH_METHOD_BEGIN;

  _testReorder("data/reorder_hex8.mesh", "fault");

  PYLITH_METHOD_END;
} // testReorderHex8Fault

// ----------------------------------------------------------------------
void
pylith::topology::TestReverseCuthillMcKee::_setupMesh(Mesh* const mesh,
						      const char* filename,
						      const char* faultGroup)
{ // _setupMesh
  PYLITH_METHOD_BEGIN;

  assert(mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.interpolate(true);

  iohandler.read(mesh);
  CPPUNIT_ASSERT(mesh->numCells() > 0);
  CPPUNIT_ASSERT(mesh->numVertices() > 0);

  // Adjust topology if necessary.
  if (faultGroup) {
    int firstLagrangeVertex = 0;
    int firstFaultCell = 0;

    faults::FaultCohesiveKin fault;
    fault.id(100);
    fault.label(faultGroup);
    const int nvertices = fault.numVerticesNoMesh(*mesh);
    firstLagrangeVertex += nvertices;
    firstFaultCell += 2*nvertices; // shadow + Lagrange vertices

    int firstFaultVertex = 0;
    fault.adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
  } // if

  PYLITH_METHOD_END;
} // _setupMesh

// ----------------------------------------------------------------------
// Test reorder().
void
pylith::topology::TestReverseCuthillMcKee::_testReorder(const char* filename,
							const char* faultGroup)
{ // _testReorder
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _setupMesh(&mesh, filename, faultGroup);
  const PetscDM dmOrig = mesh.dmMesh();
  PetscObjectReference((PetscObject) dmOrig);

  ReverseCuthillMcKee::reorder(&mesh);
  
  const PetscDM& dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices (size only)
  topology::Stratum verticesStratumE(dmOrig, topology::Stratum::DEPTH, 0);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  CPPUNIT_ASSERT_EQUAL(verticesStratumE.size(), verticesStratum.size());

  // Check cells (size only)
  topology::Stratum cellsStratumE(dmOrig, topology::Stratum::HEIGHT, 0);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();
  CPPUNIT_ASSERT_EQUAL(cellsStratumE.size(), cellsStratum.size());

  // Check groups
  PetscInt numGroupsE, numGroups, pStart, pEnd;
  PetscErrorCode err;
  err = DMPlexGetNumLabels(dmOrig, &numGroupsE);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
  { // DEBUGGING
    PetscErrorCode err;
    for (PetscInt iGroup = 0; iGroup < numGroupsE; ++iGroup) {
      const char *name = NULL;
      err = DMPlexGetLabelName(dmOrig, iGroup, &name);PYLITH_CHECK_ERROR(err);
      std::cout << "ORIG Group["<<iGroup<<"]: " << name << std::endl;
    } // for
    for (PetscInt iGroup = 0; iGroup < numGroups; ++iGroup) {
      const char *name = NULL;
      err = DMPlexGetLabelName(dmMesh, iGroup, &name);PYLITH_CHECK_ERROR(err);
      std::cout << "REORDER Group["<<iGroup<<"]: " << name << std::endl;
    } // for
  } // DEBUGGING
  CPPUNIT_ASSERT_EQUAL(numGroupsE, numGroups);

  for (PetscInt iGroup = 0; iGroup < numGroups; ++iGroup) {
    const char *name = NULL;
    err = DMPlexGetLabelName(dmMesh, iGroup, &name);PYLITH_CHECK_ERROR(err);

    PetscInt numPointsE, numPoints;
    err = DMPlexGetStratumSize(dmOrig, name, 1, &numPointsE);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetStratumSize(dmMesh, name, 1, &numPoints);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numPointsE, numPoints);
  } // for

  // Verify reduction in Jacobian bandwidth
  Field field(mesh);
  field.newSection(FieldBase::VERTICES_FIELD, mesh.dimension());
  field.allocate();
  field.zero();
  Jacobian jacobian(field);
  std::cout << "****TEST IS DISABLED. FORCING FAILURE.****" << std::endl;
  CPPUNIT_ASSERT(false);

  PYLITH_METHOD_END;
} // _testReorder


// End of file 
