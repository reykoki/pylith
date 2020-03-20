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

#include <portinfo>

#include "TestDataWriterMesh.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cppunit/extensions/HelperMacros.h>

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterMesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _mesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterMesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterMesh::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriter_Data* data = _getData();CPPUNIT_ASSERT(data);

    delete _mesh;_mesh = new topology::Mesh;CPPUNIT_ASSERT(_mesh);
    MeshIOAscii iohandler;
    iohandler.filename(data->meshFilename);
    iohandler.read(_mesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->dimension());
    _mesh->setCoordSys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_mesh, normalizer);

    if (data->faultLabel) {
	pylith::faults::FaultCohesiveStub fault;
        fault.setInterfaceId(data->faultId);
        fault.setSurfaceMarkerLabel(data->faultLabel);
        fault.adjustTopology(_mesh);
    } // if

    PYLITH_METHOD_END;
} // _initialize


// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMesh::_createVertexFields(pylith::topology::Fields* fields) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(fields);

    const TestDataWriter_Data* data = _getData();CPPUNIT_ASSERT(data);

    FieldFactory factory(*fields);
    factory.scalar(data->vertexDiscretization, data->vertexScalarValues, data->vertexNumPoints, data->vertexScalarNumComponents);
    factory.vector(data->vertexDiscretization, data->vertexVectorValues, data->vertexNumPoints, data->vertexVectorNumComponents);
    factory.tensor(data->vertexDiscretization, data->vertexTensorValues, data->vertexNumPoints, data->vertexTensorNumComponents);
    factory.other(data->vertexDiscretization, data->vertexOtherValues, data->vertexNumPoints, data->vertexOtherNumComponents);

    PYLITH_METHOD_END;
} // _createVertexFields


// ----------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterMesh::_createCellFields(pylith::topology::Fields* fields) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(fields);

    const TestDataWriter_Data* data = _getData();CPPUNIT_ASSERT(data);

    FieldFactory factory(*fields);
    factory.scalar(data->cellDiscretization, data->cellScalarValues, data->cellNumPoints, data->cellScalarNumComponents);
    factory.vector(data->cellDiscretization, data->cellVectorValues, data->cellNumPoints, data->cellVectorNumComponents);
    factory.tensor(data->cellDiscretization, data->cellTensorValues, data->cellNumPoints, data->cellTensorNumComponents);
    factory.other(data->cellDiscretization, data->cellOtherValues, data->cellNumPoints, data->cellOtherNumComponents);

    PYLITH_METHOD_END;
} // _createCellFields


// ----------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::_setDataTri(void) {
    TestDataWriter_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tri3.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------
    data->vertexNumPoints = 6;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[6*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6,
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 2;
    static const PylithScalar vertexVectorValues[6*2] = {
        1.1, 2.2,
        3.3, 4.4,
        5.5, 6.6,
        7.7, 8.8,
        9.9, 10.0,
        11.1, 12.2,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 3;
    static const PylithScalar vertexTensorValues[6*3] = {
        1.1, 1.2, 1.3,
        2.1, 2.2, 2.3,
        3.1, 3.2, 3.3,
        4.1, 4.2, 4.3,
        5.1, 5.2, 5.3,
        6.1, 6.2, 6.3,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[6*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        9.0, 10.1,
        11.2, 12.3,
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[2*1] = {
        2.1, 2.2,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 2;
    static const PylithScalar cellVectorValues[2*2] = {
        1.1, 2.2,
        3.3, 4.4,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 3;
    static const PylithScalar cellTensorValues[2*3] = {
        1.2, 2.3, 3.4,
        4.5, 5.6, 6.7,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[2*2] = {
        1.2, 2.3,
        4.5, 5.6,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataTri


// ----------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::_setDataQuad(void) {
    TestDataWriter_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    // We do not use a fault in this test case.
    data->meshFilename = "data/quad4.mesh";
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------
    data->vertexNumPoints = 6;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[6*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 2;
    static const PylithScalar vertexVectorValues[6*2] = {
        1.1, 2.2,
        3.3, 4.4,
        5.5, 6.6,
        7.7, 8.8,
        9.9, 10.1,
        11.2, 12.3,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 3;
    static const PylithScalar vertexTensorValues[6*3] = {
        1.1, 1.2, 1.3,
        2.1, 2.2, 3.3,
        3.1, 3.2, 4.3,
        4.1, 4.2, 5.3,
        5.1, 5.2, 6.3,
        6.1, 6.2, 7.3,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[6*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        9.8, 7.6,
        6.5, 5.4
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[2*1] = {
        2.1, 2.2,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 2;
    static const PylithScalar cellVectorValues[2*2] = {
        1.1, 2.2,
        3.3, 4.4,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 3;
    static const PylithScalar cellTensorValues[2*3] = {
        1.2, 2.3, 3.4,
        4.5, 5.6, 6.7,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[2*2] = {
        1.2, 2.3,
        4.5, 5.6,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataQuad


// ----------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::_setDataTet(void) {
    TestDataWriter_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tet4.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------
    data->vertexNumPoints = 8;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[8*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 3;
    static const PylithScalar vertexVectorValues[8*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
        7.7, 8.8, 9.9,
        10.0, 11.1, 12.2,
        13.3, 14.4, 15.5,
        16.6, 17.7, 18.8,
        19.9, 20.0, 21.1,
        22.2, 23.3, 24.4,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 6;
    static const PylithScalar vertexTensorValues[8*6] = {
        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
        3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
        4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
        5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
        6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
        7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
        8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[8*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        9.0, 10.1,
        11.2, 12.3,
        13.4, 14.5,
        15.6, 16.7,
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[2*1] = {
        2.1, 3.2,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 3;
    static const PylithScalar cellVectorValues[2*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 6;
    static const PylithScalar cellTensorValues[2*6] = {
        1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
        7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[2*2] = {
        1.2, 2.3,
        7.8, 8.9,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataTet


// ----------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::_setDataHex(void) {
    TestDataWriter_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/hex8.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------
    data->vertexNumPoints = 16;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[16*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
        10.0, 12.1, 11.1, 13.1, 14.1, 15.1, 16.1, 17.1,
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 3;
    static const PylithScalar vertexVectorValues[16*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
        7.7, 8.8, 9.9,
        10.1, 11.2, 12.3,
        1.2, 2.3, 3.4,
        4.5, 5.6, 6.7,
        7.8, 8.9, 9.0,
        10.2, 11.3, 12.4,
        1.3, 2.4, 3.5,
        4.6, 5.7, 6.8,
        7.9, 8.0, 9.1,
        10.2, 11.3, 12.4,
        13.5, 14.6, 15.7,
        16.8, 17.9, 18.1,
        19.2, 20.3, 21.4,
        22.5, 23.6, 24.7,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 6;
    static const PylithScalar vertexTensorValues[16*6] = {
        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
        3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
        4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
        5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
        6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
        7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
        8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
        9.1, 9.2, 9.3, 9.4, 9.5, 9.6,
        10.1, 10.2, 10.3, 10.4, 10.5, 10.6,
        11.1, 11.2, 11.3, 11.4, 11.5, 11.6,
        12.1, 12.2, 12.3, 12.4, 12.5, 12.6,
        13.1, 13.2, 13.3, 13.4, 13.5, 13.6,
        14.1, 14.2, 14.3, 14.4, 14.5, 14.6,
        15.1, 15.2, 15.3, 15.4, 15.5, 15.6,
        16.1, 16.2, 16.3, 16.4, 16.5, 16.6,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[16*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        1.3, 2.4,
        3.5, 4.6,
        5.7, 6.8,
        7.9, 8.0,
        1.3, 2.4,
        3.5, 4.6,
        5.7, 6.8,
        8.0, 1.4,
        2.5, 3.6,
        4.8, 1.5,
        2.6, 3.7,
        4.8, 5.9,
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[2*1] = {
        2.1, 3.2,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 3;
    static const PylithScalar cellVectorValues[2*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 6;
    static const PylithScalar cellTensorValues[2*6] = {
        1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
        1.1, 2.2, 3.3, 4.4, 5.5, 6.6,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[2*2] = {
        1.2, 2.3,
        1.1, 2.2,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataHex


// End of file
