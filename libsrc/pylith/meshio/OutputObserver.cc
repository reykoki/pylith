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

#include <portinfo>

#include "OutputObserver.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/FieldFilter.hh" // USES FieldFilter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputObserver::OutputObserver(void) :
    _timeScale(1.0),
    _writer(NULL),
    _fieldFilter(NULL),
    _trigger(NULL)
{}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputObserver::~OutputObserver(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputObserver::deallocate(void) {
    if (_writer) {
        _writer->close();
        _writer->deallocate();
    }
    if (_fieldFilter) { _fieldFilter->deallocate(); }

    typedef std::map<std::string, OutputSubfield*> subfield_t;
    for (subfield_t::iterator iter = _subfields.begin(); iter != _subfields.end(); ++iter) {
        delete iter->second;iter->second = NULL;
    } // for
    _subfields.clear();

    _writer = NULL; // :TODO: Use shared pointer
    _fieldFilter = NULL; // :TODO: Use shared pointer
    _trigger = NULL; // :TODO: Use shared pointer
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set trigger for how often to write output.
void
pylith::meshio::OutputObserver::setTrigger(pylith::meshio::OutputTrigger* const trigger) {
    PYLITH_COMPONENT_DEBUG("OutputObserver::setTrigger(otrigger="<<typeid(trigger).name()<<")");

    _trigger = trigger;
} // setTrigger


// ------------------------------------------------------------------------------------------------
// Get trigger for how often to write otuput.
const pylith::meshio::OutputTrigger*
pylith::meshio::OutputObserver::getTrigger(void) const {
    return _trigger;
} // getTrigger


// ------------------------------------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputObserver::setWriter(DataWriter* const writer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setWrite(datawriter="<<typeid(writer).name()<<")");

    _writer = writer; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // setWriter


// ------------------------------------------------------------------------------------------------
// Set filter for fields.
void
pylith::meshio::OutputObserver::setFieldFilter(FieldFilter* const filter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputObserver::setFieldFilter(filter="<<typeid(filter).name()<<")");

    _fieldFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // setFieldFilter


// ------------------------------------------------------------------------------------------------
// Set time scale.
void
pylith::meshio::OutputObserver::setTimeScale(const PylithReal value) {
    if (value <= 0.0) {
        std::ostringstream msg;
        msg << "Time scale ("<<value<<") for output observer is nonpositive.";
        throw std::logic_error(msg.str());
    } // if
    _timeScale = value;
} // setTimeScale


// ------------------------------------------------------------------------------------------------
// Get output subfield, creating if necessary.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputObserver::_getSubfield(const pylith::topology::Field& field,
                                             const char* name) {
    if (_subfields.count(name) == 0) {
        _subfields[name] = OutputSubfield::create(field, name, _fieldFilter);
    } // if
    return _subfields[name];
}


// ---------------------------------------------------------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputObserver::_appendField(const PylithReal t,
                                             const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_appendField(t="<<t<<", subfield="<<typeid(subfield).name()<<")");

    const int basisOrder = pylith::topology::FieldOps::getBasisOrder(subfield.getDM());
    switch (basisOrder) {
    case 0:
        _writer->writeCellField(t, subfield);
        break;

    case 1:
        _writer->writeVertexField(t, subfield);
        break;

    default:
        PYLITH_COMPONENT_ERROR(
            "Unsupported basis order ("
                << basisOrder <<") for output. Use FieldFilterProject with basis order of 0 or 1. Skipping output of '"
                << subfield.getDescription().label << "' field."
            );
    } // switch

    PYLITH_METHOD_END;
} // _appendField


// End of file
