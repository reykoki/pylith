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

#include "Fields.hh" // Implementation of class methods

#include "Field.hh" // USES Field

#include <pylith/utils/error.hh> // USES PYLITH_CHECK_ERROR

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Fields::Fields(const Mesh& mesh) :
    _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Fields::~Fields(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Fields::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    const map_type::iterator begin = _fields.begin();
    const map_type::iterator end = _fields.end();
    for (map_type::iterator iter=begin; iter != end; ++iter) {
        delete iter->second; iter->second = NULL;
    } // for
    _fields.clear();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Check if fields contains a given field.
bool
pylith::topology::Fields::hasField(const char* name) const
{ // hasField
    PYLITH_METHOD_BEGIN;

    map_type::const_iterator iter = _fields.find(name);

    PYLITH_METHOD_RETURN(iter != _fields.end());
} // hasField

// ----------------------------------------------------------------------
// Add field.
void
pylith::topology::Fields::add(const char* name,
                              const char* label)
{ // add
    PYLITH_METHOD_BEGIN;

    if (hasField(name)) {
        std::ostringstream msg;
        msg << "Could not add field '" << name << "' to fields manager, because it already exists.";
        throw std::runtime_error(msg.str());
    } // if

    _fields[name] = new Field(_mesh);
    _fields[name]->label(label);

    PYLITH_METHOD_END;
} // add

// ----------------------------------------------------------------------
// Add field.
void
pylith::topology::Fields::add(const char* name,
                              const char* label,
                              const pylith::topology::FieldBase::DomainEnum domain,
                              const int fiberDim)
{ // add
    PYLITH_METHOD_BEGIN;

    if (hasField(name)) {
        std::ostringstream msg;
        msg << "Could not add field '" << name << "' to fields manager, because it already exists.";
        throw std::runtime_error(msg.str());
    } // if

    _fields[name] = new Field(_mesh);
    _fields[name]->label(label);
    _fields[name]->newSection(domain, fiberDim);

    PYLITH_METHOD_END;
} // add

// ----------------------------------------------------------------------
// Delete field.
void
pylith::topology::Fields::del(const char* name)
{ // del
    PYLITH_METHOD_BEGIN;

    map_type::iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager to delete.";
        throw std::runtime_error(msg.str());
    } // if
    delete iter->second; iter->second = 0;
    _fields.erase(name);

    PYLITH_METHOD_END;
} // del

// ----------------------------------------------------------------------
// Get field.
const pylith::topology::Field&
pylith::topology::Fields::get(const char* name) const
{ // get
    PYLITH_METHOD_BEGIN;

    map_type::const_iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager for retrieval.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(*iter->second);
} // get

// ----------------------------------------------------------------------
// Get field.
pylith::topology::Field&
pylith::topology::Fields::get(const char* name)
{ // get
    PYLITH_METHOD_BEGIN;

    map_type::iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager for retrieval.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(*iter->second);
} // get

// ----------------------------------------------------------------------
// Copy layout to other fields.
void
pylith::topology::Fields::copyLayout(const char* name)
{ // copyLayout
    PYLITH_METHOD_BEGIN;

    map_type::const_iterator src = _fields.find(name);
    if (src == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager for retrieval.";
        throw std::runtime_error(msg.str());
    } // if

    const map_type::iterator begin = _fields.begin();
    const map_type::iterator end = _fields.end();
    for (map_type::iterator iter=begin; iter != end; ++iter) {
        if (iter != src) {
            iter->second->cloneSection(*src->second);
        } // if
    } // for

    PYLITH_METHOD_END;
} // copyLayout

// ----------------------------------------------------------------------
// Get mesh associated with fields.
const pylith::topology::Mesh&
pylith::topology::Fields::mesh(void) const
{ // mesh
    return _mesh;
} // mesh


// End of file
