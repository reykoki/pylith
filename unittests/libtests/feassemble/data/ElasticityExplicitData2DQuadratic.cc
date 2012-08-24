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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticityexplicitapp.

#include "ElasticityExplicitData2DQuadratic.hh"

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_spaceDim = 2;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_cellDim = 2;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numVertices = 6;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numBasis = 6;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numQuadPts = 6;

const char* pylith::feassemble::ElasticityExplicitData2DQuadratic::_matType = "ElasticPlaneStrain";

const char* pylith::feassemble::ElasticityExplicitData2DQuadratic::_matDBFilename = "data/elasticplanestrain.spatialdb";

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitData2DQuadratic::_matLabel = "elastic strain 2-D";

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_dtStableExplicit =   7.34858411e-05;

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_gravityVec[] = {
  0.00000000e+00, -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_vertices[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  2.00000000e-01,
 -1.50000000e+00,  5.00000000e-01,
 -2.50000000e-01,  3.50000000e-01,
 -1.25000000e+00, -2.50000000e-01,
  0.00000000e+00, -4.00000000e-01,
};

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_cells[] = {
0,1,2,3,4,5,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_quadPts[] = {
 -7.50000000e-01, -7.50000000e-01,
  7.50000000e-01, -7.50000000e-01,
 -7.50000000e-01,  7.50000000e-01,
  0.00000000e+00, -7.50000000e-01,
 -7.50000000e-01,  0.00000000e+00,
  2.50000000e-01,  2.50000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_quadWts[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_basis[] = {
  3.75000000e-01, -9.37500000e-02,
 -9.37500000e-02,  6.25000000e-02,
  3.75000000e-01,  3.75000000e-01,
  0.00000000e+00,  6.56250000e-01,
 -9.37500000e-02,  4.37500000e-01,
 -0.00000000e+00, -0.00000000e+00,
  0.00000000e+00, -9.37500000e-02,
  6.56250000e-01,  4.37500000e-01,
 -0.00000000e+00, -0.00000000e+00,
 -9.37500000e-02,  0.00000000e+00,
 -9.37500000e-02,  2.50000000e-01,
  1.87500000e-01,  7.50000000e-01,
 -9.37500000e-02, -9.37500000e-02,
  0.00000000e+00,  2.50000000e-01,
  7.50000000e-01,  1.87500000e-01,
  3.75000000e-01,  1.56250000e-01,
  1.56250000e-01,  1.56250000e+00,
 -6.25000000e-01, -6.25000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_basisDerivRef[] = {
 -1.00000000e+00, -1.00000000e+00,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  2.50000000e-01,
 -2.50000000e-01,  1.25000000e+00,
  1.25000000e+00, -2.50000000e-01,
  5.00000000e-01,  5.00000000e-01,
  1.25000000e+00,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  1.75000000e+00,
 -2.50000000e-01, -2.50000000e-01,
 -1.75000000e+00, -1.75000000e+00,
  5.00000000e-01,  5.00000000e-01,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00,  1.25000000e+00,
  1.75000000e+00,  2.50000000e-01,
 -1.75000000e+00, -1.75000000e+00,
 -2.50000000e-01, -2.50000000e-01,
 -2.50000000e-01, -2.50000000e-01,
  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  1.00000000e+00,
 -2.50000000e-01,  5.00000000e-01,
 -2.50000000e-01, -1.00000000e+00,
 -2.50000000e-01, -2.50000000e-01,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,
  1.00000000e+00,  2.50000000e-01,
 -1.00000000e+00, -2.50000000e-01,
  5.00000000e-01, -2.50000000e-01,
  1.00000000e+00,  1.00000000e+00,
  7.50000000e-01,  0.00000000e+00,
  0.00000000e+00,  7.50000000e-01,
  1.25000000e+00,  1.25000000e+00,
 -1.25000000e+00, -1.75000000e+00,
 -1.75000000e+00, -1.25000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_fieldTIncr[] = {
 -4.00000000e-01, -6.00000000e-01,
  7.00000000e-01,  8.00000000e-01,
  0.00000000e+00,  2.00000000e-01,
 -5.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  9.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_fieldT[] = {
 -3.00000000e-01, -4.00000000e-01,
  5.00000000e-01,  6.00000000e-01,
  0.00000000e+00,  1.00000000e-01,
 -2.00000000e-01, -3.00000000e-01,
  2.00000000e-01,  3.00000000e-01,
 -1.00000000e-01, -2.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_fieldTmdt[] = {
 -2.00000000e-01, -3.00000000e-01,
  3.00000000e-01,  4.00000000e-01,
  0.00000000e+00, -1.00000000e-01,
 -3.00000000e-01, -2.00000000e-01,
  1.00000000e-01,  4.00000000e-01,
 -2.00000000e-01, -6.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsResidual[] = {
  1.14469580e+10,  9.89206525e+10,
 -3.22261697e+10, -5.24164966e+09,
  5.33382294e+09,  6.79748894e+10,
  6.97040375e+10,  1.10679118e+11,
 -7.00903642e+10, -2.03315152e+11,
  1.58455201e+10, -6.90098433e+10,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsJacobian[] = {
  2.24121094e+06,  0.00000000e+00,
  2.41699219e+05,  0.00000000e+00,
  2.41699219e+05,  0.00000000e+00,
  4.21875000e+06,  0.00000000e+00,
 -1.36230469e+06,  0.00000000e+00,
 -1.36230469e+06,  0.00000000e+00,
  0.00000000e+00,  2.24121094e+06,
  0.00000000e+00,  2.41699219e+05,
  0.00000000e+00,  2.41699219e+05,
  0.00000000e+00,  4.21875000e+06,
  0.00000000e+00, -1.36230469e+06,
  0.00000000e+00, -1.36230469e+06,
  2.41699219e+05,  0.00000000e+00,
  3.61083984e+06,  0.00000000e+00,
 -6.73828125e+05,  0.00000000e+00,
  3.45703125e+06,  0.00000000e+00,
 -1.52343750e+06,  0.00000000e+00,
 -1.12792969e+06,  0.00000000e+00,
  0.00000000e+00,  2.41699219e+05,
  0.00000000e+00,  3.61083984e+06,
  0.00000000e+00, -6.73828125e+05,
  0.00000000e+00,  3.45703125e+06,
  0.00000000e+00, -1.52343750e+06,
  0.00000000e+00, -1.12792969e+06,
  2.41699219e+05,  0.00000000e+00,
 -6.73828125e+05,  0.00000000e+00,
  3.61083984e+06,  0.00000000e+00,
  3.45703125e+06,  0.00000000e+00,
 -1.12792969e+06,  0.00000000e+00,
 -1.52343750e+06,  0.00000000e+00,
  0.00000000e+00,  2.41699219e+05,
  0.00000000e+00, -6.73828125e+05,
  0.00000000e+00,  3.61083984e+06,
  0.00000000e+00,  3.45703125e+06,
  0.00000000e+00, -1.12792969e+06,
  0.00000000e+00, -1.52343750e+06,
  4.21875000e+06,  0.00000000e+00,
  3.45703125e+06,  0.00000000e+00,
  3.45703125e+06,  0.00000000e+00,
  2.21484375e+07,  0.00000000e+00,
 -5.39062500e+06,  0.00000000e+00,
 -5.39062500e+06,  0.00000000e+00,
  0.00000000e+00,  4.21875000e+06,
  0.00000000e+00,  3.45703125e+06,
  0.00000000e+00,  3.45703125e+06,
  0.00000000e+00,  2.21484375e+07,
  0.00000000e+00, -5.39062500e+06,
  0.00000000e+00, -5.39062500e+06,
 -1.36230469e+06,  0.00000000e+00,
 -1.52343750e+06,  0.00000000e+00,
 -1.12792969e+06,  0.00000000e+00,
 -5.39062500e+06,  0.00000000e+00,
  8.46679688e+06,  0.00000000e+00,
  6.09375000e+06,  0.00000000e+00,
  0.00000000e+00, -1.36230469e+06,
  0.00000000e+00, -1.52343750e+06,
  0.00000000e+00, -1.12792969e+06,
  0.00000000e+00, -5.39062500e+06,
  0.00000000e+00,  8.46679688e+06,
  0.00000000e+00,  6.09375000e+06,
 -1.36230469e+06,  0.00000000e+00,
 -1.12792969e+06,  0.00000000e+00,
 -1.52343750e+06,  0.00000000e+00,
 -5.39062500e+06,  0.00000000e+00,
  6.09375000e+06,  0.00000000e+00,
  8.46679688e+06,  0.00000000e+00,
  0.00000000e+00, -1.36230469e+06,
  0.00000000e+00, -1.12792969e+06,
  0.00000000e+00, -1.52343750e+06,
  0.00000000e+00, -5.39062500e+06,
  0.00000000e+00,  6.09375000e+06,
  0.00000000e+00,  8.46679688e+06,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsResidualLumped[] = {
  1.14454133e+10,  9.89209294e+10,
 -3.22283567e+10, -5.24308887e+09,
  5.33172309e+09,  6.79741797e+10,
  6.97057895e+10,  1.10680806e+11,
 -7.00892582e+10, -2.03318379e+11,
  1.58484937e+10, -6.90064316e+10,
};

const PylithScalar pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsJacobianLumped[] = {
  4.21875000e+06,  4.21875000e+06,
  3.98437500e+06,  3.98437500e+06,
  3.98437500e+06,  3.98437500e+06,
  2.25000000e+07,  2.25000000e+07,
  5.15625000e+06,  5.15625000e+06,
  5.15625000e+06,  5.15625000e+06,
};

pylith::feassemble::ElasticityExplicitData2DQuadratic::ElasticityExplicitData2DQuadratic(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  dtStableExplicit = _dtStableExplicit;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
  valsResidualLumped = const_cast<PylithScalar*>(_valsResidualLumped);
  valsJacobianLumped = const_cast<PylithScalar*>(_valsJacobianLumped);
} // constructor

pylith::feassemble::ElasticityExplicitData2DQuadratic::~ElasticityExplicitData2DQuadratic(void)
{}


// End of file
