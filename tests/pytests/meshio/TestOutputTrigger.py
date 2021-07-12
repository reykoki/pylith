#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
# @file tests/pytests/meshio/TestOutputTrigger.py
#
# @brief Unit testing of Python OutputTrigger object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.meshio.OutputTrigger import OutputTrigger


class TestOutputTrigger(TestAbstractComponent):
    """Unit testing of OutputTrigger object.
    """
    _class = OutputTrigger


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputTrigger))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
