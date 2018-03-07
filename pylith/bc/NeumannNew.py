#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/bc/NeumannNew.py
##
# @brief Python object for managing a Neummann (natural) boundary condition.
##
# Factory: boundary_condition

from .BoundaryConditionNew import BoundaryConditionNew
from pylith.feassemble.IntegratorPointwise import IntegratorPointwise

# NeumannNew class


class NeumannNew(BoundaryConditionNew,
                   IntegratorPointwise):
    """
    Python object for managing a Neumann (natural)
    boundary condition.

    Factory: boundary_condition
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="neumannnew"):
        """
        Constructor.
        """
        BoundaryConditionNew.__init__(self, name)
        IntegratorPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        BoundaryConditionNew.preinitialize(self, mesh)
        IntegratorPointwise.preinitialize(self, mesh)
        return

    def finalize(self):
        """
        Cleanup after running problem.
        """
        print(":TODO: @brad Implement NeumannNew.finalize() once output manager is added.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            BoundaryConditionNew._configure(self)
            IntegratorPointwise._configure(self)
        except ValueError as err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring Dirichlet boundary condition "
                             "(%s):\n%s" % (aliases, err.message))
        return

# End of file