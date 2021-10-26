from PYB11Generator3 import *
from Tessellator import *
from TessellatorCommonMethods import *

@PYB11template("int Dimension", "RealType")
class DistributedTessellator(Tessellator):
    """DistributedTessellator

Provides a parallel tessellation.

This class assumes the user provides:
1.  A serial Tessellator.
2.  The generators in parallel, distributed in any manner the user likes
    so long as the generators are not degenerate: i.e., don't repeast the
    same generator on different domains.

Based on the parallel tessellation algorithm originally implemented in 
Spheral++.
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               serialTessellator = "Tessellator<%(Dimension)s, %(RealType)s>*",
               assumeControl = ("bool", "true"),
               buildCommunicationInfo = ("bool", "false")):
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def handlesPLCs(self):
        """Override this method to return true if this Tessellator supports 
the description of a domain boundary using a PLC (as in the second 
tessellate method, above), and false if it does not. Some algorithms 
for tessellation do not naturally accommodate an explicit boundary 
description, and Tessellators using these algorithms should override 
this method to return false. A stub method for PLC-enabled
tessellation is provided for convenience.
This query mechanism prevents us from descending into the taxonomic 
hell associated with elaborate inheritance hierarchies."""
        return "bool"

    @PYB11virtual
    @PYB11const
    def name(self):
        "A unique name string per tessellation instance."
        return "std::string"

    @PYB11virtual
    @PYB11const
    def degeneracy(self):
        """Returns the accuracy to which this tessellator can distinguish coordinates.
Should be returned appropriately for normalized coordinates, i.e., if all
coordinates are in the range xi \\\\in [0,1], what is the minimum allowed 
delta in x."""
        return "%(RealType)s"

#-------------------------------------------------------------------------------
# Inject the common methods
PYB11inject(TessellatorCommonMethods, DistributedTessellator)

#-------------------------------------------------------------------------------
# Template instantiations
DistributedTessellator2d = PYB11TemplateClass(DistributedTessellator, template_parameters=("2", "double"))
DistributedTessellator3d = PYB11TemplateClass(DistributedTessellator, template_parameters=("3", "double"))
